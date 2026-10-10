# =============================================================================
# 08_partial_responses.R
# Marginal (partial) response curves for the deployed ensemble, per species.
#
# For each of the 8 predictors, hold the other 7 at their median across the
# training presences, sweep the focal predictor across its observed range,
# and predict ensemble suitability. Rug plot of actual presence values
# overlaid.
#
# Output:
#   PARTIAL_DIR/<sid>.png    2x4 panel (one per predictor)
#
# NOTE: This step is SLOW (each species runs 8 mini-projections through
# biomod2). RUN_PARTIAL is FALSE by default in RUN_ALL.R. Turn it on only
# after you've decided the ensembles are worth interpreting deeply.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

log_line("partial", "-", "starting")

orig_wd     <- getwd()
biomod_root <- normalizePath(MODELS_DIR, winslash = "/", mustWork = FALSE)

GRID_N <- 100   # points along each predictor sweep

load_biomod_out <- function(path) {
  env <- new.env()
  loaded <- tryCatch(load(path, envir = env), error = function(e) character(0))
  if (length(loaded) == 0) return(NULL)
  env[[loaded[1]]]
}

load_species_biomod <- function(sid) {
  for (d_name in unique(c(gsub("_", ".", sid), sid))) {
    sp_dir <- file.path(biomod_root, d_name)
    if (!dir.exists(sp_dir)) next
    all_out  <- list.files(sp_dir, pattern = "\\.out$", full.names = TRUE)
    em_files <- grep("ensemble\\.models\\.out$", all_out, value = TRUE)
    bm_files <- setdiff(grep("\\.models\\.out$", all_out, value = TRUE), em_files)
    if (!length(bm_files) || !length(em_files)) next
    bm_out <- load_biomod_out(bm_files[order(file.mtime(bm_files), decreasing = TRUE)][1])
    em_out <- load_biomod_out(em_files[order(file.mtime(em_files), decreasing = TRUE)][1])
    if (is.null(bm_out) || is.null(em_out)) next
    return(list(bm_out = bm_out, em_out = em_out, sp_dir = sp_dir, biomod_sid = d_name))
  }
  NULL
}

extract_deployed <- function(r_em_proj) {
  nms <- names(r_em_proj)
  hit <- grep(DEPLOYED_ENSEMBLE, nms, value = TRUE, fixed = TRUE)
  if (length(hit) == 0) {
    key <- sub("^EM", "", DEPLOYED_ENSEMBLE)
    hit <- grep(key, nms, value = TRUE, fixed = TRUE)
  }
  if (length(hit) == 0) return(NULL)
  r <- r_em_proj[[hit[1]]]
  if (max(terra::minmax(r)[2, ], na.rm = TRUE) > 1.5) r <- r / 1000
  r
}

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("partial", sid, sprintf("=== %s ===", sp$latin))

  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (!file.exists(train_path)) { log_line("partial", sid, "no training; skip"); next }
  out_png <- file.path(PARTIAL_DIR, paste0(sid, ".png"))
  if (file.exists(out_png)) { log_line("partial", sid, "png exists; skip"); next }

  bundle <- load_species_biomod(sid)
  if (is.null(bundle)) { log_line("partial", sid, "no biomod outputs; skip"); next }
  bm_out <- bundle$bm_out; em_out <- bundle$em_out
  sp_dir <- bundle$sp_dir

  train <- readRDS(train_path)$train
  pres_rows <- train[train$pa == 1L, PREDICTORS, drop = FALSE]
  med_vals  <- vapply(pres_rows, median, numeric(1), na.rm = TRUE)

  curves <- list()

  for (v in PREDICTORS) {
    rng <- range(pres_rows[[v]], na.rm = TRUE)
    grid <- as.data.frame(matrix(rep(med_vals, each = GRID_N), nrow = GRID_N))
    names(grid) <- PREDICTORS
    grid[[v]]   <- seq(rng[1], rng[2], length.out = GRID_N)

    # Tiny 1xGRID_N raster stack the ensemble can project onto
    r <- terra::rast(nrows = 1, ncols = GRID_N,
                     xmin = 0, xmax = GRID_N, ymin = 0, ymax = 1,
                     crs = "EPSG:4326", nlyrs = length(PREDICTORS))
    names(r) <- PREDICTORS
    for (p in PREDICTORS) terra::values(r[[p]]) <- grid[[p]]

    setwd(biomod_root)
    proj_name <- gsub("[^A-Za-z0-9_]", "_",
                      sprintf("partial_hab_%s_%s", sid, v))
    pf <- tryCatch(
      biomod2::BIOMOD_Projection(
        bm.mod = bm_out, proj.name = proj_name, new.env = r,
        models.chosen = "all", metric.binary = NULL,
        compress = TRUE, build.clamping.mask = FALSE,
        output.format = ".tif", nb.cpu = 1),
      error = function(e) NULL)
    pf_em <- if (!is.null(pf))
      tryCatch(biomod2::BIOMOD_EnsembleForecasting(
        bm.em = em_out, bm.proj = pf,
        output.format = ".tif", nb.cpu = 1),
        error = function(e) NULL) else NULL
    setwd(orig_wd)

    if (is.null(pf_em)) next

    em_dir <- file.path(sp_dir, paste0("proj_", proj_name))
    em_tif <- list.files(em_dir, pattern = "_ensemble\\.tif$", full.names = TRUE)
    if (!length(em_tif)) next
    r_out <- extract_deployed(terra::rast(em_tif[1]))
    if (is.null(r_out)) next

    curves[[v]] <- data.frame(
      predictor = v,
      x         = grid[[v]],
      y         = as.numeric(terra::values(r_out))[seq_len(GRID_N)]
    )
  }

  if (!length(curves)) { log_line("partial", sid, "no curves produced"); next }
  df <- dplyr::bind_rows(curves)

  # Rugs (presence values per predictor)
  rug_df <- tidyr::pivot_longer(pres_rows, cols = everything(),
                                names_to = "predictor", values_to = "x")

  g <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 0.6, colour = "#1f78b4") +
    geom_rug(data = rug_df, aes(x = x), inherit.aes = FALSE,
             sides = "b", alpha = 0.2, length = grid::unit(0.02, "npc")) +
    facet_wrap(~ predictor, ncol = 4, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste0(sid, "  —  marginal responses (", DEPLOYED_ENSEMBLE, ")"),
         x = NULL, y = "P(suitable)") +
    theme_bw(base_size = 10)

  ggsave(out_png, g, width = 11, height = 6, dpi = 130)
  log_line("partial", sid, sprintf("wrote %s.png (%d curves)", sid, length(curves)))
}

log_line("partial", "-", "done")

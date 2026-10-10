# =============================================================================
# 06_thresholds.R
# For each species, compute two binary classification thresholds:
#
#   maxTSS_thr   — threshold that maximizes True Skill Statistic (Youden's J)
#   sens90_thr   — threshold that gives ~90% sensitivity
#
# Prefers OUT-OF-FOLD ensemble predictions written by 07b_honest_scores.R
# (eval_hab/<sid>_oof_ensemble.rds). Falls back to in-sample predictions
# (ensemble raster scored at training points) if the OOF file is absent — in
# which case the thresholds and reported sens/spec/tss are resubstitution
# estimates and should be labelled as such.
#
# The maxTSS threshold is the default cutoff used by 10, 13, and
# aggregate_counties. sens90 is the recall-biased alternative for management
# applications where missing an early invader is worse than a false alarm.
#
# Output:
#   EVAL_DIR/thresholds_summary.csv  species, maxTSS_thr, sens, spec, tss,
#                                    sens90_thr, source ("oof" or "insample")
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
  library(pROC)
  library(dplyr)
})

log_line("thresholds", "-", "starting")

orig_wd     <- getwd()
biomod_root <- normalizePath(MODELS_DIR, winslash = "/", mustWork = TRUE)

# biomod2 uses save() (not saveRDS()) for its .out files
load_biomod_out <- function(path) {
  env <- new.env()
  loaded <- tryCatch(load(path, envir = env), error = function(e) character(0))
  if (length(loaded) == 0) return(NULL)
  env[[loaded[1]]]
}

load_species_em <- function(sid) {
  # biomod2 replaces underscores with dots in the resp.name folder
  for (d_name in unique(c(gsub("_", ".", sid), sid))) {
    sp_dir <- file.path(biomod_root, d_name)
    if (!dir.exists(sp_dir)) next
    em_files <- list.files(sp_dir, pattern = "ensemble\\.models\\.out$",
                           full.names = TRUE)
    if (length(em_files) == 0) next
    em_file <- em_files[order(file.mtime(em_files), decreasing = TRUE)][1]
    em_out  <- load_biomod_out(em_file)
    if (!is.null(em_out)) return(list(em_out = em_out, sp_dir = sp_dir))
  }
  NULL
}

# Locate the ensemble prediction raster for the current climate
load_current_ensemble_raster <- function(sid) {
  f <- file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid))
  if (!file.exists(f)) return(NULL)
  terra::rast(f)
}

find_maxTSS <- function(obs, pred) {
  r <- pROC::roc(obs, pred, quiet = TRUE, direction = "<")
  j <- pROC::coords(r, "best", best.method = "youden",
                    ret = c("threshold", "sensitivity", "specificity"))
  if (is.data.frame(j) && nrow(j) > 1) j <- j[1, ]
  list(thr  = as.numeric(j$threshold),
       sens = as.numeric(j$sensitivity),
       spec = as.numeric(j$specificity),
       tss  = as.numeric(j$sensitivity + j$specificity - 1))
}

find_sens_thr <- function(obs, pred, target = 0.90) {
  r <- pROC::roc(obs, pred, quiet = TRUE, direction = "<")
  cc <- pROC::coords(r, "all",
                     ret = c("threshold", "sensitivity", "specificity"))
  cc <- cc[order(abs(cc$sensitivity - target)), ]
  b  <- cc[1, ]
  list(thr  = as.numeric(b$threshold),
       sens = as.numeric(b$sensitivity),
       spec = as.numeric(b$specificity))
}

rows <- list()

load_predictions <- function(sid) {
  oof_path <- file.path(EVAL_DIR, sprintf("%s_oof_ensemble.rds", sid))
  if (file.exists(oof_path)) {
    oof <- readRDS(oof_path)
    ok  <- is.finite(oof$em_mean) & is.finite(oof$obs)
    if (sum(ok) >= 30 && length(unique(oof$obs[ok])) >= 2) {
      return(list(obs = oof$obs[ok], pred = oof$em_mean[ok],
                  source = "oof", n = sum(ok)))
    }
  }
  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (!file.exists(train_path)) return(NULL)
  cur_r <- load_current_ensemble_raster(sid)
  if (is.null(cur_r)) return(NULL)
  if (max(terra::minmax(cur_r)[2, ], na.rm = TRUE) > 1.5) cur_r <- cur_r / 1000
  train <- readRDS(train_path)$train
  pts   <- terra::vect(train[, c("x", "y")], geom = c("x", "y"), crs = "EPSG:4326")
  pred  <- terra::extract(cur_r, pts)[, 2]
  obs   <- train$pa
  ok    <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 30 || length(unique(obs[ok])) < 2) return(NULL)
  list(obs = obs[ok], pred = pred[ok], source = "insample", n = sum(ok))
}

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("thresholds", sid, sprintf("=== %s ===", sp$latin))

  dat <- load_predictions(sid)
  if (is.null(dat)) {
    log_line("thresholds", sid, "no usable predictions; skipping"); next
  }

  mt <- find_maxTSS(dat$obs, dat$pred)
  s9 <- find_sens_thr(dat$obs, dat$pred, target = 0.90)

  rows[[sid]] <- data.frame(
    species     = sid,
    latin       = sp$latin,
    maxTSS_thr  = round(mt$thr, 4),
    sens        = round(mt$sens, 3),
    spec        = round(mt$spec, 3),
    tss         = round(mt$tss, 3),
    sens90_thr  = round(s9$thr, 4),
    n_train     = dat$n,
    source      = dat$source,
    stringsAsFactors = FALSE
  )

  log_line("thresholds", sid,
           sprintf("[%s] maxTSS thr=%.3f (tss=%.3f)  sens90 thr=%.3f",
                   dat$source, mt$thr, mt$tss, s9$thr))
}

if (length(rows)) {
  out <- dplyr::bind_rows(rows)
  out_path <- file.path(EVAL_DIR, "thresholds_summary.csv")
  write.csv(out, out_path, row.names = FALSE)
  log_line("thresholds", "-",
           sprintf("wrote %s (%d species)", basename(out_path), nrow(out)))
}

log_line("thresholds", "-", "done")
cat("\nThresholds written to", file.path(EVAL_DIR, "thresholds_summary.csv"), "\n")

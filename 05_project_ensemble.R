# =============================================================================
# 05_project_ensemble.R  —  habitat-augmented (final_hab_v1)
#
# WHAT CHANGED vs. the climate-only version
# -------------------------------------------------------------------------
#   1. Projections now use the full env stack (climate + habitat) via
#      load_env_current() and load_env_future(gcm, ssp, label). Habitat is
#      static across futures; only climate varies.
#   2. Pre-flight climate check now expects 4 layers (CLIM_PREDICTORS), not
#      8 (PREDICTORS) — the future files still only carry bioclim vars.
#   3. Added pre-flight check that hab_static.tif exists at the right extent.
#   4. proj_name changed from "current_v3" → "current_hab" so biomod2 writes
#      into a fresh proj_current_hab/ subdir per species. Future proj names
#      also carry a "_hab" flag. That way even inside models_hab/<sid>/ the
#      per-run projection folders don't collide with anything.
#
# All the resilience machinery (extent-aware skip, honest success detection,
# per-worker logging, serial-only driver, base-R median aggregation) is
# preserved. This file remains SAFE to re-run — it skips anything that's
# already on disk at the right extent.
#
# Usage (from the pipeline root):
#   FORCE_REPROJECT <- FALSE
#   source("05_project_ensemble.R")
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
})

# ── EXTENT-AWARE RESUME ──────────────────────────────────────────────────────
extent_matches <- function(path, target_ext) {
  if (!file.exists(path)) return(FALSE)
  r <- tryCatch(terra::rast(path), error = function(e) NULL)
  if (is.null(r)) return(FALSE)
  ex  <- terra::ext(r)
  tol <- max(terra::res(r)) * 0.51
  all(abs(c(ex$xmin - target_ext$xmin,
            ex$xmax - target_ext$xmax,
            ex$ymin - target_ext$ymin,
            ex$ymax - target_ext$ymax)) <= tol)
}

if (!exists("FORCE_REPROJECT")) FORCE_REPROJECT <- FALSE
if (!exists("N_CPU"))           N_CPU <- max(1, parallel::detectCores() - 1)
n_cpu_bio <- N_CPU

# ── LOGGING ──────────────────────────────────────────────────────────────────
dir.create("logs", recursive = TRUE, showWarnings = FALSE)
LOG_DIR <- normalizePath("logs", winslash = "/", mustWork = TRUE)

checkpoint <- function(sid, msg) {
  tryCatch(log_line("project", sid, msg), error = function(e) NULL)
  tryCatch({
    f  <- file.path(LOG_DIR, sprintf("05_hab_worker_%s.log", sid))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    con <- file(f, open = "at")
    writeLines(sprintf("%s  %s", ts, msg), con)
    close(con)
  }, error = function(e) NULL)
  cat(sprintf("[%s] %s: %s\n",
              format(Sys.time(), "%H:%M:%S"), sid, msg))
}

mem_mb <- function() round(as.numeric(gc(verbose = FALSE)[2, 2]), 1)

log_line("project", "-", "starting (SERIAL, honest, hab)")
log_line("project", "-",
         sprintf("PROJ_EXTENT: xmin=%g xmax=%g ymin=%g ymax=%g%s",
                 PROJ_EXTENT$xmin, PROJ_EXTENT$xmax,
                 PROJ_EXTENT$ymin, PROJ_EXTENT$ymax,
                 if (FORCE_REPROJECT) "  [FORCE_REPROJECT=TRUE]" else ""))
log_line("project", "-", sprintf("predictors (%d): %s",
                                 length(PREDICTORS),
                                 paste(PREDICTORS, collapse = ", ")))

ok_at_extent <- function(path) {
  !FORCE_REPROJECT && extent_matches(path, PROJ_EXTENT)
}

orig_wd         <- getwd()
biomod_root     <- normalizePath(MODELS_DIR,       winslash = "/", mustWork = TRUE)
cur_dir_abs     <- normalizePath(CURRENT_DIR,      winslash = "/", mustWork = FALSE)
fut_dir_abs     <- normalizePath(FUTURE_DIR,       winslash = "/", mustWork = FALSE)
fut_dir_climate <- normalizePath(file.path(CLIMATE_DIR, "future"),
                                 winslash = "/", mustWork = FALSE)
hab_static_path <- normalizePath(file.path(HAB_DIR, "hab_static.tif"),
                                 winslash = "/", mustWork = FALSE)
dir.create(cur_dir_abs, recursive = TRUE, showWarnings = FALSE)
dir.create(fut_dir_abs, recursive = TRUE, showWarnings = FALSE)

# ── HELPERS ──────────────────────────────────────────────────────────────────
rescale_biomod_raster <- function(r) {
  if (max(terra::minmax(r)[2, ], na.rm = TRUE) > 1.5) r <- r / 1000
  r
}

load_biomod_out <- function(path) {
  env <- new.env()
  loaded <- tryCatch(load(path, envir = env), error = function(e) character(0))
  if (length(loaded) == 0) return(NULL)
  env[[loaded[1]]]
}

load_species_biomod <- function(sid) {
  biomod_sid <- gsub("_", ".", sid)
  for (d_name in unique(c(biomod_sid, sid))) {
    sp_dir <- file.path(biomod_root, d_name)
    if (!dir.exists(sp_dir)) next
    all_out  <- list.files(sp_dir, pattern = "\\.out$", full.names = TRUE)
    em_files <- grep("ensemble\\.models\\.out$", all_out, value = TRUE)
    bm_files <- grep("\\.models\\.out$",        all_out, value = TRUE)
    bm_files <- setdiff(bm_files, em_files)
    if (length(bm_files) == 0 || length(em_files) == 0) next
    bm_file <- bm_files[order(file.mtime(bm_files), decreasing = TRUE)][1]
    em_file <- em_files[order(file.mtime(em_files), decreasing = TRUE)][1]
    bm_out <- load_biomod_out(bm_file)
    em_out <- load_biomod_out(em_file)
    if (is.null(bm_out) || is.null(em_out)) next
    return(list(bm_out = bm_out, em_out = em_out,
                biomod_sid = d_name, sp_dir = sp_dir))
  }
  return(NULL)
}

read_biomod_proj_rasters <- function(sp_dir, proj_name) {
  proj_dir <- file.path(sp_dir, paste0("proj_", proj_name))
  if (!dir.exists(proj_dir)) return(list(alg = NULL, em = NULL))
  tifs <- list.files(proj_dir, pattern = "\\.tif$", full.names = TRUE)
  em_tifs  <- grep("_ensemble\\.tif$", tifs, value = TRUE)
  alg_tifs <- setdiff(tifs, em_tifs)
  list(
    alg = if (length(alg_tifs)) terra::rast(alg_tifs[1]) else NULL,
    em  = if (length(em_tifs))  terra::rast(em_tifs[1])  else NULL
  )
}

extract_ensemble_layer <- function(r_em_proj, deployed_tag = DEPLOYED_ENSEMBLE) {
  if (is.null(r_em_proj)) return(NULL)
  nms <- names(r_em_proj)
  hit <- grep(deployed_tag, nms, value = TRUE, fixed = TRUE)
  if (length(hit) == 0) {
    key <- sub("^EM", "", deployed_tag)
    hit <- grep(key, nms, value = TRUE, fixed = TRUE)
  }
  if (length(hit) == 0) return(NULL)
  rescale_biomod_raster(r_em_proj[[hit[1]]])
}

extract_algo_layer <- function(r_proj, algo) {
  if (is.null(r_proj)) return(NULL)
  nms <- names(r_proj)
  hit <- grep(paste0("_", algo, "$"), nms, value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0)
    hit <- grep(algo, nms, value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) return(NULL)
  r_sub <- r_proj[[hit]]
  if (terra::nlyr(r_sub) > 1) r_sub <- terra::mean(r_sub, na.rm = TRUE)
  rescale_biomod_raster(r_sub)
}

species_is_done <- function(sid) {
  cur <- file.path(cur_dir_abs, sprintf("%s_ensemble_current.tif", sid))
  if (!ok_at_extent(cur)) return(FALSE)
  for (ssp in SSPS) {
    for (label in PERIOD_LABELS) {
      f <- file.path(fut_dir_abs,
                     sprintf("%s_ensemble_%s_ssp%s_median.tif", sid, label, ssp))
      if (!ok_at_extent(f)) return(FALSE)
    }
  }
  TRUE
}

# ── SAFE AGGREGATION — NO terra::app(). Base R matrix median. ────────────────
aggregate_gcms <- function(rlist, out_median, out_agreement, out_cv) {
  if (length(rlist) == 0) stop("aggregate_gcms: empty rlist")
  template <- rlist[[1]]
  vals <- matrix(NA_real_, nrow = terra::ncell(template), ncol = length(rlist))
  for (i in seq_along(rlist)) {
    vals[, i] <- terra::values(rlist[[i]])[, 1]
  }
  med_vals <- apply(vals, 1, function(v) median(v, na.rm = TRUE))
  med_r <- terra::setValues(terra::rast(template), med_vals)
  terra::writeRaster(med_r, out_median, overwrite = TRUE, gdal = GDAL_OPTS)

  agree_vals <- apply(vals, 1, function(v) mean(v >= 0.5, na.rm = TRUE))
  agree_r <- terra::setValues(terra::rast(template), agree_vals)
  terra::writeRaster(agree_r, out_agreement, overwrite = TRUE, gdal = GDAL_OPTS)

  cv_vals <- apply(vals, 1, function(v) {
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu) || mu == 0) return(NA_real_)
    sd(v, na.rm = TRUE) / mu
  })
  cv_r <- terra::setValues(terra::rast(template), cv_vals)
  terra::writeRaster(cv_r, out_cv, overwrite = TRUE, gdal = GDAL_OPTS)

  if (!file.exists(out_median))    stop(sprintf("aggregate_gcms: median NOT written: %s", out_median))
  if (!file.exists(out_agreement)) stop(sprintf("aggregate_gcms: agreement NOT written: %s", out_agreement))
  if (!file.exists(out_cv))        stop(sprintf("aggregate_gcms: cv NOT written: %s", out_cv))
  invisible(TRUE)
}

# ── PER-SPECIES PROJECTION ───────────────────────────────────────────────────
project_one_species <- function(sp, n_cpu_bio = 1L) {
  sid <- sp$short_id

  if (species_is_done(sid)) {
    checkpoint(sid, "already complete — skipping")
    return(invisible(TRUE))
  }

  checkpoint(sid, sprintf("=== %s === START  (mem=%.1f MB, nb.cpu=%d)",
                          sp$latin, mem_mb(), n_cpu_bio))

  worker_wd <- getwd()
  on.exit(setwd(worker_wd), add = TRUE)

  bundle <- load_species_biomod(sid)
  if (is.null(bundle)) {
    checkpoint(sid, "no biomod workspace; skipping")
    return(invisible(FALSE))
  }
  bm_out     <- bundle$bm_out
  em_out     <- bundle$em_out
  biomod_sid <- bundle$biomod_sid
  sp_dir     <- bundle$sp_dir
  checkpoint(sid, sprintf("loaded biomod workspace  (mem=%.1f MB)", mem_mb()))

  setwd(biomod_root)

  # ── 1. CURRENT ─────────────────────────────────────────────────────────
  cur_out_path <- file.path(cur_dir_abs, sprintf("%s_ensemble_current.tif", sid))
  cur_ok <- ok_at_extent(cur_out_path)
  if (!cur_ok) {
    checkpoint(sid, sprintf("current: projecting  (mem=%.1f MB)", mem_mb()))
    # NEW: combined climate + habitat stack (8 layers, order = PREDICTORS)
    env_current <- load_env_current()

    proj_current <- biomod2::BIOMOD_Projection(
      bm.mod = bm_out, proj.name = "current_hab", new.env = env_current,
      models.chosen = "all", metric.binary = NULL, metric.filter = NULL,
      compress = TRUE, build.clamping.mask = FALSE,
      output.format = ".tif", nb.cpu = n_cpu_bio)

    proj_current_em <- biomod2::BIOMOD_EnsembleForecasting(
      bm.em = em_out, bm.proj = proj_current, models.chosen = "all",
      metric.binary = NULL, compress = TRUE, output.format = ".tif",
      nb.cpu = n_cpu_bio)

    setwd(worker_wd)
    rr <- read_biomod_proj_rasters(sp_dir, "current_hab")

    for (algo in BIOMOD_MODELS) {
      r <- extract_algo_layer(rr$alg, algo)
      if (is.null(r)) next
      algo_tag <- tolower(algo)
      if (algo_tag == "maxnet") algo_tag <- "maxent"
      terra::writeRaster(r,
                         file.path(cur_dir_abs, sprintf("%s_%s_current.tif", sid, algo_tag)),
                         overwrite = TRUE, gdal = GDAL_OPTS)
    }

    r_deploy <- extract_ensemble_layer(rr$em, DEPLOYED_ENSEMBLE)
    if (!is.null(r_deploy)) {
      terra::writeRaster(r_deploy, cur_out_path,
                         overwrite = TRUE, gdal = GDAL_OPTS)
      checkpoint(sid, sprintf("wrote current  (mem=%.1f MB)", mem_mb()))
    }

    cv_hit <- grep("EMcv", names(rr$em), value = TRUE, fixed = TRUE)
    if (length(cv_hit)) {
      cv_r <- rescale_biomod_raster(rr$em[[cv_hit[1]]])
      terra::writeRaster(cv_r,
                         file.path(cur_dir_abs, sprintf("%s_ensemble_current_cv.tif", sid)),
                         overwrite = TRUE, gdal = GDAL_OPTS)
    }
    rm(rr, r_deploy, proj_current, proj_current_em, env_current)
    gc(verbose = FALSE)
  } else {
    checkpoint(sid, "current already exists at correct extent — skipping")
    setwd(worker_wd)
  }

  # ── 2. FUTURES ─────────────────────────────────────────────────────────
  scenarios_written <- 0
  for (ssp in SSPS) {
    for (label in PERIOD_LABELS) {
      key <- sprintf("ssp%s_%s", ssp, label)
      out_med_path   <- file.path(fut_dir_abs, sprintf("%s_ensemble_%s_ssp%s_median.tif",    sid, label, ssp))
      out_agree_path <- file.path(fut_dir_abs, sprintf("%s_ensemble_%s_ssp%s_agreement.tif", sid, label, ssp))
      out_cv_path    <- file.path(fut_dir_abs, sprintf("%s_ensemble_%s_ssp%s_cv.tif",        sid, label, ssp))

      if (ok_at_extent(out_med_path)) {
        checkpoint(sid, sprintf("%s already exists — skipping", basename(out_med_path)))
        scenarios_written <- scenarios_written + 1
        next
      }

      checkpoint(sid, sprintf(">>> scenario %s: starting %d GCMs  (mem=%.1f MB)",
                              key, length(GCMS), mem_mb()))

      rlist <- list()
      for (gcm in GCMS) {
        fut_tif <- file.path(fut_dir_climate, sprintf("%s_ssp%s_%s.tif", gcm, ssp, label))
        if (!file.exists(fut_tif)) {
          checkpoint(sid, sprintf("  %s MISSING climate tif, skipping", gcm))
          next
        }

        checkpoint(sid, sprintf("  %s: START  (mem=%.1f MB)", gcm, mem_mb()))
        gcm_t0 <- Sys.time()

        # NEW: combined future-climate + static-habitat stack
        env_fut <- load_env_future(gcm, ssp, label)

        proj_name <- gsub("[^A-Za-z0-9_]", "_",
                          sprintf("fut_hab_%s_ssp%s_%s", gcm, ssp, label))

        setwd(biomod_root)
        pf <- biomod2::BIOMOD_Projection(
          bm.mod = bm_out, proj.name = proj_name, new.env = env_fut,
          models.chosen = "all", metric.binary = NULL,
          compress = TRUE, build.clamping.mask = FALSE,
          output.format = ".tif", nb.cpu = n_cpu_bio)

        pf_em <- biomod2::BIOMOD_EnsembleForecasting(
          bm.em = em_out, bm.proj = pf,
          output.format = ".tif", nb.cpu = n_cpu_bio)
        setwd(worker_wd)

        rr_fut <- read_biomod_proj_rasters(sp_dir, proj_name)
        r_deploy_fut <- extract_ensemble_layer(rr_fut$em, DEPLOYED_ENSEMBLE)
        if (!is.null(r_deploy_fut)) {
          rlist[[length(rlist) + 1]] <- r_deploy_fut
        }
        rm(rr_fut, r_deploy_fut, env_fut, pf, pf_em)
        gc(verbose = FALSE)

        gcm_elapsed <- round(as.numeric(difftime(Sys.time(), gcm_t0, units = "mins")), 2)
        checkpoint(sid, sprintf("  %s: DONE in %g min  (mem=%.1f MB)",
                                gcm, gcm_elapsed, mem_mb()))
      }

      if (length(rlist) == 0) {
        stop(sprintf("scenario %s produced no GCM outputs — cannot aggregate", key))
      }

      checkpoint(sid, sprintf("<<< scenario %s: aggregating %d GCMs (base R median)",
                              key, length(rlist)))
      aggregate_gcms(rlist, out_med_path, out_agree_path, out_cv_path)
      checkpoint(sid, sprintf("wrote %s (%d GCMs)  (mem=%.1f MB)",
                              basename(out_med_path), length(rlist), mem_mb()))
      scenarios_written <- scenarios_written + 1
      rm(rlist)
      gc(verbose = FALSE)
    }
  }

  # HONEST success check — did we actually write all 8 futures?
  expected_scenarios <- length(SSPS) * length(PERIOD_LABELS)
  if (scenarios_written < expected_scenarios) {
    stop(sprintf("%s finished with only %d/%d scenarios written",
                 sid, scenarios_written, expected_scenarios))
  }

  rm(bm_out, em_out, bundle)
  gc(verbose = FALSE)
  checkpoint(sid, sprintf("=== DONE (%d/%d scenarios written)  (mem=%.1f MB)",
                          scenarios_written, expected_scenarios, mem_mb()))
  invisible(TRUE)
}

# ── PRE-FLIGHT ───────────────────────────────────────────────────────────────
# Climate futures: expect 4 layers (CLIM_PREDICTORS only)
log_line("project", "-", "pre-flight: checking future climate rasters")
bad_fut <- character(0)
for (gcm in GCMS) for (ssp in SSPS) for (label in PERIOD_LABELS) {
  f <- file.path(fut_dir_climate, sprintf("%s_ssp%s_%s.tif", gcm, ssp, label))
  if (!file.exists(f)) { bad_fut <- c(bad_fut, paste("MISSING", basename(f))); next }
  r <- tryCatch(terra::rast(f), error = function(e) NULL)
  if (is.null(r))               { bad_fut <- c(bad_fut, paste("UNREADABLE", basename(f))); next }
  if (terra::nlyr(r) != length(CLIM_PREDICTORS))
    bad_fut <- c(bad_fut, sprintf("WRONG LAYERS (%d): %s", terra::nlyr(r), basename(f)))
  ex <- terra::ext(r); tol <- max(terra::res(r)) * 0.51
  if (any(abs(c(ex$xmin - PROJ_EXTENT$xmin, ex$xmax - PROJ_EXTENT$xmax,
                ex$ymin - PROJ_EXTENT$ymin, ex$ymax - PROJ_EXTENT$ymax)) > tol))
    bad_fut <- c(bad_fut, paste("WRONG EXTENT", basename(f)))
}
if (length(bad_fut)) {
  log_line("project", "-", "pre-flight FAILED:")
  for (m in bad_fut) log_line("project", "-", paste("  ", m))
  stop("pre-flight climate check failed — see log")
}
log_line("project", "-",
         sprintf("pre-flight OK (climate): %d future rasters",
                 length(GCMS) * length(SSPS) * length(PERIOD_LABELS)))

# NEW: pre-flight for hab_static.tif — must exist, have 4 layers, right extent
log_line("project", "-", "pre-flight: checking hab_static.tif")
if (!file.exists(hab_static_path))
  stop("Missing hab_static.tif — run 01b_download_habitat.R")
hs <- terra::rast(hab_static_path)
if (terra::nlyr(hs) != length(HAB_PREDICTORS))
  stop(sprintf("hab_static.tif has %d layers, expected %d (%s)",
               terra::nlyr(hs), length(HAB_PREDICTORS),
               paste(HAB_PREDICTORS, collapse = ", ")))
if (!extent_matches(hab_static_path, PROJ_EXTENT))
  stop("hab_static.tif extent doesn't match PROJ_EXTENT — re-run 01b_download_habitat.R")
log_line("project", "-", "pre-flight OK (habitat)")

# ── SERIAL DRIVER — one species at a time, honest reporting ─────────────────
overall_t0 <- Sys.time()
outcomes <- character(length(SPECIES))
names(outcomes) <- vapply(SPECIES, function(sp) sp$short_id, character(1))

for (i in seq_along(SPECIES)) {
  sp <- SPECIES[[i]]
  sp_t0 <- Sys.time()
  message(sprintf("\n>>> species %d/%d: %s", i, length(SPECIES), sp$short_id))
  outcomes[i] <- tryCatch({
    project_one_species(sp, n_cpu_bio = n_cpu_bio)
    "OK"
  }, error = function(e) {
    msg <- conditionMessage(e)
    log_line("project", sp$short_id, sprintf("ERROR: %s", msg))
    message(sprintf("!!! %s failed: %s", sp$short_id, msg))
    sprintf("ERROR: %s", msg)
  })
  sp_elapsed <- round(as.numeric(difftime(Sys.time(), sp_t0, units = "mins")), 1)
  message(sprintf("<<< species %d/%d done in %g min  (outcome: %s)",
                  i, length(SPECIES), sp_elapsed, outcomes[i]))
}

# ── FINAL SUMMARY ────────────────────────────────────────────────────────────
total_min <- round(as.numeric(difftime(Sys.time(), overall_t0, units = "mins")), 1)
n_ok  <- sum(outcomes == "OK")
n_err <- sum(outcomes != "OK")

cat("\n=====================================================\n")
cat(sprintf("Projection complete.  Total time: %g min\n", total_min))
cat(sprintf("  %d OK, %d errored (of %d species)\n", n_ok, n_err, length(SPECIES)))
if (n_err > 0) {
  cat("\nFailed species:\n")
  for (sid in names(outcomes)[outcomes != "OK"]) {
    cat(sprintf("  %-20s  %s\n", sid, outcomes[sid]))
  }
}

n_cur_disk <- sum(vapply(SPECIES, function(sp)
  file.exists(file.path(cur_dir_abs, sprintf("%s_ensemble_current.tif", sp$short_id))),
  logical(1)))
n_fut_disk <- length(list.files(fut_dir_abs, pattern = "_ensemble_.*_median\\.tif$"))
cat(sprintf("\nDisk state: %d/%d currents, %d/96 future medians\n",
            n_cur_disk, length(SPECIES), n_fut_disk))
cat("=====================================================\n")

setwd(orig_wd)

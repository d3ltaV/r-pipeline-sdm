# =============================================================================
# 01_download_climate.R
# Download WorldClim 2.1 current + CMIP6 futures for the 4 bioclim predictors,
# crop to PROJ_EXTENT, save GeoTIFFs.
#
# NOTE This script produces the CLIMATE-ONLY rasters (4 layers = CLIM_PREDICTORS).
#      Habitat layers are produced separately by 01b_download_habitat.R.
#      The two are merged into the full 8-predictor stack at read-time by
#      load_env_current() / load_env_future() / load_env_global() in 00_config.R.
#      Therefore layer names here MUST be CLIM_PREDICTORS, NOT PREDICTORS
#      (assigning 8 names to a 4-layer stack throws
#      "[names<-] incorrect number of names").
#
# YOU DON'T NEED TO RUN THIS if climate/clim_current.tif and the future stacks
# already exist from your v2a run — they're identical between pipeline versions
# (same predictors, same source). This script is here for reproducibility.
#
# Runtime on a fresh machine: ~1-2 hours (dominated by CMIP6 downloads;
# 5 GCMs x 2 SSPs x 4 periods = 40 raster stacks, each ~15 MB compressed).
#
# Outputs
#   CLIMATE_DIR/clim_current_global.tif       global; needed by 03 for niche extraction
#   CLIMATE_DIR/clim_current.tif              cropped to PROJ_EXTENT; used for prediction
#   CLIMATE_DIR/future/<gcm>_ssp<ssp>_<label>.tif   one per future scenario
#
# Uses `geodata` package — the canonical R front-end to WorldClim + CMIP6.
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(geodata)
  library(terra)
})

# Defensive sanity check — if someone edits 00_config.R and gets these
# out of sync, fail loudly with a clear message instead of the cryptic
# "[names<-] incorrect number of names" three lines down.
stopifnot(
  "BIO_IDX and CLIM_PREDICTORS must have the same length" =
    length(BIO_IDX) == length(CLIM_PREDICTORS)
)

# ── EXTENT-AWARE RESUME ──────────────────────────────────────────────────────
# If the user changes PROJ_EXTENT between runs, plain file.exists() would
# silently keep stale cropped rasters at the old extent. These helpers let
# 01/05 detect that and regenerate only what's affected. Backward-compatible:
# first-time runs (no files) behave identically; repeat runs at the same
# PROJ_EXTENT still fast-skip; repeat runs at a *different* extent detect the
# mismatch and rebuild.
#
# Tolerance is half a pixel — terra::crop snaps to cell edges so requested vs
# actual extent can differ by <1 cell width even at "matching" resolution.
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

# Escape hatch: `FORCE_REPROJECT <- TRUE` at the top of a session forces
# regeneration of every raster regardless of what's on disk. Undefined => FALSE.
if (!exists("FORCE_REPROJECT")) FORCE_REPROJECT <- FALSE

log_line("download_climate", "-", "starting")
log_line("download_climate", "-",
         sprintf("PROJ_EXTENT: xmin=%g xmax=%g ymin=%g ymax=%g%s",
                 PROJ_EXTENT$xmin, PROJ_EXTENT$xmax,
                 PROJ_EXTENT$ymin, PROJ_EXTENT$ymax,
                 if (FORCE_REPROJECT) "  [FORCE_REPROJECT=TRUE]" else ""))

# ── CURRENT CLIMATE ──────────────────────────────────────────────────────────
out_global <- file.path(CLIMATE_DIR, "clim_current_global.tif")
out_proj   <- file.path(CLIMATE_DIR, "clim_current.tif")

# The GLOBAL stack is extent-independent — either it's on disk or it isn't.
# The PROJ stack is a crop, so it has to match PROJ_EXTENT to be reused.
proj_ok <- !FORCE_REPROJECT && extent_matches(out_proj, PROJ_EXTENT)

if (file.exists(out_global) && proj_ok) {
  log_line("download_climate", "-",
           "current climate cached and extent matches; skipping")
} else if (file.exists(out_global) && !proj_ok) {
  # Common path when the user widens PROJ_EXTENT: raw global data is fine,
  # we just need to re-crop. No re-download.
  log_line("download_climate", "-",
           if (FORCE_REPROJECT)
             "FORCE: re-cropping current from cached global stack"
           else
             "PROJ_EXTENT changed since last run - re-cropping current from cached global stack")
  clim <- terra::rast(out_global)
  names(clim) <- CLIM_PREDICTORS                     # 4-layer climate-only stack
  writeRaster(terra::crop(clim, PROJ_EXTENT), out_proj,
              overwrite = TRUE, gdal = GDAL_OPTS)
  log_line("download_climate", "-", "wrote current (from cache, no download)")
  rm(clim); gc(verbose = FALSE)
} else {
  # First-time run: download WorldClim + write both global and cropped stacks.
  log_line("download_climate", "-",
           sprintf("downloading WorldClim current (bio, res=%d')", WORLDCLIM_RES))
  bio <- geodata::worldclim_global(var = "bio", res = WORLDCLIM_RES,
                                   path = CACHE_DIR)
  # geodata returns 19 bioclim layers as "wc2.1_5m_bio_1" ... "wc2.1_5m_bio_19"
  names(bio) <- paste0("bio", seq_len(nlyr(bio)))
  clim <- bio[[ paste0("bio", BIO_IDX) ]]            # subset to our 4 bioclims
  names(clim) <- CLIM_PREDICTORS                     # 4 names for 4 layers
  
  # Save the global stack for niche extraction (script 03 samples background
  # from within a buffered calibration region that may extend beyond
  # PROJ_EXTENT, so the extraction raster has to be global).
  writeRaster(clim, out_global, overwrite = TRUE, gdal = GDAL_OPTS)
  
  # Cropped version for prediction rasters.
  writeRaster(terra::crop(clim, PROJ_EXTENT), out_proj,
              overwrite = TRUE, gdal = GDAL_OPTS)
  
  log_line("download_climate", "-", "wrote current + global stacks")
  rm(bio, clim); gc(verbose = FALSE)
}

# ── FUTURE CLIMATE (CMIP6) ───────────────────────────────────────────────────
fut_dir <- file.path(CLIMATE_DIR, "future")
dir.create(fut_dir, recursive = TRUE, showWarnings = FALSE)

for (gcm in GCMS) {
  for (ssp in SSPS) {
    for (i in seq_along(PERIODS)) {
      period <- PERIODS[i]        # "2021-2040" etc.
      label  <- PERIOD_LABELS[i]  # "2030"      etc.
      out    <- file.path(fut_dir, sprintf("%s_ssp%s_%s.tif", gcm, ssp, label))
      
      if (!FORCE_REPROJECT && extent_matches(out, PROJ_EXTENT)) {
        log_line("download_climate", "-", sprintf("cached: %s", basename(out)))
        next
      }
      if (file.exists(out) && !FORCE_REPROJECT) {
        log_line("download_climate", "-",
                 sprintf("stale extent: %s - re-cropping", basename(out)))
      }
      log_line("download_climate", "-",
               sprintf("fetching %s ssp%s %s", gcm, ssp, period))
      # cmip6_world uses CACHE_DIR too — so on a re-run with widened
      # PROJ_EXTENT this hits the cache and skips the network, then we
      # just re-crop and re-write.
      fb <- tryCatch(
        geodata::cmip6_world(model = gcm, ssp = ssp, time = period,
                             var = "bioc", res = WORLDCLIM_RES,
                             path = CACHE_DIR),
        error = function(e) {
          log_line("download_climate", "-",
                   sprintf("FAILED %s ssp%s %s : %s", gcm, ssp, period,
                           conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fb)) next
      
      names(fb) <- paste0("bio", seq_len(nlyr(fb)))
      fb <- fb[[ paste0("bio", BIO_IDX) ]]           # subset to our 4 bioclims
      names(fb) <- CLIM_PREDICTORS                   # 4 names for 4 layers
      fb <- terra::crop(fb, PROJ_EXTENT)
      writeRaster(fb, out, overwrite = TRUE, gdal = GDAL_OPTS)
      log_line("download_climate", "-", sprintf("wrote %s", basename(out)))
      rm(fb); gc(verbose = FALSE)
    }
  }
}

log_line("download_climate", "-", "done")
cat("\nClimate prep complete.\n")
cat("  Current:  ", out_proj, "\n")
cat("  Global:   ", out_global, "\n")
cat("  Futures:  ", length(list.files(fut_dir, pattern = "\\.tif$")),
    " scenarios in ", fut_dir, "\n", sep = "")
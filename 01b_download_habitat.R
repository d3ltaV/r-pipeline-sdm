# =============================================================================
# 01b_download_habitat.R  —  new for final_hab_v1
#
# Download and align the 4 habitat predictors to the WorldClim 2.1 5' grid,
# then save them as one 4-layer stack (hab_static.tif for PROJ_EXTENT,
# hab_static_global.tif for niche extraction at native-range presences).
#
# Predictors (in canonical layer order):
#   elev     elevation (m)                — geodata::elevation_global (SRTM/GMTED)
#   hfp      Human Footprint Index (0-50) — geodata::footprint (Venter et al. 2018)
#   treecov  tree cover fraction (0-100)  — geodata::landcover("trees") (Copernicus)
#   soilph   soil pH * 10, top 5 cm       — geodata::soil_world("phh2o")
#
# Everything is aggregated/resampled to WorldClim's 5' grid so the raster
# geometries line up cell-for-cell — that's what makes downstream stacking
# in 03/04/05 fast and lossless.
#
# YOU DON'T NEED TO RUN THIS again if hab_static.tif already exists at the
# correct extent and layer count. Re-runs are cached.
#
# Runtime on a fresh machine: ~15-45 min (dominated by SoilGrids and Copernicus
# tree cover tiles). geodata caches downloads under CACHE_DIR so subsequent
# extent changes only trigger re-cropping, not re-downloads.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(geodata)
  library(terra)
})

if (!exists("FORCE_REPROJECT")) FORCE_REPROJECT <- FALSE

log_line("download_hab", "-", "starting")

# ── REFERENCE GRID (WorldClim global from script 01) ────────────────────────
clim_global_path <- file.path(CLIMATE_DIR, "clim_current_global.tif")
if (!file.exists(clim_global_path))
  stop("Need clim_current_global.tif from script 01 first.\n",
       "  Expected at: ", clim_global_path)

ref_global <- terra::rast(clim_global_path)[[1]]   # any single layer is fine
log_line("download_hab", "-",
         sprintf("reference grid: %d cols × %d rows, res=%s deg",
                 terra::ncol(ref_global), terra::nrow(ref_global),
                 paste(signif(terra::res(ref_global), 3), collapse = " × ")))

out_global <- file.path(HAB_DIR, "hab_static_global.tif")
out_proj   <- file.path(HAB_DIR, "hab_static.tif")

# ── EXTENT-AWARE RESUME ─────────────────────────────────────────────────────
extent_matches <- function(path, target_ext) {
  if (!file.exists(path)) return(FALSE)
  r <- tryCatch(terra::rast(path), error = function(e) NULL)
  if (is.null(r)) return(FALSE)
  ex  <- terra::ext(r)
  tol <- max(terra::res(r)) * 0.51
  all(abs(c(ex$xmin - target_ext$xmin, ex$xmax - target_ext$xmax,
            ex$ymin - target_ext$ymin, ex$ymax - target_ext$ymax)) <= tol)
}

already_ok <- function() {
  if (!file.exists(out_global) || !file.exists(out_proj)) return(FALSE)
  rg <- tryCatch(terra::rast(out_global), error = function(e) NULL)
  rp <- tryCatch(terra::rast(out_proj),   error = function(e) NULL)
  if (is.null(rg) || is.null(rp)) return(FALSE)
  if (terra::nlyr(rg) != length(HAB_PREDICTORS)) return(FALSE)
  if (terra::nlyr(rp) != length(HAB_PREDICTORS)) return(FALSE)
  extent_matches(out_proj, PROJ_EXTENT)
}

if (!FORCE_REPROJECT && already_ok()) {
  log_line("download_hab", "-",
           "hab_static already at correct extent + layer count; skipping")
  cat("Habitat stack cached and OK.\n")
  cat("  Global: ", out_global, "\n")
  cat("  Proj:   ", out_proj,   "\n")
} else {

  # If globals exist but proj is stale, we can re-crop without re-downloading.
  if (!FORCE_REPROJECT && file.exists(out_global) &&
      terra::nlyr(terra::rast(out_global)) == length(HAB_PREDICTORS) &&
      !extent_matches(out_proj, PROJ_EXTENT)) {
    log_line("download_hab", "-",
             "PROJ_EXTENT changed since last run — re-cropping from cached global stack")
    hab_g <- terra::rast(out_global)
    names(hab_g) <- HAB_PREDICTORS
    terra::writeRaster(terra::crop(hab_g, PROJ_EXTENT), out_proj,
                       overwrite = TRUE, gdal = GDAL_OPTS)
    log_line("download_hab", "-", "wrote hab_static (from cache, no downloads)")
  } else {

    # ── FULL DOWNLOAD + ALIGN + WRITE ────────────────────────────────────

    # Helper: resample a raw layer onto the WorldClim global grid.
    #   - elevation, HFP, tree cover, soil pH are all continuous.
    #   - "average" is the right aggregation when going from fine → coarse
    #     (which is the case here — all sources are finer than 5').
    #   - "bilinear" is used for elevation because SRTM has void-filled areas
    #     where averaging vs. interpolation both work, but bilinear preserves
    #     terrain gradients better where the source is smoother than target.
    align_to_ref <- function(r, layer_name, method = "average") {
      # Some geodata rasters have odd/multiple layers — collapse first
      if (terra::nlyr(r) > 1) r <- terra::mean(r, na.rm = TRUE)
      r <- terra::resample(r, ref_global, method = method,
                           threads = TRUE)
      names(r) <- layer_name
      r
    }

    layers <- list()

    # 1. ELEVATION ────────────────────────────────────────────────────────
    log_line("download_hab", "-",
             sprintf("(1/4) elevation — geodata::elevation_global(res=%d)",
                     WORLDCLIM_RES))
    elev <- tryCatch(
      geodata::elevation_global(res = WORLDCLIM_RES, path = CACHE_DIR),
      error = function(e) {
        log_line("download_hab", "-",
                 paste("elevation_global failed:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(elev)) {
      layers$elev <- align_to_ref(elev, "elev", method = "bilinear")
      log_line("download_hab", "-",
               sprintf("  elev OK — range [%s, %s] m",
                       signif(min(terra::values(layers$elev), na.rm = TRUE), 4),
                       signif(max(terra::values(layers$elev), na.rm = TRUE), 4)))
      rm(elev); gc(verbose = FALSE)
    }

    # 2. HUMAN FOOTPRINT ─────────────────────────────────────────────────
    # Venter et al. 2018 — 8 pressures combined (built environments, crop
    # lands, pasture lands, human population density, night lights, railways,
    # roads, waterways). Range 0-50. This one predictor typically gives the
    # largest single boost to invasive SDMs (Bradley et al. 2010,
    # Ibáñez et al. 2020).
    log_line("download_hab", "-",
             "(2/4) Human Footprint — geodata::footprint(year=2009)")
    hfp <- tryCatch(
      geodata::footprint(year = 2009, path = CACHE_DIR),
      error = function(e) {
        log_line("download_hab", "-",
                 paste("footprint failed:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(hfp)) {
      layers$hfp <- align_to_ref(hfp, "hfp", method = "average")
      log_line("download_hab", "-",
               sprintf("  hfp  OK — range [%s, %s]",
                       signif(min(terra::values(layers$hfp), na.rm = TRUE), 3),
                       signif(max(terra::values(layers$hfp), na.rm = TRUE), 3)))
      rm(hfp); gc(verbose = FALSE)
    }

    # 3. TREE COVER FRACTION ─────────────────────────────────────────────
    # Copernicus Global Land Cover, 100 m, tree fractional cover 0-100.
    # Downloaded as 30°×30° tiles by geodata, then mosaicked.
    log_line("download_hab", "-",
             "(3/4) tree cover — geodata::landcover(var='trees') (~15 min)")
    treecov <- tryCatch(
      geodata::landcover(var = "trees", path = CACHE_DIR),
      error = function(e) {
        log_line("download_hab", "-",
                 paste("landcover(trees) failed:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(treecov)) {
      layers$treecov <- align_to_ref(treecov, "treecov", method = "average")
      log_line("download_hab", "-",
               sprintf("  treecov OK — range [%s, %s] %%",
                       signif(min(terra::values(layers$treecov), na.rm = TRUE), 3),
                       signif(max(terra::values(layers$treecov), na.rm = TRUE), 3)))
      rm(treecov); gc(verbose = FALSE)
    }

    # 4. SOIL pH ─────────────────────────────────────────────────────────
    # SoilGrids v2, top 0-5 cm layer, pH in H2O * 10. Reasonable for
    # temperate plants; drop or swap for `cec` or `soc` if pH turns out
    # to be uninformative in the VIF report.
    log_line("download_hab", "-",
             "(4/4) soil pH — geodata::soil_world('phh2o', depth=5)")
    soilph <- tryCatch(
      geodata::soil_world(var = "phh2o", depth = 5, path = CACHE_DIR),
      error = function(e) {
        log_line("download_hab", "-",
                 paste("soil_world failed:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(soilph)) {
      layers$soilph <- align_to_ref(soilph, "soilph", method = "average")
      log_line("download_hab", "-",
               sprintf("  soilph OK — range [%s, %s]",
                       signif(min(terra::values(layers$soilph), na.rm = TRUE), 3),
                       signif(max(terra::values(layers$soilph), na.rm = TRUE), 3)))
      rm(soilph); gc(verbose = FALSE)
    }

    # ── SANITY CHECK — all 4 present? ────────────────────────────────────
    missing <- setdiff(HAB_PREDICTORS, names(layers))
    if (length(missing)) {
      stop("Habitat downloads failed for: ", paste(missing, collapse = ", "),
           "\n  See ", LOG_FILE, " for details.",
           "\n  Fix connectivity / geodata version and re-run 01b.")
    }

    # ── WRITE COMBINED STACKS ────────────────────────────────────────────
    # Layer order MUST match HAB_PREDICTORS so downstream stacking with
    # climate produces the canonical PREDICTORS order.
    hab_stack <- terra::rast(layers[HAB_PREDICTORS])
    names(hab_stack) <- HAB_PREDICTORS

    terra::writeRaster(hab_stack, out_global,
                       overwrite = TRUE, gdal = GDAL_OPTS)
    log_line("download_hab", "-", "wrote hab_static_global.tif (global)")

    terra::writeRaster(terra::crop(hab_stack, PROJ_EXTENT), out_proj,
                       overwrite = TRUE, gdal = GDAL_OPTS)
    log_line("download_hab", "-",
             sprintf("wrote hab_static.tif (%d cells at PROJ_EXTENT)",
                     terra::ncell(terra::crop(hab_stack, PROJ_EXTENT))))
  }
}

log_line("download_hab", "-", "done")
cat("\nHabitat prep complete.\n")
cat("  Global: ", out_global, "\n")
cat("  Proj:   ", out_proj,   "\n")
cat("\nNext: source(\"03_build_training.R\") to rebuild training tables\n")
cat("      with the expanded predictor set.\n")

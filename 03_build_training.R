# =============================================================================
# 03_build_training.R  —  habitat-augmented (final_hab_v1)
#
# WHAT CHANGED vs. the climate-only version
# -------------------------------------------------------------------------
# 1. env stack: load_env_global() (8 predictors) instead of the 4-layer
#    climate-only raster.
#
# 2. *** BACKGROUND / CALIBRATION REGION FIX (2026-08-22) ***
#    build_calibration_region() now unions per-presence buffers (an accessible
#    area) instead of taking a single convex hull. A convex hull over a species
#    with disjunct native + invaded ranges (e.g. Europe + North America) bridges
#    the ocean between them -- and in lon/lat can wrap the long way around the
#    globe -- so background got scattered across open ocean (dropped as NA later
#    -> n_bg far below N_BACKGROUND) and across far-flung, trivially unsuitable
#    regions (which inflate AUC/TSS). Buffered presences keep background near
#    where the species actually occurs, in BOTH ranges. This restores the
#    accessible-area approach.
#
#    NOTE: BUFFER_KM in 00_config.R is now the per-presence accessible RADIUS,
#    not a pad around a hull. Set it to ~100-200 km (100 = tightest / hardest
#    negatives; 200 = a bit more environmental coverage). 500 is too loose now.
#
# The training table schema (pa, x, y, range_type, <PREDICTORS>) is unchanged.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(rgbif)
  library(CoordinateCleaner)
  library(usdm)
  library(dplyr)
  library(sf)
})

# ── ACCESSIBLE-AREA RADIUS (the background fix) ───────────────────────────────
# Background is drawn within this many km of a presence (union of per-presence
# buffers -> build_calibration_region below). This REPLACES the old convex hull.
# Sourced from BUFFER_KM in 00_config.R so the config file is the single source
# of truth. 100 = tightest/hardest negatives (default); 200 = a little more
# environmental coverage. Smaller = stricter test. Both defensible.
ACCESS_KM <- as.integer(BUFFER_KM)

log_line("build_training", "-", "starting")
log_line("build_training", "-", sprintf("accessible-area radius: %d km", ACCESS_KM))
log_line("build_training", "-",
         sprintf("predictors (%d): %s",
                 length(PREDICTORS), paste(PREDICTORS, collapse = ", ")))

# ── LOAD GLOBAL ENV STACK (climate + habitat) ────────────────────────────────
env_global <- load_env_global()

# ── HELPERS ──────────────────────────────────────────────────────────────────
find_cache <- function(sid, tag) {
  cand <- c(
    file.path(GBIF_DIR,        paste0(sid, "_gbif_", tag, ".rds")),
    file.path(GBIF_DIR_LEGACY, paste0(sid, "_gbif_", tag, ".rds"))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) return(NA_character_) else hit[1]
}

load_gbif <- function(path, range_type) {
  d <- readRDS(path)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("decimalLongitude" %in% names(d))
    d <- dplyr::rename(d, longitude = decimalLongitude,
                       latitude = decimalLatitude)
  d$range_type <- range_type
  d
}

# ACCESSIBLE AREA = union of per-presence buffers, computed in an equal-area
# CRS. Follows the species' ACTUAL distribution -- including disjunct native and
# invaded ranges -- instead of a single convex hull. A hull over presences on
# two continents bridges the ocean between them (and, in lon/lat, can wrap the
# long way around the globe), flooding background with ocean cells (dropped
# later as NA -> the n_bg shortfall) and far-flung, trivially unsuitable
# negatives (which inflate AUC/TSS). Buffered points don't.
build_calibration_region <- function(occ, buffer_km) {
  pts     <- terra::vect(occ, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  pts_ea  <- terra::project(pts, "ESRI:54034")          # equal-area (metres)
  buff_ea <- terra::buffer(pts_ea, width = buffer_km * 1000)
  buff_ea <- terra::aggregate(buff_ea)                  # dissolve to one (multi)polygon
  bv      <- terra::project(buff_ea, "EPSG:4326")
  bv$id   <- 1L
  bv
}

# ── OPTIONAL: PRE-LOAD TGB POINTS ────────────────────────────────────────────
tgb_points <- NULL
if (USE_TGB) {
  tgb_v3     <- file.path(GBIF_DIR,        "tracheophyta_tgb.rds")
  tgb_legacy <- file.path(GBIF_DIR_LEGACY, "tracheophyta_tgb.rds")
  tgb_path <- if (file.exists(tgb_v3)) tgb_v3
  else if (file.exists(tgb_legacy)) tgb_legacy
  else NA_character_
  if (!is.na(tgb_path)) {
    tgb_points <- readRDS(tgb_path)
    log_line("build_training", "-",
             sprintf("TGB pool loaded: %d locations", nrow(tgb_points)))
  } else {
    log_line("build_training", "-",
             "WARN: USE_TGB=TRUE but no TGB cache found; falling back to random")
  }
}

# ── MAIN LOOP ────────────────────────────────────────────────────────────────
for (sp in SPECIES) {
  sid     <- sp$short_id
  out_rds <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (file.exists(out_rds)) {
    log_line("build_training", sid, "training table cached; skipping")
    next
  }
  
  log_line("build_training", sid, sprintf("=== %s ===", sp$latin))
  
  # 1. LOAD ---------------------------------------------------------------
  native_cache  <- find_cache(sid, "native")
  invaded_cache <- find_cache(sid, "invaded")
  if (is.na(native_cache)) {
    log_line("build_training", sid, "NO native cache; skipping"); next
  }
  occ_native  <- load_gbif(native_cache, "native")
  occ_invaded <- if (!is.na(invaded_cache)) load_gbif(invaded_cache, "invaded")
  else NULL
  
  if (is.null(occ_native)) {
    log_line("build_training", sid, "empty native cache; skipping"); next
  }
  log_line("build_training", sid,
           sprintf("native=%d  invaded(USA+CAN+…)=%d",
                   nrow(occ_native),
                   if (is.null(occ_invaded)) 0 else nrow(occ_invaded)))
  
  # 2. BASIC FILTERS ------------------------------------------------------
  occ <- dplyr::bind_rows(occ_native, occ_invaded) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude),
                  longitude >= -180, longitude <= 180,
                  latitude  >=  -90, latitude  <=  90) %>%
    dplyr::distinct(longitude, latitude, .keep_all = TRUE)
  
  if ("basisOfRecord" %in% names(occ))
    occ <- dplyr::filter(occ,
                         !basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))
  if ("coordinateUncertaintyInMeters" %in% names(occ))
    occ <- dplyr::filter(occ,
                         is.na(coordinateUncertaintyInMeters) |
                           coordinateUncertaintyInMeters <= 10000)
  
  log_line("build_training", sid,
           sprintf("after basic filter: %d", nrow(occ)))
  if (nrow(occ) < MIN_PRESENCES) {
    log_line("build_training", sid, "TOO FEW; skipping"); next
  }
  
  # 3. COORDINATECLEANER --------------------------------------------------
  cc_cols <- c("longitude", "latitude")
  if ("countryCode" %in% names(occ)) cc_cols <- c(cc_cols, "countryCode")
  occ_cc <- as.data.frame(occ[, cc_cols, drop = FALSE])
  flags <- tryCatch(
    CoordinateCleaner::clean_coordinates(
      x = occ_cc, lon = "longitude", lat = "latitude",
      countries = if ("countryCode" %in% names(occ_cc)) "countryCode" else NULL,
      tests = c("capitals", "centroids", "equal", "gbif",
                "institutions", "seas", "zeros"),
      value = "flagged", verbose = FALSE
    ),
    error = function(e) {
      log_line("build_training", sid,
               sprintf("cleaner warning: %s", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(flags)) occ <- occ[flags, , drop = FALSE]
  log_line("build_training", sid,
           sprintf("after CoordinateCleaner: %d", nrow(occ)))
  if (nrow(occ) < MIN_PRESENCES) {
    log_line("build_training", sid, "TOO FEW; skipping"); next
  }
  
  # 4. CALIBRATION REGION (accessible area = union of buffered presences) --
  #    ACCESS_KM (set at top of this script) is the per-presence radius.
  cal_region <- tryCatch(
    build_calibration_region(occ, ACCESS_KM),
    error = function(e) {
      log_line("build_training", sid,
               sprintf("buffer failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(cal_region)) next
  
  pts    <- terra::vect(occ, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  inside <- !is.na(terra::extract(cal_region, pts)[, 2])
  occ    <- occ[inside, , drop = FALSE]
  pts    <- pts[inside, ]
  log_line("build_training", sid,
           sprintf("inside %dkm accessible area: %d", ACCESS_KM, nrow(occ)))
  if (nrow(occ) < MIN_PRESENCES) {
    log_line("build_training", sid, "TOO FEW; skipping"); next
  }
  
  # 5. SPATIAL THINNING ---------------------------------------------------
  deg <- THIN_KM / 111
  occ$cellx <- floor(occ$longitude / deg)
  occ$celly <- floor(occ$latitude  / deg)
  occ <- occ %>%
    dplyr::distinct(cellx, celly, .keep_all = TRUE) %>%
    dplyr::select(-cellx, -celly)
  log_line("build_training", sid,
           sprintf("after ~%dkm thinning: %d", THIN_KM, nrow(occ)))
  if (nrow(occ) < MIN_PRESENCES) {
    log_line("build_training", sid, "TOO FEW; skipping"); next
  }
  pres <- terra::vect(occ, geom = c("longitude", "latitude"),
                      crs = "EPSG:4326")
  
  # 6. BACKGROUND ---------------------------------------------------------
  if (USE_TGB && !is.null(tgb_points) && nrow(tgb_points) > 0) {
    tgb_v      <- terra::vect(tgb_points, geom = c("longitude", "latitude"),
                              crs = "EPSG:4326")
    tgb_inside <- !is.na(terra::extract(cal_region, tgb_v)[, 2])
    tgb_cal    <- tgb_v[tgb_inside, ]
    n_avail    <- length(tgb_cal)
    if (n_avail < N_BACKGROUND) {
      bg <- tgb_cal
      log_line("build_training", sid,
               sprintf("background(TGB): using all %d available", n_avail))
    } else {
      bg <- tgb_cal[sample.int(n_avail, N_BACKGROUND), ]
      log_line("build_training", sid,
               sprintf("background(TGB): sampled %d of %d", N_BACKGROUND, n_avail))
    }
  } else {
    # Random background within the accessible area. Cropping the env template
    # to the (now near-presence) region keeps this fast and off the ocean.
    poly_env <- terra::crop(env_global[[1]], cal_region)
    poly_r   <- terra::rasterize(cal_region, poly_env, field = 1)
    bg <- terra::spatSample(poly_r, size = N_BACKGROUND, method = "random",
                            as.points = TRUE, na.rm = TRUE, exhaustive = TRUE)
    log_line("build_training", sid,
             sprintf("background(random): %d points", length(bg)))
  }
  
  # 7. EXTRACT ENV + ASSEMBLE TRAINING TABLE ------------------------------
  ex_pres <- terra::extract(env_global, pres)[, PREDICTORS]
  ex_bg   <- terra::extract(env_global, bg  )[, PREDICTORS]
  xy_pres <- terra::crds(pres); xy_bg <- terra::crds(bg)
  
  train <- rbind(
    data.frame(pa = 1L, x = xy_pres[, 1], y = xy_pres[, 2],
               range_type = if ("range_type" %in% names(occ))
                 occ$range_type else "invaded",
               ex_pres),
    data.frame(pa = 0L, x = xy_bg[, 1],   y = xy_bg[, 2],
               range_type = "background",
               ex_bg)
  )
  # Drop rows with any NA in predictors — mostly ocean/coastal cells where
  # SoilGrids returns NA. With the accessible-area fix this should be a small
  # fraction now, not a big chunk of the background.
  n_before <- nrow(train)
  train <- train[complete.cases(train[, PREDICTORS]), ]
  n_dropped <- n_before - nrow(train)
  log_line("build_training", sid,
           sprintf("training rows: %d (presences=%d, background=%d)  [%d dropped for NA]",
                   nrow(train), sum(train$pa == 1L), sum(train$pa == 0L), n_dropped))
  
  # 8. VIF / CORRELATION REPORT -------------------------------------------
  vif_tbl <- tryCatch(usdm::vif(train[, PREDICTORS]),
                      error = function(e) NULL)
  cor_mat <- cor(train[, PREDICTORS])
  writeLines(capture.output({
    cat("Species:  ", sid, " (", sp$latin, ")\n", sep = "")
    cat("Predictors: ", paste(PREDICTORS, collapse = ", "), "\n", sep = "")
    cat("Presences: ", sum(train$pa == 1L),
        "   Background: ", sum(train$pa == 0L), "\n\n", sep = "")
    if (!is.null(vif_tbl)) {
      cat("VIF (target <", VIF_THRESHOLD, "):\n"); print(vif_tbl)
    } else {
      cat("VIF failed (usually rank-deficient — inspect correlations)\n")
    }
    cat("\nPairwise correlation (warn if |r| >", COR_THRESHOLD, "):\n")
    print(round(cor_mat, 2))
    cm <- cor_mat; diag(cm) <- 0
    idx <- which(abs(cm) == max(abs(cm), na.rm = TRUE), arr.ind = TRUE)[1, ]
    cat(sprintf("\nHighest |r|: %s <-> %s = %.2f\n",
                rownames(cor_mat)[idx[1]], colnames(cor_mat)[idx[2]],
                cor_mat[idx[1], idx[2]]))
  }), file.path(DIAG_DIR, paste0(sid, "_vif.txt")))
  
  if (!is.null(vif_tbl) && any(vif_tbl$VIF > VIF_THRESHOLD, na.rm = TRUE)) {
    hi <- vif_tbl$Variables[vif_tbl$VIF > VIF_THRESHOLD]
    log_line("build_training", sid,
             sprintf("WARN: VIF > %d for %s (see %s_vif.txt)",
                     VIF_THRESHOLD, paste(hi, collapse = ", "), sid))
  }
  
  # 9. SAVE ---------------------------------------------------------------
  pred_range <- apply(train[train$pa == 1L, PREDICTORS], 2,
                      function(v) range(v, na.rm = TRUE))
  rownames(pred_range) <- c("min", "max")
  
  saveRDS(list(
    train        = train,
    meta         = sp,
    pred_range   = pred_range,
    n_pres       = sum(train$pa == 1L),
    n_bg         = sum(train$pa == 0L),
    predictors   = PREDICTORS,
    clim_preds   = CLIM_PREDICTORS,
    hab_preds    = HAB_PREDICTORS,
    buffer_km    = ACCESS_KM,
    thin_km      = THIN_KM,
    use_tgb      = USE_TGB,
    pipeline_tag = PIPELINE_TAG
  ), out_rds)
  log_line("build_training", sid,
           sprintf("wrote %s (%d rows)", basename(out_rds), nrow(train)))
}

log_line("build_training", "-", "done")
cat("\nTraining tables complete. Ready for 04_fit_ensemble.R\n")
# =============================================================================
# 00_config.R  —  HABITAT-AUGMENTED (tag: final_hab_v1)
#
# WHAT'S NEW vs. your previous 00_config.R
# ---------------------------------------------------------------------------
# 1. Predictor set doubled from 4 (climate) to 8 (climate + habitat):
#      CLIMATE (dynamic, changes across CMIP6 futures):
#        mat     BIO1   annual mean temperature
#        trng    BIO7   temperature annual range
#        ap      BIO12  annual precipitation
#        pseas   BIO15  precipitation seasonality (CV%)
#      HABITAT (static across futures — see note below):
#        elev    elevation (m)                    — SRTM/GMTED via geodata
#        hfp     Human Footprint Index (0–50)     — Venter et al. 2018
#        treecov tree cover fraction (0–100)      — Copernicus GLC 100 m
#        soilph  soil pH * 10, top 5 cm           — SoilGrids v2
#
# 2. All output directories get a "_hab" suffix (MODELS_DIR = "models_hab",
#    CURRENT_DIR = "output/current_hab", etc.) so your existing v2a/final
#    outputs stay put. Nothing under models/, output/, eval/, figures/,
#    webapp_output/ gets touched by this rebuild.
#
# 3. New helper functions load_env_global() / load_env_current() /
#    load_env_future(gcm, ssp, label) build the combined climate+habitat
#    SpatRaster stack that downstream scripts (03, 04, 05) consume.
#
# ---------------------------------------------------------------------------
# WHY THE HABITAT PREDICTORS ARE STATIC ACROSS FUTURES
# ---------------------------------------------------------------------------
# Only WorldClim futures carry a time signal. Soil pH, elevation, tree cover,
# and Human Footprint are held at their present-day values for every future
# scenario. This is standard practice in invasive SDMs that add non-climate
# predictors (Ibáñez et al. 2020 Frontiers, Jarnevich et al. 2024 Ecology
# Letters). Rationale: at ~9 km resolution over 60-year horizons, climate
# change dwarfs plausible change in soil pH, elevation, or aggregated land
# cover. HFP does have SSP-based projections available (Kennedy et al. 2019)
# if you ever want to relax this — that would need a per-scenario HFP grid
# and adjustments to load_env_future().
# =============================================================================

# ── PROJECT ROOTS ────────────────────────────────────────────────────────────
PROJECT_ROOT <- "C:/Users/Joelle/OneDrive/Desktop/INVASIVE-SDM-PREDICTION-FINAL-PIPELINE"
WEBAPP_ROOT  <- "C:/Users/Joelle/OneDrive/Desktop/webapp"
WEBAPP_DATA  <- file.path(WEBAPP_ROOT, "data")

PIPELINE_TAG <- "final_hab_v1"
HAB_SUFFIX   <- "_hab"     # applied to output dirs so we don't clobber v2a/final

# ── DATA DIRS ────────────────────────────────────────────────────────────────
DATA_DIR       <- "data"
CLIMATE_DIR    <- "data/climate"
HAB_DIR        <- "data/habitat"       # NEW — populated by 01b_download_habitat.R
GBIF_DIR       <- "data/gbif"
GBIF_DIR_LEGACY <- GBIF_DIR

# ── PIPELINE OUTPUT DIRS (all suffixed with _hab) ────────────────────────────
CACHE_DIR         <- "cache"                            # geodata/GADM caches (shared)
MODELS_DIR        <- paste0("models",         HAB_SUFFIX)
CURRENT_DIR       <- paste0("output/current", HAB_SUFFIX)
FUTURE_DIR        <- paste0("output/future",  HAB_SUFFIX)
EVAL_DIR          <- paste0("eval",           HAB_SUFFIX)
DIAG_DIR          <- paste0("diag",           HAB_SUFFIX)
PARTIAL_DIR       <- paste0("figures/partial",HAB_SUFFIX)
MAPS_DIR          <- paste0("figures/maps",   HAB_SUFFIX)
SUMMARY_DIR       <- paste0("summary",        HAB_SUFFIX)
LOG_DIR           <- "logs"
LOG_FILE          <- normalizePath(file.path(LOG_DIR, "pipeline_hab.log"),
                                   winslash = "/", mustWork = FALSE)

WEBAPP_OUTPUT_DIR <- paste0("webapp_output",  HAB_SUFFIX)

for (d in c(DATA_DIR, CLIMATE_DIR, HAB_DIR, GBIF_DIR,
            CACHE_DIR, MODELS_DIR,
            CURRENT_DIR, FUTURE_DIR, EVAL_DIR, DIAG_DIR,
            PARTIAL_DIR, MAPS_DIR, SUMMARY_DIR, LOG_DIR,
            WEBAPP_OUTPUT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── CLIMATE PREDICTORS (WorldClim 2.1 bioclim) ──────────────────────────────
# BIO1  = annual mean temperature          → mat
# BIO7  = temperature annual range         → trng
# BIO12 = annual precipitation             → ap
# BIO15 = precipitation seasonality (CV%)  → pseas
BIO_IDX          <- c(1, 7, 12, 15)
CLIM_PREDICTORS  <- c("mat", "trng", "ap", "pseas")
WORLDCLIM_RES    <- 5

# ── HABITAT PREDICTORS (all static across futures) ──────────────────────────
HAB_PREDICTORS   <- c("elev", "hfp", "treecov", "soilph")

# The combined predictor set used everywhere downstream. ORDER MATTERS —
# raster layer names in hab_static.tif and clim_current.tif must match this
# order after they're stacked.
PREDICTORS       <- c(CLIM_PREDICTORS, HAB_PREDICTORS)

# ── PROJECTION EXTENT ────────────────────────────────────────────────────────
PROJ_EXTENT      <- terra::ext(-141, -52, 24, 75)   # US + Canada

# ── SPECIES LIST ─────────────────────────────────────────────────────────────
# Invaded ranges cover USA + CAN for every species; knotweed and
# tree_of_heaven also carry European invaded records.
SPECIES <- list(
  list(short_id = "autumn_olive",       latin = "Elaeagnus umbellata",
       common_name = "Autumn olive",
       native_iso  = c("CHN", "JPN", "KOR"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "bittersweet",        latin = "Celastrus orbiculatus",
       common_name = "Oriental bittersweet",
       native_iso  = c("CHN", "JPN", "KOR"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "burning_bush",       latin = "Euonymus alatus",
       common_name = "Winged burning bush",
       native_iso  = c("CHN", "JPN", "KOR"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "common_buckthorn",   latin = "Rhamnus cathartica",
       common_name = "Common buckthorn",
       native_iso  = c("GBR", "FRA", "DEU", "POL", "ITA", "ESP",
                       "AUT", "CZE", "BEL", "NLD", "SWE", "DNK",
                       "CHE", "HUN", "ROU", "BGR"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "garlic_mustard",     latin = "Alliaria petiolata",
       common_name = "Garlic mustard",
       native_iso  = c("GBR", "FRA", "DEU", "POL", "ITA", "ESP",
                       "AUT", "CZE", "BEL", "NLD", "SWE", "DNK"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "japanese_barberry",  latin = "Berberis thunbergii",
       common_name = "Japanese barberry",
       native_iso  = c("JPN", "KOR"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "knotweed",           latin = "Reynoutria japonica",
       common_name = "Japanese knotweed",
       native_iso  = c("JPN", "KOR", "CHN"),
       invaded_iso = c("USA", "CAN", "GBR", "FRA", "DEU")),
  
  list(short_id = "multiflora_rose",    latin = "Rosa multiflora",
       common_name = "Multiflora rose",
       native_iso  = c("JPN", "KOR", "CHN"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "norway_maple",       latin = "Acer platanoides",
       common_name = "Norway maple",
       native_iso  = c("GBR", "FRA", "DEU", "POL", "ITA", "ESP",
                       "AUT", "CZE", "BEL", "NLD", "SWE", "DNK",
                       "NOR", "FIN"),
       invaded_iso = c("USA", "CAN")),
  
  list(short_id = "stiltgrass",         latin = "Microstegium vimineum",
       common_name = "Japanese stiltgrass",
       native_iso  = c("CHN", "JPN", "KOR", "IND", "NPL", "VNM"),
       invaded_iso = c("USA")),
  
  list(short_id = "tree_of_heaven",     latin = "Ailanthus altissima",
       common_name = "Tree of heaven",
       native_iso  = c("CHN", "TWN"),
       invaded_iso = c("USA", "CAN", "ITA", "FRA", "ESP", "DEU"))
)

SPECIES_IDS <- vapply(SPECIES, `[[`, character(1), "short_id")

# ── CLIMATE SCENARIOS (CMIP6) ────────────────────────────────────────────────
GCMS <- c("ACCESS-CM2", "CNRM-CM6-1", "MPI-ESM1-2-HR",
          "MRI-ESM2-0",  "IPSL-CM6A-LR")
SSPS          <- c("245", "585")
PERIODS       <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
PERIOD_LABELS <- c("2030",      "2050",      "2070",      "2090")

KEEP_PER_GCM_RASTERS <- FALSE

# ── TRAINING-DATA KNOBS ──────────────────────────────────────────────────────
# BUFFER_KM is the per-presence accessible-area RADIUS used by
# 03_build_training.R::build_calibration_region (union of per-presence buffers,
# equal-area). 100 = tightest/hardest negatives (default); 200 = a little more
# environmental coverage. The old convex-hull-plus-pad meaning is gone.
BUFFER_KM     <- 100
USE_TGB       <- FALSE
TGB_LIMIT     <- 100000
N_BACKGROUND  <- 10000
THIN_KM       <- 5
MIN_PRESENCES <- 30

VIF_THRESHOLD <- 5      # with 8 predictors this WILL fire on some species;
COR_THRESHOLD <- 0.7    # the diag report tells you which pair(s) to worry about

# ── biomod2 MODELING OPTIONS ─────────────────────────────────────────────────
BIOMOD_MODELS     <- c("GAM", "GLM", "GBM", "RF", "MAXNET")
# biomod2's "block" strategy wraps ENMeval::get.block, which always produces
# 4 data-driven lon/lat quadrants. CV.k and a block-size argument are ignored.
# Keep CV_STRATEGY here as documentation; a true k-fold spatial scheme would
# require CV.strategy="user.defined" with blockCV::cv_spatial().
CV_STRATEGY       <- "block"
METRIC_EVAL       <- c("TSS", "AUCroc")
METRIC_SELECT     <- "TSS"
METRIC_THRESH     <- 0.5
BIOMOD_EM_ALGO    <- c("EMmean", "EMmedian", "EMwmean", "EMca", "EMcv")
DEPLOYED_ENSEMBLE <- "EMwmeanByTSS"

# ── DEFAULT BINARY THRESHOLD (fallback if 06_thresholds hasn't run) ─────────
THRESH <- 0.5

# ── GAM legacy options ──────────────────────────────────────────────────────
GAM_K       <- 5
GAM_METHOD  <- "REML"
GAM_SELECT  <- TRUE
GAM_GAMMA   <- 1.4

# ── MANAGEMENT-UNIT SETTINGS (script 13) ────────────────────────────────────
DROP_STATES <- c("Alaska", "Hawaii", "Puerto Rico",
                 "United States Virgin Islands",
                 "Northern Mariana Islands", "Guam", "American Samoa")

# ── RASTER I/O ───────────────────────────────────────────────────────────────
GDAL_OPTS <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")

# ── LOGGING ──────────────────────────────────────────────────────────────────
log_line <- function(stage, sid, msg = "") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("%s  [%s]  %s  %s\n", ts, stage, sid, msg)
  cat(line)
  cat(line, file = LOG_FILE, append = TRUE)
}

# =============================================================================
# ENV-STACK HELPERS
# ---------------------------------------------------------------------------
# Everything that consumes predictors (03 extraction, 04 fit-time raster ref,
# 05 projection) goes through these three functions. That way if you change
# the predictor set, add a new source, or swap in dynamic HFP futures later,
# there's exactly one place to edit.
#
# Each returns a SpatRaster with layer order = PREDICTORS.
# =============================================================================

# terra needs to be loaded before these are called — usually the calling
# script does that itself, but require it here defensively so sourcing
# 00_config.R in a fresh session doesn't break the helper definitions.
if (requireNamespace("terra", quietly = TRUE)) {
  suppressPackageStartupMessages(library(terra))
}

# Resolve pipeline-relative paths to absolute paths anchored at PROJECT_ROOT.
# This makes the load_env_* helpers safe to call from ANY working directory.
#
# WHY THIS MATTERS: script 05 calls `setwd(MODELS_DIR)` inside its per-species
# loop to keep biomod2 happy. If load_env_current() then used relative paths
# like `data/climate/clim_current.tif`, they would resolve to
# `models_hab/data/climate/clim_current.tif` — which doesn't exist — and 05
# would fail on every species with a "Missing … run 01_download_climate.R
# first" error even though the file is fine. This bug hit on 2026-08-20 after
# a full 8-hour step-04 run wasted its projections.
#
# Fix: build absolute paths from PROJECT_ROOT (which IS absolute, line 42).
# Kept as a tiny helper so all three load_env_* functions share the same
# resolver and future maintainers don't have to remember the pattern.
.abs <- function(...) file.path(PROJECT_ROOT, ...)

# Load the full env stack, GLOBAL scope. Used by 03 for niche extraction at
# presence + background points, which live in native + invaded ranges anywhere
# on Earth.
load_env_global <- function() {
  clim_path <- .abs(CLIMATE_DIR, "clim_current_global.tif")
  hab_path  <- .abs(HAB_DIR,     "hab_static_global.tif")
  if (!file.exists(clim_path))
    stop("Missing ", clim_path, " — run 01_download_climate.R first.")
  if (!file.exists(hab_path))
    stop("Missing ", hab_path,  " — run 01b_download_habitat.R first.")
  
  clim <- terra::rast(clim_path); names(clim) <- CLIM_PREDICTORS
  hab  <- terra::rast(hab_path);  names(hab)  <- HAB_PREDICTORS
  
  # 01b resamples habitat onto the global climate grid, so this stack is aligned.
  # Defensive check in case someone regenerates one but not the other:
  if (!terra::compareGeom(clim, hab, stopOnError = FALSE)) {
    warning("hab_static_global geometry disagrees with clim_current_global — resampling")
    hab <- terra::resample(hab, clim, method = "average")
  }
  s <- c(clim, hab)
  names(s) <- PREDICTORS
  s
}

# Load the full env stack for the CURRENT scenario, cropped to PROJ_EXTENT.
# Used by 04 (for biomod2's model-options niche-context) and 05 (current proj).
load_env_current <- function() {
  clim_path <- .abs(CLIMATE_DIR, "clim_current.tif")
  hab_path  <- .abs(HAB_DIR,     "hab_static.tif")
  if (!file.exists(clim_path))
    stop("Missing ", clim_path, " — run 01_download_climate.R first.")
  if (!file.exists(hab_path))
    stop("Missing ", hab_path,  " — run 01b_download_habitat.R first.")
  
  clim <- terra::rast(clim_path); names(clim) <- CLIM_PREDICTORS
  hab  <- terra::rast(hab_path);  names(hab)  <- HAB_PREDICTORS
  
  if (!terra::compareGeom(clim, hab, stopOnError = FALSE)) {
    warning("hab_static geometry disagrees with clim_current — resampling")
    hab <- terra::resample(hab, clim, method = "average")
  }
  s <- c(clim, hab)
  names(s) <- PREDICTORS
  s
}

# Load the full env stack for a FUTURE scenario (climate varies, habitat static).
# Used by 05.
load_env_future <- function(gcm, ssp, label) {
  clim_path <- .abs(CLIMATE_DIR, "future",
                    sprintf("%s_ssp%s_%s.tif", gcm, ssp, label))
  hab_path  <- .abs(HAB_DIR, "hab_static.tif")
  if (!file.exists(clim_path))
    stop("Missing ", clim_path, " — run 01_download_climate.R first.")
  if (!file.exists(hab_path))
    stop("Missing ", hab_path,  " — run 01b_download_habitat.R first.")
  
  clim <- terra::rast(clim_path); names(clim) <- CLIM_PREDICTORS
  hab  <- terra::rast(hab_path);  names(hab)  <- HAB_PREDICTORS
  
  if (!terra::compareGeom(clim, hab, stopOnError = FALSE)) {
    hab <- terra::resample(hab, clim, method = "average")
  }
  s <- c(clim, hab)
  names(s) <- PREDICTORS
  s
}

# ── SANITY REPORT ────────────────────────────────────────────────────────────
cat(sprintf("Pipeline config loaded (tag = %s)\n", PIPELINE_TAG))
cat(sprintf("  project root:    %s\n", PROJECT_ROOT))
cat(sprintf("  climate dir:     %s\n", CLIMATE_DIR))
cat(sprintf("  habitat dir:     %s\n", HAB_DIR))
cat(sprintf("  gbif dir:        %s\n", GBIF_DIR))
cat(sprintf("  models dir:      %s\n", MODELS_DIR))
cat(sprintf("  output current:  %s\n", CURRENT_DIR))
cat(sprintf("  output future:   %s\n", FUTURE_DIR))
cat(sprintf("  webapp mirror:   %s\n", WEBAPP_OUTPUT_DIR))
cat(sprintf("  webapp target:   %s\n", WEBAPP_DATA))
cat(sprintf("  species:         %d\n", length(SPECIES)))
cat(sprintf("  predictors:      %s\n", paste(PREDICTORS, collapse = ", ")))
cat(sprintf("  ensemble algos:  %s\n", paste(BIOMOD_MODELS, collapse = ", ")))

CONFIG_V3_LOADED <- TRUE
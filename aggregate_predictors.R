# ---- aggregate_predictors.R --------------------------------------------------
# Aggregate the model's INPUT rasters (climate + habitat) to per-county means so
# the atlas can show the climate/habitat conditions at a clicked county (in the
# result panel, under the predicted-suitability readout). All 8 predictors.
# Writes ONE shared, species-independent file: data/predictors.json, id-keyed to
# match map.js. Run after aggregate_counties.R. Paste into the RStudio Console.
# -----------------------------------------------------------------------------
if (!file.exists("00_config.R"))
  stop("Run from the pipeline root (where 00_config.R lives).")
if (!exists("CONFIG_V3_LOADED")) source("00_config.R")
suppressPackageStartupMessages({
  library(terra); library(sf); library(exactextractr); library(jsonlite)
})

counties_path <- file.path(WEBAPP_DATA, "counties.geojson")
if (!file.exists(counties_path))
  stop("Missing ", counties_path, " — run aggregate_counties.R first.")
counties <- sf::st_read(counties_path, quiet = TRUE)
ids <- as.character(counties$id)
cat("counties:", length(ids), "\n")

# Input rasters (same ones the SDM was trained/projected on), current climate.
clim <- terra::rast(file.path(CLIMATE_DIR, "clim_current.tif")); names(clim) <- CLIM_PREDICTORS
hab  <- terra::rast(file.path(HAB_DIR,     "hab_static.tif"));    names(hab)  <- HAB_PREDICTORS
env  <- c(clim, hab)

# Predictors shown in the click panel (must match PREDICTOR_LAYERS in config.js).
want <- c("mat", "ap", "trng", "pseas", "hfp", "treecov", "soilph", "elev")

out <- list()
for (p in want) {
  if (!(p %in% names(env))) { cat("skip (layer not in rasters):", p, "\n"); next }
  m <- exactextractr::exact_extract(env[[p]], counties, "mean", progress = FALSE)
  m <- round(m, 3); m[is.nan(m)] <- NA
  out[[p]] <- setNames(as.list(m), ids)          # id -> value (matches map.js)
  cat("aggregated:", p, "  (range",
      round(min(m, na.rm = TRUE), 2), "..", round(max(m, na.rm = TRUE), 2), ")\n")
}

writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, na = "null", digits = 3),
           file.path(WEBAPP_DATA, "predictors.json"))
cat("\nWrote", file.path(WEBAPP_DATA, "predictors.json"), "\n")
cat("Ranges above are your real data — if they differ a lot from the min/max in\n",
    "config.js PREDICTOR_LAYERS, tweak those so the color ramp uses the full scale.\n", sep = "")
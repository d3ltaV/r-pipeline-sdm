# =============================================================================
# aggregate_counties.R
# Turn per-pixel suitability rasters into per-county means so the webapp can
# render EDDMapS-style choropleths. Also pre-rasterizes the county layer
# onto the sim grid so the simulator can do fast pixel→county lookups.
#
# Outputs (into WEBAPP_DATA):
#   counties.geojson        shared polygon layer (~3-5 MB, simplified)
#   counties_grid.bin       UInt16 raster: pixel value = 1-based county index
#   counties_grid.json      { width, height, xmin, xmax, ymin, ymax, n_counties }
#   <sid>/counties.json     per-species map of scenario_key -> county-mean array
#
# One-shot: run once after 12_copy_to_webapp.R. Idempotent — reuses cached
# counties files across reruns.
#
# Usage:  source("aggregate_counties.R")
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

need <- c("terra", "sf", "exactextractr", "jsonlite", "geodata", "rmapshaper")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(terra); library(sf); library(exactextractr)
  library(jsonlite); library(geodata); library(rmapshaper)
})

log_line("counties", "-", "starting")

counties_out      <- file.path(WEBAPP_DATA, "counties.geojson")
counties_grid_out <- file.path(WEBAPP_DATA, "counties_grid.bin")
counties_info_out <- file.path(WEBAPP_DATA, "counties_grid.json")

dir.create(WEBAPP_DATA, recursive = TRUE, showWarnings = FALSE)

# Grid template: any species' current.tif shares the grid
species_dirs <- list.dirs(WEBAPP_DATA, recursive = FALSE, full.names = TRUE)
species_dirs <- species_dirs[!grepl("(_gadm_cache|/data$)", species_dirs)]
template_tif <- NULL
for (sp_dir in species_dirs) {
  cand <- file.path(sp_dir, "current.tif")
  if (file.exists(cand)) { template_tif <- cand; break }
}
if (is.null(template_tif)) {
  # Fall back to pipeline output
  cand <- list.files(CURRENT_DIR, pattern = "_ensemble_current\\.tif$",
                     full.names = TRUE)
  if (length(cand)) template_tif <- cand[1]
}
if (is.null(template_tif)) stop("No current.tif found for grid template.")
ref <- terra::rast(template_tif)
cat(sprintf("Grid template: %s  (%d x %d)\n", template_tif,
            terra::ncol(ref), terra::nrow(ref)))

# ── STEP 1: shared counties layer + rasterized grid (one-time) ───────────────
need_rebuild <- !file.exists(counties_out) ||
  !file.exists(counties_grid_out) ||
  !file.exists(counties_info_out)

if (need_rebuild) {
  cat("Building US + Canada county layer (~2 min first time)...\n")
  
  gadm_cache <- file.path(dirname(WEBAPP_DATA), "_gadm_cache")
  dir.create(gadm_cache, showWarnings = FALSE, recursive = TRUE)
  
  us <- sf::st_as_sf(geodata::gadm("USA", level = 2, path = gadm_cache))
  ca <- sf::st_as_sf(geodata::gadm("CAN", level = 2, path = gadm_cache))
  
  slim <- function(x, cc) {
    sf::st_sf(
      id       = paste0(cc, "-", x$GID_2),
      name     = x$NAME_2,
      state    = x$NAME_1,
      country  = cc,
      geometry = sf::st_geometry(x)
    )
  }
  counties <- rbind(slim(us, "US"), slim(ca, "CA"))
  
  if (sf::st_crs(counties)$epsg != 4326)
    counties <- sf::st_transform(counties, 4326)
  
  # 1-based positional index; rasterize BEFORE simplifying (better pixel accuracy)
  counties$county_idx <- seq_len(nrow(counties))
  
  cat("Rasterizing county lookup grid...\n")
  counties_v <- terra::vect(counties)
  grid_r <- terra::rasterize(counties_v, ref, field = "county_idx",
                             background = 0)
  
  # Write as UInt16 binary blob
  vals <- terra::values(grid_r)[, 1]
  vals[is.na(vals)] <- 0
  vals <- as.integer(vals)
  writeBin(as.integer(vals), counties_grid_out, size = 2)
  
  info <- list(
    width      = terra::ncol(ref),
    height     = terra::nrow(ref),
    xmin       = terra::xmin(ref),
    xmax       = terra::xmax(ref),
    ymin       = terra::ymin(ref),
    ymax       = terra::ymax(ref),
    n_counties = nrow(counties)
  )
  writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE),
             counties_info_out)
  
  cat("Simplifying + writing counties.geojson ...\n")
  counties_slim <- rmapshaper::ms_simplify(counties, keep = 0.1,
                                           keep_shapes = TRUE)
  sf::st_write(counties_slim, counties_out, driver = "GeoJSON",
               delete_dsn = TRUE, quiet = TRUE)
  
  log_line("counties", "-",
           sprintf("built layer: %d counties (%.1f MB geojson)",
                   nrow(counties_slim),
                   file.info(counties_out)$size / 1e6))
} else {
  counties <- sf::st_read(counties_out, quiet = TRUE)
  counties$county_idx <- seq_len(nrow(counties))
  log_line("counties", "-",
           sprintf("using cached layer: %d counties", nrow(counties)))
}

# ── STEP 2: per-species per-scenario county means ────────────────────────────
scenarios <- c(list(list(key = "current", path_fn = function(sid)
  file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid)))))

for (ssp in SSPS) for (label in PERIOD_LABELS) {
  scenarios[[length(scenarios) + 1]] <- list(
    key = sprintf("%s_ssp%s", label, ssp),
    path_fn = function(sid, lb = label, sp = ssp)
      file.path(FUTURE_DIR,
                sprintf("%s_ensemble_%s_ssp%s_median.tif", sid, lb, sp))
  )
}

for (sp in SPECIES) {
  sid   <- sp$short_id
  spdir <- file.path(WEBAPP_DATA, sid)
  dir.create(spdir, recursive = TRUE, showWarnings = FALSE)
  out_json <- file.path(spdir, "counties.json")
  
  county_means <- list()
  for (sc in scenarios) {
    p <- sc$path_fn(sid)
    if (!file.exists(p)) next
    r <- terra::rast(p)
    if (max(terra::minmax(r)[2, ], na.rm = TRUE) > 1.5) r <- r / 1000
    means <- exactextractr::exact_extract(r, counties, "mean",
                                          progress = FALSE)
    means <- round(means, 4)
    means[is.nan(means)] <- NA
    names(means) <- counties$id           # key by county id (matches map.js activeValues[id])
    county_means[[sc$key]] <- as.list(means)
  }
  
  if (length(county_means)) {
    writeLines(jsonlite::toJSON(county_means,
                                auto_unbox = TRUE, na = "null", digits = 4),
               out_json)
    log_line("counties", sid,
             sprintf("wrote counties.json (%d scenarios, %d counties)",
                     length(county_means), nrow(counties)))
  }
}

log_line("counties", "-", "done")
cat("\nCounties output complete. Next: 12_copy_to_webapp.R to sync.\n")
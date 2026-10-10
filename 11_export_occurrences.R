# =============================================================================
# 11_export_occurrences.R
# For each species, write webapp-consumable artefacts into project root:
#
#   webapp_output_hab/<sid>/presence.tif           binary presence raster (0/1)
#   webapp_output_hab/<sid>/occurrences.geojson    invaded-range GBIF points
#   webapp_output_hab/<sid>/occurrence_meta.json   {n_native, n_invaded, ranges}
#
# 12_copy_to_webapp.R (which runs AFTER this) syncs webapp_output_hab/ into
# WEBAPP_DATA. This script writes to the project root only — never to the
# live webapp directory.
#
# Sourcing: reads <sid>_train_hab.rds (built by 03) for thinned presences;
# reads raw <sid>_gbif_invaded.rds for the un-thinned GeoJSON so the map
# shows all citizen-science points, not just the thinned training subset.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

log_line("export_occurrences", "-", "starting")

GEOJSON_DECIMALS <- 4

# ── HELPERS ──────────────────────────────────────────────────────────────────
write_presence_raster <- function(pres_xy, template, out_path) {
  if (nrow(pres_xy) == 0) {
    log_line("export_occurrences", basename(dirname(out_path)),
             "no presences in extent; skipping presence.tif")
    return(invisible(FALSE))
  }
  pv <- terra::vect(pres_xy, geom = c("x", "y"), crs = "EPSG:4326")
  pr <- terra::rasterize(pv, template, field = 1, background = 0)
  names(pr) <- "presence"
  terra::writeRaster(pr, out_path, overwrite = TRUE, gdal = GDAL_OPTS)
  n_cells <- terra::global(pr, "sum", na.rm = TRUE)[1, 1]
  log_line("export_occurrences", basename(dirname(out_path)),
           sprintf("wrote presence.tif (%d cells)", n_cells))
  invisible(TRUE)
}

# Manual GeoJSON — no geojsonio dep, control decimal precision
write_occurrences_geojson <- function(df, out_path) {
  if (nrow(df) == 0) {
    writeLines('{"type":"FeatureCollection","features":[]}', out_path)
    return(invisible(0))
  }
  lon <- round(df$x, GEOJSON_DECIMALS)
  lat <- round(df$y, GEOJSON_DECIMALS)
  rng <- if ("range" %in% names(df)) df$range else rep("invaded", nrow(df))
  feats <- vapply(seq_len(nrow(df)), function(i) {
    sprintf(
      '{"type":"Feature","geometry":{"type":"Point","coordinates":[%s,%s]},"properties":{"range":"%s"}}',
      format(lon[i], trim = TRUE, scientific = FALSE),
      format(lat[i], trim = TRUE, scientific = FALSE),
      rng[i]
    )
  }, character(1))
  writeLines(paste0('{"type":"FeatureCollection","features":[',
                    paste(feats, collapse = ","), "]}"), out_path)
  invisible(nrow(df))
}

# Load raw invaded-range GBIF (for un-thinned webapp display)
load_raw_invaded <- function(sid) {
  cand <- c(file.path(GBIF_DIR,        paste0(sid, "_gbif_invaded.rds")),
            file.path(GBIF_DIR_LEGACY, paste0(sid, "_gbif_invaded.rds")))
  hit  <- cand[file.exists(cand)]
  if (length(hit) == 0) return(NULL)
  d <- readRDS(hit[1])
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("decimalLongitude" %in% names(d)) {
    d$x <- d$decimalLongitude; d$y <- d$decimalLatitude
  } else {
    d$x <- d$longitude; d$y <- d$latitude
  }
  d <- d[!is.na(d$x) & !is.na(d$y), , drop = FALSE]
  d
}

# ── MAIN ─────────────────────────────────────────────────────────────────────
summary_rows <- list()

for (sp in SPECIES) {
  sid <- sp$short_id
  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  cur_path   <- file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid))

  if (!file.exists(train_path) || !file.exists(cur_path)) {
    log_line("export_occurrences", sid, "missing training or current; skip")
    next
  }

  spdir <- file.path(WEBAPP_OUTPUT_DIR, sid)
  dir.create(spdir, recursive = TRUE, showWarnings = FALSE)

  template <- terra::rast(cur_path)

  # Thinned presences from training (used for the raster)
  train <- readRDS(train_path)$train
  pres  <- train[train$pa == 1L, c("x", "y")]

  # Raw invaded-range points for GeoJSON (all points, not just thinned)
  raw   <- load_raw_invaded(sid)

  # Restrict to invaded range for map display
  pres_invaded <- if ("range_type" %in% names(train))
    train[train$pa == 1L & train$range_type == "invaded", c("x", "y")]
  else pres

  # Clip to raster extent
  ex <- terra::ext(template)
  in_extent <- function(df) {
    df[df$x >= ex$xmin & df$x <= ex$xmax &
       df$y >= ex$ymin & df$y <= ex$ymax, , drop = FALSE]
  }
  pres_inv_in <- in_extent(pres_invaded)

  # presence.tif from thinned invaded presences
  write_presence_raster(pres_inv_in, template,
                        file.path(spdir, "presence.tif"))

  # occurrences.geojson from raw invaded points (for busier map display)
  geojson_df <- if (!is.null(raw)) {
    d <- in_extent(raw[, c("x", "y")])
    d$range <- "invaded"
    d
  } else {
    d <- pres_inv_in; d$range <- "invaded"; d
  }
  n_geo <- write_occurrences_geojson(geojson_df, file.path(spdir, "occurrences.geojson"))
  log_line("export_occurrences", sid,
           sprintf("wrote occurrences.geojson (%d points)", n_geo))

  # Meta
  meta <- list(
    species     = sid,
    latin       = sp$latin,
    common_name = sp$common_name,
    n_train_presence = as.integer(nrow(pres)),
    n_invaded_raw    = as.integer(if (is.null(raw)) 0 else nrow(raw)),
    n_geojson_points = as.integer(n_geo),
    native_iso   = sp$native_iso,
    invaded_iso  = sp$invaded_iso,
    pipeline_tag = PIPELINE_TAG,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
             file.path(spdir, "occurrence_meta.json"))

  summary_rows[[sid]] <- data.frame(
    species = sid,
    n_train_pres = nrow(pres),
    n_geojson    = n_geo,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows)) {
  out <- do.call(rbind, summary_rows)
  write.csv(out, file.path(SUMMARY_DIR, "export_occurrences_summary.csv"),
            row.names = FALSE)
}

log_line("export_occurrences", "-", "done")
cat("\nOccurrences exported to", WEBAPP_OUTPUT_DIR, "\n")

# =============================================================================
# 13_management_units.R
# Aggregate ensemble suitability to management units (US states + Canadian
# provinces) following INHABIT methodology (Jarnevich et al. 2024 Ecol Lett).
#
# For each unit x species x time point, computes:
#   suitable_fraction    fraction of unit area above per-species maxTSS threshold
#   suitable_area_km2    suitable_fraction * unit area
#   status               "currently_observed" | "at_doorstep" | "distant"
#
# "at_doorstep" = presence within DOORSTEP_BUFFER_KM of the unit boundary
# but not yet inside — the actionable watchlist category.
#
# Output:
#   EVAL_DIR/mgmt_units_suitability.csv    long: unit x species x time_point
#   EVAL_DIR/mgmt_units_doorstep.csv       2050 watchlist rollup
#   EVAL_DIR/mgmt_units_summary.csv        counts by unit x time_point
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(geodata)
  library(dplyr)
})

log_line("mgmt_units", "-", "starting")

UNIT_SUITABLE_THRESHOLD <- 0.05    # 5% suitable area to flag a unit
DOORSTEP_BUFFER_KM      <- 50

# ── THRESHOLDS ───────────────────────────────────────────────────────────────
thr_file <- file.path(EVAL_DIR, "thresholds_summary.csv")
thr_df   <- if (file.exists(thr_file)) read.csv(thr_file) else NULL
get_threshold <- function(sid) {
  if (!is.null(thr_df) && sid %in% thr_df$species)
    thr_df[thr_df$species == sid, "maxTSS_thr"][1]
  else THRESH
}

# ── UNITS ────────────────────────────────────────────────────────────────────
log_line("mgmt_units", "-", "loading GADM states + provinces")

states <- terra::vect(geodata::gadm("USA", level = 1, path = CACHE_DIR))
states <- states[!(states$NAME_1 %in% DROP_STATES), ]
states$unit_id   <- paste0("STATE_", states$NAME_1)
states$unit_type <- "state"
states$unit_name <- states$NAME_1

provinces <- terra::vect(geodata::gadm("CAN", level = 1, path = CACHE_DIR))
provinces$unit_id   <- paste0("PROV_", provinces$NAME_1)
provinces$unit_type <- "province"
provinces$unit_name <- provinces$NAME_1

mgmt <- rbind(
  states[,    c("unit_id", "unit_type", "unit_name")],
  provinces[, c("unit_id", "unit_type", "unit_name")]
)
mgmt$area_km2 <- terra::expanse(mgmt, unit = "km")

log_line("mgmt_units", "-",
         sprintf("units: %d states + %d provinces = %d",
                 nrow(states), nrow(provinces), nrow(mgmt)))

# ── HELPERS ──────────────────────────────────────────────────────────────────
aggregate_suit <- function(suit_r, thr, units) {
  bn <- suit_r >= thr
  ex <- terra::extract(bn, units, fun = mean, na.rm = TRUE, ID = TRUE)
  frac <- ex[, 2]
  data.frame(
    unit_id           = units$unit_id,
    unit_type         = units$unit_type,
    unit_name         = units$unit_name,
    area_km2          = round(units$area_km2, 1),
    suitable_fraction = round(frac, 4),
    suitable_area_km2 = round(frac * units$area_km2, 1),
    stringsAsFactors  = FALSE
  )
}

# Doorstep: which units are within BUFFER_KM of any known invaded-range presence
compute_doorstep <- function(sid, units, buffer_km) {
  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (!file.exists(train_path)) return(character(0))
  train <- readRDS(train_path)$train
  if ("range_type" %in% names(train))
    pres <- train[train$pa == 1L & train$range_type == "invaded", c("x", "y")]
  else
    pres <- train[train$pa == 1L, c("x", "y")]
  if (nrow(pres) == 0) return(character(0))
  pv <- terra::vect(pres, geom = c("x", "y"), crs = "EPSG:4326")
  # Buffer in equal-area CRS
  units_sf <- sf::st_as_sf(units)
  pv_sf    <- sf::st_as_sf(pv)
  units_ea <- sf::st_transform(units_sf, "ESRI:54034")
  pv_ea    <- sf::st_transform(pv_sf,    "ESRI:54034")
  pv_buf   <- sf::st_buffer(pv_ea, dist = buffer_km * 1000)
  # Unit is "observed" if it contains any presence, "at_doorstep" if it
  # intersects any buffer but doesn't contain a presence
  observed_ix   <- unique(unlist(sf::st_intersects(units_ea, pv_ea)))
  observed_id   <- units$unit_id[observed_ix]
  buf_hit_ix    <- unique(unlist(sf::st_intersects(units_ea, pv_buf)))
  doorstep_id   <- setdiff(units$unit_id[buf_hit_ix], observed_id)
  list(observed = observed_id, doorstep = doorstep_id)
}

# ── MAIN LOOP ────────────────────────────────────────────────────────────────
time_points <- c(
  list(list(label = "current", ssp = NA, path_fn = function(sid)
    file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid)))),
  lapply(PERIOD_LABELS, function(lb) list(label = lb, ssp = "585",
    path_fn = function(sid) file.path(FUTURE_DIR,
      sprintf("%s_ensemble_%s_ssp585_median.tif", sid, lb))))
)

all_rows <- list()

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("mgmt_units", sid, sprintf("=== %s ===", sp$latin))
  thr <- get_threshold(sid)
  door <- compute_doorstep(sid, mgmt, DOORSTEP_BUFFER_KM)

  for (tp in time_points) {
    p <- tp$path_fn(sid)
    if (!file.exists(p)) next
    r <- terra::rast(p)
    if (max(terra::minmax(r)[2, ], na.rm = TRUE) > 1.5) r <- r / 1000
    df <- aggregate_suit(r, thr, mgmt)
    df$species    <- sid
    df$time_point <- tp$label
    df$ssp        <- ifelse(is.na(tp$ssp), "n/a", tp$ssp)
    df$status <- dplyr::case_when(
      df$unit_id %in% door$observed ~ "currently_observed",
      df$unit_id %in% door$doorstep ~ "at_doorstep",
      df$suitable_fraction >= UNIT_SUITABLE_THRESHOLD ~ "suitable_distant",
      TRUE ~ "unsuitable"
    )
    all_rows[[paste(sid, tp$label, tp$ssp, sep = "_")]] <- df
  }
}

if (length(all_rows)) {
  out <- dplyr::bind_rows(all_rows)
  write.csv(out, file.path(EVAL_DIR, "mgmt_units_suitability.csv"),
            row.names = FALSE)

  # Doorstep-only watchlist for 2050 SSP585
  door_df <- out[out$time_point == "2050" & out$status == "at_doorstep", ]
  door_df <- door_df[order(-door_df$suitable_fraction), ]
  write.csv(door_df, file.path(EVAL_DIR, "mgmt_units_doorstep.csv"),
            row.names = FALSE)

  # Summary counts per unit x time_point
  summ <- out %>%
    dplyr::group_by(unit_id, unit_type, unit_name, time_point) %>%
    dplyr::summarise(
      n_currently_observed = sum(status == "currently_observed"),
      n_at_doorstep        = sum(status == "at_doorstep"),
      n_suitable_distant   = sum(status == "suitable_distant"),
      .groups = "drop"
    )
  write.csv(summ, file.path(EVAL_DIR, "mgmt_units_summary.csv"),
            row.names = FALSE)

  log_line("mgmt_units", "-",
           sprintf("wrote %d rows to mgmt_units_suitability.csv", nrow(out)))
}

log_line("mgmt_units", "-", "done")

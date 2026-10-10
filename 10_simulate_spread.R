# =============================================================================
# 10_simulate_spread.R
# Cellular automaton (CA) dispersal simulation on top of ensemble suitability.
#
# For each species x SSP, step year-by-year from START_YEAR to END_YEAR:
#   1. Cells within an exponential dispersal kernel of currently-occupied cells
#      colonize with probability 1 - exp(-pressure * amp), gated by whether
#      the target cell is climatically suitable at that year.
#   2. Rare long-distance dispersal (LDD) events fire at LDD_RATE_PER_YEAR
#      per occupied cell, jumping LDD_MIN..LDD_MAX km.
#   3. Stochastic: N_REPLICATES replicates; probabilities averaged across.
#
# Between decadal snapshots, suitability is linearly interpolated between
# available time slices (current, 2030, 2050, 2070, 2090).
#
# NOTE: This is expensive (hours per species) and RUN_SPREAD is FALSE by
# default. Turn on only for a paper / final analysis, not routine iteration.
#
# Output per species x SSP:
#   FUTURE_DIR/<sid>_spread_<year>_ssp<ssp>.tif   colonization P per decade
#   FUTURE_DIR/<sid>_priority_ssp<ssp>.tif        4=now, 3=2030, 2=2050, 1=2090, 0=never
#   EVAL_DIR/spread_trajectories.csv              annual area by species/ssp
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
})

log_line("simulate_spread", "-", "starting")

# ── PARAMETERS ───────────────────────────────────────────────────────────────
TIMESTEP_YEARS    <- 5
START_YEAR        <- 1985
END_YEAR          <- 2090
N_REPLICATES      <- 3
REPORT_YEARS      <- c(2030, 2050, 2070, 2090)

LDD_RATE_PER_YEAR <- 0.02
LDD_MIN_KM        <- 50
LDD_MAX_KM        <- 500

# Per-species mean dispersal rates (km/yr) — literature values
DISPERSAL_RATES <- list(
  autumn_olive = 4, bittersweet = 4, burning_bush = 3,
  common_buckthorn = 4, garlic_mustard = 6,
  japanese_barberry = 3, knotweed = 5, multiflora_rose = 5,
  norway_maple = 3, stiltgrass = 5, tree_of_heaven = 10
)
DEFAULT_DISPERSAL <- 5

# ── THRESHOLDS ───────────────────────────────────────────────────────────────
thr_file <- file.path(EVAL_DIR, "thresholds_summary.csv")
thr_df   <- if (file.exists(thr_file)) read.csv(thr_file) else NULL

get_threshold <- function(sid) {
  if (!is.null(thr_df) && sid %in% thr_df$species)
    thr_df[thr_df$species == sid, "maxTSS_thr"][1]
  else THRESH
}
get_dispersal <- function(sid) {
  if (sid %in% names(DISPERSAL_RATES)) DISPERSAL_RATES[[sid]] else DEFAULT_DISPERSAL
}

# ── KERNEL ───────────────────────────────────────────────────────────────────
build_kernel <- function(mean_km_per_step, cell_km) {
  radius_cells <- max(1, ceiling(3 * mean_km_per_step / cell_km))
  size <- 2 * radius_cells + 1
  ker  <- matrix(0, size, size)
  for (dy in seq(-radius_cells, radius_cells))
    for (dx in seq(-radius_cells, radius_cells)) {
      d <- sqrt(dx^2 + dy^2) * cell_km
      ker[dy + radius_cells + 1, dx + radius_cells + 1] <-
        if (d == 0) 0 else exp(-d / mean_km_per_step)
    }
  ker / sum(ker)
}

# ── SUITABILITY TIMELINE HELPER ──────────────────────────────────────────────
load_suitability_timeline <- function(sid, ssp) {
  cur <- file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid))
  if (!file.exists(cur)) return(NULL)
  timeline <- list("1985" = terra::rast(cur))
  for (label in PERIOD_LABELS) {
    f <- file.path(FUTURE_DIR,
                   sprintf("%s_ensemble_%s_ssp%s_median.tif", sid, label, ssp))
    if (file.exists(f)) timeline[[label]] <- terra::rast(f)
  }
  timeline
}

# Linear interp between adjacent decadal snapshots for a given year
suit_at_year <- function(timeline, year) {
  yrs <- as.numeric(names(timeline))
  if (year <= min(yrs)) return(timeline[[as.character(min(yrs))]])
  if (year >= max(yrs)) return(timeline[[as.character(max(yrs))]])
  lo <- max(yrs[yrs <= year]); hi <- min(yrs[yrs > year])
  w  <- (year - lo) / (hi - lo)
  (1 - w) * timeline[[as.character(lo)]] + w * timeline[[as.character(hi)]]
}

# ── MAIN LOOP ────────────────────────────────────────────────────────────────
traj_rows <- list()

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("simulate_spread", sid, sprintf("=== %s ===", sp$latin))

  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (!file.exists(train_path)) next
  train <- readRDS(train_path)$train
  seeds <- train[train$pa == 1L, c("x", "y")]
  if (nrow(seeds) < 5) next

  thr    <- get_threshold(sid)
  disp   <- get_dispersal(sid)

  for (ssp in SSPS) {
    log_line("simulate_spread", sid, sprintf("SSP %s (thr=%.3f, disp=%d km/yr)",
                                             ssp, thr, disp))
    timeline <- load_suitability_timeline(sid, ssp)
    if (is.null(timeline)) next
    template <- timeline[[1]]
    cell_km  <- mean(terra::res(template)) * 111

    ker <- build_kernel(disp * TIMESTEP_YEARS, cell_km)

    # Seed presence grid
    seed_v <- terra::vect(seeds, geom = c("x", "y"), crs = "EPSG:4326")
    occ0   <- terra::rasterize(seed_v, template, field = 1, background = 0)

    # Accumulate colonization probability
    p_col_sum <- rep(0, terra::ncell(template))
    yr_first  <- rep(NA_integer_, terra::ncell(template))
    yr_first[terra::values(occ0)[, 1] == 1] <- START_YEAR

    for (rep_i in seq_len(N_REPLICATES)) {
      occ <- occ0
      for (yr in seq(START_YEAR + TIMESTEP_YEARS, END_YEAR, by = TIMESTEP_YEARS)) {
        suit_r <- suit_at_year(timeline, yr)
        # Neighborhood pressure
        pressure <- terra::focal(occ, w = ker, fun = "sum", na.rm = TRUE,
                                 fillvalue = 0)
        p_col   <- 1 - exp(-3 * pressure)
        # Gate by suitability > threshold
        suit_ok <- terra::values(suit_r)[, 1] >= thr
        p_vals  <- terra::values(p_col)[, 1] * suit_ok
        # Stochastic colonization
        occ_v <- terra::values(occ)[, 1]
        new_col <- which(occ_v == 0 & runif(length(p_vals)) < p_vals)
        occ_v[new_col] <- 1
        # LDD
        n_occ <- sum(occ_v == 1, na.rm = TRUE)
        n_ldd <- rpois(1, n_occ * LDD_RATE_PER_YEAR * TIMESTEP_YEARS)
        if (n_ldd > 0 && n_occ > 0) {
          from <- sample(which(occ_v == 1), min(n_ldd, n_occ), replace = TRUE)
          xy_f <- terra::xyFromCell(template, from)
          for (k in seq_len(nrow(xy_f))) {
            d_km <- runif(1, LDD_MIN_KM, LDD_MAX_KM)
            ang  <- runif(1, 0, 2 * pi)
            dx   <- d_km * cos(ang) / (111 * cos(xy_f[k, 2] * pi / 180))
            dy   <- d_km * sin(ang) / 111
            newc <- terra::cellFromXY(template,
                                      matrix(c(xy_f[k, 1] + dx,
                                               xy_f[k, 2] + dy), ncol = 2))
            if (!is.na(newc) && length(suit_ok) >= newc && suit_ok[newc])
              occ_v[newc] <- 1
          }
        }
        terra::values(occ) <- occ_v
        # First-arrival tracking (only on rep 1 for a deterministic priority map)
        if (rep_i == 1) {
          just_arrived <- which(is.na(yr_first) & occ_v == 1)
          yr_first[just_arrived] <- yr
        }
      }
      p_col_sum <- p_col_sum + terra::values(occ)[, 1]
    }
    p_col_avg <- p_col_sum / N_REPLICATES

    # Snapshot rasters per REPORT_YEAR (approximated as p_col_avg for now)
    out_r <- template
    terra::values(out_r) <- p_col_avg
    names(out_r) <- sprintf("%s_spread_end_ssp%s", sid, ssp)
    terra::writeRaster(out_r,
                       file.path(FUTURE_DIR,
                                 sprintf("%s_spread_end_ssp%s.tif", sid, ssp)),
                       overwrite = TRUE, gdal = GDAL_OPTS)

    # Priority raster: 4=now, 3=2030, 2=2050, 1=2090, 0=never
    pri <- integer(length(yr_first))
    pri[!is.na(yr_first) & yr_first == START_YEAR]           <- 4
    pri[!is.na(yr_first) & yr_first > START_YEAR & yr_first <= 2030] <- 3
    pri[!is.na(yr_first) & yr_first > 2030 & yr_first <= 2050]      <- 2
    pri[!is.na(yr_first) & yr_first > 2050]                          <- 1
    pri_r <- template
    terra::values(pri_r) <- pri
    terra::writeRaster(pri_r,
                       file.path(FUTURE_DIR,
                                 sprintf("%s_priority_ssp%s.tif", sid, ssp)),
                       overwrite = TRUE, gdal = GDAL_OPTS)

    traj_rows[[paste0(sid, "_", ssp)]] <- data.frame(
      species = sid, ssp = ssp,
      colonized_km2 = round(sum(p_col_avg > 0.5, na.rm = TRUE) * cell_km^2, 1)
    )
    log_line("simulate_spread", sid,
             sprintf("SSP %s: %.0f km2 colonized (>0.5 P)",
                     ssp, tail(traj_rows, 1)[[1]]$colonized_km2))
  }
}

if (length(traj_rows)) {
  out <- do.call(rbind, traj_rows)
  write.csv(out, file.path(EVAL_DIR, "spread_trajectories.csv"),
            row.names = FALSE)
}

log_line("simulate_spread", "-", "done")

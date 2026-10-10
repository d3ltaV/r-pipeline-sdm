# =============================================================================
# 02_download_gbif.R
# Download GBIF occurrence records for the 11 species — both native and
# invaded ranges — and cache to disk.
#
# YOU DON'T NEED TO RUN THIS if your v2a caches already exist. Script 03
# (build_training) checks GBIF_DIR first, then falls back to GBIF_DIR_LEGACY
# (gbif_v2a/). No new downloads happen if either cache is populated.
#
# Runtime on a fresh machine: ~1-3 hours (GBIF rate-limited to ~300 records
# per request; large-range species like Norway maple can pull >50k records).
#
# Requires a GBIF account. Credentials via env vars OR ~/.Renviron:
#   GBIF_USER=<your gbif username>
#   GBIF_PWD=<your gbif password>
#   GBIF_EMAIL=<your email>
#
# Outputs (per species):
#   GBIF_DIR/<sid>_gbif_native.rds     raw native-range records
#   GBIF_DIR/<sid>_gbif_invaded.rds    raw invaded-range records
#
# The training tables (<sid>_train.rds) are built later by 03_build_training.R
# — this script only handles the download/cache step so it can be resumed.
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(rgbif)
  library(dplyr)
})

log_line("download_gbif", "-", "starting")

# ── ISO3 → ISO2 country-code map for rgbif::occ_search() ─────────────────────
ISO3_TO_ISO2 <- c(
  JPN="JP", KOR="KR", CHN="CN", TWN="TW", VNM="VN", IND="IN", NPL="NP",
  USA="US", CAN="CA", MEX="MX",
  GBR="GB", FRA="FR", DEU="DE", POL="PL", ITA="IT", ESP="ES",
  AUT="AT", CZE="CZ", BEL="BE", NLD="NL", SWE="SE", NOR="NO",
  FIN="FI", DNK="DK", CHE="CH", HUN="HU", ROU="RO", BGR="BG",
  RUS="RU", UKR="UA"
)

# ── PER-COUNTRY FETCH (paginated) ────────────────────────────────────────────
# rgbif::occ_search caps at 300 records per call. Loop until exhausted or
# `cap` reached. Slow but reliable; you can also switch to occ_download for
# large pulls (async, but requires post-processing of DwC-A zips).
fetch_one_country <- function(latin, iso2, chunk = 300, cap = 30000) {
  recs <- list(); start <- 0
  repeat {
    res <- tryCatch(
      rgbif::occ_search(scientificName = latin, country = iso2,
                        hasCoordinate = TRUE, limit = chunk, start = start),
      error = function(e) { Sys.sleep(5); NULL }
    )
    if (is.null(res) || is.null(res$data) || nrow(res$data) == 0) break
    keep <- intersect(
      c("decimalLongitude", "decimalLatitude", "year", "countryCode",
        "basisOfRecord", "coordinateUncertaintyInMeters",
        "institutionCode", "datasetKey"),
      names(res$data)
    )
    recs[[length(recs) + 1]] <- res$data[, keep, drop = FALSE]
    start <- start + chunk
    if (nrow(res$data) < chunk || start >= cap) break
  }
  if (!length(recs)) return(NULL)
  out <- dplyr::bind_rows(recs)
  if ("decimalLongitude" %in% names(out))
    out <- dplyr::rename(out, longitude = decimalLongitude,
                              latitude  = decimalLatitude)
  out
}

# ── PER-SPECIES FETCH (native + invaded, both cached) ────────────────────────
fetch_range <- function(sp, iso_codes, cache_file, label) {
  if (file.exists(cache_file)) {
    d <- readRDS(cache_file)
    if (is.data.frame(d) && nrow(d) > 0) {
      log_line("download_gbif", sp$short_id,
               sprintf("%s cached (%d rows)", label, nrow(d)))
      return(invisible(d))
    }
  }
  all_recs <- list()
  for (iso3 in iso_codes) {
    iso2 <- ISO3_TO_ISO2[[iso3]]
    if (is.null(iso2) || is.na(iso2)) {
      log_line("download_gbif", sp$short_id,
               sprintf("no ISO2 mapping for %s (skipping)", iso3))
      next
    }
    d <- fetch_one_country(sp$latin, iso2)
    n <- if (is.null(d)) 0 else nrow(d)
    log_line("download_gbif", sp$short_id,
             sprintf("%s %s (%s): %d records", label, iso3, iso2, n))
    if (!is.null(d) && nrow(d) > 0) all_recs[[iso3]] <- d
    Sys.sleep(1)   # be polite to GBIF
  }
  if (!length(all_recs)) {
    log_line("download_gbif", sp$short_id,
             sprintf("%s: NO RECORDS ANYWHERE", label))
    return(invisible(NULL))
  }
  out <- dplyr::bind_rows(all_recs, .id = "iso3")
  saveRDS(out, cache_file)
  log_line("download_gbif", sp$short_id,
           sprintf("%s: cached %d rows to %s", label, nrow(out),
                   basename(cache_file)))
  invisible(out)
}

# ── MAIN LOOP ────────────────────────────────────────────────────────────────
for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("download_gbif", sid, sprintf("=== %s ===", sp$latin))

  # native — check v3 dir, then v2a legacy, then download
  native_v3     <- file.path(GBIF_DIR,        paste0(sid, "_gbif_native.rds"))
  native_legacy <- file.path(GBIF_DIR_LEGACY, paste0(sid, "_gbif_native.rds"))
  native_target <- if (file.exists(native_v3)) native_v3
                   else if (file.exists(native_legacy)) native_legacy
                   else native_v3
  fetch_range(sp, sp$native_iso, native_target, "native")

  # invaded — same logic
  invaded_v3     <- file.path(GBIF_DIR,        paste0(sid, "_gbif_invaded.rds"))
  invaded_legacy <- file.path(GBIF_DIR_LEGACY, paste0(sid, "_gbif_invaded.rds"))
  invaded_target <- if (file.exists(invaded_v3)) invaded_v3
                    else if (file.exists(invaded_legacy)) invaded_legacy
                    else invaded_v3
  fetch_range(sp, sp$invaded_iso, invaded_target, "invaded")
}

# ── OPTIONAL: TARGET-GROUP BACKGROUND (Tracheophyta) ─────────────────────────
# Only run if USE_TGB = TRUE in config. Pulls up to TGB_LIMIT records per
# country for all vascular plants. Used as pseudo-absence pool that reflects
# recorder effort rather than random geography (Phillips et al. 2009).
if (USE_TGB) {
  tgb_cache <- file.path(GBIF_DIR, "tracheophyta_tgb.rds")
  if (file.exists(tgb_cache)) {
    log_line("download_gbif", "TGB", "cached; skipping")
  } else {
    all_isos <- unique(unlist(lapply(SPECIES,
                                     function(s) c(s$native_iso, s$invaded_iso))))
    tgb_country_dir <- file.path(GBIF_DIR, "tgb_countries")
    dir.create(tgb_country_dir, recursive = TRUE, showWarnings = FALSE)
    all_recs <- list()
    for (iso3 in all_isos) {
      iso2 <- ISO3_TO_ISO2[[iso3]]
      if (is.null(iso2) || is.na(iso2)) next
      per_country <- file.path(tgb_country_dir, paste0("tgb_", iso3, ".rds"))
      if (file.exists(per_country)) {
        all_recs[[iso3]] <- readRDS(per_country)
        log_line("download_gbif", "TGB",
                 sprintf("%s cached (%d)", iso3, nrow(all_recs[[iso3]])))
        next
      }
      log_line("download_gbif", "TGB", sprintf("fetching %s", iso3))
      d <- tryCatch(
        rgbif::occ_search(taxonKey = 7707728,   # Tracheophyta
                          country = iso2, hasCoordinate = TRUE,
                          limit = TGB_LIMIT)$data,
        error = function(e) NULL
      )
      if (!is.null(d) && nrow(d) > 0) {
        d <- dplyr::rename(d[, c("decimalLongitude", "decimalLatitude")],
                           longitude = decimalLongitude,
                           latitude  = decimalLatitude)
        saveRDS(d, per_country)
        all_recs[[iso3]] <- d
      } else {
        saveRDS(data.frame(longitude = numeric(), latitude = numeric()),
                per_country)
      }
      Sys.sleep(1)
    }
    tgb <- dplyr::bind_rows(all_recs) |>
      dplyr::filter(!is.na(longitude), !is.na(latitude)) |>
      dplyr::distinct(longitude, latitude)
    saveRDS(tgb, tgb_cache)
    log_line("download_gbif", "TGB",
             sprintf("wrote %d unique TGB locations", nrow(tgb)))
  }
}

log_line("download_gbif", "-", "done")
cat("\nGBIF download complete. Caches in ", GBIF_DIR, "\n", sep = "")

# =============================================================================
# 12_copy_to_webapp.R
# Stages ensemble rasters into webapp_output_hab/<sid>/ and (optionally) syncs
# that folder to WEBAPP_DATA/<sid>/ for live deployment.
#
# What lands in each species folder:
#   current.tif                              (from 05, renamed)
#   <label>_ssp<ssp>.tif                     (from 05, renamed)
#   <label>_ssp<ssp>_agreement.tif           (from 05, GCM agreement)
#   <label>_ssp<ssp>_cv.tif                  (from 05, GCM uncertainty)
#   presence.tif                             (from 11)
#   occurrences.geojson                      (from 11)
#   occurrence_meta.json                     (from 11)
#   counties.json (or *_counties.json)       (from aggregate_counties, if present)
#
# Toggles at the top control DRY_RUN, cleaning behavior, and whether Stage 2
# (sync to WEBAPP_DATA) actually runs. Set DRY_RUN <- TRUE to preview.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

# ── KNOBS ────────────────────────────────────────────────────────────────────
MODEL_TAG      <- "ensemble"
DRY_RUN        <- FALSE   # TRUE = print plan, don't touch files
CLEAN_FIRST    <- TRUE    # wipe stale *.tif in stage before restaging
COPY_META      <- FALSE   # also copy agreement + CV rasters (webapp doesn't read them; flip on if you add an uncertainty overlay)
SYNC_TO_WEBAPP <- TRUE    # push staged mirror → WEBAPP_DATA

log_line("copy_to_webapp", "-",
         sprintf("MODEL_TAG=%s  DRY_RUN=%s  CLEAN=%s  SYNC=%s",
                 MODEL_TAG, DRY_RUN, CLEAN_FIRST, SYNC_TO_WEBAPP))

# ── DIRS ─────────────────────────────────────────────────────────────────────
src_current_dir <- normalizePath(CURRENT_DIR,       mustWork = FALSE)
src_future_dir  <- normalizePath(FUTURE_DIR,        mustWork = FALSE)
stage_dir       <- normalizePath(WEBAPP_OUTPUT_DIR, mustWork = FALSE)
webapp_data     <- normalizePath(WEBAPP_DATA,       mustWork = FALSE)

if (!dir.exists(src_current_dir)) stop("Missing: ", src_current_dir)
if (!dir.exists(src_future_dir))  stop("Missing: ", src_future_dir)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

cat("SRC current: ", src_current_dir, "\n")
cat("SRC future:  ", src_future_dir,  "\n")
cat("STAGE:       ", stage_dir,       "\n")
cat("WEBAPP:      ", webapp_data,     "\n\n")

# ── HELPERS ──────────────────────────────────────────────────────────────────
do_copy <- function(src, dst) {
  if (!file.exists(src)) {
    cat("  x MISSING SRC:", basename(src), "\n"); return(FALSE)
  }
  if (DRY_RUN) {
    cat("  [dry]", basename(src), "->", basename(dst), "\n"); return(TRUE)
  }
  ok <- file.copy(src, dst, overwrite = TRUE)
  cat(if (ok) "  v " else "  x ", basename(src), " -> ", basename(dst), "\n", sep = "")
  ok
}

# ── STAGE 1: PER-SPECIES COPY ────────────────────────────────────────────────
n_ok <- 0
for (sp in SPECIES) {
  sid   <- sp$short_id
  cat("\n== ", sid, " ==\n", sep = "")
  spdir <- file.path(stage_dir, sid)
  dir.create(spdir, recursive = TRUE, showWarnings = FALSE)

  if (CLEAN_FIRST && !DRY_RUN) {
    # Keep presence.tif from 11; only clean this stage's rasters.
    stale <- setdiff(
      list.files(spdir, pattern = "\\.tif$", full.names = TRUE),
      file.path(spdir, "presence.tif")
    )
    if (length(stale)) file.remove(stale)
  }

  # Current ensemble
  src_cur <- file.path(src_current_dir, sprintf("%s_%s_current.tif", sid, MODEL_TAG))
  do_copy(src_cur, file.path(spdir, "current.tif"))

  # Futures (median + agreement + cv per label x ssp)
  for (ssp in SSPS) for (label in PERIOD_LABELS) {
    src_med   <- file.path(src_future_dir,
                           sprintf("%s_%s_%s_ssp%s_median.tif",    sid, MODEL_TAG, label, ssp))
    dst_med   <- file.path(spdir, sprintf("%s_ssp%s.tif", label, ssp))
    do_copy(src_med, dst_med)

    if (COPY_META) {
      src_agr <- file.path(src_future_dir,
                           sprintf("%s_%s_%s_ssp%s_agreement.tif", sid, MODEL_TAG, label, ssp))
      dst_agr <- file.path(spdir, sprintf("%s_ssp%s_agreement.tif", label, ssp))
      do_copy(src_agr, dst_agr)
      src_cv  <- file.path(src_future_dir,
                           sprintf("%s_%s_%s_ssp%s_cv.tif",        sid, MODEL_TAG, label, ssp))
      dst_cv  <- file.path(spdir, sprintf("%s_ssp%s_cv.tif", label, ssp))
      do_copy(src_cv, dst_cv)
    }
  }
  n_ok <- n_ok + 1
}

# ── STAGE 2: SYNC STAGE → WEBAPP_DATA ────────────────────────────────────────
if (SYNC_TO_WEBAPP) {
  cat("\n\n== SYNC: ", stage_dir, " -> ", webapp_data, " ==\n", sep = "")
  if (!dir.exists(webapp_data)) {
    if (DRY_RUN) cat("  [dry] would create ", webapp_data, "\n", sep = "")
    else dir.create(webapp_data, recursive = TRUE, showWarnings = FALSE)
  }
  for (sp in SPECIES) {
    sid   <- sp$short_id
    src_d <- file.path(stage_dir, sid)
    dst_d <- file.path(webapp_data, sid)
    if (!dir.exists(src_d)) next
    if (DRY_RUN) {
      cat("  [dry] would sync ", sid, "/\n", sep = "")
      next
    }
    dir.create(dst_d, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(src_d, full.names = TRUE, recursive = FALSE)
    for (f in files) file.copy(f, file.path(dst_d, basename(f)), overwrite = TRUE)
    cat("  synced ", sid, " (", length(files), " files)\n", sep = "")
  }

  # Shared root files (counties.topo.json, counties_grid.bin, etc.)
  for (fn in c("counties.topo.json", "counties.geojson", "counties_grid.bin",
               "counties_grid.json")) {
    src <- file.path(getwd(), fn)
    if (file.exists(src)) {
      dst <- file.path(webapp_data, fn)
      if (DRY_RUN) cat("  [dry] shared: ", fn, "\n", sep = "")
      else { file.copy(src, dst, overwrite = TRUE); cat("  shared: ", fn, "\n", sep = "") }
    }
  }
} else {
  cat("\nSYNC skipped (SYNC_TO_WEBAPP = FALSE)\n")
}

log_line("copy_to_webapp", "-",
         sprintf("done (%d species staged, sync=%s, dry=%s)",
                 n_ok, SYNC_TO_WEBAPP, DRY_RUN))

cat("\nDONE.\n")
if (DRY_RUN) cat("This was a DRY RUN. Set DRY_RUN <- FALSE and re-source to actually copy.\n")

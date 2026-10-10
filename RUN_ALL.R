# =============================================================================
# RUN_ALL.R  —  habitat-augmented pipeline (final_hab_v1)
# End-to-end driver. Runs from scratch or resumes an interrupted run.
#
# NEW USER
# --------
# 1. Read QUICKSTART.md first — it walks through R install, package install,
#    directory layout, credentials, and where to put existing data.
# 2. Set your working directory to the project root, then source this file:
#      setwd("/path/to/invasive-sdm-pipeline")
#      source("RUN_ALL.R")
# 3. First run downloads all raw data (climate + GBIF + habitat) into data/,
#    then fits models and generates outputs. Everything caches — you can Ctrl-C
#    and re-run source("RUN_ALL.R") to resume from where it stopped.
#
# EXISTING USER (already have climate + GBIF from an earlier v2a/final run)
# ------------------------------------------------------------------------
# Same command. The download scripts (01, 02) look on disk first and skip if
# data is already there at the right extent. Only 01b (habitat) will actually
# download anything the first time. Toggle any RUN_* flag to FALSE below if
# you want to skip a step explicitly.
#
# What's different from your climate-only RUN_ALL.R:
#   1. New step 01b downloads elev/hfp/treecov/soilph and aligns them to the
#      WorldClim grid. First time: 15-45 min. After that: cached.
#   2. Output paths carry a _hab suffix (models_hab/, output/current_hab/…)
#      so this rebuild leaves your existing v2a/final results alone.
#   3. Training tables land at data/gbif/<sid>_train_hab.rds — coexists with
#      the climate-only <sid>_train.rds files.
#   4. Scripts 06/07/08/09/10/11/12/13 and aggregate_counties.R DID NOT need
#      structural changes — they read via CURRENT_DIR / FUTURE_DIR / EVAL_DIR
#      (now suffixed) and iterate PREDICTORS (now 8 items).
#
# Runtime estimate on a laptop (fresh, everything from scratch):
#   01  climate download:     1-2 hours (dominated by 40 CMIP6 fetches)
#   01b habitat download:     15-45 min ONCE
#   02  GBIF download:        1-3 hours (rate-limited)
#   03  training tables:      ~6 min for 12 species
#   04  biomod2 fits:         2-8 hours (biggest single cost)
#   05  projections:          3-8 hours (12 species × 41 rasters)
#   06-13:                    ~30 min combined
# Total cold-start:           8-18 hours. Every step is cached and resumable.
# =============================================================================

# ── STEP FLAGS ───────────────────────────────────────────────────────────────
# All ON by default. Each script checks its own cache before doing work, so
# leaving downloads ON is safe even when data is already present. Flip to
# FALSE only if you want to explicitly skip a step for debugging.

RUN_DL_CLIMATE   <- TRUE    # 01   WorldClim current + 40 CMIP6 futures
RUN_DL_HABITAT   <- TRUE    # 01b  elev, HFP, tree cover, soil pH
RUN_DL_GBIF      <- TRUE    # 02   GBIF occurrence records (needs credentials)

RUN_BUILD_TRAIN  <- TRUE    # 03   thin, buffer, background, extract
RUN_FIT          <- TRUE    # 04   biomod2 ensemble on 8 predictors     ← slow
RUN_HONEST_EVAL  <- TRUE    # 07b  spatial-block OOF ensemble scores + OOF rds
RUN_PROJECT      <- TRUE    # 05   current + 40 futures on 8 predictors ← slow
RUN_THRESHOLDS   <- TRUE    # 06   per-species maxTSS (uses 07b OOF predictions)
RUN_EVALUATE     <- TRUE    # 07   in-sample AUC/TSS/Boyce (diagnostic only)
RUN_PARTIAL      <- FALSE   # 08   response curves (opt-in, slow)
RUN_MAPS         <- TRUE    # 09   QC 4-panel figures

RUN_SPREAD       <- FALSE   # 10   cellular automaton spread (opt-in, slow)

RUN_OCCURRENCES  <- TRUE    # 11   presence + geojson → webapp_output_hab/
RUN_COPY_WEBAPP  <- TRUE    # 12   stage + sync webapp (DRY_RUN off by default)

RUN_MGMT_UNITS   <- FALSE   # 13   state/province rollup (opt-in, slow)

# ── LOAD CONFIG ──────────────────────────────────────────────────────────────
source("00_config.R")

# ── STEP RUNNER ──────────────────────────────────────────────────────────────
# Each step runs in its own tryCatch so a failure in (e.g.) 08 doesn't kill
# the pipeline. Failures are logged; execution continues with the next step.
run_step <- function(flag, path, label) {
  if (!isTRUE(flag)) {
    log_line("run_all", "-", sprintf("SKIP  %s", label))
    return(invisible(NULL))
  }
  if (!file.exists(path)) {
    log_line("run_all", "-", sprintf("MISSING FILE  %s (skipping)", path))
    message(sprintf("!!! Missing script: %s — skipping %s", path, label))
    return(invisible(NULL))
  }
  log_line("run_all", "-", sprintf(">>> %s", label))
  t0 <- Sys.time()
  ok <- tryCatch({
    source(path, local = FALSE)
    TRUE
  }, error = function(e) {
    log_line("run_all", "-",
             sprintf("FAILED %s : %s", label, conditionMessage(e)))
    message("\n!!!  ", path, " failed: ", conditionMessage(e))
    message("     Continuing with next step.\n")
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_line("run_all", "-",
           sprintf("<<< %s  (%.1f min, ok=%s)", label, dt, ok))
  invisible(ok)
}

# ── EXECUTE IN ORDER ─────────────────────────────────────────────────────────
overall_t0 <- Sys.time()
log_line("run_all", "-", "===== PIPELINE START (final_hab_v1) =====")

run_step(RUN_DL_CLIMATE,  "01_download_climate.R",   "01 download climate")
run_step(RUN_DL_HABITAT,  "01b_download_habitat.R",  "01b download habitat")
run_step(RUN_DL_GBIF,     "02_download_gbif.R",      "02 download GBIF")
run_step(RUN_BUILD_TRAIN, "03_build_training.R",     "03 build training tables (8 preds)")
run_step(RUN_FIT,         "04_fit_ensemble.R",       "04 fit biomod2 ensemble")
run_step(RUN_HONEST_EVAL, "07b_honest_scores.R",     "07b honest OOF ensemble scores")
run_step(RUN_PROJECT,     "05_project_ensemble.R",   "05 project onto climates")
run_step(RUN_THRESHOLDS,  "06_thresholds.R",         "06 species thresholds (OOF-based)")
run_step(RUN_EVALUATE,    "07_evaluate.R",           "07 in-sample eval (diagnostic)")
run_step(RUN_PARTIAL,     "08_partial_responses.R",  "08 partial responses")
run_step(RUN_MAPS,        "09_map_panels.R",         "09 QC map panels")
run_step(RUN_SPREAD,      "10_simulate_spread.R",    "10 spread simulation")
run_step(RUN_OCCURRENCES, "11_export_occurrences.R", "11 export occurrences")
run_step(RUN_COPY_WEBAPP, "12_copy_to_webapp.R",     "12 copy to webapp")
run_step(RUN_MGMT_UNITS,  "13_management_units.R",   "13 management units")

# ── SUMMARY ──────────────────────────────────────────────────────────────────
dt_all <- as.numeric(difftime(Sys.time(), overall_t0, units = "mins"))
log_line("run_all", "-",
         sprintf("===== PIPELINE DONE (%.1f min) =====", dt_all))

cat("\n\n")
cat("========================================\n")
cat("  Pipeline finished in ", round(dt_all, 1), " minutes\n", sep = "")
cat("========================================\n\n")
cat("Check:\n")
cat("  eval_hab/ensemble_blockcv_heldout.csv — honest OOF AUC/TSS/Boyce (report these)\n")
cat("  eval_hab/eval_insample.csv         — in-sample (diagnostic only, ~0.99)\n")
cat("  eval_hab/biomod_var_importance.csv — which predictors mattered\n")
cat("  diag_hab/<sid>_vif.txt             — collinearity per species\n")
cat("  figures/maps_hab/*.png             — QC visuals\n")
cat("  webapp_output_hab/<sid>/           — webapp mirror in project root\n")
cat("  logs/pipeline_hab.log              — full run log\n\n")
cat("If 12 was in DRY_RUN mode, open 12_copy_to_webapp.R, flip DRY_RUN <-\n")
cat("FALSE, and source it again to actually sync the webapp.\n\n")

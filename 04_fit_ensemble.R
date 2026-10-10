# =============================================================================
# 04_fit_ensemble.R  —  habitat-augmented (final_hab_v1)
#
# WHAT CHANGED vs. the climate-only version
# -------------------------------------------------------------------------
#   1. Reads <sid>_train_hab.rds instead of <sid>_train.rds (so both pipelines
#      can coexist).
#   2. Loads the env stack via load_env_current() (climate + habitat) instead
#      of just the climate raster. Fixes the names() bug you'd hit if you
#      pointed the old code at the new 8-predictor training tables.
#   3. Writes biomod paths into the suffixed MODELS_DIR ("models_hab/") and
#      eval CSVs into the suffixed EVAL_DIR ("eval_hab/").
#
# biomod2's fit doesn't actually need the raster for training — it uses the
# data.frame of expl.var + resp.var. The env raster is only referenced for
# response-curve machinery. So the change is purely cosmetic for the fit
# itself, but keeps everything consistent for 05/08 that DO need the raster.
#
# The rest of the script (BIOMOD_FormatingData, BIOMOD_Modeling,
# BIOMOD_EnsembleModeling, tidy-score aggregation) is unchanged.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
  library(dplyr)
})

# Base R < 4.4 doesn't have %||%; define our own so this file is self-contained
`%||%` <- function(a, b) if (is.null(a)) b else a

log_line("fit_ensemble", "-",
         sprintf("starting (biomod2 v%s), %d predictors",
                 packageVersion("biomod2"), length(PREDICTORS)))

orig_wd     <- getwd()
biomod_root <- normalizePath(MODELS_DIR, winslash = "/", mustWork = FALSE)
dir.create(biomod_root, recursive = TRUE, showWarnings = FALSE)

# Env raster (climate + habitat) — kept for consistency with 05, not strictly
# needed by BIOMOD_Modeling for the fit itself.
env_current <- load_env_current()

# ── PER-SPECIES FIT ──────────────────────────────────────────────────────────
score_rows      <- list()
importance_rows <- list()

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("fit_ensemble", sid, sprintf("=== %s ===", sp$latin))

  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  if (!file.exists(train_path)) {
    log_line("fit_ensemble", sid, "no training data; skipping"); next
  }
  paths_out <- normalizePath(
    file.path(MODELS_DIR, sid, paste0(sid, "_biomod_paths.rds")),
    winslash = "/", mustWork = FALSE)
  dir.create(dirname(paths_out), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(paths_out)) {
    log_line("fit_ensemble", sid, "already fit; skipping"); next
  }

  obj  <- readRDS(train_path)
  d    <- obj$train

  # Sanity: does the training table's predictor set match our config?
  if (!identical(sort(obj$predictors %||% PREDICTORS), sort(PREDICTORS))) {
    log_line("fit_ensemble", sid,
             sprintf("WARN: training predictors (%s) don't match config (%s); using config",
                     paste(obj$predictors, collapse = ","),
                     paste(PREDICTORS,     collapse = ",")))
  }

  n_p  <- sum(d$pa == 1L); n_b <- sum(d$pa == 0L)
  log_line("fit_ensemble", sid, sprintf("n_pres=%d  n_bg=%d", n_p, n_b))

  setwd(biomod_root)
  on.exit(setwd(orig_wd), add = TRUE, after = FALSE)

  # 1. BIOMOD_FormatingData -------------------------------------------------
  bm_data <- tryCatch(
    biomod2::BIOMOD_FormatingData(
      resp.var  = d$pa,
      expl.var  = d[, PREDICTORS, drop = FALSE],
      resp.xy   = as.matrix(d[, c("x", "y")]),
      resp.name = sid,
      filter.raster = FALSE
    ),
    error = function(e) {
      log_line("fit_ensemble", sid,
               sprintf("FormatingData failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(bm_data)) { setwd(orig_wd); next }

  # 2. Modeling options -----------------------------------------------------
  bm_opts <- biomod2::bm_ModelingOptions(
    data.type = "binary",
    models    = BIOMOD_MODELS,
    strategy  = "default"
  )

  # 3. BIOMOD_Modeling — the big fit call ------------------------------------
  bm_out <- tryCatch(
    biomod2::BIOMOD_Modeling(
      bm.format         = bm_data,
      modeling.id       = paste0("hab_", format(Sys.Date(), "%Y%m%d")),
      models            = BIOMOD_MODELS,
      OPT.strategy      = "default",
      # biomod2's "block" strategy ignores CV.k and CV.perc — it always
      # produces 4 ENMeval lon/lat quadrants (RUN1..RUN4). Left here only
      # because biomod2 requires CV.nb.rep.
      CV.strategy       = CV_STRATEGY,
      CV.nb.rep         = 1,
      metric.eval       = METRIC_EVAL,
      var.import        = 3,
      scale.models      = TRUE,
      nb.cpu            = 1,
      seed.val          = 42
    ),
    error = function(e) {
      log_line("fit_ensemble", sid,
               sprintf("BIOMOD_Modeling failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(bm_out)) { setwd(orig_wd); next }

  # 4. BIOMOD_EnsembleModeling ---------------------------------------------
  em_out <- tryCatch(
    biomod2::BIOMOD_EnsembleModeling(
      bm.mod                = bm_out,
      em.by                 = "all",
      em.algo               = BIOMOD_EM_ALGO,
      metric.select         = METRIC_SELECT,
      metric.select.thresh  = METRIC_THRESH,
      metric.eval           = METRIC_EVAL,
      var.import            = 3,
      EMci.alpha            = 0.05,
      EMwmean.decay         = "proportional",
      seed.val              = 42
    ),
    error = function(e) {
      log_line("fit_ensemble", sid,
               sprintf("EnsembleModeling failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(em_out)) { setwd(orig_wd); next }

  # 5. Persist paths so downstream scripts can reload -----------------------
  saveRDS(list(
    bm_out_link   = bm_out@link,
    em_out_link   = em_out@link,
    modeling_id   = bm_out@modeling.id,
    resp_name     = sid,
    workspace     = biomod_root,
    predictors    = PREDICTORS,
    pipeline_tag  = PIPELINE_TAG
  ), paths_out)

  # 6. Tidy evaluation scores + variable importance ------------------------
  scores <- biomod2::get_evaluations(bm_out)
  if (!is.null(scores) && nrow(scores) > 0) {
    scores$species <- sid
    score_rows[[sid]] <- scores
  }
  em_scores <- tryCatch(biomod2::get_evaluations(em_out),
                        error = function(e) NULL)
  if (!is.null(em_scores) && nrow(em_scores) > 0) {
    em_scores$species <- sid
    score_rows[[paste0(sid, "_ensemble")]] <- em_scores
  }

  vip <- tryCatch(biomod2::get_variables_importance(bm_out),
                  error = function(e) NULL)
  if (!is.null(vip) && nrow(vip) > 0) {
    vip$species <- sid
    importance_rows[[sid]] <- vip
  }

  log_line("fit_ensemble", sid,
           sprintf("wrote biomod paths → %s", paths_out))
  setwd(orig_wd)
}

setwd(orig_wd)

# ── AGGREGATE TIDY TABLES ────────────────────────────────────────────────────
if (length(score_rows)) {
  scores_df <- dplyr::bind_rows(score_rows)
  write.csv(scores_df, file.path(EVAL_DIR, "biomod_scores.csv"),
            row.names = FALSE)
  log_line("fit_ensemble", "-",
           sprintf("wrote %d rows to biomod_scores.csv", nrow(scores_df)))
}
if (length(importance_rows)) {
  imp_df <- dplyr::bind_rows(importance_rows)
  write.csv(imp_df, file.path(EVAL_DIR, "biomod_var_importance.csv"),
            row.names = FALSE)
  log_line("fit_ensemble", "-",
           sprintf("wrote %d rows to biomod_var_importance.csv", nrow(imp_df)))
  # Also print a quick top-3 per species so you see impact of new predictors
  # right away without opening the CSV:
  cat("\nTop 3 predictors by mean var.importance per species:\n")
  suppressMessages({
    top <- imp_df %>%
      dplyr::group_by(species, expl.var) %>%
      dplyr::summarise(imp = mean(var.imp, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(species) %>%
      dplyr::slice_max(imp, n = 3, with_ties = FALSE) %>%
      dplyr::mutate(imp = round(imp, 3))
  })
  print(top, n = 100)
}

log_line("fit_ensemble", "-", "done")
cat("\nEnsemble fits complete.\n")
cat("  Workspace root: ", biomod_root, "\n", sep = "")
cat("  Next step:      source(\"05_project_ensemble.R\")\n")



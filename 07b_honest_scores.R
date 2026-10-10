# =============================================================================
# 07c_ensemble_heldout_eval.R
# TRUE spatially-blocked, out-of-fold evaluation of the DEPLOYED ENSEMBLE.
#
# WHAT THIS IS AND WHY IT'S DIFFERENT FROM 07 / 07b
# ---------------------------------------------------------------------------
# 07_evaluate.R  -> scores the final ensemble raster at the SAME points it was
#                   trained on. In-sample. ~0.99 for everything. Not honest.
# 07b            -> surfaces biomod2's per-SINGLE-MODEL held-out CV scores.
#                   Honest, but it's the base models, not the ensemble you ship.
# 07c (this)     -> reconstructs a genuine OUT-OF-FOLD prediction for the
#                   ENSEMBLE and scores that. This is the number a reviewer
#                   wants for the deployed suitability surface.
#
# HOW IT WORKS (no refitting, no raster work)
# ---------------------------------------------------------------------------
# During block CV, biomod2 fit each base model on 3 of the 4 ENMeval spatial
# quadrants and stored its prediction at EVERY training point -- including the
# held-out quadrant. get_predictions() gives those stored predictions;
# get_calib_lines() flags, per run, which points were held out (calibration ==
# FALSE). Each point is held out in exactly one run, so for every point we can
# take the base-model predictions made by the models that did NOT see it, then
# combine them into an ensemble. Pool across all points -> a fully out-of-fold
# ensemble prediction vector. Score it.
#
# This is faithful to your pipeline: it uses biomod2's own fitted models and
# its own fold assignment. It also works unchanged for any species that fell
# back to k-fold CV -- it just uses whatever folds biomod2 stored.
#
# HEADLINE = EMmean (unweighted committee average of the selected base models).
# It's the clean number: no weight leakage from the held-out fold. An
# EMwmean-by-TSS column is also produced for comparison with what you deploy;
# it weights by each base model's CV TSS, so treat it as indicative, not as the
# primary honest figure.
#
# PREREQUISITES: run AFTER 04 (needs models_hab/<sid>/ intact and the
# <sid>_biomod_paths.rds files). biomod2 must load the objects, so don't move
# the models_hab folder after fitting (the stored links are absolute).
#
#   setwd("/path/to/pipeline"); source("07c_ensemble_heldout_eval.R")
#
# NOTE: reloading biomod2 objects touches version-sensitive internals. The
# script prints the structure of the first species' prediction/calib tables
# once, and detects columns tolerantly. If something mismatches, it stops with
# that printout -- send it to me and it's a one-line fix.
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(biomod2)
  library(dplyr)
  library(pROC)
})

log_line("heldout_eval", "-", "starting")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ── METRIC HELPERS (identical logic to 07 so numbers are comparable) ─────────
boyce_index <- function(fit_pres, fit_bg) {
  pres_ok <- fit_pres[is.finite(fit_pres)]
  bg_ok   <- fit_bg[is.finite(fit_bg)]
  if (length(pres_ok) < 10 || length(bg_ok) < 10) return(NA_real_)
  if (requireNamespace("ecospat", quietly = TRUE)) {
    out <- tryCatch(
      ecospat::ecospat.boyce(fit = c(pres_ok, bg_ok), obs = pres_ok,
                             nclass = 0, PEplot = FALSE),
      error = function(e) NULL)
    if (!is.null(out) && is.finite(out$cor)) return(as.numeric(out$cor))
  }
  nclass <- 20
  fit_all <- c(pres_ok, bg_ok)
  breaks <- seq(min(fit_all), max(fit_all), length.out = nclass + 1)
  mids   <- head(breaks, -1) + diff(breaks) / 2
  Fi <- vapply(seq_along(mids), function(i) {
    lo <- breaks[i]; hi <- breaks[i + 1]
    p_pres <- mean(pres_ok >= lo & pres_ok < hi)
    p_bg   <- mean(bg_ok   >= lo & bg_ok   < hi)
    if (p_bg == 0) NA_real_ else p_pres / p_bg
  }, numeric(1))
  ok <- is.finite(Fi)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(mids[ok], Fi[ok], method = "spearman"))
}

calibration_slope <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred) & pred > 0 & pred < 1
  if (sum(ok) < 20 || length(unique(obs[ok])) < 2) return(NA_real_)
  x <- log(pred[ok] / (1 - pred[ok]))
  m <- tryCatch(glm(obs[ok] ~ x, family = binomial()), error = function(e) NULL)
  if (is.null(m)) NA_real_ else as.numeric(coef(m)[2])
}

metrics_from <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (length(unique(obs)) < 2 || sum(obs == 1) < 10)
    return(c(AUC = NA, TSS = NA, Boyce = NA, calib = NA))
  roc <- tryCatch(pROC::roc(obs, pred, quiet = TRUE, direction = "<"),
                  error = function(e) NULL)
  auc <- if (is.null(roc)) NA_real_ else as.numeric(pROC::auc(roc))
  tss <- NA_real_
  if (!is.null(roc)) {
    cc  <- pROC::coords(roc, "all", ret = c("sensitivity", "specificity"))
    tss <- max(cc$sensitivity + cc$specificity - 1, na.rm = TRUE)
  }
  c(AUC   = auc,
    TSS   = tss,
    Boyce = boyce_index(pred[obs == 1], pred[obs == 0]),
    calib = calibration_slope(obs, pred))
}

# ── ONE-TIME STRUCTURE PRINT so column drift is visible, not silent ──────────
.printed_structure <- FALSE
show_structure_once <- function(preds, calib) {
  if (.printed_structure) return(invisible())
  cat("\n--- get_predictions() columns ---\n");  print(names(preds))
  cat("--- head ---\n"); print(utils::head(preds, 3))
  cat("\n--- get_calib_lines(as.data.frame=TRUE) columns ---\n"); print(names(calib))
  cat("--- head ---\n"); print(utils::head(calib, 3))
  cat("\n")
  .printed_structure <<- TRUE
}

pick <- function(nms, cands) { hit <- intersect(cands, nms); if (length(hit)) hit[1] else NA_character_ }
run_int <- function(x) suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))

# ── PER-SPECIES OUT-OF-FOLD ENSEMBLE EVALUATION ──────────────────────────────
eval_one <- function(sp) {
  sid <- sp$short_id
  paths_file <- file.path(MODELS_DIR, sid, paste0(sid, "_biomod_paths.rds"))
  if (!file.exists(paths_file)) { log_line("heldout_eval", sid, "no paths; skip"); return(NULL) }
  paths <- readRDS(paths_file)
  
  bm_out <- tryCatch({ nm <- load(paths$bm_out_link); get(nm[1]) },
                     error = function(e) { log_line("heldout_eval", sid,
                                                    sprintf("couldn't load models.out (%s)", conditionMessage(e))); NULL })
  if (is.null(bm_out)) return(NULL)
  
  preds <- tryCatch(biomod2::get_predictions(bm_out), error = function(e) NULL)
  calib <- tryCatch(biomod2::get_calib_lines(bm_out, as.data.frame = TRUE),
                    error = function(e) NULL)
  if (is.null(preds) || is.null(calib)) {
    log_line("heldout_eval", sid, "no stored predictions/calib lines; skip"); return(NULL)
  }
  show_structure_once(preds, calib)
  
  # tolerant column mapping
  c_pt_p  <- pick(names(preds), c("points", "point"))
  c_val   <- pick(names(preds), c("pred", "prediction", "value"))
  c_run_p <- pick(names(preds), c("run", "Run", "PA.run"))
  c_algo  <- pick(names(preds), c("algo", "Model", "model"))
  c_pt_c  <- pick(names(calib), c("points", "point"))
  c_run_c <- pick(names(calib), c("run", "Run", "PA.run"))
  c_flag  <- pick(names(calib), c("calib.lines", "value", "calib"))
  if (any(is.na(c(c_pt_p, c_val, c_run_p, c_algo, c_pt_c, c_run_c, c_flag))))
    stop("Column detection failed for ", sid,
         " -- see the structure printout above and report the column names.")
  
  P <- data.frame(points = preds[[c_pt_p]],
                  run_i  = run_int(preds[[c_run_p]]),
                  algo   = as.character(preds[[c_algo]]),
                  pred   = suppressWarnings(as.numeric(preds[[c_val]])),
                  stringsAsFactors = FALSE)
  C <- data.frame(points = calib[[c_pt_c]],
                  run_i  = run_int(calib[[c_run_c]]),
                  held   = !as.logical(calib[[c_flag]]),   # calib==FALSE -> held out
                  stringsAsFactors = FALSE)
  
  if (max(P$pred, na.rm = TRUE) > 1.5) P$pred <- P$pred / 1000  # biomod stores 0..1000
  
  # keep only held-out (point, run) combinations, then average per point x algo
  oof <- P %>%
    dplyr::inner_join(C[C$held, c("points", "run_i")], by = c("points", "run_i")) %>%
    dplyr::group_by(points, algo) %>%
    dplyr::summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")
  if (!nrow(oof)) { log_line("heldout_eval", sid, "no held-out rows; skip"); return(NULL) }
  
  # observed response, aligned to biomod's point index
  resp <- tryCatch(biomod2::get_formal_data(bm_out, subinfo = "resp.var"),
                   error = function(e) NULL)
  if (is.null(resp)) {
    tr <- readRDS(file.path(GBIF_DIR, paste0(sid, "_train_hab.rds")))$train
    resp <- tr$pa
  }
  resp[is.na(resp)] <- 0L
  obs_df <- data.frame(points = seq_along(resp), obs = as.integer(resp))
  
  # per-base-model CV TSS -> selection + optional weights
  ev <- tryCatch(biomod2::get_evaluations(bm_out), error = function(e) NULL)
  algo_tss <- NULL
  if (!is.null(ev)) {
    m_metric <- pick(names(ev), c("metric.eval", "metric"))
    m_valid  <- pick(names(ev), c("validation", "Testing.data", "valid"))
    m_algo   <- pick(names(ev), c("algo", "Model"))
    if (!any(is.na(c(m_metric, m_valid, m_algo)))) {
      e2 <- ev[toupper(as.character(ev[[m_metric]])) == "TSS", ]
      e2$.v <- suppressWarnings(as.numeric(e2[[m_valid]]))
      algo_tss <- tapply(e2$.v, as.character(e2[[m_algo]]), mean, na.rm = TRUE)
    }
  }
  
  algos_all <- sort(unique(oof$algo))
  selected  <- if (!is.null(algo_tss)) {
    s <- names(algo_tss)[is.finite(algo_tss) & algo_tss >= METRIC_THRESH]
    if (length(s)) intersect(algos_all, s) else algos_all
  } else algos_all
  
  # Combine base-model OOF predictions into the ensemble, per point, with
  # tapply over a single grouping vector. (base reshape() silently mangles a
  # grouped dplyr table and returned an all-NA matrix -> all metrics NA.)
  sel <- as.data.frame(oof[oof$algo %in% selected, , drop = FALSE])
  em_mean_v <- tapply(sel$pred, sel$points, mean, na.rm = TRUE)
  
  if (!is.null(algo_tss)) {
    wv  <- as.numeric(algo_tss[sel$algo]); wv[!is.finite(wv)] <- 0
    num <- tapply(sel$pred * wv, sel$points, sum, na.rm = TRUE)
    den <- tapply(wv,            sel$points, sum, na.rm = TRUE)
    em_wmean_v <- num / den
  } else {
    em_wmean_v <- em_mean_v
  }
  
  pts <- as.integer(names(em_mean_v))
  ens <- data.frame(points   = pts,
                    em_mean  = as.numeric(em_mean_v),
                    em_wmean = as.numeric(em_wmean_v[as.character(pts)]),
                    stringsAsFactors = FALSE)
  ens <- merge(ens, obs_df, by = "points")
  
  if (!any(is.finite(ens$em_mean)))
    log_line("heldout_eval", sid, "WARN: ensemble predictions all non-finite")
  
  mm <- metrics_from(ens$obs, ens$em_mean)
  mw <- metrics_from(ens$obs, ens$em_wmean)

  # Persist per-species OOF ensemble vector so 06_thresholds.R can derive
  # thresholds from honest predictions instead of training-set predictions.
  oof_path <- file.path(EVAL_DIR, sprintf("%s_oof_ensemble.rds", sid))
  saveRDS(list(species = sid, points = ens$points, obs = ens$obs,
               em_mean = ens$em_mean, em_wmean = ens$em_wmean),
          oof_path)
  
  # per-base-model out-of-fold scores (diagnostic)
  per_algo <- lapply(algos_all, function(a) {
    d <- merge(oof[oof$algo == a, c("points", "pred")], obs_df, by = "points")
    m <- metrics_from(d$obs, d$pred)
    data.frame(algo = a, oof_AUC = round(m["AUC"], 3), oof_TSS = round(m["TSS"], 3))
  })
  per_algo_df <- dplyr::bind_rows(per_algo)
  per_algo_df$species <- sid
  
  log_line("heldout_eval", sid,
           sprintf("EMmean OOF  AUC=%.3f TSS=%.3f Boyce=%.3f calib=%.2f  (n_pres=%d n_bg=%d, algos=%s)",
                   mm["AUC"], mm["TSS"], mm["Boyce"], mm["calib"],
                   sum(ens$obs == 1), sum(ens$obs == 0), paste(selected, collapse = "/")))
  
  list(
    summary = data.frame(
      species        = sid,
      latin          = sp$latin,
      n_pres         = sum(ens$obs == 1),
      n_bg           = sum(ens$obs == 0),
      algos_used     = paste(selected, collapse = "/"),
      EMmean_AUC     = round(mm["AUC"], 3),
      EMmean_TSS     = round(mm["TSS"], 3),
      EMmean_Boyce   = round(mm["Boyce"], 3),
      EMmean_calib   = round(mm["calib"], 3),
      EMwmean_AUC    = round(mw["AUC"], 3),
      EMwmean_TSS    = round(mw["TSS"], 3),
      row.names = NULL, stringsAsFactors = FALSE),
    per_algo = per_algo_df)
}

# ── RUN ──────────────────────────────────────────────────────────────────────
sum_rows <- list(); algo_rows <- list()
for (sp in SPECIES) {
  res <- tryCatch(eval_one(sp), error = function(e) {
    log_line("heldout_eval", sp$short_id, sprintf("ERROR: %s", conditionMessage(e))); NULL })
  if (!is.null(res)) { sum_rows[[sp$short_id]] <- res$summary; algo_rows[[sp$short_id]] <- res$per_algo }
}

if (!length(sum_rows)) stop("No species evaluated. Check the structure printout / logs above.")

out <- dplyr::bind_rows(sum_rows)
out_path <- file.path(EVAL_DIR, "ensemble_blockcv_heldout.csv")
write.csv(out, out_path, row.names = FALSE)

algo_out <- dplyr::bind_rows(algo_rows)
write.csv(algo_out, file.path(EVAL_DIR, "ensemble_blockcv_heldout_by_algo.csv"),
          row.names = FALSE)

# ── SIDE-BY-SIDE WITH THE IN-SAMPLE NUMBERS ──────────────────────────────────
app_path <- file.path(EVAL_DIR, "eval_insample.csv")
cat("\n===================================================================\n")
cat("DEPLOYED ENSEMBLE — spatial block, OUT-OF-FOLD (honest)\n")
cat("===================================================================\n\n")
if (file.exists(app_path)) {
  app <- read.csv(app_path, stringsAsFactors = FALSE)[, c("species", "AUC_insample", "TSS_insample")]
  names(app) <- c("species", "inSample_AUC", "inSample_TSS")
  cmp <- merge(out[, c("species", "EMmean_AUC", "EMmean_TSS", "EMmean_Boyce", "EMmean_calib")],
               app, by = "species")
  cmp$AUC_drop <- round(cmp$inSample_AUC - cmp$EMmean_AUC, 3)
  cmp$TSS_drop <- round(cmp$inSample_TSS - cmp$EMmean_TSS, 3)
  print(cmp[, c("species", "EMmean_TSS", "inSample_TSS", "TSS_drop",
                "EMmean_AUC", "inSample_AUC", "AUC_drop", "EMmean_Boyce", "EMmean_calib")],
        row.names = FALSE)
} else {
  print(out[, c("species", "EMmean_AUC", "EMmean_TSS", "EMmean_Boyce", "EMmean_calib")],
        row.names = FALSE)
}

cat("\nWrote:\n  ", out_path, "\n  ",
    file.path(EVAL_DIR, "ensemble_blockcv_heldout_by_algo.csv"), "\n", sep = "")
cat("\nEMmean_* are your honest deployed-surface numbers. Report these.\n")
cat("EMmean_Boyce and EMmean_calib on held-out data tell you whether the\n")
cat("suitability values mean anything out-of-sample -- the calibration slope\n")
cat("especially, since that's what your webapp color ramp is showing.\n")

log_line("heldout_eval", "-", "done")
# =============================================================================
# 07_evaluate.R
# IN-SAMPLE evaluation of the deployed ensemble, for diagnostic use only.
#
# Scores the final ensemble raster at the SAME points that were used to fit
# it (training presences + background). The resulting AUC/TSS/Boyce/calib are
# resubstitution numbers, not validation. They are typically ~0.99 and should
# NOT be reported as model performance. For honest held-out scores, see
# 07b_honest_scores.R -> eval_hab/ensemble_blockcv_heldout.csv.
#
# Output:
#   EVAL_DIR/eval_insample.csv   one row per species, columns suffixed _insample
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(pROC)
})

log_line("evaluate", "-", "starting")

# base R < 4.4 doesn't have %||%; define our own defensively
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ── METRIC HELPERS ───────────────────────────────────────────────────────────
# Boyce index (Hirzel et al. 2006). Prefers ecospat::ecospat.boyce (moving
# window, correct formulation). Falls back to a 20-bin hand-rolled version if
# ecospat isn't available — the old 10-bin version saturated at 1.00 for most
# species and was effectively useless.
boyce_index <- function(fit_at_presences, fit_background) {
  pres_ok <- fit_at_presences[is.finite(fit_at_presences)]
  bg_ok   <- fit_background[is.finite(fit_background)]
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
  if (sum(ok) < 20 || length(unique(obs[ok])) < 2)
    return(list(slope = NA, intercept = NA))
  x <- log(pred[ok] / (1 - pred[ok]))
  m <- tryCatch(glm(obs[ok] ~ x, family = binomial()), error = function(e) NULL)
  if (is.null(m)) return(list(slope = NA, intercept = NA))
  list(slope     = as.numeric(coef(m)[2]),
       intercept = as.numeric(coef(m)[1]))
}

# ── MAIN LOOP ────────────────────────────────────────────────────────────────
rows <- list()

for (sp in SPECIES) {
  sid <- sp$short_id
  log_line("evaluate", sid, sprintf("=== %s ===", sp$latin))

  train_path <- file.path(GBIF_DIR, paste0(sid, "_train_hab.rds"))
  cur_path   <- file.path(CURRENT_DIR, sprintf("%s_ensemble_current.tif", sid))
  if (!file.exists(train_path) || !file.exists(cur_path)) {
    log_line("evaluate", sid, "missing training or current raster; skipping")
    next
  }

  cur_r <- terra::rast(cur_path)
  if (max(terra::minmax(cur_r)[2, ], na.rm = TRUE) > 1.5) cur_r <- cur_r / 1000

  train <- readRDS(train_path)$train
  pts   <- terra::vect(train[, c("x", "y")],
                       geom = c("x", "y"), crs = "EPSG:4326")
  pred  <- terra::extract(cur_r, pts)[, 2]
  obs   <- train$pa

  ok <- is.finite(pred) & is.finite(obs)
  pred_ok <- pred[ok]
  obs_ok  <- obs[ok]

  if (length(unique(obs_ok)) < 2 || sum(obs_ok == 1) < 10) {
    log_line("evaluate", sid, "too few presences with valid predictions")
    next
  }

  # AUC
  auc_val <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(obs_ok, pred_ok, quiet = TRUE, direction = "<"))),
    error = function(e) NA_real_)

  # TSS (max across all thresholds)
  roc_obj <- tryCatch(pROC::roc(obs_ok, pred_ok, quiet = TRUE, direction = "<"),
                      error = function(e) NULL)
  tss_val <- NA_real_
  if (!is.null(roc_obj)) {
    cc <- pROC::coords(roc_obj, "all",
                       ret = c("sensitivity", "specificity"))
    tss_val <- max(cc$sensitivity + cc$specificity - 1, na.rm = TRUE)
  }

  # Boyce
  boy_val <- boyce_index(pred_ok[obs_ok == 1], pred_ok[obs_ok == 0])

  # Calibration
  calib <- calibration_slope(obs_ok, pred_ok)

  rows[[sid]] <- data.frame(
    species             = sid,
    latin               = sp$latin,
    n_pres              = sum(obs_ok == 1),
    n_bg                = sum(obs_ok == 0),
    AUC_insample        = round(auc_val, 3),
    TSS_insample        = round(tss_val, 3),
    Boyce_insample      = round(boy_val, 3),
    calib_slope_insample     = round(calib$slope, 3),
    calib_intercept_insample = round(calib$intercept, 3),
    stringsAsFactors = FALSE
  )

  log_line("evaluate", sid,
           sprintf("in-sample AUC=%.3f TSS=%.3f Boyce=%.3f calib=%.2f",
                   auc_val, tss_val, boy_val, calib$slope %||% NA))
}

if (length(rows)) {
  out <- dplyr::bind_rows(rows)
  out_path <- file.path(EVAL_DIR, "eval_insample.csv")
  write.csv(out, out_path, row.names = FALSE)
  # Clean up the misleadingly-named legacy file so no reader confuses them.
  legacy <- file.path(EVAL_DIR, "eval_summary.csv")
  if (file.exists(legacy)) file.remove(legacy)
  log_line("evaluate", "-",
           sprintf("wrote %s (%d species, in-sample)", basename(out_path), nrow(out)))
  cat("\nIn-sample evaluation (resubstitution, not honest performance):\n")
  print(out[, c("species", "n_pres", "AUC_insample", "TSS_insample",
                "Boyce_insample", "calib_slope_insample")],
        row.names = FALSE)
  cat("\nReport ensemble_blockcv_heldout.csv (from 07b) as the honest numbers.\n")
}

log_line("evaluate", "-", "done")

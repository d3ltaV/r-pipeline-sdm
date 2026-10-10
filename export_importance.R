# ---- export_importance.R -----------------------------------------------------
# Turn the biomod variable-importance table into data/importance.json so the
# atlas click panel can show a per-species "what drives this species" breakdown.
# Structure: { "<species>": { "<predictor>": importance, ... }, ... }
# Defensive: it auto-detects columns and PRINTS what it found — if it can't,
# paste the "columns:" line back to me. Paste into the RStudio Console.
# -----------------------------------------------------------------------------
if (!file.exists("00_config.R"))
  setwd("C:/Users/Joelle/OneDrive/Desktop/INVASIVE-SDM-PREDICTION-FINAL-PIPELINE")
if (!exists("CONFIG_V3_LOADED")) source("00_config.R")
suppressPackageStartupMessages(library(jsonlite))

# 1. Locate the CSV (searches recursively if not in eval_hab/).
cand <- unique(c(
  file.path("eval_hab", "biomod_var_importance.csv"),
  "biomod_var_importance.csv",
  list.files(".", pattern = "biomod_var_importance\\.csv$", recursive = TRUE, full.names = TRUE)
))
f <- cand[file.exists(cand)][1]
if (is.na(f))
  stop("Couldn't find biomod_var_importance.csv under ", getwd(), " — tell me its path.")
cat("reading:", f, "\n")

df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
cat("columns:", paste(names(df), collapse = ", "), "\n")
print(utils::head(df, 4))

# 2. Auto-detect the species / predictor / importance columns.
lc   <- tolower(trimws(names(df)))
pick <- function(...) { for (x in c(...)) { i <- which(lc == x); if (length(i)) return(names(df)[i[1]]) }; NA_character_ }
sp_c <- pick("species", "short_id", "sid", "sp")
pv_c <- pick("expl.var", "predictor", "variable", "var", "expl_var", "explvar")
im_c <- pick("var.imp", "importance", "value", "imp", "var_imp", "varimp", "mean")

if (any(is.na(c(sp_c, pv_c, im_c))))
  stop("Couldn't auto-detect columns.\n",
       "  Detected -> species: ", sp_c, " | predictor: ", pv_c, " | importance: ", im_c, "\n",
       "  Paste the 'columns:' line above to me and I'll adjust the script.")
cat(sprintf("using   -> species: %s | predictor: %s | importance: %s\n", sp_c, pv_c, im_c))

# 3. Average importance over algorithms / CV runs, per (species, predictor).
d <- data.frame(species   = as.character(df[[sp_c]]),
                predictor = as.character(df[[pv_c]]),
                imp       = suppressWarnings(as.numeric(df[[im_c]])),
                stringsAsFactors = FALSE)
d <- d[!is.na(d$imp), ]
agg <- aggregate(imp ~ species + predictor, data = d, FUN = function(x) mean(x, na.rm = TRUE))

# 4. Emit JSON.
out <- list()
for (s in unique(agg$species)) {
  sub <- agg[agg$species == s, ]
  out[[s]] <- setNames(as.list(round(sub$imp, 4)), sub$predictor)
}
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, na = "null", digits = 4),
           file.path(WEBAPP_DATA, "importance.json"))

cat("\nWrote", file.path(WEBAPP_DATA, "importance.json"), "for", length(out), "species.\n")
cat("species keys:  ", paste(names(out), collapse = ", "), "\n")
cat("predictor keys:", paste(sort(unique(agg$predictor)), collapse = ", "), "\n")
cat("The keys above must match the webapp's ids (species: autumn_olive, burning_bush, ...;\n",
    "predictors: mat, ap, trng, pseas, hfp, treecov, soilph, elev). If they differ,\n",
    "tell me and I'll add a name map.\n", sep = "")

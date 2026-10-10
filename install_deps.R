# =============================================================================
# install_deps.R
# One-shot dependency installer for the habitat-augmented SDM pipeline.
#
# Run once per machine, before RUN_ALL.R:
#   source("install_deps.R")
#
# Skips anything already installed. Uses the RStudio CRAN mirror by default —
# override with options(repos = c(CRAN = "https://your.mirror/")) if you have
# an institutional preference.
#
# System-level dependencies (install these OUTSIDE R):
#   Windows: RTools 4.3+ (matching your R version) — https://cran.r-project.org/bin/windows/Rtools/
#   macOS:   Xcode Command Line Tools (`xcode-select --install`)
#            + GDAL/GEOS/PROJ (`brew install gdal geos proj`)
#   Linux:   libgdal-dev, libgeos-dev, libproj-dev, libudunits2-dev
#            (e.g. `sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev`)
#
# After install, a GBIF account is required for script 02 (free at gbif.org).
# Add to ~/.Renviron (or Sys.setenv() in .Rprofile) and restart R:
#   GBIF_USER=<your username>
#   GBIF_PWD=<your password>
#   GBIF_EMAIL=<your email>
# =============================================================================

if (getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# ── PACKAGE LIST ─────────────────────────────────────────────────────────────
# Grouped by pipeline stage so it's obvious which package each script needs.

pkg_groups <- list(
  # Raster / spatial core (all scripts)
  "spatial core"     = c("terra", "sf"),

  # Data download (01, 01b, 02)
  "downloads"        = c("geodata", "rgbif"),

  # Training data prep (03)
  "training prep"    = c("CoordinateCleaner", "usdm", "dplyr", "tidyr"),

  # Modeling (04, 05, 06, 07)
  "modeling"         = c("biomod2", "blockCV", "pROC", "ecospat"),

  # Visuals (09)
  "figures"          = c("ggplot2", "tidyterra", "patchwork"),

  # Export + web (11, 12, 13, aggregate_counties)
  "export / webapp"  = c("jsonlite", "exactextractr", "rmapshaper")
)

# Optional: for parallel projection in 05 (nice-to-have, not required)
optional_pkgs <- c("future.apply", "future")

# ── INSTALLER ────────────────────────────────────────────────────────────────
installed <- rownames(installed.packages())

report <- data.frame(
  group   = character(),
  package = character(),
  status  = character(),
  stringsAsFactors = FALSE
)

install_group <- function(group_name, pkgs) {
  cat(sprintf("\n== %s (%d packages) ==\n", group_name, length(pkgs)))
  for (p in pkgs) {
    if (p %in% installed) {
      cat(sprintf("  [ok]      %s\n", p))
      report[nrow(report) + 1, ] <<- c(group_name, p, "already installed")
      next
    }
    cat(sprintf("  [install] %s ...\n", p))
    ok <- tryCatch({
      install.packages(p, dependencies = TRUE, quiet = TRUE)
      p %in% rownames(installed.packages())
    }, error = function(e) {
      cat(sprintf("    !! %s\n", conditionMessage(e)))
      FALSE
    })
    report[nrow(report) + 1, ] <<- c(
      group_name, p,
      if (ok) "installed" else "FAILED"
    )
  }
}

for (g in names(pkg_groups)) install_group(g, pkg_groups[[g]])
install_group("parallel (optional)", optional_pkgs)

# ── VERIFY ───────────────────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("  SUMMARY\n")
cat("============================================================\n\n")
print(report, row.names = FALSE)

failed <- report[report$status == "FAILED", "package"]
if (length(failed)) {
  cat("\n!!! FAILED PACKAGES:", paste(failed, collapse = ", "), "\n")
  cat("Common fixes:\n")
  cat("  - Windows: install RTools matching your R version and restart R\n")
  cat("  - macOS:   brew install gdal geos proj udunits, then retry\n")
  cat("  - Linux:   sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev\n")
  cat("  - biomod2 sometimes needs 'Rcpp' rebuilt first: install.packages('Rcpp')\n\n")
} else {
  cat("\nAll good. Next: read QUICKSTART.md, then source('RUN_ALL.R').\n\n")
}

# ── GBIF CREDENTIAL CHECK ────────────────────────────────────────────────────
cat("--- GBIF credential check ---\n")
u <- Sys.getenv("GBIF_USER")
p <- Sys.getenv("GBIF_PWD")
e <- Sys.getenv("GBIF_EMAIL")
if (nzchar(u) && nzchar(p) && nzchar(e)) {
  cat(sprintf("  GBIF_USER=%s  GBIF_EMAIL=%s  (password set)\n", u, e))
  cat("  If script 02 fails with auth errors, verify these at gbif.org.\n")
} else {
  cat("  GBIF credentials not set. Script 02 will fail.\n")
  cat("  Add to ~/.Renviron and restart R:\n")
  cat("    GBIF_USER=<your username>\n")
  cat("    GBIF_PWD=<your password>\n")
  cat("    GBIF_EMAIL=<your email>\n\n")
}

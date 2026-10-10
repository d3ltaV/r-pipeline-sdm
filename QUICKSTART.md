# QUICKSTART — Invasive-Species SDM Pipeline (habitat-augmented, from scratch)

End-to-end species distribution modeling for 12 northeastern invasive plants across the US and Canada. Uses biomod2 (GAM / GLM / GBM / RF / MAXNET ensemble) with 8 predictors: 4 WorldClim bioclim variables plus 4 static habitat layers (elevation, Human Footprint, tree cover, soil pH). Projects onto WorldClim current + 40 CMIP6 futures (5 GCMs × 2 SSPs × 4 time periods).

This bundle contains **every file needed to run the pipeline from scratch on a fresh machine**. No prior data required — all downloads are automatic.

---

## The 4-command version

```r
setwd("/path/to/your/project/folder")   # after you copy this bundle there
source("install_deps.R")                # first time only, ~10 min
# Set GBIF credentials in ~/.Renviron and restart R (see Step 3)
source("RUN_ALL.R")                     # 8–18 hours cold start
```

That's it. Read on if you want to know what's happening.

---

## Step 1 — System prerequisites

### R
- R 4.3 or newer. Download from https://cran.r-project.org/.
- RStudio recommended (https://posit.co/download/rstudio-desktop/) but not required.

### Native libraries (needed to compile `terra`, `sf`, `biomod2`)

**Windows.** Install RTools matching your R version from https://cran.r-project.org/bin/windows/Rtools/. That's it.

**macOS.**
```
xcode-select --install
brew install gdal geos proj udunits
```

**Linux (Ubuntu/Debian).**
```
sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev
```

### Disk space
About 35–45 GB total when everything's downloaded and fit:
- Climate rasters: ~15 GB (WorldClim current + 40 CMIP6 futures)
- Habitat rasters: ~500 MB
- GBIF caches: ~500 MB
- biomod2 workspaces + projection outputs: ~10–15 GB

---

## Step 2 — Copy this bundle to your project folder

Pick a project root. Windows example: `C:/Users/<you>/OneDrive/Desktop/TRAILS-hab-pipeline/`. macOS/Linux: `~/projects/trails-hab-pipeline/`. Absolute path, no spaces if possible.

Copy every file from this bundle into that directory. When you're done it should look like this:

```
TRAILS-hab-pipeline/
├── install_deps.R              ← Run first (once)
├── QUICKSTART.md               ← This file
├── README_hab.md               ← Deep-dive on the habitat predictors
├── RUN_ALL.R                   ← Main driver
├── 00_config.R                 ← Paths, species list, model settings
├── 01_download_climate.R       ← WorldClim + CMIP6
├── 01b_download_habitat.R      ← elev, HFP, tree cover, soil pH
├── 02_download_gbif.R          ← GBIF occurrences (needs credentials)
├── 03_build_training.R         ← Clean, thin, background, extract
├── 04_fit_ensemble.R           ← biomod2 fitting (slow)
├── 05_project_ensemble.R       ← Current + 40 futures (slow)
├── 06_thresholds.R             ← per-species maxTSS
├── 07_evaluate.R               ← AUC / TSS / Boyce / calibration
├── 08_partial_responses.R      ← Response curves (optional, OFF by default)
├── 09_map_panels.R             ← QC 4-panel figures
├── 10_simulate_spread.R        ← CA dispersal (optional, OFF by default)
├── 11_export_occurrences.R     ← Presence + GeoJSON for webapp
├── 12_copy_to_webapp.R         ← Sync webapp mirror
├── 13_management_units.R       ← State/province rollup (optional, OFF)
└── aggregate_counties.R        ← County choropleth (run manually after pipeline)
```

Do NOT create `data/`, `cache/`, `models_hab/`, `output/`, `eval_hab/`, or any other output folder — the scripts create them automatically.

---

## Step 3 — Credentials (GBIF only)

Script 02 needs a free GBIF account. Register at https://www.gbif.org/user/profile.

Add these three lines to `~/.Renviron` (create the file if it doesn't exist) and **restart R**:

```
GBIF_USER=<your GBIF username>
GBIF_PWD=<your GBIF password>
GBIF_EMAIL=<your GBIF email>
```

On Windows, `~/.Renviron` lives at `C:/Users/<you>/.Renviron` (or `C:/Users/<you>/Documents/.Renviron`).

Verify after restarting R:
```r
Sys.getenv("GBIF_USER")     # should print your username
```

If this is blank, the file isn't where R is looking. Try `usethis::edit_r_environ()` — it opens the right file for your R install.

---

## Step 4 — Install R packages

```r
setwd("/path/to/TRAILS-hab-pipeline")
source("install_deps.R")
```

Takes 5–20 min on a fresh machine (biomod2 has many transitive dependencies). The script prints a summary at the end. If any package fails, it prints platform-specific fixes.

---

## Step 5 — Run the pipeline

```r
setwd("/path/to/TRAILS-hab-pipeline")
source("RUN_ALL.R")
```

RUN_ALL.R executes 13 scripts in order. Every download step checks its cache first, so re-running after an interrupt resumes rather than restarting.

### Runtime expectations (laptop, 8 GB RAM, no parallelism)

| Step | What it does | Cold start | Warm (cached) |
| --- | --- | --- | --- |
| 01 | WorldClim + 40 CMIP6 futures | 1–2 h | seconds |
| 01b | 4 habitat layers | 15–45 min | seconds |
| 02 | GBIF occurrences | 1–3 h | seconds |
| 03 | Training tables (12 species) | ~6 min | seconds |
| 04 | biomod2 fits | 2–8 h | seconds per already-fit species |
| 05 | Projections (current + 40 futures × 12 species) | 3–8 h | seconds per already-projected scenario |
| 06 | Species thresholds | ~2 min | ~2 min |
| 07 | Evaluation report | ~5 min | ~5 min |
| 09 | QC maps | ~3 min | ~3 min |
| 11 | Occurrences export | ~2 min | ~2 min |
| 12 | Webapp sync | ~1 min | ~1 min |

Total cold start: **8–18 hours** depending on CPU and network. Kick it off before bed; check in the morning.

If your laptop sleeps or you Ctrl-C, just `source("RUN_ALL.R")` again — everything is cached and resumable at species-level granularity for 04 and species×scenario granularity for 05.

---

## Step 6 — Optional post-pipeline: county aggregation

If you're using the webapp, run this after everything else finishes:

```r
setwd("/path/to/TRAILS-hab-pipeline")
source("aggregate_counties.R")
```

Downloads US + Canada county polygons (~2 min first time), builds a shared pixel→county lookup grid, and writes per-species per-scenario county-mean JSON files into `WEBAPP_DATA/<species>/counties.json`. Only needed if your webapp renders EDDMapS-style choropleths.

---

## Step 7 — What outputs to look at

After a successful run:

- **`eval_hab/biomod_var_importance.csv`** — for each species, which of the 8 predictors mattered most. Top 3 per species also printed to the R console when 04 finishes.
- **`eval_hab/ensemble_blockcv_heldout.csv`** — honest per-species AUC / TSS / Boyce / calibration from spatial-block OOF. Report these numbers.
- **`eval_hab/eval_insample.csv`** — resubstitution scores on training data. Diagnostic only; not model performance.
- **`eval_hab/thresholds_summary.csv`** — per-species maxTSS threshold (used by 10, 13, aggregate_counties).
- **`diag_hab/<species>_vif.txt`** — collinearity per species. Watch for VIF > 5 or |r| > 0.7. In the 8-predictor setup, `hfp ↔ treecov` and `elev ↔ soilph` (in mountainous regions) are the most likely to fight.
- **`figures/maps_hab/<species>.png`** — 4-panel QC map (current | 2090 SSP585 | GCM agreement | GCM CV). First place to eyeball whether projections look sensible.
- **`webapp_output_hab/<species>/`** — mirror of what would sync to the webapp: `current.tif`, `<year>_ssp<ssp>.tif`, `presence.tif`, `occurrences.geojson`, plus per-scenario agreement + CV rasters.
- **`logs/pipeline_hab.log`** — full run log with timestamps per step.

---

## Directory layout after everything runs

```
TRAILS-hab-pipeline/
├── (all .R files + docs)
├── data/
│   ├── climate/
│   │   ├── clim_current.tif
│   │   ├── clim_current_global.tif
│   │   └── future/
│   │       └── <GCM>_ssp<ssp>_<year>.tif   × 40
│   ├── gbif/
│   │   ├── <species>_gbif_native.rds       × 12
│   │   ├── <species>_gbif_invaded.rds      × 12
│   │   └── <species>_train_hab.rds         × 12  ← built by 03
│   └── habitat/
│       ├── hab_static.tif
│       └── hab_static_global.tif
├── cache/                    ← geodata + GADM caches (auto-managed)
├── models_hab/               ← biomod2 workspaces per species
├── output/
│   ├── current_hab/          ← <species>_ensemble_current.tif
│   └── future_hab/           ← <species>_ensemble_<year>_ssp<ssp>_median.tif etc
├── eval_hab/                 ← scores + thresholds
├── diag_hab/                 ← VIF reports
├── figures/
│   ├── maps_hab/             ← 4-panel QC PNGs
│   └── partial_hab/          ← response curves (if 08 ran)
├── summary_hab/
├── webapp_output_hab/        ← staged webapp mirror
└── logs/
    └── pipeline_hab.log
```

---

## Common issues

- **`Error: unable to load shared object 'terra.so'`** — GDAL/PROJ/GEOS missing at the OS level. Re-do the native-libraries part of Step 1, then `install.packages("terra")` again.
- **`Error in cmip6_world: HTTP 500`** — WorldClim's CMIP6 mirror is intermittently flaky. The script logs the failure and continues; re-run RUN_ALL.R later and it'll retry only the missing files.
- **`Error in occ_search: authentication failed`** — GBIF credentials wrong or not set. Recheck Step 3, verify `Sys.getenv("GBIF_USER")` returns your username.
- **`Warning: pre-flight climate check failed`** in script 05 — one or more `future/*.tif` files are missing, unreadable, at the wrong extent, or have the wrong number of layers. The log lists which ones. Re-run 01; it'll download just the missing ones from the geodata cache.
- **biomod2 blockCV fails on a species with few presences** — flip `CV_STRATEGY` in `00_config.R` from `"block"` to `"kfold"` and re-run 04.
- **Out of memory during step 05** — this pipeline runs serial by design. If you're still hitting limits, reduce `PROJ_EXTENT` to a smaller region in `00_config.R` and delete `output/current_hab/` + `output/future_hab/` to force a rebuild at the smaller extent.
- **You want to change the species list** — edit the `SPECIES <- list(...)` block in `00_config.R`. Native + invaded ISO3 codes drive which countries GBIF pulls from. After changing, delete the old `data/gbif/*_train_hab.rds` for affected species so 03 rebuilds them.

---

## Deeper reading

- **README_hab.md** — why these 4 habitat predictors, why habitat is static across futures, citations.
- **Comments at the top of each script** — every script documents its inputs, outputs, and design decisions in a `# =====` header block.
- **`logs/pipeline_hab.log` + `logs/05_hab_worker_<species>.log`** — timestamped, per-species checkpoints.

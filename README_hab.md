# README — Habitat-Augmented SDM Pipeline (final_hab_v1)

This is the "why" companion to QUICKSTART.md. Read that first if you just want to run the pipeline.

## What changed vs. the climate-only pipeline

The original pipeline used 4 WorldClim bioclim variables. This one adds 4 static habitat predictors and refits all 11 species ensembles on 8 total predictors.

| # | Var | Source | Rationale |
|---|-----|--------|-----------|
| 1 | mat | WorldClim BIO1 | Annual mean temperature — first-order climate control |
| 2 | trng | WorldClim BIO7 | Temperature annual range — cold-tolerance limit |
| 3 | ap | WorldClim BIO12 | Annual precipitation — water availability |
| 4 | pseas | WorldClim BIO15 | Precipitation seasonality — Mediterranean vs monsoonal signal |
| 5 | **elev** | SRTM/GMTED via `geodata::elevation_global` | Cold-air pooling, ecotone edges, mountain barriers |
| 6 | **hfp** | Venter et al. 2018 (Human Footprint) | Propagule pressure proxy — highest-expected-impact addition for invasive plants |
| 7 | **treecov** | Copernicus GLC 100m via `geodata::landcover` | Shade tolerance / canopy openness for understory invaders (garlic mustard, stiltgrass) |
| 8 | **soilph** | SoilGrids v2 top 5 cm via `geodata::soil_world` | Substrate constraint for calcicoles vs calcifuges (buckthorns strongly prefer calcareous soils) |

### Why these 4 and not others

- **hfp** first because for invasive plants specifically, the single strongest predictor after macroclimate is propagule pressure, and HFP is the best available proxy for it at continental scale (Bradley et al. 2010 Biol Invasions).
- **elev** second because it captures cold-air-pooling and topographic microclimate effects that BIO variables at 5-arcmin resolution miss.
- **treecov** third because our 11 species split cleanly into understory (need shade) and edge/canopy (need light) — climate alone can't separate common_buckthorn from norway_maple at the same latitude.
- **soilph** fourth because the buckthorns and knotweed have well-documented substrate preferences that pH captures well.

### What I rejected

- **Impervious surface** — highly correlated with `hfp` (r > 0.8), no independent signal.
- **Distance to roads** — same problem, plus `hfp` explicitly incorporates it.
- **NDVI** — collinear with `treecov`, and NDVI is very seasonal (WorldClim already handles that indirectly via precipitation seasonality).
- **Slope / aspect** — at 5-arcmin (~9 km) these are almost meaningless. Would only matter at 30 m or finer.
- **Categorical land cover** — biomod2's handling of factor variables in ensemble modeling is fragile. Continuous tree cover fraction gives the same signal without the headaches.

## Why habitat is static across CMIP6 futures

Only WorldClim futures carry a time signal. Soil pH, elevation, tree cover, and Human Footprint are held at their present-day values for every future scenario. Standard practice in invasive SDMs that add non-climate predictors (Ibáñez et al. 2020 Frontiers, Jarnevich et al. 2024 Ecol Lett). Rationale: at ~9 km resolution over 60-year horizons, climate change dwarfs plausible change in soil pH, elevation, or aggregated land cover. HFP does have SSP-based projections available (Kennedy et al. 2019 Nature Comm) — if you ever want to relax this, you'd add per-scenario HFP grids and expand `load_env_future()` in 00_config.R.

## What the `_hab` suffix does

All new outputs land in `_hab`-suffixed directories so this rebuild doesn't touch your existing v2a/final results:

| Old (climate-only) | New (with habitat) |
|-------------------|--------------------|
| `models/` | `models_hab/` |
| `output/current/` | `output/current_hab/` |
| `output/future/` | `output/future_hab/` |
| `eval/` | `eval_hab/` |
| `diag/` | `diag_hab/` |
| `figures/maps/` | `figures/maps_hab/` |
| `figures/partial/` | `figures/partial_hab/` |
| `webapp_output/` | `webapp_output_hab/` |
| `data/gbif/<sid>_train.rds` | `data/gbif/<sid>_train_hab.rds` |

So you can A/B compare climate-only vs habitat-augmented ensembles side by side. Point your webapp at whichever gives better validation.

## What to look at when it finishes

1. **`eval_hab/biomod_var_importance.csv`** — this is the payoff. For each species, which of the 4 new predictors actually did work? Sorted per algorithm, so you'll see if GBM cares about hfp but MAXNET doesn't. Script 04 also prints top-3 per species to the console.

2. **`diag_hab/<sid>_vif.txt`** — flags collinearity. Predictor pairs to watch:
   - `hfp` ↔ `treecov` (dense forest tends to be low-HFP)
   - `elev` ↔ `mat` (temperature falls with altitude — expected, usually manageable)
   - `elev` ↔ `soilph` (mountainous soils are often acidic)

   If any pair hits |r| > 0.7 or VIF > 5, biomod2's random forest and GBM will still handle it fine but the response curves become harder to interpret.

3. **`eval_hab/ensemble_blockcv_heldout.csv`** — honest spatial-block OOF AUC/TSS/Boyce/calibration for the deployed ensemble. These are the numbers to report. `eval_insample.csv` is resubstitution on training data and only useful as a sanity check (it will be ~0.99 for everything).

## Design decisions worth knowing

- **8 predictors, not more.** Adding predictors past ~8 in a presence-background SDM starts causing more overfitting than signal-gain, especially when several predictors share variance. If you want to try more, edit `PREDICTORS` in 00_config.R and expect worse Boyce scores.
- **Static habitat across futures.** See above. If you want dynamic HFP, you'd need to (a) source future HFP grids from Kennedy et al. 2019, (b) rewrite `load_env_future()` to swap the habitat block by scenario, (c) reproject — hours of work, marginal science gain.
- **Same species list, same extent, same GBIF caches.** This rebuild is intentionally minimal-change so we can attribute any AUC delta to the predictor set, not to confounding data changes.
- **Serial-only step 05.** biomod2's projection code has memory issues under future.apply parallelism on wide extents. Serial is slower but reliable. Comments at the top of `05_project_ensemble.R` explain the 2026-08-14 debugging history.

## Property-level analysis (asked about — not built here)

Continental-scale parcel data doesn't exist as a unified dataset. US parcels are fragmented across ~3000 county GIS offices; Canadian provinces vary widely (some free, some paywalled). Two realistic paths if you want sub-county resolution:

1. **MassGIS parcel demo** — MA has free statewide parcels (https://www.mass.gov/info-details/massgis-data-property-tax-parcels). Could support a "TRAILS regional case study" section that shows habitat suitability at parcel level for a slice of your service area.
2. **OSM building-footprint proximity** — for each county, compute "% suitable within 500m of any building footprint" from OpenStreetMap. Continental coverage, ~1-2 hour compute, useful proxy for "developed edge" habitat.

Neither is built in this bundle. Ask if you want a script for either.

## Citations for the design choices above

- Bradley BA et al. 2010. Predicting plant invasions in an era of global change. *Trends in Ecology & Evolution* 25(6): 310–318.
- Hirzel AH et al. 2006. Evaluating the ability of habitat suitability models to predict species presences. *Ecological Modelling* 199: 142–152.
- Ibáñez I et al. 2020. Anthropogenic drivers of forest invasions across scales. *Frontiers in Ecology and Evolution* 8: 573325.
- Jarnevich CS et al. 2024. Invasion informatics with the INHABIT model. *Ecology Letters* 27: e14403.
- Kennedy CM et al. 2019. Managing the middle: a shift in conservation priorities based on the global human modification gradient. *Nature Communications* 10: 4844.
- Valavi R et al. 2021. Predictive performance of presence-only species distribution models. *Ecological Monographs* 92(1): e01486.
- Venter O et al. 2018. Global Terrestrial Human Footprint maps for 1993 and 2009. *Scientific Data* 5: 180067.

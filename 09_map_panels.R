# =============================================================================
# 09_map_panels.R
# QC visualization: per-species 4-panel PNG showing
#   [current suitability] | [end-century ssp585 median]
#   [GCM agreement]       | [across-GCM CV (uncertainty)]
#
# Fast eyeball check that projections look sensible before webapp deployment.
# Not intended for publication.
#
# Output: MAPS_DIR/<sid>.png  (11 x 7.5 in, 140 dpi)
# =============================================================================

if (!exists("CONFIG_V3_LOADED")) source("00_config.R")

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(patchwork)
})

log_line("map_panels", "-", "starting")

DISPLAY_LABEL <- tail(PERIOD_LABELS, 1)     # last period (e.g. "2090")
DISPLAY_SSP   <- "585"

suit_scale  <- scale_fill_distiller(
  palette = "YlOrRd", direction = 1, na.value = NA,
  limits = c(0, 1), name = "P(suit)", trans = "sqrt")
agree_scale <- scale_fill_viridis_c(
  option = "viridis", na.value = NA, limits = c(0, 1), name = "Agreement")
cv_scale <- scale_fill_viridis_c(
  option = "plasma", na.value = NA, name = "CV")

map_theme <- theme_void(base_size = 9) +
  theme(legend.key.width  = grid::unit(0.4, "cm"),
        legend.key.height = grid::unit(0.6, "cm"),
        plot.title = element_text(size = 10, face = "plain", hjust = 0))

for (sp in SPECIES) {
  sid <- sp$short_id
  f_cur   <- file.path(CURRENT_DIR, paste0(sid, "_ensemble_current.tif"))
  f_med   <- file.path(FUTURE_DIR,
                       sprintf("%s_ensemble_%s_ssp%s_median.tif",
                               sid, DISPLAY_LABEL, DISPLAY_SSP))
  f_agree <- file.path(FUTURE_DIR,
                       sprintf("%s_ensemble_%s_ssp%s_agreement.tif",
                               sid, DISPLAY_LABEL, DISPLAY_SSP))
  f_cv    <- file.path(FUTURE_DIR,
                       sprintf("%s_ensemble_%s_ssp%s_cv.tif",
                               sid, DISPLAY_LABEL, DISPLAY_SSP))

  if (!all(file.exists(c(f_cur, f_med)))) {
    log_line("map_panels", sid, "missing current or future; skipping"); next
  }

  p1 <- ggplot() + geom_spatraster(data = rast(f_cur)) + suit_scale +
        coord_sf() + labs(title = "Current") + map_theme
  p2 <- ggplot() + geom_spatraster(data = rast(f_med)) + suit_scale +
        coord_sf() +
        labs(title = sprintf("%s ssp%s median", DISPLAY_LABEL, DISPLAY_SSP)) +
        map_theme
  p3 <- if (file.exists(f_agree)) {
    ggplot() + geom_spatraster(data = rast(f_agree)) + agree_scale +
      coord_sf() + labs(title = "GCM agreement (>=0.5)") + map_theme
  } else patchwork::plot_spacer()
  p4 <- if (file.exists(f_cv)) {
    ggplot() + geom_spatraster(data = rast(f_cv)) + cv_scale +
      coord_sf() + labs(title = "Across-GCM CV (uncertainty)") + map_theme
  } else patchwork::plot_spacer()

  panel <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = paste0(sid, "  (", DEPLOYED_ENSEMBLE, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = 12)))
  out_png <- file.path(MAPS_DIR, paste0(sid, ".png"))
  ggsave(out_png, panel, width = 11, height = 7.5, dpi = 140)
  log_line("map_panels", sid, sprintf("wrote %s.png", sid))
}

log_line("map_panels", "-", "done")

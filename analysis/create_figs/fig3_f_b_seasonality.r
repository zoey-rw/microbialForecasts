# Seasonal amplitude and forecast accuracy relationships
# Central question: does seasonal variability predict how well a taxon is forecast,
# and is this driven by mean abundance (a simpler explanation)?
#
# Panel A: amplitude ~ mean_abun — are abundant taxa more seasonal?
# Panel B: nRMSE ~ amplitude    — do more seasonal taxa have worse (or better) forecasts?
# Panel C: nRMSE ~ mean_abun    — is forecast error simply driven by abundance?

library(tidyverse)
library(ggpubr)

source("source.R")

# ── Data loading ──────────────────────────────────────────────────────────────
seasonal_amplitude_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_vals_scores      <- seasonal_amplitude_in[[6]]
cycl_vals_scores$model_id <- gsub("_beta_regression$", "", cycl_vals_scores$model_id)

scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged   <- gsub("_beta_regression$", "", scores_list$converged_list)

# Hindcast performance (env_cycl), averaged per taxon rank
hindcast_rsq <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged,
         model_name == "env_cycl",
         grepl("observed", site_prediction)) %>%
  group_by(rank_name) %>%
  summarise(
    RMSE.norm      = mean(RMSE.norm,       na.rm = TRUE),
    RSQ            = mean(RSQ,             na.rm = TRUE),
    CRPS_truncated = mean(CRPS_truncated,  na.rm = TRUE),
    .groups = "drop"
  )

# Mean calibration-period abundance per model
# Hindcasts via the package loader (reads + unions the per-model parquet files)
hindcast_data <- load_hindcasts()

mean_abun <- hindcast_data %>%
  filter(model_id %in% converged, !is.na(truth), fcast_period != "hindcast") %>%
  group_by(model_id) %>%
  summarise(mean_abun = mean(truth, na.rm = TRUE), .groups = "drop")

# ── Build analysis dataset ────────────────────────────────────────────────────
# Use cycl_only amplitude; hindcast nRMSE from env_cycl
to_plot <- cycl_vals_scores %>%
  filter(model_name == "cycl_only") %>%
  left_join(hindcast_rsq, by = c("rank" = "rank_name")) %>%
  left_join(mean_abun,    by = "model_id") %>%
  filter(!is.na(RMSE.norm), !is.na(mean_abun), !is.na(amplitude), amplitude > 0)

# ── Shared aesthetics ─────────────────────────────────────────────────────────
# kingdom_colors comes from source.R

base_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# ── Panel A: seasonal amplitude ~ mean abundance ──────────────────────────────
pA <- ggplot(to_plot,
             aes(x = mean_abun, y = amplitude, color = pretty_group)) +
  geom_point(alpha = 0.35, size = 2) +
  geom_smooth(method = "glm", formula = y ~ x,
              method.args = list(family = Gamma(link = "log")),
              se = FALSE, linewidth = 1.2) +
  # Data cloud is dense upper-left; place labels in the empty upper-right
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           label.x.npc = 0.45, label.y.npc = 0.97, size = 3.5) +
  scale_x_sqrt() +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(
    x = "Mean observed abundance (sqrt scale)",
    y = "Seasonal amplitude (cycl_only)"
  ) +
  base_theme

# ── Panel B: forecast error ~ seasonal amplitude ──────────────────────────────
pB <- ggplot(to_plot,
             aes(x = amplitude, y = RMSE.norm, color = pretty_group)) +
  geom_point(alpha = 0.35, size = 2) +
  geom_smooth(method = "glm", formula = y ~ x,
              method.args = list(family = Gamma(link = "log")),
              se = FALSE, linewidth = 1.2) +
  # Bacteria (~0.8) and Fungi (~1.7) clusters leave a clear mid gap; place labels there
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           label.x.npc = 0.05, label.y.npc = 0.60, size = 3.5) +
  scale_x_sqrt() +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(
    x = "Seasonal amplitude (sqrt scale)",
    y = "Relative forecast error (nRMSE)"
  ) +
  base_theme

# ── Panel C: forecast error ~ mean abundance ──────────────────────────────────
pC <- ggplot(to_plot,
             aes(x = mean_abun, y = RMSE.norm, color = pretty_group)) +
  geom_point(alpha = 0.35, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 1.2) +
  # Both clouds decline, opening a wide empty band through the center; place labels there
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           label.x.npc = 0.38, label.y.npc = 0.60, size = 3.5) +
  scale_x_sqrt() +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(
    x = "Mean observed abundance (sqrt scale)",
    y = "Relative forecast error (nRMSE)"
  ) +
  base_theme

# ── Combine and save ──────────────────────────────────────────────────────────
fig <- ggarrange(
  pA, pB, pC,
  ncol          = 3,
  labels        = c("A", "B", "C"),
  common.legend = TRUE,
  legend        = "bottom"
)

print(fig)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "f_b_seasonality.png"), fig, width = 13, height = 5, dpi = 200)
cat("Saved: figures/f_b_seasonality.png\n")

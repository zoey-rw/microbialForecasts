# Functional vs taxonomic group forecast comparison
# Central question: do functional or taxonomic groupings produce better forecasts,
# and are environmental drivers more or less informative for each?

library(tidyverse)
library(ggpubr)

source("source.R")

# ── Data loading ──────────────────────────────────────────────────────────────
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged   <- gsub("_beta_regression$", "", scores_list$converged_list)
fg_names    <- microbialForecast:::keep_fg_names

sum.all <- readRDS(here("data", "summary/predictor_effects.rds"))

# ── Hindcast accuracy (env_cycl, observed-site predictions) ──────────────────
hindcast_rsq <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged,
         model_name == "env_cycl",
         grepl("observed", site_prediction)) %>%
  distinct() %>%
  mutate(
    fcast_type = ifelse(species %in% fg_names, "Functional group", "Taxonomic group"),
    RMSE.norm  = pmin(RMSE.norm, 5)
  )

# ── Predictor effect sizes (calibration period, env_cycl) ────────────────────
df_cal <- sum.all %>%
  filter(time_period == "20130601_20180101", model_name == "env_cycl") %>%
  mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id_base %in% converged,
         !beta %in% c("sin", "cos", "rho")) %>%
  mutate(fcast_type = ifelse(taxon %in% fg_names, "Functional group", "Taxonomic group"))

# ── Shared aesthetics (Okabe-Ito color-blind safe) ──────────────────────────
type_colors <- c("Functional group" = "#0072B2", "Taxonomic group" = "#D55E00")
type_shapes <- c("Functional group" = 16, "Taxonomic group" = 17)
comparisons <- list(c("Functional group", "Taxonomic group"))

base_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# ── Panel A: Hindcast nRMSE by forecast type ─────────────────────────────────
pA <- ggplot(hindcast_rsq,
             aes(x = fcast_type, y = RMSE.norm, fill = fcast_type)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, show.legend = FALSE) +
  geom_point(aes(color = fcast_type, shape = fcast_type),
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.3, size = 1.2, show.legend = FALSE) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  facet_wrap(~pretty_group, scales = "fixed") +
  scale_fill_manual(values  = type_colors) +
  scale_color_manual(values = type_colors) +
  scale_shape_manual(values = type_shapes) +
  labs(x = NULL, y = "Relative forecast error (nRMSE)", tag = "A") +
  base_theme

# ── Panel B: Hindcast CRPS by forecast type ──────────────────────────────────
pB <- ggplot(hindcast_rsq,
             aes(x = fcast_type, y = CRPS_truncated, fill = fcast_type)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, show.legend = FALSE) +
  geom_point(aes(color = fcast_type, shape = fcast_type),
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.3, size = 1.2, show.legend = FALSE) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  facet_wrap(~pretty_group, scales = "fixed") +
  scale_fill_manual(values  = type_colors) +
  scale_color_manual(values = type_colors) +
  scale_shape_manual(values = type_shapes) +
  labs(x = NULL, y = "Probabilistic forecast error (CRPS)", tag = "B") +
  base_theme

# ── Panel C: Predictor effect sizes by forecast type (forest plot) ────────────
eff_summary <- df_cal %>%
  group_by(pretty_group, fcast_type, beta) %>%
  summarise(
    mean_eff = mean(effSize, na.rm = TRUE),
    se_eff   = sd(effSize, na.rm = TRUE) / sqrt(sum(!is.na(effSize))),
    .groups  = "drop"
  )

pC <- ggplot(eff_summary,
             aes(x = mean_eff, y = beta, color = fcast_type, shape = fcast_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_linerange(aes(xmin = mean_eff - 1.96 * se_eff,
                     xmax = mean_eff + 1.96 * se_eff),
                 position = position_dodge(width = 0.55),
                 linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) +
  facet_wrap(~pretty_group) +
  scale_color_manual(values = type_colors, name = "Forecast type") +
  scale_shape_manual(values = type_shapes, name = "Forecast type") +
  labs(x = "Mean absolute effect size (\u00b195% CI)", y = NULL, tag = "C") +
  base_theme +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom"
  )

# ── Combine and save ──────────────────────────────────────────────────────────
fig1 <- ggarrange(
  ggarrange(pA, pB, ncol = 2,
            common.legend = TRUE, legend = "bottom"),
  pC,
  nrow = 2,
  heights = c(1, 0.8)
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig1_forecast_map_examples.png"), fig1,
       width = 12, height = 10, dpi = 300)
cat("Saved: figures/fig1_forecast_map_examples.png\n")

# ── Diagnostic: why are nRMSE and CRPS inverted? ─────────────────────────────
# nRMSE = RMSE / mean_abundance; CRPS_truncated uses full predictive distribution.
# Load raw hindcasts to check whether the inversion is driven by:
#   (a) abundance normalization inflating/deflating nRMSE, or
#   (b) forecast spread (sd) differing between functional vs taxonomic groups.
cat("\n── Metric inversion diagnostic ──\n")
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")
if (file.exists(parquet_file) && requireNamespace("nanoparquet", quietly = TRUE)) {
  hindcast_data <- nanoparquet::read_parquet(parquet_file)
} else if (file.exists(parquet_file) && requireNamespace("arrow", quietly = TRUE)) {
  hindcast_data <- arrow::read_parquet(parquet_file)
} else {
  hindcast_data <- readRDS(rds_file)
}
hindcast_diag <- hindcast_data %>%
  as_tibble() %>%
  filter(model_name == "env_cycl",
         grepl("observed", site_prediction),
         !is.na(truth)) %>%
  mutate(fcast_type = ifelse(species %in% fg_names, "Functional group", "Taxonomic group"))

diag_summary <- hindcast_diag %>%
  group_by(fcast_type) %>%
  summarise(
    n_obs            = n(),
    mean_abundance   = mean(truth, na.rm = TRUE),
    median_abundance = median(truth, na.rm = TRUE),
    mean_sd          = mean(sd, na.rm = TRUE),
    median_sd        = median(sd, na.rm = TRUE),
    mean_abs_error   = mean(abs(mean - truth), na.rm = TRUE),
    mean_rel_error   = mean(abs(mean - truth) / pmax(truth, 0.005), na.rm = TRUE),
    .groups = "drop"
  )
cat("Per-observation summary by forecast type:\n")
print(as.data.frame(diag_summary), digits = 4)

# Also break down by pretty_group to see if the pattern holds within kingdoms
diag_by_group <- hindcast_diag %>%
  group_by(pretty_group, fcast_type) %>%
  summarise(
    n_obs          = n(),
    mean_abundance = mean(truth, na.rm = TRUE),
    mean_sd        = mean(sd, na.rm = TRUE),
    mean_abs_error = mean(abs(mean - truth), na.rm = TRUE),
    .groups = "drop"
  )
cat("\nBroken down by kingdom:\n")
print(as.data.frame(diag_by_group), digits = 4)
rm(hindcast_data, hindcast_diag)

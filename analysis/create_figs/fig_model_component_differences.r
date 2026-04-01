# Explore how model predictions differ based on their components
# The R² scatterplot shows similar overall fit, but models may differ in:
# 1. Forecast horizon (how far ahead predictions remain skillful)
# 2. Degradation of accuracy with forecast lead time
# 3. Scoring metrics (CRPS, RMSE) that capture calibration, not just correlation
# 4. Performance for taxa with strong vs weak seasonality

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(data.table)
source("source.R")

# ============================================================
# Load data
# ============================================================

# Forecast horizon results (model-level)
fcast_horizon_in <- readRDS(here("data", "summary/fcast_horizon_df.rds"))
fcast_horizon_results <- fcast_horizon_in[[1]]
fcast_horizon_long <- fcast_horizon_in[[3]]

# Scoring metrics
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged <- scores_list$converged_strict_list
converged_base <- gsub("_beta_regression$", "", converged)

# Site-level scoring data for degradation analysis
fi <- readRDS(here("data/summary/fcast_horizon_input.rds"))
model_mean <- as.data.table(fi[[4]])  # per model x site x months_since_obs
null_site <- as.data.table(fi[[3]])
# Strip null_ prefix and filter to site_mean null type for backward compat
if ("null_type" %in% names(null_site)) null_site <- null_site[null_type == "site_mean"]
nc <- grep("^null_", names(null_site), value = TRUE)
if (length(nc) > 0) setnames(null_site, nc, gsub("^null_", "", nc))

# Phenology data for seasonality stratification
pheno_data <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))[[1]] %>%
  mutate(model_id = gsub("_beta_regression$", "", model_id))

# Assign functional vs taxonomic
fg_names <- microbialForecast:::keep_fg_names

# ── Shared aesthetics ────────────────────────────────────────
model_colors <- c(cycl_only = "#2166AC", env_cov = "#B2182B", env_cycl = "#4DAF4A")
kingdom_colors <- c(Bacteria = "#F8766D", Fungi = "#00BFC4")

base_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )


# ============================================================
# 1. FORECAST HORIZON BY MODEL TYPE
#    The key question: do env_cycl models stay skillful longer?
# ============================================================

hz_data <- fcast_horizon_results %>%
  filter(model_id %in% converged_base) %>%
  mutate(
    forecast_horizon = coalesce(rsq_fcast_horizon, crps_fcast_horizon, rmse_fcast_horizon),
    fcast_type = ifelse(species %in% fg_names, "Functional groups", "Taxonomic groups")
  ) %>%
  filter(is.finite(forecast_horizon) & forecast_horizon > 0 & !is.na(pretty_group))

cat("Forecast horizon data:", nrow(hz_data), "model fits\n")
cat("By model type:\n")
print(hz_data %>% group_by(model_name) %>% summarize(
  n = n(), median_hz = median(forecast_horizon), mean_hz = mean(forecast_horizon)))

# Panel 1a: Paired comparison - for each taxon, compare horizons across models
hz_wide <- hz_data %>%
  select(species, model_name, pretty_group, forecast_horizon) %>%
  group_by(species, model_name, pretty_group) %>%
  summarize(forecast_horizon = mean(forecast_horizon, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = model_name, values_from = forecast_horizon)

panel_1a <- ggplot(hz_data,
                   aes(x = model_name, y = forecast_horizon,
                       fill = model_name, color = model_name)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(size = 1.5, alpha = 0.3,
             position = position_jitter(width = 0.15, height = 0)) +
  facet_wrap(~ pretty_group) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model.labs) +
  stat_compare_means(comparisons = list(
    c("cycl_only", "env_cycl"), c("env_cov", "env_cycl"), c("cycl_only", "env_cov")),
    method = "wilcox.test", label = "p.signif", hide.ns = FALSE,
    step.increase = 0.08) +
  labs(y = "Forecast horizon (months)", x = NULL) +
  base_theme + theme(legend.position = "none")

# Panel 1b: Same but split by functional vs taxonomic
panel_1b <- ggplot(hz_data,
                   aes(x = model_name, y = forecast_horizon,
                       fill = model_name, color = model_name)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(size = 1.2, alpha = 0.3,
             position = position_jitter(width = 0.15, height = 0)) +
  facet_grid(fcast_type ~ pretty_group) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model.labs) +
  labs(y = "Forecast horizon (months)", x = NULL) +
  base_theme + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 20, hjust = 1))


# ============================================================
# 2. ACCURACY DEGRADATION OVER LEAD TIME
#    Do models diverge as forecast lead time increases?
# ============================================================

# Get per-model-type degradation curves
degrade_data <- model_mean[model_name != "all_covariates" &
                             !is.na(pretty_group) &
                             is.finite(RSQ.1)]

# Mean across taxa, per model type x kingdom x lead time
degrade_summary <- degrade_data[,
  .(mean_rsq = mean(RSQ.1, na.rm = TRUE),
    mean_rmse = mean(RMSE, na.rm = TRUE),
    mean_crps = mean(mean_crps, na.rm = TRUE),
    n_taxa = .N),
  by = .(model_name, pretty_group, months_since_obs)]

# Also compute null baseline
null_rsq_overall <- null_site[, .(null_RSQ = mean(RSQ.1, na.rm = TRUE)),
                               by = .(model_name, pretty_group, months_since_obs)]

degrade_merged <- merge(degrade_summary, null_rsq_overall,
                        by = c("model_name", "pretty_group", "months_since_obs"),
                        all.x = TRUE)

panel_2a <- ggplot(degrade_summary,
                   aes(x = months_since_obs, y = mean_rsq,
                       color = model_name, linetype = model_name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_wrap(~ pretty_group) +
  scale_color_manual(values = model_colors, labels = model.labs, name = "Model") +
  scale_linetype_manual(values = c(cycl_only = "dashed", env_cov = "dotted", env_cycl = "solid"),
                        labels = model.labs, name = "Model") +
  labs(x = "Months since last observation",
       y = expression(Mean~R^2~(1:1))) +
  coord_cartesian(ylim = c(0, NA)) +
  base_theme

panel_2b <- ggplot(degrade_summary,
                   aes(x = months_since_obs, y = mean_crps,
                       color = model_name, linetype = model_name)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_wrap(~ pretty_group) +
  scale_color_manual(values = model_colors, labels = model.labs, name = "Model") +
  scale_linetype_manual(values = c(cycl_only = "dashed", env_cov = "dotted", env_cycl = "solid"),
                        labels = model.labs, name = "Model") +
  labs(x = "Months since last observation",
       y = "Mean CRPS") +
  base_theme


# ============================================================
# 3. SCORING METRICS THAT CAPTURE CALIBRATION
#    R² only captures correlation; CRPS and RMSE capture bias + spread
# ============================================================

# Debug: check available values
cat("\nscoring_metrics_long site_prediction values:\n")
print(unique(scores_list$scoring_metrics_long$site_prediction))
cat("scoring_metrics_long model_name values:\n")
print(unique(scores_list$scoring_metrics_long$model_name))
cat("N converged IDs matching scoring_metrics_long:\n")
cat(sum(converged %in% scores_list$scoring_metrics_long$model_id), "of", length(converged), "\n")
cat("pretty_group values:\n")
print(unique(scores_list$scoring_metrics_long$pretty_group))

# Use greedy filter: pick the site_prediction value that contains "observed" or take all
avail_site_pred <- unique(scores_list$scoring_metrics_long$site_prediction)
obs_site_pred <- grep("bserved|new time", avail_site_pred, value = TRUE, ignore.case = TRUE)
if (length(obs_site_pred) == 0) obs_site_pred <- avail_site_pred

scoring_long <- scores_list$scoring_metrics_long %>%
  filter(model_name != "all_covariates" &
           site_prediction %in% obs_site_pred &
           !is.na(pretty_group))

# If converged filter removes everything, skip it
if (sum(scoring_long$model_id %in% converged) > 10) {
  scoring_long <- scoring_long %>% filter(model_id %in% converged)
}
cat("scoring_long rows after filter:", nrow(scoring_long), "\n")
cat("metrics available:", paste(unique(scoring_long$metric), collapse = ", "), "\n")

# CRPS (calibration-sensitive) by model type
# Try CRPS_truncated first, fall back to CRPS
crps_metric <- if ("CRPS_truncated" %in% scoring_long$metric) "CRPS_truncated" else "CRPS"
crps_data <- scoring_long %>% filter(metric == crps_metric)
cat("CRPS data rows:", nrow(crps_data), "using metric:", crps_metric, "\n")

# RMSE (normalized) by model type - try RMSE.norm, fall back to RMSE
rmse_metric <- if ("RMSE.norm" %in% scoring_long$metric) "RMSE.norm" else "RMSE"
rmse_data <- scoring_long %>% filter(metric == rmse_metric)
cat("RMSE data rows:", nrow(rmse_data), "using metric:", rmse_metric, "\n")

panel_3a <- panel_3b <- NULL

if (nrow(crps_data) > 0 && length(unique(crps_data$pretty_group)) > 0) {
  panel_3a <- ggplot(crps_data,
                     aes(x = model_name, y = score,
                         fill = model_name, color = model_name)) +
    geom_violin(alpha = 0.3, quantiles = 0.5, show.legend = FALSE) +
    geom_point(size = 1.5, alpha = 0.3,
               position = position_jitter(width = 0.15, height = 0),
               show.legend = FALSE) +
    facet_wrap(~ pretty_group, scales = "free_y") +
    scale_fill_manual(values = model_colors) +
    scale_color_manual(values = model_colors) +
    scale_x_discrete(labels = model.labs) +
    stat_compare_means(comparisons = list(
      c("cycl_only", "env_cycl"), c("env_cov", "env_cycl")),
      method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
    labs(y = "CRPS (lower = better)", x = NULL, title = "Hindcast CRPS by model type") +
    base_theme + theme(legend.position = "none")
}

if (nrow(rmse_data) > 0 && length(unique(rmse_data$pretty_group)) > 0) {
  panel_3b <- ggplot(rmse_data,
                     aes(x = model_name, y = score,
                         fill = model_name, color = model_name)) +
    geom_violin(alpha = 0.3, quantiles = 0.5, show.legend = FALSE) +
    geom_point(size = 1.5, alpha = 0.3,
               position = position_jitter(width = 0.15, height = 0),
               show.legend = FALSE) +
    facet_wrap(~ pretty_group, scales = "free_y") +
    scale_fill_manual(values = model_colors) +
    scale_color_manual(values = model_colors) +
    scale_x_discrete(labels = model.labs) +
    stat_compare_means(comparisons = list(
      c("cycl_only", "env_cycl"), c("env_cov", "env_cycl")),
      method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
    labs(y = "Normalized RMSE (lower = better)", x = NULL, title = "Hindcast RMSE by model type") +
    base_theme + theme(legend.position = "none")
}


# ============================================================
# 4. PERFORMANCE STRATIFIED BY SEASONAL AMPLITUDE
#    env_cycl should shine for taxa with strong seasonality
# ============================================================

# Get amplitude per taxon from the cycl_only model (baseline seasonal signal)
cat("\npheno_data columns:", paste(names(pheno_data), collapse = ", "), "\n")
cat("pheno_data taxon examples:", paste(head(unique(pheno_data$taxon)), collapse = ", "), "\n")
cat("scoring_long species examples:", paste(head(unique(scoring_long$species)), collapse = ", "), "\n")

amplitude_by_taxon <- pheno_data %>%
  filter(model_name == "cycl_only") %>%
  select(taxon, amplitude) %>%
  distinct() %>%
  filter(!is.na(amplitude)) %>%
  mutate(amp_tercile = cut(amplitude,
    quantile(amplitude, c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"), include.lowest = TRUE))

cat("Amplitude data:", nrow(amplitude_by_taxon), "taxa\n")

# Merge amplitude with scoring metrics - try both species and model_id-based join
panel_4 <- NULL
scoring_with_amp <- scoring_long %>%
  filter(metric == crps_metric) %>%
  left_join(amplitude_by_taxon, by = c("species" = "taxon"))

# If no matches via species, try extracting taxon from model_id
if (all(is.na(scoring_with_amp$amp_tercile))) {
  cat("species-taxon join failed, trying model_id-based join\n")
  amp_by_model <- pheno_data %>%
    filter(!is.na(amplitude)) %>%
    mutate(model_id = gsub("_beta_regression$", "", model_id)) %>%
    select(model_id, amplitude) %>%
    distinct() %>%
    mutate(amp_tercile = cut(amplitude,
      quantile(amplitude, c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("Low", "Medium", "High"), include.lowest = TRUE))
  scoring_with_amp <- scoring_long %>%
    filter(metric == crps_metric) %>%
    left_join(amp_by_model, by = "model_id")
}

scoring_with_amp <- scoring_with_amp %>% filter(!is.na(amp_tercile))
cat("Scoring x amplitude rows:", nrow(scoring_with_amp), "\n")

if (nrow(scoring_with_amp) > 0 &&
    length(unique(scoring_with_amp$pretty_group)) > 0 &&
    length(unique(scoring_with_amp$amp_tercile)) > 1) {
  panel_4 <- ggplot(scoring_with_amp,
                    aes(x = model_name, y = score,
                        fill = model_name, color = model_name)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.15, height = 0),
               show.legend = FALSE) +
    facet_grid(amp_tercile ~ pretty_group,
               labeller = labeller(amp_tercile = c(
                 Low = "Low seasonality", Medium = "Med. seasonality", High = "High seasonality"))) +
    scale_fill_manual(values = model_colors, labels = model.labs, name = "Model") +
    scale_color_manual(values = model_colors, labels = model.labs, name = "Model") +
    scale_x_discrete(labels = model.labs) +
    labs(y = "CRPS (lower = better)", x = NULL,
         title = "Does model advantage depend on seasonal amplitude?") +
    base_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}


# ============================================================
# 5. PAIRED DIFFERENCE IN HORIZON: env_cycl - cycl_only per taxon
#    Shows whether adding environment extends forecasts
# ============================================================

hz_paired <- hz_data %>%
  select(species, model_name, pretty_group, forecast_horizon) %>%
  filter(model_name %in% c("cycl_only", "env_cycl")) %>%
  group_by(species, model_name, pretty_group) %>%
  summarize(forecast_horizon = mean(forecast_horizon, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = model_name, values_from = forecast_horizon) %>%
  filter(!is.na(cycl_only) & !is.na(env_cycl)) %>%
  mutate(horizon_diff = env_cycl - cycl_only)

cat("\nPaired horizon difference (env_cycl - cycl_only):\n")
print(hz_paired %>% group_by(pretty_group) %>% summarize(
  n = n(),
  mean_diff = mean(horizon_diff),
  median_diff = median(horizon_diff),
  pct_env_wins = mean(horizon_diff > 0) * 100,
  wilcox_p = wilcox.test(horizon_diff)$p.value))

panel_5 <- ggplot(hz_paired,
                  aes(x = cycl_only, y = env_cycl, color = pretty_group)) +
  geom_point(size = 2.5, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(x = "Forecast horizon: Seasonality only (months)",
       y = "Forecast horizon: Env. + Seasonality (months)",
       title = "Paired comparison of forecast horizons") +
  coord_equal() +
  base_theme


# ============================================================
# Compose figures
# ============================================================

# Figure 1: Forecast horizons by model type
fig1 <- ggarrange(
  panel_1a, panel_5,
  labels = c("A", "B"), nrow = 1, widths = c(1.2, 1))

ggsave(here("figures", "model_component_horizons.png"), fig1,
       width = 14, height = 6, dpi = 200)
cat("Saved: figures/model_component_horizons.png\n")

# Figure 2: Degradation over time
fig2 <- ggarrange(panel_2a, panel_2b,
                  labels = c("A", "B"), nrow = 2, common.legend = TRUE)

ggsave(here("figures", "model_component_degradation.png"), fig2,
       width = 10, height = 8, dpi = 200)
cat("Saved: figures/model_component_degradation.png\n")

# Figure 3: Scoring metrics + seasonality interaction
scoring_panels <- Filter(Negate(is.null), list(panel_3a, panel_3b))
if (length(scoring_panels) > 0) {
  scoring_row <- ggarrange(plotlist = scoring_panels,
                           labels = LETTERS[seq_along(scoring_panels)], nrow = 1)
  if (!is.null(panel_4)) {
    fig3 <- ggarrange(scoring_row, panel_4,
                      labels = c(NA, LETTERS[length(scoring_panels) + 1]),
                      nrow = 2, heights = c(0.8, 1))
  } else {
    fig3 <- scoring_row
  }
  ggsave(here("figures", "model_component_scoring.png"), fig3,
         width = 12, height = if (!is.null(panel_4)) 10 else 5, dpi = 200)
  cat("Saved: figures/model_component_scoring.png\n")
}

# Combined overview - only include non-NULL panels
overview_panels <- Filter(Negate(is.null), list(
  panel_1a,
  if (!is.null(panel_2a) && !is.null(panel_2b))
    ggarrange(panel_2a, panel_2b, nrow = 1, common.legend = TRUE, legend = "top"),
  if (!is.null(panel_3a))
    ggarrange(plotlist = Filter(Negate(is.null), list(panel_3a, panel_5)),
              nrow = 1, widths = c(1, 1)),
  panel_4
))

if (length(overview_panels) > 1) {
  fig_all <- ggarrange(plotlist = overview_panels,
                       nrow = length(overview_panels),
                       heights = rep(1, length(overview_panels)))
  ggsave(here("figures", "model_component_differences_overview.png"), fig_all,
         width = 14, height = 5 * length(overview_panels), dpi = 200)
  cat("Saved: figures/model_component_differences_overview.png\n")
}

# ============================================================
# Summary statistics
# ============================================================
cat("\n======= SUMMARY =======\n")

cat("\n-- Forecast horizon by model type --\n")
hz_data %>% group_by(model_name, pretty_group) %>%
  summarize(n = n(), median = median(forecast_horizon),
            mean = round(mean(forecast_horizon), 1),
            sd = round(sd(forecast_horizon), 1), .groups = "drop") %>%
  print(n = 20)

if (nrow(crps_data) > 0) {
  cat("\n-- CRPS by model type --\n")
  crps_data %>% group_by(model_name, pretty_group) %>%
    summarize(n = n(), median = round(median(score, na.rm = TRUE), 4),
              mean = round(mean(score, na.rm = TRUE), 4), .groups = "drop") %>%
    print(n = 20)
}

if (nrow(scoring_with_amp) > 0) {
  cat("\n-- CRPS by model type x seasonality tercile --\n")
  scoring_with_amp %>% group_by(model_name, amp_tercile, pretty_group) %>%
    summarize(n = n(), mean_crps = round(mean(score, na.rm = TRUE), 4), .groups = "drop") %>%
    print(n = 30)
}


# ============================================================
# 6. CALIBRATION VS HINDCAST: where do models actually differ?
#    In-sample fit should favor env_cov/env_cycl (more parameters),
#    but out-of-sample may penalize overfitting.
# ============================================================

cat("\n\n====== CALIBRATION VS HINDCAST ANALYSIS ======\n")

# Taxon-level calibration metrics
cal_metrics <- scores_list$calibration_metrics %>%
  filter(model_name != "all_covariates" & !is.na(pretty_group))

# Taxon-level hindcast metrics
hind_metrics <- scores_list$scoring_metrics %>%
  filter(model_name != "all_covariates" &
           site_prediction == "New time (observed site)" &
           !is.na(pretty_group))

# Combine calibration and hindcast for direct comparison
cal_long <- cal_metrics %>%
  select(model_id, model_name, pretty_group, species, RSQ, CRPS = CRPS_truncated) %>%
  mutate(period = "Calibration")
hind_long <- hind_metrics %>%
  select(model_id, model_name, pretty_group, species, RSQ, CRPS = CRPS_truncated) %>%
  mutate(period = "Hindcast")
both_periods <- bind_rows(cal_long, hind_long) %>%
  mutate(period = factor(period, levels = c("Calibration", "Hindcast")))

cat("\n-- RSQ: Calibration vs Hindcast --\n")
both_periods %>% group_by(model_name, pretty_group, period) %>%
  summarize(n = n(), med_RSQ = round(median(RSQ, na.rm = TRUE), 3),
            med_CRPS = round(median(CRPS, na.rm = TRUE), 4), .groups = "drop") %>%
  print(n = 20)

# Panel 6a: RSQ drop from calibration to hindcast
panel_6a <- ggplot(both_periods,
                   aes(x = model_name, y = RSQ,
                       fill = period, color = period)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.7,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ pretty_group) +
  scale_x_discrete(labels = model.labs) +
  scale_fill_manual(values = c(Calibration = "grey60", Hindcast = "steelblue"), name = NULL) +
  scale_color_manual(values = c(Calibration = "grey40", Hindcast = "steelblue4"), name = NULL) +
  labs(y = expression(R^2), x = NULL,
       title = "In-sample vs out-of-sample fit") +
  base_theme

# Panel 6b: Paired RSQ drop per taxon - which model loses the most going out-of-sample?
rsq_drop <- both_periods %>%
  select(model_id, model_name, pretty_group, species, RSQ, period) %>%
  group_by(model_id, model_name, pretty_group, species, period) %>%
  summarize(RSQ = mean(RSQ, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = RSQ) %>%
  filter(!is.na(Calibration) & !is.na(Hindcast)) %>%
  mutate(rsq_drop = Calibration - Hindcast,
         rsq_retention = Hindcast / Calibration)

cat("\n-- RSQ drop (Calibration - Hindcast) by model type --\n")
rsq_drop %>% group_by(model_name, pretty_group) %>%
  summarize(n = n(),
            med_cal = round(median(Calibration, na.rm = TRUE), 3),
            med_hind = round(median(Hindcast, na.rm = TRUE), 3),
            med_drop = round(median(rsq_drop, na.rm = TRUE), 3),
            med_retention = round(median(rsq_retention, na.rm = TRUE), 2),
            .groups = "drop") %>%
  print(n = 20)

panel_6b <- ggplot(rsq_drop,
                   aes(x = model_name, y = rsq_drop,
                       fill = model_name, color = model_name)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(size = 1.5, alpha = 0.3,
             position = position_jitter(width = 0.15, height = 0),
             show.legend = FALSE) +
  facet_wrap(~ pretty_group) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model.labs) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_compare_means(comparisons = list(
    c("cycl_only", "env_cycl"), c("env_cov", "env_cycl"), c("cycl_only", "env_cov")),
    method = "wilcox.test", label = "p.signif", hide.ns = FALSE,
    step.increase = 0.08) +
  labs(y = expression(R^2~drop~(calibration - hindcast)),
       x = NULL,
       title = "Overfitting: how much accuracy is lost out-of-sample?") +
  base_theme + theme(legend.position = "none")


# ============================================================
# 7. SITE-LEVEL: observed vs new sites
#    Models differ most when predicting at new sites
# ============================================================

site_scoring <- scores_list$scoring_metrics_site_long %>%
  filter(model_name != "all_covariates" &
           metric == "CRPS_truncated" &
           !is.na(pretty_group))

# Simplify site_prediction labels
site_scoring <- site_scoring %>%
  mutate(site_type = case_when(
    site_prediction == "New time (observed site)" ~ "Observed site",
    site_prediction == "New time x site (random effect)" ~ "New site\n(random effect)",
    site_prediction == "New time x site (modeled effect)" ~ "New site\n(PLSR predicted)",
    TRUE ~ site_prediction))

cat("\n-- Site-level CRPS by model type x site type --\n")
site_scoring %>% group_by(model_name, pretty_group, site_type) %>%
  summarize(n = n(), med_crps = round(median(score, na.rm = TRUE), 4), .groups = "drop") %>%
  print(n = 30)

panel_7a <- ggplot(site_scoring,
                   aes(x = model_name, y = score,
                       fill = model_name, color = model_name)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  facet_grid(site_type ~ pretty_group, scales = "free_y") +
  scale_fill_manual(values = model_colors, labels = model.labs, name = "Model") +
  scale_color_manual(values = model_colors, labels = model.labs, name = "Model") +
  scale_x_discrete(labels = model.labs) +
  scale_y_log10() +
  labs(y = "Site-level CRPS (log scale, lower = better)", x = NULL,
       title = "Model performance by site type") +
  base_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                     legend.position = "none")

# Panel 7b: Difference in CRPS between site types, per model
# How much worse is each model at new sites vs observed?
site_wide <- site_scoring %>%
  filter(site_type %in% c("Observed site", "New site\n(random effect)")) %>%
  mutate(site_cat = ifelse(grepl("Observed", site_type), "observed", "new")) %>%
  group_by(model_name, pretty_group, species, site_cat) %>%
  summarize(med_crps = median(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site_cat, values_from = med_crps) %>%
  filter(!is.na(observed) & !is.na(new)) %>%
  mutate(crps_ratio = new / observed)

cat("\n-- Transferability: new-site / observed-site CRPS ratio --\n")
site_wide %>% group_by(model_name, pretty_group) %>%
  summarize(n = n(),
            med_ratio = round(median(crps_ratio, na.rm = TRUE), 2),
            mean_ratio = round(mean(crps_ratio, na.rm = TRUE), 2),
            .groups = "drop") %>%
  print(n = 20)

panel_7b <- ggplot(site_wide,
                   aes(x = model_name, y = crps_ratio,
                       fill = model_name, color = model_name)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(size = 1.5, alpha = 0.3,
             position = position_jitter(width = 0.15, height = 0),
             show.legend = FALSE) +
  facet_wrap(~ pretty_group) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  scale_x_discrete(labels = model.labs) +
  stat_compare_means(comparisons = list(
    c("cycl_only", "env_cycl"), c("env_cov", "env_cycl"), c("cycl_only", "env_cov")),
    method = "wilcox.test", label = "p.signif", hide.ns = FALSE,
    step.increase = 0.08) +
  labs(y = "CRPS ratio (new site / observed site)",
       x = NULL,
       title = "Transferability penalty by model type") +
  base_theme + theme(legend.position = "none")


# ============================================================
# Compose calibration + site figures
# ============================================================

fig_cal <- ggarrange(panel_6a, panel_6b,
                     labels = c("A", "B"), nrow = 1, widths = c(1, 1))
ggsave(here("figures", "model_component_calibration_vs_hindcast.png"), fig_cal,
       width = 14, height = 6, dpi = 200)
cat("\nSaved: figures/model_component_calibration_vs_hindcast.png\n")

fig_site <- ggarrange(panel_7a, panel_7b,
                      labels = c("A", "B"), nrow = 1, widths = c(1.2, 1))
ggsave(here("figures", "model_component_site_transferability.png"), fig_site,
       width = 14, height = 8, dpi = 200)
cat("Saved: figures/model_component_site_transferability.png\n")

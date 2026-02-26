# Forecast horizon by seasonality, kingdom, forecast type, and latitude
# Uses pre-computed phenology from 09_assignPeakPhenophase_optimized.r
# and site-level scoring data for latitude analysis

library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
source("source.R")

# ============================================================
# Load data
# ============================================================

# Model-level forecast horizons
fcast_horizon_df <- readRDS(here("data", "summary/fcast_horizon_df.rds"))[[1]]
fcast_horizon_df$forecast_horizon <- fcast_horizon_df$rsq_fcast_horizon
fcast_horizon_df$forecast_horizon <- ifelse(
  is.na(fcast_horizon_df$forecast_horizon) | !is.finite(fcast_horizon_df$forecast_horizon),
  ifelse(is.finite(fcast_horizon_df$crps_fcast_horizon) & fcast_horizon_df$crps_fcast_horizon > 0,
         fcast_horizon_df$crps_fcast_horizon,
         fcast_horizon_df$rmse_fcast_horizon),
  fcast_horizon_df$forecast_horizon)

# Assign functional vs taxonomic
fg_names <- microbialForecast:::keep_fg_names
fcast_horizon_df$fcast_type <- ifelse(fcast_horizon_df$species %in% fg_names,
                                       "Functional Groups", "Taxonomic Groups")

# Phenology amplitude + seasonality significance flags
pheno_sig <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))[[1]] %>%
  mutate(model_id = gsub("_beta_regression$", "", model_id)) %>%
  select(model_id, amplitude, any_of(c("significant_sin", "significant_cos"))) %>%
  distinct()

# Merge and filter to env_cycl
plot_data <- merge(fcast_horizon_df, pheno_sig, by = "model_id", all.x = TRUE) %>%
  filter(model_name == "env_cycl" & is.finite(forecast_horizon) & forecast_horizon > 0)

if ("pretty_group.y" %in% names(plot_data)) {
  plot_data$pretty_group <- ifelse(is.na(plot_data$pretty_group.y),
                                   plot_data$pretty_group.x,
                                   plot_data$pretty_group.y)
}

cat("Model-level data:", nrow(plot_data), "env_cycl models\n")

# ── Shared aesthetics ──────────────────────────────────────────────────────────
kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

# ============================================================
# Site-level horizons for latitude analysis
# ============================================================

fi <- readRDS(here("data/summary/fcast_horizon_input.rds"))
model_mean <- as.data.table(fi[[4]])
null_site <- as.data.table(fi[[3]])

# Null RSQ per model
null_rsq <- null_site[model_name == "env_cycl",
                       .(null_RSQ = mean(RSQ.1, na.rm = TRUE)), by = model_id]

# Per-site scoring data
site_data <- model_mean[model_name == "env_cycl"]
site_data <- merge(site_data, null_rsq, by = "model_id", all.x = TRUE)

# Site-level horizon: max months_since_obs where model beats null
site_hz <- site_data[RSQ.1 > null_RSQ & is.finite(RSQ.1),
                      .(site_horizon = max(months_since_obs, na.rm = TRUE)),
                      by = .(model_id, species, siteID, pretty_group)]

# Include model x site combos that never beat null (horizon = 0)
all_combos <- unique(site_data[, .(model_id, species, siteID, pretty_group)])
site_hz <- merge(all_combos, site_hz,
                 by = c("model_id", "species", "siteID", "pretty_group"), all.x = TRUE)
site_hz[is.na(site_horizon), site_horizon := 0]

# Merge latitude
lat_df <- readRDS(here("data/clean/site_effect_predictors.rds"))[, c("siteID", "latitude")]
site_hz <- merge(site_hz, lat_df, by = "siteID", all.x = TRUE)

# Aggregate: mean horizon per site x kingdom
site_kingdom <- site_hz[, .(mean_horizon = mean(site_horizon, na.rm = TRUE),
                             n_models = .N),
                          by = .(siteID, latitude, pretty_group)]

# Aggregate: mean horizon per site (all taxa)
site_avg <- site_hz[, .(mean_horizon = mean(site_horizon, na.rm = TRUE),
                          n_models = .N,
                          prop_skilled = mean(site_horizon > 0)),
                      by = .(siteID, latitude)]

cat("Site-level data:", nrow(site_hz), "model x site combos across",
    length(unique(site_hz$siteID)), "sites\n")

# ============================================================
# Panel A: Amplitude terciles by kingdom
# ============================================================
plot_data_amp <- plot_data %>% filter(!is.na(amplitude))
plot_data_amp$seasonality <- cut(
  plot_data_amp$amplitude,
  quantile(plot_data_amp$amplitude, c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE)

panel_a <- ggplot(plot_data_amp,
                  aes(x = seasonality, y = forecast_horizon, color = pretty_group)) +
  geom_boxplot(aes(fill = pretty_group), alpha = 0.3, outlier.shape = NA) +
  geom_point(size = 1.5, alpha = 0.3,
             position = position_jitter(height = 0, width = 0.15)) +
  facet_wrap(~ pretty_group) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  labs(y = "Forecast horizon (months)", x = "Seasonal amplitude (terciles)") +
  base_theme

# ============================================================
# Panel B: Continuous amplitude scatter per kingdom
# ============================================================
panel_b <- ggplot(plot_data_amp,
                  aes(x = amplitude, y = forecast_horizon, color = pretty_group)) +
  geom_point(size = 1.8, alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = kingdom_colors, name = "Kingdom") +
  labs(y = "Forecast horizon (months)", x = "Seasonal amplitude") +
  base_theme + theme(legend.position = "top")

tryCatch({
  cor_bac <- with(plot_data_amp %>% filter(pretty_group == "Bacteria"),
                  cor.test(amplitude, forecast_horizon))
  cor_fun <- with(plot_data_amp %>% filter(pretty_group == "Fungi"),
                  cor.test(amplitude, forecast_horizon))
  cor_label <- paste0("Bacteria: r = ", round(cor_bac$estimate, 2),
                      ", p = ", round(cor_bac$p.value, 3),
                      "\nFungi: r = ", round(cor_fun$estimate, 2),
                      ", p = ", round(cor_fun$p.value, 3))
  panel_b <- panel_b +
    annotate("text", x = max(plot_data_amp$amplitude) * 0.6,
             y = min(plot_data_amp$forecast_horizon) + 0.5,
             label = cor_label, size = 4, hjust = 0)
}, error = function(e) cat("Correlation annotation failed:", e$message, "\n"))

# ============================================================
# Panel C: Kingdom comparison (horizon horizontal)
# ============================================================
panel_c <- ggplot(plot_data,
                  aes(x = forecast_horizon, y = pretty_group,
                      color = pretty_group, fill = pretty_group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.4) +
  geom_point(size = 1.5, alpha = 0.2,
             position = position_jitter(height = 0.1, width = 0)) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  labs(x = "Forecast horizon (months)", y = NULL) +
  base_theme

tryCatch({
  t_res <- plot_data %>% rstatix::t_test(forecast_horizon ~ pretty_group, detailed = TRUE)
  p_label <- ifelse(t_res$p < 0.001, "p < 0.001", paste0("p = ", round(t_res$p, 3)))
  panel_c <- panel_c + labs(subtitle = p_label)
}, error = function(e) cat("Kingdom t-test failed:", e$message, "\n"))

# ============================================================
# Panel D: Functional vs Taxonomic (horizon horizontal)
# ============================================================
panel_d <- ggplot(plot_data,
                  aes(x = forecast_horizon, y = fcast_type, color = pretty_group)) +
  geom_boxplot(aes(fill = pretty_group), alpha = 0.3, outlier.shape = NA, width = 0.4) +
  geom_point(size = 1.5, alpha = 0.2,
             position = position_jitterdodge(jitter.height = 0.05, dodge.width = 0)) +
  scale_color_manual(values = kingdom_colors, name = "Kingdom") +
  scale_fill_manual(values  = kingdom_colors, guide = "none") +
  labs(x = "Forecast horizon (months)", y = NULL) +
  base_theme + theme(legend.position = "top")

tryCatch({
  t_res2 <- plot_data %>% rstatix::t_test(forecast_horizon ~ fcast_type, detailed = TRUE)
  p_label2 <- ifelse(t_res2$p < 0.001, "p < 0.001", paste0("p = ", round(t_res2$p, 3)))
  panel_d <- panel_d + labs(subtitle = p_label2)
}, error = function(e) cat("Fcast type t-test failed:", e$message, "\n"))

# ============================================================
# Panel E: Site-level horizon vs latitude (per kingdom)
# ============================================================
panel_e <- ggplot(site_kingdom,
                  aes(x = latitude, y = mean_horizon, color = pretty_group)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = kingdom_colors, name = "Kingdom") +
  labs(x = "Latitude (\u00b0N)", y = "Mean site-level\nforecast horizon (months)") +
  base_theme + theme(legend.position = "top")

tryCatch({
  cor_bac_lat <- with(site_kingdom[pretty_group == "Bacteria"],
                      cor.test(latitude, mean_horizon))
  cor_fun_lat <- with(site_kingdom[pretty_group == "Fungi"],
                      cor.test(latitude, mean_horizon))
  lat_label <- paste0("Bacteria: r = ", round(cor_bac_lat$estimate, 2),
                       ", p = ", round(cor_bac_lat$p.value, 3),
                       "\nFungi: r = ", round(cor_fun_lat$estimate, 2),
                       ", p = ", round(cor_fun_lat$p.value, 3))
  panel_e <- panel_e +
    annotate("text", x = 55, y = max(site_kingdom$mean_horizon) * 0.9,
             label = lat_label, size = 4, hjust = 0)
}, error = function(e) cat("Latitude correlation failed:", e$message, "\n"))

# ============================================================
# Panel F: Proportion of skilled models per site vs latitude
# ============================================================
panel_f <- ggplot(site_avg,
                  aes(x = latitude, y = prop_skilled)) +
  geom_point(size = 2.5, alpha = 0.6, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "grey30") +
  ggrepel::geom_text_repel(aes(label = siteID), size = 2.5, alpha = 0.7,
                           max.overlaps = 10) +
  labs(x = "Latitude (\u00b0N)", y = "Proportion of taxa\nwith skilled forecasts") +
  ylim(c(0, 1)) +
  base_theme

tryCatch({
  cor_prop <- cor.test(site_avg$latitude, site_avg$prop_skilled)
  prop_label <- paste0("r = ", round(cor_prop$estimate, 2),
                        ", p = ", round(cor_prop$p.value, 3))
  panel_f <- panel_f +
    annotate("text", x = 55, y = 0.9, label = prop_label, size = 4, hjust = 0)
}, error = function(e) cat("Proportion correlation failed:", e$message, "\n"))

# ============================================================
# Panel G: Horizon by seasonality significance (all model types)
# Uses significant_sin / significant_cos from phenology data.
# Compares forecast horizon for taxa with vs without significant seasonality.
# ============================================================
has_sig_flags <- all(c("significant_sin", "significant_cos") %in% names(pheno_sig))

if (has_sig_flags) {
  pheno_sig_flag <- pheno_sig %>%
    mutate(seasonal = (significant_sin == 1 | significant_cos == 1),
           seasonal_label = ifelse(seasonal, "Significant\nseasonality", "No significant\nseasonality")) %>%
    select(model_id, seasonal_label)

  horizon_seas <- fcast_horizon_df %>%
    left_join(pheno_sig_flag, by = "model_id") %>%
    filter(is.finite(forecast_horizon), forecast_horizon > 0,
           !is.na(seasonal_label), !is.na(pretty_group))

  cat("Horizon-by-seasonality data:", nrow(horizon_seas), "rows\n")

  comparisons_g <- list(c("No significant\nseasonality", "Significant\nseasonality"))

  panel_g <- ggplot(horizon_seas,
                    aes(x = forecast_horizon, y = seasonal_label,
                        fill = pretty_group, color = pretty_group)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, width = 0.5,
                 position = position_dodge(width = 0.7)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                               jitter.height = 0,
                                               dodge.width  = 0.7),
               alpha = 0.25, size = 1) +
    facet_wrap(~ model_name,
               labeller = labeller(model_name = c(cycl_only = "Cycl. only",
                                                   env_cov   = "Env. only",
                                                   env_cycl  = "Env. + Cycl."))) +
    scale_color_manual(values = kingdom_colors, name = "Kingdom") +
    scale_fill_manual(values  = kingdom_colors, name = "Kingdom") +
    labs(x = "Forecast horizon (months)", y = NULL) +
    base_theme + theme(legend.position = "top")
} else {
  cat("significant_sin/cos columns not found in phenology data; skipping Panel G\n")
  panel_g <- NULL
}

# ============================================================
# Compose and save
# ============================================================
fig_top <- ggarrange(panel_a, panel_b, labels = c("A", "B"), widths = c(1, 1))
fig_mid <- ggarrange(panel_c, panel_d, labels = c("C", "D"), widths = c(0.8, 1.2))
fig_bot <- ggarrange(panel_e, panel_f, labels = c("E", "F"), widths = c(1, 1))

if (!is.null(panel_g)) {
  fig_all <- ggarrange(fig_top, fig_mid, fig_bot,
                       ggarrange(panel_g, labels = "G"),
                       nrow = 4, heights = c(1, 0.7, 1, 0.7))
} else {
  fig_all <- ggarrange(fig_top, fig_mid, fig_bot, nrow = 3, heights = c(1, 0.7, 1))
}

out_dir <- here("data", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_horizon_by_seasonality.pdf"), fig_all,
       width = 12, height = if (!is.null(panel_g)) 16 else 12)
ggsave(file.path(out_dir, "fig_horizon_by_seasonality.png"), fig_all,
       width = 12, height = if (!is.null(panel_g)) 16 else 12, dpi = 200)

cat("Saved: data/figures/fig_horizon_by_seasonality.pdf / .png\n")

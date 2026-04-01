# Spatial heterogeneity vs seasonal amplitude
# Shows the relationship between spatial autocorrelation (Moran's I),
# seasonal amplitude, and forecast accuracy across taxa.
# Panel A: Moran's I vs seasonal amplitude scatter with density margins
# Panel B: Multi-predictor AICc model selection (dredge) for forecast accuracy

source("source.R")
library(cowplot)
library(MuMIn)
options(na.action = "na.fail")

# ── Colorblind-friendly palette ─────────────────────────────────────────────
kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")

# ── Data loading ─────────────────────────────────────────────────────────────
scores_list    <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_strict <- scores_list$converged_strict_list

moran_df <- readRDS("data/clean/moran_stat.rds")

rho_core_in <- readRDS(here("data/summary/rho_core_sd_effects.rds")) %>%
  filter(model_name != "all_covariates",
         model_id %in% converged_strict) %>%
  select(-pretty_name) %>%
  mutate(model_id = gsub("_beta_regression$", "", model_id))

# Seasonal amplitude
seasonal_amplitude_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_only_vals <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "cycl_only") %>%
  rename(cycl_amplitude = amplitude, cycl_max = max) %>%
  select(taxon, time_period, cycl_amplitude, cycl_max) %>%
  distinct(taxon, time_period, .keep_all = TRUE)

env_cycl_vals <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "env_cycl") %>%
  rename(env_cycl_amplitude = amplitude, env_cycl_max = max) %>%
  select(taxon, time_period, env_cycl_amplitude, env_cycl_max) %>%
  distinct(taxon, time_period, .keep_all = TRUE)

# ── Merge scoring metrics with rho/precision ────────────────────────────────
# Use observed-site predictions if available, otherwise all
scoring_metrics <- scores_list$scoring_metrics
if ("site_prediction" %in% colnames(scoring_metrics)) {
  obs_rows <- grepl("observed", scoring_metrics$site_prediction, ignore.case = TRUE)
  if (sum(obs_rows) > 0) scoring_metrics <- scoring_metrics[obs_rows, ]
}

rho_core_unique <- rho_core_in %>% distinct(model_id, rowname, .keep_all = TRUE)
model_df <- merge(scoring_metrics, rho_core_unique, by = "model_id", allow.cartesian = TRUE)

# Clean duplicate columns from merge
if ("pretty_group.x" %in% colnames(model_df)) {
  model_df$pretty_group <- coalesce(model_df$pretty_group.x, model_df$pretty_group.y)
  model_df$pretty_group.x <- NULL
  model_df$pretty_group.y <- NULL
}
if ("model_name.x" %in% colnames(model_df)) {
  model_df$model_name <- model_df$model_name.x
  model_df$model_name.x <- NULL
  model_df$model_name.y <- NULL
}

# Merge Moran's I
if ("species" %in% colnames(model_df)) {
  model_df <- merge(model_df, moran_df, by.x = "species", by.y = "taxon", all.x = TRUE)
} else {
  model_df <- merge(model_df, moran_df, by = "taxon", all.x = TRUE)
}

# Merge seasonal amplitude
model_df <- merge(model_df, cycl_only_vals, by = c("taxon", "time_period"), all.x = TRUE)
model_df <- merge(model_df, env_cycl_vals, by = c("taxon", "time_period"), all.x = TRUE)

# Pivot rho/precision to wide
if ("rowname" %in% colnames(model_df) && "Mean" %in% colnames(model_df)) {
  cols_to_drop <- intersect(c("SD", "Naive SE", "Time-series SE"), colnames(model_df))
  if (length(cols_to_drop) > 0) model_df <- model_df %>% select(-all_of(cols_to_drop))
  model_df <- model_df %>%
    pivot_wider(names_from = "rowname", values_from = "Mean")
}

# ── Panel A: Moran's I vs seasonal amplitude ────────────────────────────────
plot_data <- model_df %>%
  filter(is.finite(mean_morans), is.finite(cycl_amplitude),
         !is.na(pretty_group))

panel_a <- ggplot(plot_data,
                  aes(x = mean_morans, y = cycl_amplitude, color = pretty_group)) +
  geom_jitter(alpha = 0.3, size = 2.5, height = 0.01, width = 0.01) +
  geom_smooth(method = "lm", linewidth = 1.5, alpha = 0.2, se = FALSE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           size = 4) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  scale_x_sqrt() +
  labs(x = "Spatial autocorrelation (Moran's I)",
       y = "Seasonal amplitude (cycl_only)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

# ── Panel B: Core-level variability vs seasonal amplitude ───────────────────
if ("core_sd" %in% colnames(model_df) || "precision" %in% colnames(model_df)) {
  spatial_var <- if ("core_sd" %in% colnames(model_df)) "core_sd" else "precision"

  plot_data_b <- model_df %>%
    filter(is.finite(.data[[spatial_var]]), is.finite(cycl_amplitude),
           !is.na(pretty_group))

  panel_b <- ggplot(plot_data_b,
                    aes(x = .data[[spatial_var]], y = cycl_amplitude,
                        color = pretty_group)) +
    geom_jitter(alpha = 0.3, size = 2.5, height = 0.01, width = 0.01) +
    geom_smooth(method = "lm", linewidth = 1.5, alpha = 0.2, se = FALSE) +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
             size = 4) +
    scale_color_manual(values = kingdom_colors, name = NULL) +
    labs(x = paste0("Core-level variation (", spatial_var, ")"),
         y = "Seasonal amplitude (cycl_only)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank())
} else {
  panel_b <- NULL
}

# ── Combine and save ────────────────────────────────────────────────────────
if (!is.null(panel_b)) {
  fig <- plot_grid(panel_a, panel_b, labels = c("A", "B"),
                   ncol = 2, align = "h", label_size = 16)
} else {
  fig <- panel_a
}

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS6_moran_vs_rsq.png"), fig,
       width = 14, height = 7, dpi = 200)

cat("Saved: figures/figS6_moran_vs_rsq.png\n")

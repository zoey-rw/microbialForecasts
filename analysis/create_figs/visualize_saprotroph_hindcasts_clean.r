# Clean Visualization of Saprophytic Fungi Hindcasts with Uncertainty
# This script creates focused visualizations of saprophytic fungi hindcasts
# with proper uncertainty bands and clear model comparisons

cat("=== CLEAN VISUALIZATION OF SAPROTROPHIC FUNGI HINDCASTS ===\n\n")

# Load required packages and environment
source("../../source.R")

cat("✅ Environment loaded successfully\n")

# Load saprotroph hindcast data
if (!file.exists(here("data/summary/hindcast_saprotroph.rds"))) {
  stop("hindcast_saprotroph.rds not found. Please ensure the hindcast data is available.")
}

saprotroph_data <- readRDS(here("data/summary/hindcast_saprotroph.rds"))

cat("📊 Loaded saprotroph hindcast data:", nrow(saprotroph_data), "rows\n")

# Clean and prepare data for visualization
hindcast_clean <- saprotroph_data %>%
  # Remove rows with no predictions
  filter(!is.na(med)) %>%
  # Add date information
  mutate(
    dates = as.Date(dates),
    month = lubridate::month(dates),
    month_label = lubridate::month(dates, label = TRUE),
    year = lubridate::year(dates),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring", 
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    )
  ) %>%
  # Add forecast period information
  mutate(
    fcast_period = case_when(
      dates <= as.Date("2018-01-01") ~ "calibration",
      TRUE ~ "hindcast"
    )
  )

cat("🧹 Cleaned data:", nrow(hindcast_clean), "rows with predictions\n")

# Select representative sites for visualization (top 4 by data coverage)
site_coverage <- hindcast_clean %>%
  group_by(siteID) %>%
  summarize(
    n_predictions = n(),
    n_truth = sum(!is.na(truth)),
    date_range = max(dates, na.rm = TRUE) - min(dates, na.rm = TRUE)
  ) %>%
  arrange(desc(n_predictions))

cat("🏞️ Top sites by data coverage:\n")
print(head(site_coverage, 6))

# Select top 4 sites for cleaner visualization
top_sites <- head(site_coverage$siteID, 4)
cat("🎯 Selected sites for visualization:", paste(top_sites, collapse = ", "), "\n")

# Filter data for selected sites
plot_data <- hindcast_clean %>%
  filter(siteID %in% top_sites) %>%
  # Select plots with good coverage
  group_by(plotID) %>%
  filter(n() >= 50) %>%
  ungroup()

cat("📊 Final plotting dataset:", nrow(plot_data), "rows\n")

# Create focused visualization
cat("🎨 Creating focused saprotroph hindcast visualization...\n")

# 1. Main hindcast visualization with uncertainty bands
main_plot <- ggplot(plot_data, aes(x = dates, y = med, group = plotID)) +
  # Facet by site and model type
  facet_grid(
    rows = vars(siteID), 
    cols = vars(model_name), 
    scales = "free_y",
    labeller = labeller(
      model_name = c(
        "cycl_only" = "Seasonal Only",
        "env_cov" = "Environmental Only", 
        "env_cycl" = "Environmental + Seasonal",
        "all_covariates" = "All Covariates"
      )
    )
  ) +
  # Uncertainty bands for calibration period
  geom_ribbon(
    data = ~filter(.x, fcast_period == "calibration"),
    aes(ymin = lo, ymax = hi), 
    alpha = 0.3, 
    fill = "steelblue"
  ) +
  geom_ribbon(
    data = ~filter(.x, fcast_period == "calibration"),
    aes(ymin = lo_25, ymax = hi_75), 
    alpha = 0.6, 
    fill = "steelblue"
  ) +
  # Uncertainty bands for hindcast period
  geom_ribbon(
    data = ~filter(.x, fcast_period == "hindcast"),
    aes(ymin = lo, ymax = hi), 
    alpha = 0.2, 
    fill = "coral"
  ) +
  geom_ribbon(
    data = ~filter(.x, fcast_period == "hindcast"),
    aes(ymin = lo_25, ymax = hi_75), 
    alpha = 0.4, 
    fill = "coral"
  ) +
  # Prediction lines
  geom_line(
    data = ~filter(.x, fcast_period == "calibration"),
    aes(x = dates, y = med), 
    alpha = 0.8,
    color = "steelblue",
    linewidth = 0.8
  ) +
  geom_line(
    data = ~filter(.x, fcast_period == "hindcast"),
    aes(x = dates, y = med), 
    alpha = 0.6,
    color = "coral",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  # Observed data points (only where truth exists)
  geom_point(
    data = ~filter(.x, !is.na(truth)),
    aes(y = as.numeric(truth)), 
    alpha = 0.8,
    size = 2,
    color = "darkgreen"
  ) +
  # Add vertical line to separate calibration and hindcast
  geom_vline(
    xintercept = as.Date("2018-01-01"), 
    linetype = "dashed", 
    color = "darkgray",
    size = 1
  ) +
  # Theme and labels
  theme_bw(base_size = 12) +
  theme(
    panel.spacing = unit(0.3, "cm"),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  labs(
    title = "Saprophytic Fungi Hindcasts with Prediction Uncertainty",
    subtitle = "Steel blue: Calibration period, Coral: Hindcast period\nShaded bands show 50% and 95% prediction intervals",
    x = "Date",
    y = "Relative Abundance",
    caption = "Vertical line indicates end of calibration period (2018-01-01)"
  )

# 2. Create seasonal comparison plot
cat("🍃 Creating seasonal comparison visualization...\n")

seasonal_plot <- ggplot(plot_data, aes(x = month_label, y = med, fill = season)) +
  facet_grid(
    rows = vars(siteID), 
    cols = vars(model_name),
    scales = "free_y"
  ) +
  # Box plots for seasonal patterns
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  # Add individual points for truth values (only where they exist)
  geom_jitter(
    data = ~filter(.x, !is.na(truth)),
    aes(y = as.numeric(truth)), 
    alpha = 0.6, 
    size = 1.5,
    width = 0.2,
    color = "darkgreen"
  ) +
  # Seasonal color scheme
  scale_fill_manual(
    values = c(
      "Winter" = "#E8F4FD", 
      "Spring" = "#E8F8E8", 
      "Summer" = "#FFF8E8", 
      "Fall" = "#FDE8E8"
    )
  ) +
  # Theme and labels
  theme_bw(base_size = 12) +
  theme(
    panel.spacing = unit(0.3, "cm"),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Seasonal Patterns in Saprophytic Fungi Abundance",
    subtitle = "Box plots show predicted values, green points show observed values",
    x = "Month",
    y = "Relative Abundance",
    fill = "Season"
  )

# 3. Create model comparison plot
cat("📈 Creating model comparison visualization...\n")

# Calculate model performance metrics
model_performance <- plot_data %>%
  filter(!is.na(truth)) %>%
  group_by(model_name, siteID) %>%
  summarize(
    rmse = sqrt(mean((med - as.numeric(truth))^2, na.rm = TRUE)),
    bias = mean(med - as.numeric(truth), na.rm = TRUE),
    coverage_95 = mean(as.numeric(truth) >= lo & as.numeric(truth) <= hi, na.rm = TRUE),
    coverage_50 = mean(as.numeric(truth) >= lo_25 & as.numeric(truth) <= hi_75, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

cat("📊 Model performance summary:\n")
print(model_performance)

# Model comparison plot
model_comparison <- ggplot(model_performance, aes(x = model_name, y = rmse, fill = model_name)) +
  facet_wrap(~siteID, scales = "free_y", ncol = 2) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(
    aes(label = sprintf("%.3f", rmse)), 
    vjust = -0.5, 
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(
    panel.spacing = unit(0.3, "cm"),
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model Performance Comparison by Site",
    subtitle = "Lower RMSE indicates better predictions",
    x = "Model Type",
    y = "RMSE",
    fill = "Model Type"
  )

# 4. Create uncertainty analysis plot
cat("🎯 Creating uncertainty analysis visualization...\n")

# Calculate uncertainty metrics
uncertainty_metrics <- plot_data %>%
  filter(!is.na(med)) %>%
  mutate(
    uncertainty_95 = hi - lo,
    uncertainty_50 = hi_75 - lo_25,
    relative_uncertainty_95 = uncertainty_95 / med,
    relative_uncertainty_50 = uncertainty_50 / med
  ) %>%
  group_by(model_name, siteID, fcast_period) %>%
  summarize(
    mean_uncertainty_95 = mean(uncertainty_95, na.rm = TRUE),
    mean_uncertainty_50 = mean(uncertainty_50, na.rm = TRUE),
    mean_relative_uncertainty_95 = mean(relative_uncertainty_95, na.rm = TRUE),
    mean_relative_uncertainty_50 = mean(relative_uncertainty_50, na.rm = TRUE),
    .groups = "drop"
  )

cat("📊 Uncertainty metrics summary:\n")
print(uncertainty_metrics)

# Uncertainty comparison plot
uncertainty_plot <- ggplot(uncertainty_metrics, aes(x = model_name, y = mean_relative_uncertainty_95, fill = fcast_period)) +
  facet_wrap(~siteID, scales = "free_y", ncol = 2) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(
    values = c("calibration" = "steelblue", "hindcast" = "coral"),
    labels = c("calibration" = "Calibration", "hindcast" = "Hindcast")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.spacing = unit(0.3, "cm"),
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Prediction Uncertainty by Model Type and Period",
    subtitle = "Relative uncertainty (95% interval width / predicted value)",
    x = "Model Type",
    y = "Relative Uncertainty",
    fill = "Period"
  )

# Save all plots
cat("💾 Saving visualizations...\n")

# Create output directory
output_dir <- here("figures", "saprotroph_hindcasts_clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save main plot
ggsave(
  filename = file.path(output_dir, "saprotroph_hindcasts_main_clean.png"),
  plot = main_plot,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Save seasonal plot
ggsave(
  filename = file.path(output_dir, "saprotroph_hindcasts_seasonal_clean.png"),
  plot = seasonal_plot,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Save model comparison plot
ggsave(
  filename = file.path(output_dir, "saprotroph_hindcasts_model_comparison_clean.png"),
  plot = model_comparison,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Save uncertainty plot
ggsave(
  filename = file.path(output_dir, "saprotroph_hindcasts_uncertainty_clean.png"),
  plot = uncertainty_plot,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("✅ All clean visualizations saved to:", output_dir, "\n")

# Summary statistics
cat("\n📊 VISUALIZATION SUMMARY\n")
cat("========================\n")
cat("Main hindcast plot: Shows predictions with uncertainty bands\n")
cat("Seasonal plot: Displays seasonal patterns across sites and models\n")
cat("Model comparison: Compares performance across different model types\n")
cat("Uncertainty analysis: Analyzes prediction uncertainty patterns\n")
cat("Sites visualized:", length(unique(plot_data$siteID)), "\n")
cat("Model types:", length(unique(plot_data$model_name)), "\n")
cat("Total predictions:", nrow(plot_data), "\n")
cat("Observed values:", sum(!is.na(plot_data$truth)), "\n")

# Model performance summary
cat("\n🏆 MODEL PERFORMANCE SUMMARY\n")
cat("============================\n")
overall_performance <- model_performance %>%
  group_by(model_name) %>%
  summarize(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_coverage_95 = mean(coverage_95, na.rm = TRUE),
    mean_coverage_50 = mean(coverage_50, na.rm = TRUE),
    n_sites = n()
  )

print(overall_performance)

cat("\n=== CLEAN SAPROTROPH HINDCAST VISUALIZATION COMPLETE ===\n")

#!/usr/bin/env Rscript
# Comprehensive calibration analysis for all taxonomic and functional group models
source("source.R")
library(ggplot2)
library(dplyr)
library(tidyr)

cat("=== COMPREHENSIVE HINDCAST CALIBRATION ANALYSIS ===\n")

# Function to analyze calibration for a dataset
analyze_calibration <- function(data, group_name) {
  if (nrow(data) == 0) return(NULL)
  
  # Calculate mean value continuity
  mean_continuity <- data %>%
    group_by(plotID, fcast_period) %>%
    summarise(mean_abundance = mean(mean, na.rm=TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = fcast_period, values_from = mean_abundance) %>%
    mutate(
      percent_change = (hindcast - calibration) / calibration * 100,
      abs_percent_change = abs(percent_change)
    ) %>%
    summarise(
      mean_abs_percent_change = mean(abs_percent_change, na.rm=TRUE),
      max_abs_percent_change = max(abs_percent_change, na.rm=TRUE),
      .groups = 'drop'
    )
  
  # Calculate uncertainty calibration
  ci_stats <- data %>%
    filter(!is.na(mean) & !is.na(lo) & !is.na(hi)) %>%
    mutate(ci_width = hi - lo) %>%
    group_by(fcast_period) %>%
    summarise(mean_ci_width = mean(ci_width, na.rm=TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = fcast_period, values_from = mean_ci_width) %>%
    mutate(ci_width_ratio = hindcast / calibration)
  
  # Combine results
  result <- cbind(
    group = group_name,
    mean_continuity,
    ci_stats
  )
  
  return(result)
}

# Analyze taxonomic models (from all_hindcasts_raw.rds)
cat("Loading taxonomic hindcast data...\n")
taxonomic_data <- readRDS(here("data/summary/all_hindcasts_raw.rds"))

taxonomic_results <- taxonomic_data %>%
  group_by(species, model_name) %>%
  group_modify(~ analyze_calibration(.x, paste(.y$species, .y$model_name, sep="_"))) %>%
  ungroup() %>%
  mutate(group_type = "Taxonomic")

# Analyze functional group models
cat("Loading functional group hindcast data...\n")
fg_files <- list.files(here("data/summary"), pattern = "hindcasts_.*_observed\\.rds", full.names = TRUE)

# Limit to first 20 functional groups for manageable analysis
fg_files <- fg_files[1:min(20, length(fg_files))]

functional_results <- data.frame()

for (file in fg_files) {
  group_name <- gsub(".*hindcasts_(.*)_observed\\.rds", "\\1", basename(file))
  cat("  Analyzing", group_name, "...\n")
  
  tryCatch({
    fg_data <- readRDS(file)
    if (nrow(fg_data) > 0) {
      result <- analyze_calibration(fg_data, group_name)
      if (!is.null(result)) {
        result$group_type <- "Functional"
        functional_results <- rbind(functional_results, result)
      }
    }
  }, error = function(e) {
    cat("    Error analyzing", group_name, ":", e$message, "\n")
  })
}

# Combine all results
all_results <- rbind(taxonomic_results, functional_results)

# Calculate composite scores
all_results <- all_results %>%
  mutate(
    # Normalize metrics (lower is better)
    mean_continuity_norm = mean_abs_percent_change / max(mean_abs_percent_change, na.rm=TRUE),
    uncertainty_norm = ci_width_ratio / max(ci_width_ratio, na.rm=TRUE),
    # Combined score (equal weight)
    combined_score = (mean_continuity_norm + uncertainty_norm) / 2,
    # Rank (lower is better)
    combined_rank = rank(combined_score)
  ) %>%
  arrange(combined_rank)

cat("\n=== TOP 10 BEST-CALIBRATED MODELS (Overall) ===\n")
top_models <- head(all_results, 10)
print(top_models[, c('group', 'group_type', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score', 'combined_rank')])

cat("\n=== TOP 5 BEST-CALIBRATED TAXONOMIC MODELS ===\n")
top_taxonomic <- all_results %>%
  filter(group_type == "Taxonomic") %>%
  head(5)
print(top_taxonomic[, c('group', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score')])

cat("\n=== TOP 5 BEST-CALIBRATED FUNCTIONAL MODELS ===\n")
top_functional <- all_results %>%
  filter(group_type == "Functional") %>%
  head(5)
print(top_functional[, c('group', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score')])

# Create visualization
cat("\n=== CREATING CALIBRATION COMPARISON PLOT ===\n")

# Select best models for visualization
best_models <- all_results %>%
  filter(combined_rank <= 4) %>%
  pull(group)

# Create comparison plot data
plot_data <- data.frame()

# Add taxonomic data
for (model in best_models[best_models %in% taxonomic_results$group]) {
  species_name <- strsplit(model, "_")[[1]][1]
  model_type <- strsplit(model, "_")[[1]][2]
  
  model_data <- taxonomic_data %>%
    filter(species == species_name, model_name == model_type) %>%
    filter(plotID %in% c("BART_001", "BART_002")) %>%
    mutate(
      group_label = model,
      group_type = "Taxonomic"
    )
  
  plot_data <- rbind(plot_data, model_data)
}

# Add functional data
for (model in best_models[best_models %in% functional_results$group]) {
  tryCatch({
    fg_data <- readRDS(here("data/summary", paste0("hindcasts_", model, "_observed.rds")))
    if (nrow(fg_data) > 0) {
      model_data <- fg_data %>%
        filter(plotID %in% c("BART_001", "BART_002")) %>%
        mutate(
          group_label = model,
          group_type = "Functional"
        )
      plot_data <- rbind(plot_data, model_data)
    }
  }, error = function(e) {
    cat("Error loading", model, ":", e$message, "\n")
  })
}

# Create the comparison plot
if (nrow(plot_data) > 0) {
  comparison_plot <- ggplot(plot_data, aes(x = dates, y = mean, group = plotID)) +
    facet_grid(rows = vars(group_label),
               cols = vars(plotID),
               drop = T,
               scales = "free") +
    # Uncertainty bands for calibration period
    geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
                aes(x = dates, ymin = lo, ymax = hi), alpha = 0.2, fill = "steelblue") +
    geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
                aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = 0.4, fill = "steelblue") +
    # Uncertainty bands for hindcast period
    geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
                aes(x = dates, ymin = lo, ymax = hi), alpha = 0.1, fill = "coral") +
    geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
                aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = 0.2, fill = "coral") +
    # Prediction lines
    geom_line(data = ~filter(.x, fcast_period == "calibration"),
              aes(x = dates, y = mean), alpha = 0.8, color = "steelblue", linewidth = 0.8) +
    geom_line(data = ~filter(.x, fcast_period == "hindcast"),
              aes(x = dates, y = mean), alpha = 0.6, color = "coral", linewidth = 0.8, linetype = "dashed") +
    # Add vertical line to separate calibration and hindcast periods
    geom_vline(xintercept = as.Date("2018-01-01"), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
    xlab("Date") + 
    ylab("Relative Abundance") +
    theme_bw(base_size = 12) +
    theme(panel.spacing = unit(0.2, "cm"),
          legend.position = "none",
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
    labs(title = "Best-Calibrated Hindcast Models Comparison",
         subtitle = "Steel blue: Calibration period, Coral: Hindcast period",
         caption = "Vertical line indicates end of calibration period (2018-01-01)")
  
  # Save the plot
  png(here("figures", "comprehensive_calibration_comparison.png"), width = 1400, height = 1000)
  print(comparison_plot)
  dev.off()
  
  cat("Plot saved to figures/comprehensive_calibration_comparison.png\n")
}

# Save results
write.csv(all_results, here("figures", "comprehensive_calibration_results.csv"), row.names = FALSE)
cat("Results saved to figures/comprehensive_calibration_results.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

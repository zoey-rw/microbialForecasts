# Create Seasonal Trend Plots from Best Converged Models
# This script generates actual seasonal abundance curves using the converged parameters

cat("=== CREATING SEASONAL TREND PLOTS ===\n\n")

# Load required packages and environment
source("source.R")

cat("✅ Environment loaded successfully\n")

# Load the best converged models data
best_models_file <- here("figures/phenology_comparison/best_converged_phenology_models.csv")
best_models <- read.csv(best_models_file)

cat("📊 Loaded", nrow(best_models), "best converged models\n")

# Function to generate seasonal trend for a single model
generate_seasonal_trend <- function(sin_coef, cos_coef, months = 1:12) {
  # Generate seasonal trend using sin/cos coefficients
  # Note: sin/cos parameters are already scaled in the model fitting
  seasonal_effect <- sin_coef * sin(2 * pi * months / 12) + cos_coef * cos(2 * pi * months / 12)
  return(seasonal_effect)
}

# Function to create seasonal abundance plot
create_seasonal_plot <- function(model_data, months = 1:12) {
  # Generate seasonal trends for each model
  seasonal_trends <- list()
  
  for (i in 1:nrow(model_data)) {
    row <- model_data[i, ]
    trend <- generate_seasonal_trend(row$sin_coef, row$cos_coef, months)
    
    seasonal_trends[[i]] <- data.frame(
      month = months,
      month_name = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[months],
      abundance = trend,
      species = row$species,
      model_type = row$model_type,
      amplitude = row$amplitude,
      peak_month = row$peak_month,
      convergence_score = row$convergence_score
    )
  }
  
  # Combine all trends
  all_trends <- do.call(rbind, seasonal_trends)
  
  # Create the plot
  library(ggplot2)
  library(viridis)
  library(dplyr)
  
  p <- ggplot(all_trends, aes(x = month, y = abundance, color = species, 
                              linetype = model_type, group = species)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis(discrete = TRUE) +
    scale_linetype_manual(values = c("Logit" = "solid", "CLR" = "dashed")) +
    scale_x_continuous(breaks = 1:12, 
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
    labs(title = "Seasonal Abundance Trends from Best Converged Models",
         subtitle = paste("Based on", nrow(model_data), "models with convergence_score < 10"),
         x = "Month",
         y = "Seasonal Effect (logit scale)",
         color = "Species",
         linetype = "Model Type") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE),
           linetype = guide_legend(nrow = 1))
  
  return(p)
}

# Function to create amplitude comparison plot
create_amplitude_plot <- function(model_data) {
  library(ggplot2)
  library(viridis)
  
  # Add month labels for better interpretation
  model_data$month_label <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[round(model_data$peak_month)]
  
  p <- ggplot(model_data, aes(x = reorder(species, -amplitude), y = amplitude, 
                              fill = model_type, color = model_type)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.7) +
    geom_text(aes(label = month_label, y = amplitude + max(amplitude) * 0.05), 
              size = 3, angle = 45, hjust = 0) +
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE) +
    labs(title = "Seasonal Amplitude Comparison - Best Converged Models",
         subtitle = "Peak month shown above each bar",
         x = "Species",
         y = "Seasonal Amplitude",
         fill = "Model Type",
         color = "Model Type") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

# Function to create peak timing distribution plot
create_timing_plot <- function(model_data) {
  library(ggplot2)
  library(viridis)
  
  # Create month labels
  model_data$month_label <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[round(model_data$peak_month)]
  
  p <- ggplot(model_data, aes(x = peak_month, y = amplitude, 
                              color = model_type, size = amplitude)) +
    geom_point(alpha = 0.8) +
    geom_text(aes(label = species), vjust = -0.8, size = 3, 
              color = "black", fontface = "bold") +
    scale_color_viridis(discrete = TRUE) +
    scale_size_continuous(range = c(3, 8)) +
    scale_x_continuous(breaks = 1:12, 
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
    labs(title = "Seasonal Peak Timing vs Amplitude - Best Converged Models",
         subtitle = "Point size indicates amplitude strength",
         x = "Peak Month",
         y = "Seasonal Amplitude",
         color = "Model Type",
         size = "Amplitude") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

# Main execution
cat("🎨 Creating seasonal trend visualizations...\n")

# 1. Seasonal abundance trends over time
cat("   📈 Creating seasonal abundance trends...\n")
seasonal_plot <- create_seasonal_plot(best_models)

# 2. Amplitude comparison
cat("   📊 Creating amplitude comparison...\n")
amplitude_plot <- create_amplitude_plot(best_models)

# 3. Peak timing distribution
cat("   🕐 Creating peak timing distribution...\n")
timing_plot <- create_timing_plot(best_models)

# Save all plots
output_dir <- here("figures/phenology_comparison")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("   💾 Saving plots...\n")

ggsave(file.path(output_dir, "seasonal_abundance_trends.png"), seasonal_plot, 
       width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "seasonal_amplitude_comparison.png"), amplitude_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "seasonal_peak_timing.png"), timing_plot, 
       width = 10, height = 8, dpi = 300)

cat("   ✅ All plots saved successfully!\n")

# Summary of best converged models
cat("\n📋 SUMMARY OF BEST CONVERGED MODELS\n")
cat("=", paste(rep("=", 50), collapse = ""), "\n")

for (i in 1:nrow(best_models)) {
  row <- best_models[i, ]
  month_name <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[round(row$peak_month)]
  
  cat(sprintf("%d. %s (%s)\n", i, row$species, row$model_type))
  cat("   Amplitude: %.3f\n", row$amplitude)
  cat("   Peak Month: %s (%.1f)\n", month_name, row$peak_month)
  cat("   Convergence Score: %.2f (lower = better)\n", row$convergence_score)
  cat("   Effective Sample Size: %.0f\n", row$mean_ess)
  cat("\n")
}

# Seasonal pattern analysis
cat("🌍 SEASONAL PATTERN ANALYSIS\n")
cat("-", paste(rep("-", 40), collapse = ""), "\n")

# Calculate seasonal statistics
mean_amplitude <- mean(best_models$amplitude)
median_amplitude <- median(best_models$amplitude)
mean_peak_timing <- mean(best_models$peak_timing)

cat("Mean seasonal amplitude:", round(mean_amplitude, 3), "\n")
cat("Median seasonal amplitude:", round(median_amplitude, 3), "\n")
cat("Mean peak timing:", round(mean_peak_timing, 1), "months\n")

# Seasonal distribution
peak_months <- round(best_models$peak_timing)
month_counts <- table(peak_months)
cat("\nSeasonal distribution:\n")
for (month in sort(unique(peak_months))) {
  month_name <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")[month]
  count <- month_counts[as.character(month)]
  cat("  %s (%d): %d models\n", month_name, month, count)
}

# Model type comparison
logit_models <- best_models[best_models$model_type == "Logit", ]
if (nrow(logit_models) > 0) {
  cat("\nLogit Beta (Fungal Functional Groups) Summary:\n")
  cat("  Count:", nrow(logit_models), "\n")
  cat("  Mean amplitude:", round(mean(logit_models$amplitude), 3), "\n")
  cat("  Mean peak timing:", round(mean(logit_models$peak_timing), 1), "months\n")
}

clr_models <- best_models[best_models$model_type == "CLR", ]
if (nrow(clr_models) > 0) {
  cat("\nCLR (Bacterial Phyla) Summary:\n")
  cat("  Count:", nrow(clr_models), "\n")
  cat("  Mean amplitude:", round(mean(clr_models$amplitude), 3), "\n")
  cat("  Mean peak timing:", round(mean(clr_models$peak_timing), 1), "months\n")
}

cat("\n=== SEASONAL TREND PLOTS COMPLETE ===\n")
cat("Plots saved to:", output_dir, "\n")

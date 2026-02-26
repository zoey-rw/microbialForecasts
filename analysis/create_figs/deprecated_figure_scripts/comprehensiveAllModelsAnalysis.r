#!/usr/bin/env Rscript
# Comprehensive analysis of ALL available hindcast models (taxonomic + functional)
source("source.R")
library(ggplot2)
library(dplyr)
library(tidyr)

cat("=== COMPREHENSIVE ANALYSIS OF ALL HINDCAST MODELS ===\n")

# Function to convert list columns to numeric
convert_list_columns <- function(data) {
  # Convert list columns to numeric
  list_cols <- c("lo", "lo_25", "med", "hi_75", "hi", "mean")
  
  for (col in list_cols) {
    if (col %in% names(data)) {
      data[[col]] <- sapply(data[[col]], function(x) {
        if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)
      })
    }
  }
  
  return(data)
}

# Function to analyze calibration for a dataset
analyze_calibration <- function(data, group_name) {
  if (nrow(data) == 0) return(NULL)
  
  # Convert list columns
  data <- convert_list_columns(data)
  
  # Filter out rows with no predictions
  data <- data %>% filter(!is.na(mean) & !is.na(lo) & !is.na(hi))
  
  if (nrow(data) == 0) return(NULL)
  
  # Calculate mean value continuity
  mean_continuity <- data %>%
    group_by(plotID, fcast_period) %>%
    summarise(mean_abundance = mean(mean, na.rm=TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = fcast_period, values_from = mean_abundance) %>%
    filter(!is.na(calibration) & !is.na(hindcast)) %>%
    mutate(
      percent_change = (hindcast - calibration) / calibration * 100,
      abs_percent_change = abs(percent_change)
    ) %>%
    summarise(
      mean_abs_percent_change = mean(abs_percent_change, na.rm=TRUE),
      max_abs_percent_change = max(abs_percent_change, na.rm=TRUE),
      n_plots = n(),
      .groups = 'drop'
    )
  
  if (nrow(mean_continuity) == 0) return(NULL)
  
  # Calculate uncertainty calibration
  ci_stats <- data %>%
    group_by(fcast_period) %>%
    summarise(mean_ci_width = mean(hi - lo, na.rm=TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = fcast_period, values_from = mean_ci_width) %>%
    filter(!is.na(calibration) & !is.na(hindcast)) %>%
    mutate(ci_width_ratio = hindcast / calibration)
  
  if (nrow(ci_stats) == 0) return(NULL)
  
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

cat("Taxonomic models analyzed:", nrow(taxonomic_results), "\n")

# Analyze functional group models
cat("Loading functional group hindcast data...\n")
fg_files <- list.files(here("data/summary"), pattern = "hindcasts_.*_observed\\.rds", full.names = TRUE)

cat("Found", length(fg_files), "functional group files\n")

functional_results <- data.frame()
successful_analyses <- 0

for (i in seq_along(fg_files)) {
  file <- fg_files[i]
  group_name <- gsub(".*hindcasts_(.*)_observed\\.rds", "\\1", basename(file))
  
  if (i %% 50 == 0) cat("  Processing file", i, "of", length(fg_files), ":", group_name, "\n")
  
  tryCatch({
    fg_data <- readRDS(file)
    if (nrow(fg_data) > 0) {
      result <- analyze_calibration(fg_data, group_name)
      if (!is.null(result)) {
        result$group_type <- "Functional"
        functional_results <- rbind(functional_results, result)
        successful_analyses <- successful_analyses + 1
      }
    }
  }, error = function(e) {
    # Silently skip problematic files
  })
}

cat("Successfully analyzed", successful_analyses, "functional group models\n")

# Combine all results
all_results <- rbind(taxonomic_results, functional_results)

cat("Total models analyzed:", nrow(all_results), "\n")
cat("  Taxonomic:", sum(all_results$group_type == "Taxonomic"), "\n")
cat("  Functional:", sum(all_results$group_type == "Functional"), "\n")

# Calculate composite scores
all_results <- all_results %>%
  filter(!is.na(mean_abs_percent_change) & !is.na(ci_width_ratio)) %>%
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

cat("\n=== TOP 20 BEST-CALIBRATED MODELS (Overall) ===\n")
top_models <- head(all_results, 20)
print(top_models[, c('group', 'group_type', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score', 'combined_rank')])

cat("\n=== TOP 10 BEST-CALIBRATED TAXONOMIC MODELS ===\n")
top_taxonomic <- all_results %>%
  filter(group_type == "Taxonomic") %>%
  head(10)
print(top_taxonomic[, c('group', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score')])

cat("\n=== TOP 10 BEST-CALIBRATED FUNCTIONAL MODELS ===\n")
top_functional <- all_results %>%
  filter(group_type == "Functional") %>%
  head(10)
print(top_functional[, c('group', 'mean_abs_percent_change', 'ci_width_ratio', 'combined_score')])

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
summary_stats <- all_results %>%
  group_by(group_type) %>%
  summarise(
    n_models = n(),
    mean_continuity_mean = mean(mean_abs_percent_change, na.rm=TRUE),
    mean_continuity_median = median(mean_abs_percent_change, na.rm=TRUE),
    uncertainty_mean = mean(ci_width_ratio, na.rm=TRUE),
    uncertainty_median = median(ci_width_ratio, na.rm=TRUE),
    .groups = 'drop'
  )
print(summary_stats)

# Save results
write.csv(all_results, here("figures", "comprehensive_all_models_calibration_results.csv"), row.names = FALSE)
cat("\nResults saved to figures/comprehensive_all_models_calibration_results.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

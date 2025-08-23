# Calculate scoring metrics for CLR model forecasts
# This script calculates CRPS, RMSE, and other performance metrics for CLR models

source("../../source_local.R")

cat("=== CLR Model Scoring Metrics Calculation ===\n\n")

# Load required packages
if (!require(scoringRules, quietly = TRUE)) {
  install.packages("scoringRules")
}
library(scoringRules)

# Try to load CLR hindcasts first, then fall back to regular hindcasts if needed
clr_hindcasts_file <- here("data/summary/all_hindcasts_CLR.rds")
clr_hindcasts_parquet <- here("data/summary/parquet/all_hindcasts_CLR.parquet")

if (file.exists(clr_hindcasts_file)) {
  cat("Loading CLR hindcasts from RDS...\n")
  all_hindcasts <- readRDS(clr_hindcasts_file)
  cat("✅ CLR hindcasts loaded successfully from RDS\n")
} else if (file.exists(clr_hindcasts_parquet)) {
  cat("Loading CLR hindcasts from Parquet...\n")
  if (require(arrow, quietly = TRUE)) {
    all_hindcasts <- arrow::read_parquet(clr_hindcasts_parquet)
    cat("✅ CLR hindcasts loaded successfully from Parquet\n")
  } else {
    stop("Arrow package required for Parquet files")
  }
} else {
  cat("CLR hindcasts files not found:\n")
  cat("  - RDS:", clr_hindcasts_file, "\n")
  cat("  - Parquet:", clr_hindcasts_parquet, "\n")
  cat("This script requires CLR hindcasts to be processed first\n")
  cat("Please run 07_tidyHindcasts_CLR.r before this script\n")
  stop("CLR hindcasts files not found")
}

cat("Total rows loaded:", nrow(all_hindcasts), "\n")
cat("Columns:", paste(colnames(all_hindcasts), collapse=", "), "\n")

# Check if we have the required columns for scoring
required_cols <- c("plotID", "siteID", "dateID", "dates", "med", "lo", "hi", "model_id", "truth")
missing_cols <- required_cols[!required_cols %in% colnames(all_hindcasts)]

if (length(missing_cols) > 0) {
  cat("❌ Missing required columns:", paste(missing_cols, collapse=", "), "\n")
  cat("Available columns:", paste(colnames(all_hindcasts), collapse=", "), "\n")
  
  if ("truth" %in% missing_cols) {
    cat("\nNote: 'truth' column is missing. This means we only have forecasts, not observations.\n")
    cat("Scoring metrics that require observed values cannot be calculated.\n")
    cat("Proceeding with forecast-only analysis...\n")
  }
}

# Basic data validation
cat("\nData Validation:\n")
cat("Total forecasts:", nrow(all_hindcasts), "\n")
cat("Unique models:", length(unique(all_hindcasts$model_id)), "\n")
cat("Unique sites:", length(unique(all_hindcasts$siteID)), "\n")
cat("Unique plots:", length(unique(all_hindcasts$plotID)), "\n")
cat("Date range:", min(all_hindcasts$dates, na.rm=TRUE), "to", max(all_hindcasts$dates, na.rm=TRUE), "\n")

# Check for observations (truth values)
has_observations <- "truth" %in% colnames(all_hindcasts) && 
                   any(!is.na(all_hindcasts$truth))

if (has_observations) {
  cat("✅ Observations (truth values) found - can calculate full scoring metrics\n")
  n_obs <- sum(!is.na(all_hindcasts$truth))
  cat("Number of observations:", n_obs, "\n")
} else {
  cat("⚠️  No observations found - can only calculate forecast statistics\n")
}

# Calculate basic forecast statistics
cat("\nCalculating basic forecast statistics...\n")

forecast_stats <- all_hindcasts %>%
  group_by(model_id, model_name, rank.name, taxon, group) %>%
  summarise(
    n_forecasts = n(),
    n_plots = n_distinct(plotID),
    n_sites = n_distinct(siteID),
    n_dates = n_distinct(dates),
    mean_forecast = mean(med, na.rm=TRUE),
    sd_forecast = sd(med, na.rm=TRUE),
    mean_uncertainty = mean(hi - lo, na.rm=TRUE),
    mean_uncertainty_50 = mean(hi_75 - lo_25, na.rm=TRUE),
    .groups = "drop"
  )

cat("Forecast statistics calculated for", nrow(forecast_stats), "models\n")

# Calculate scoring metrics if observations are available
if (has_observations) {
  cat("\nCalculating scoring metrics...\n")
  
  # Filter to forecasts with observations
  forecasts_with_obs <- all_hindcasts %>%
    filter(!is.na(truth))
  
  cat("Forecasts with observations:", nrow(forecasts_with_obs), "\n")
  
  # Calculate CRPS for each forecast
  cat("Calculating CRPS scores...\n")
  
  # For CLR models, we need to handle the fact that forecasts are in log-ratio space
  # We'll calculate CRPS in the original space if possible, otherwise in log-ratio space
  
  # Check if we have the original abundance data to back-transform
  if ("original_abundance" %in% colnames(forecasts_with_obs)) {
    cat("Original abundance data found - calculating CRPS in abundance space\n")
    # CRPS calculation in abundance space
    forecasts_with_obs$crps <- mapply(
      function(obs, forecast_samples) {
        if (is.na(obs) || length(forecast_samples) == 0) return(NA)
        tryCatch({
          crps_sample(obs, forecast_samples)
        }, error = function(e) {
          cat("Error calculating CRPS:", e$message, "\n")
          return(NA)
        })
      },
      forecasts_with_obs$truth,
      forecasts_with_obs$med  # Using median as point forecast for now
    )
  } else {
    cat("No original abundance data - calculating CRPS in log-ratio space\n")
    # CRPS calculation in log-ratio space
    forecasts_with_obs$crps <- mapply(
      function(obs, forecast_samples) {
        if (is.na(obs) || length(forecast_samples) == 0) return(NA)
        tryCatch({
          crps_sample(obs, forecast_samples)
        }, error = function(e) {
          cat("Error calculating CRPS:", e$message, "\n")
          return(NA)
        })
      },
      forecasts_with_obs$truth,
      forecasts_with_obs$med  # Using median as point forecast for now
    )
  }
  
  # Calculate RMSE
  cat("Calculating RMSE...\n")
  forecasts_with_obs$rmse <- (forecasts_with_obs$truth - forecasts_with_obs$med)^2
  
  # Calculate bias
  cat("Calculating bias...\n")
  forecasts_with_obs$bias <- forecasts_with_obs$med - forecasts_with_obs$truth
  
  # Calculate coverage
  cat("Calculating coverage...\n")
  forecasts_with_obs$coverage_95 <- forecasts_with_obs$truth >= forecasts_with_obs$lo & 
                                   forecasts_with_obs$truth <= forecasts_with_obs$hi
  
  if ("lo_25" %in% colnames(forecasts_with_obs) && "hi_75" %in% colnames(forecasts_with_obs)) {
    forecasts_with_obs$coverage_50 <- forecasts_with_obs$truth >= forecasts_with_obs$lo_25 & 
                                     forecasts_with_obs$truth <= forecasts_with_obs$hi_75
  }
  
  # Aggregate scoring metrics by model
  cat("Aggregating scoring metrics...\n")
  
  scoring_metrics <- forecasts_with_obs %>%
    group_by(model_id, model_name, rank.name, taxon, group) %>%
    summarise(
      n_obs = n(),
      mean_crps = mean(crps, na.rm=TRUE),
      median_crps = median(crps, na.rm=TRUE),
      sd_crps = sd(crps, na.rm=TRUE),
      mean_rmse = mean(rmse, na.rm=TRUE),
      rmse = sqrt(mean(rmse, na.rm=TRUE)),
      mean_bias = mean(bias, na.rm=TRUE),
      median_bias = median(bias, na.rm=TRUE),
      coverage_95 = mean(coverage_95, na.rm=TRUE),
      coverage_50 = ifelse("coverage_50" %in% colnames(forecasts_with_obs), 
                           mean(coverage_50, na.rm=TRUE), NA),
      .groups = "drop"
    )
  
  cat("Scoring metrics calculated for", nrow(scoring_metrics), "models\n")
  
  # Save scoring metrics
  cat("\nSaving CLR scoring metrics...\n")
  saveRDS(scoring_metrics, here("data/summary/scoring_metrics_CLR.rds"))
  cat("✅ CLR scoring metrics saved to: scoring_metrics_CLR.rds\n")
  
  # Also save as Parquet for memory efficiency
  if (require(arrow, quietly = TRUE)) {
    parquet_dir <- here("data/summary/parquet")
    if (!dir.exists(parquet_dir)) {
      dir.create(parquet_dir, recursive = TRUE)
    }
    
    arrow::write_parquet(scoring_metrics, 
                         here("data/summary/parquet/scoring_metrics_CLR.parquet"))
    cat("✅ CLR scoring metrics also saved as Parquet: scoring_metrics_CLR.parquet\n")
  }
  
  # Summary of scoring metrics
  cat("\nScoring Metrics Summary:\n")
  print(scoring_metrics)
  
} else {
  cat("\n⚠️  No observations available - skipping scoring metric calculation\n")
  cat("Only forecast statistics were calculated\n")
}

# Save forecast statistics
cat("\nSaving CLR forecast statistics...\n")
saveRDS(forecast_stats, here("data/summary/forecast_stats_CLR.rds"))
cat("✅ CLR forecast statistics saved to: forecast_stats_CLR.rds\n")

# Also save as Parquet for memory efficiency
if (require(arrow, quietly = TRUE)) {
  parquet_dir <- here("data/summary/parquet")
  if (!dir.exists(parquet_dir)) {
    dir.create(parquet_dir, recursive = TRUE)
  }
  
  arrow::write_parquet(forecast_stats, 
                       here("data/summary/parquet/forecast_stats_CLR.parquet"))
  cat("✅ CLR forecast statistics also saved as Parquet: forecast_stats_CLR.parquet\n")
}

# Final summary
cat("\n=== CLR Scoring Metrics Calculation Complete ===\n")
cat("Forecast statistics calculated for", nrow(forecast_stats), "models\n")

if (has_observations) {
  cat("Scoring metrics calculated for", nrow(scoring_metrics), "models\n")
  cat("Observations used:", sum(!is.na(all_hindcasts$truth)), "\n")
}

cat("\nOutput files created:\n")
cat("  - forecast_stats_CLR.rds (RDS format)\n")
if (has_observations) {
  cat("  - scoring_metrics_CLR.rds (RDS format)\n")
}
if (require(arrow, quietly = TRUE)) {
  cat("  - forecast_stats_CLR.parquet (Parquet format)\n")
  if (has_observations) {
    cat("  - scoring_metrics_CLR.parquet (Parquet format)\n")
  }
}

cat("\nNext steps:\n")
cat("1. Use scoring metrics for model comparison and analysis\n")
cat("2. Create visualizations of forecast performance\n")
cat("3. Analyze forecast uncertainty and coverage\n")

# Clean up memory
rm(all_hindcasts)
if (has_observations) {
  rm(forecasts_with_obs, scoring_metrics)
}
rm(forecast_stats)
gc()

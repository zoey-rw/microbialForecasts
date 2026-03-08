# Process and tidy CLR hindcast outputs for downstream analysis
# This script processes the raw CLR hindcast data and prepares it for scoring metrics

source("../../source.R")
library(data.table)
library(dplyr)

cat("=== CLR Hindcast Processing and Tidying ===\n\n")

# Try to load CLR hindcasts first, then fall back to regular hindcasts if needed
clr_hindcasts_file <- here("data/model_outputs/CLR_hindcasts.rds")
if (file.exists(clr_hindcasts_file)) {
  cat("Loading CLR hindcasts...\n")
  all_hindcasts <- readRDS(clr_hindcasts_file)
  cat("✅ CLR hindcasts loaded successfully\n")
  cat("Total rows:", nrow(all_hindcasts), "\n")
  cat("Columns:", paste(colnames(all_hindcasts), collapse=", "), "\n")
} else {
  cat("CLR hindcasts file not found:", clr_hindcasts_file, "\n")
  cat("This script requires CLR hindcasts to be generated first\n")
  cat("Please run 06_createHindcasts_CLR.r before this script\n")
  stop("CLR hindcasts file not found")
}

# Check if we have the required columns for processing
required_cols <- c("plotID", "siteID", "dateID", "dates", "med", "lo", "hi", "model_id")
missing_cols <- required_cols[!required_cols %in% colnames(all_hindcasts)]

if (length(missing_cols) > 0) {
  cat("❌ Missing required columns:", paste(missing_cols, collapse=", "), "\n")
  cat("Available columns:", paste(colnames(all_hindcasts), collapse=", "), "\n")
  stop("Missing required columns for hindcast processing")
}

# Convert to data.table for better performance
all_hindcasts <- as.data.table(all_hindcasts)

# Basic data validation
cat("\nData Validation:\n")
cat("Total forecasts:", nrow(all_hindcasts), "\n")
cat("Unique models:", length(unique(all_hindcasts$model_id)), "\n")
cat("Unique sites:", length(unique(all_hindcasts$siteID)), "\n")
cat("Unique plots:", length(unique(all_hindcasts$plotID)), "\n")
cat("Date range:", min(all_hindcasts$dates, na.rm=TRUE), "to", max(all_hindcasts$dates, na.rm=TRUE), "\n")

# Check for missing values in critical columns
missing_med <- sum(is.na(all_hindcasts$med))
missing_lo <- sum(is.na(all_hindcasts$lo))
missing_hi <- sum(is.na(all_hindcasts$hi))

cat("Missing values:\n")
cat("  med:", missing_med, "\n")
cat("  lo:", missing_lo, "\n")
cat("  hi:", missing_hi, "\n")

# Clean data - remove rows with missing critical values
cat("\nCleaning data...\n")
all_hindcasts_clean <- all_hindcasts[!is.na(med) & !is.na(lo) & !is.na(hi)]
cat("After cleaning:", nrow(all_hindcasts_clean), "rows\n")

# Add derived columns for analysis
cat("\nAdding derived columns...\n")

# Add uncertainty measures
all_hindcasts_clean[, uncertainty := hi - lo]
all_hindcasts_clean[, uncertainty_50 := hi_75 - lo_25]

# Add model metadata if missing
if (!"model_name" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, model_name := "CLR"]
}

if (!"fcast_type" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, fcast_type := "CLR"]
}

if (!"pretty_group" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, pretty_group := "CLR"]
}

# Add site prediction if missing
if (!"site_prediction" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, site_prediction := "CLR Model"]
}

# Add timepoint if missing
if (!"timepoint" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, timepoint := date_num]
}

# Add start date processing for script 8 compatibility
cat("\nAdding start date information...\n")

# Add required columns for script 8 compatibility
# Now that hindcast functions include plot_start and site_start, use those values
all_hindcasts_clean[, is_site_start_date := FALSE]
all_hindcasts_clean[, is_plot_start_date := FALSE]
all_hindcasts_clean[, is_any_start_date := FALSE]

# Use the actual start date values from the hindcast functions
if ("site_start" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, site_start_date := site_start]
} else {
  all_hindcasts_clean[, site_start_date := NA]
}

if ("plot_start" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, plot_start_date := plot_start]
} else {
  all_hindcasts_clean[, plot_start_date := NA]
}

# Determine if dates are start dates using the actual start date values
if ("site_start_date" %in% colnames(all_hindcasts_clean) && "plot_start_date" %in% colnames(all_hindcasts_clean)) {
  all_hindcasts_clean[, is_site_start_date := (date_num == site_start_date)]
  all_hindcasts_clean[, is_plot_start_date := (date_num == plot_start_date)]
  all_hindcasts_clean[, is_any_start_date := (is_site_start_date | is_plot_start_date)]
} else {
  # Fallback if start date columns don't exist
  all_hindcasts_clean[, is_site_start_date := FALSE]
  all_hindcasts_clean[, is_plot_start_date := FALSE]
  all_hindcasts_clean[, is_any_start_date := FALSE]
}

# Model summary
cat("\nModel Summary:\n")
model_summary <- all_hindcasts_clean %>%
  group_by(model_id, model_name) %>%
  summarise(
    n_forecasts = n(),
    n_plots = n_distinct(plotID),
    n_sites = n_distinct(siteID),
    n_dates = n_distinct(dates),
    mean_uncertainty = mean(uncertainty, na.rm=TRUE),
    mean_uncertainty_50 = mean(uncertainty_50, na.rm=TRUE),
    .groups = "drop"
  )

print(model_summary)

# Save processed hindcasts
cat("\nSaving processed CLR hindcasts...\n")
saveRDS(all_hindcasts_clean, here("data/summary/all_hindcasts_CLR.rds"))
cat("✅ Processed CLR hindcasts saved to: all_hindcasts_CLR.rds\n")

# Also save as Parquet for memory efficiency
if (require(arrow, quietly = TRUE)) {
  cat("Saving as Parquet for memory efficiency...\n")
  parquet_dir <- here("data/summary/parquet")
  if (!dir.exists(parquet_dir)) {
    dir.create(parquet_dir, recursive = TRUE)
  }
  
  arrow::write_parquet(all_hindcasts_clean, 
                       here("data/summary/parquet/all_hindcasts_CLR.parquet"))
  cat("✅ CLR hindcasts also saved as Parquet: all_hindcasts_CLR.parquet\n")
} else {
  cat("Arrow package not available - skipping Parquet export\n")
}

# Final summary
cat("\n=== CLR Hindcast Processing Complete ===\n")
cat("Final dataset size:", nrow(all_hindcasts_clean), "rows\n")
cat("Models processed:", length(unique(all_hindcasts_clean$model_id)), "\n")
cat("Sites covered:", length(unique(all_hindcasts_clean$siteID)), "\n")
cat("Plots covered:", length(unique(all_hindcasts_clean$plotID)), "\n")
cat("Date range:", min(all_hindcasts_clean$dates, na.rm=TRUE), "to", max(all_hindcasts_clean$dates, na.rm=TRUE), "\n")

cat("\nOutput files created:\n")
cat("  - all_hindcasts_CLR.rds (RDS format)\n")
if (require(arrow, quietly = TRUE)) {
  cat("  - all_hindcasts_CLR.parquet (Parquet format)\n")
}

cat("\nNext steps:\n")
cat("1. Run 08_calculateScoringMetrics_CLR.r to calculate forecast performance metrics\n")
cat("2. Use processed hindcasts for downstream analysis and visualization\n")
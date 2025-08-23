# Process and tidy CLR hindcast outputs for downstream analysis
# This script processes the raw CLR hindcast data and prepares it for scoring metrics

source("../../source_local.R")

cat("=== CLR Hindcast Processing and Tidying ===\n\n")

# Load data.table for optimization
if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

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
cat("  Median:", missing_med, "(", round(100*missing_med/nrow(all_hindcasts), 1), "%)\n")
cat("  Lower bound:", missing_lo, "(", round(100*missing_lo/nrow(all_hindcasts), 1), "%)\n")
cat("  Upper bound:", missing_hi, "(", round(100*missing_hi/nrow(all_hindcasts), 1), "%)\n")

# Remove rows with missing critical values
cat("\nCleaning data...\n")
initial_rows <- nrow(all_hindcasts)
all_hindcasts_clean <- all_hindcasts %>%
  filter(!is.na(med), !is.na(lo), !is.na(hi))

cat("Rows after cleaning:", nrow(all_hindcasts_clean), "\n")
cat("Rows removed:", initial_rows - nrow(all_hindcasts_clean), "\n")

# Add additional metadata columns if they don't exist
if (!"rank.name" %in% colnames(all_hindcasts_clean)) {
  cat("Adding rank.name column...\n")
  all_hindcasts_clean$rank.name <- sapply(strsplit(all_hindcasts_clean$model_id, "_"), function(x) x[2])
}

if (!"taxon" %in% colnames(all_hindcasts_clean)) {
  cat("Adding taxon column...\n")
  all_hindcasts_clean$taxon <- sapply(strsplit(all_hindcasts_clean$model_id, "_"), function(x) x[3])
}

if (!"time_period" %in% colnames(all_hindcasts_clean)) {
  cat("Adding time_period column...\n")
  all_hindcasts_clean$time_period <- sapply(strsplit(all_hindcasts_clean$model_id, "_"), function(x) paste(x[4:5], collapse="_"))
}

if (!"model_name" %in% colnames(all_hindcasts_clean)) {
  cat("Adding model_name column...\n")
  all_hindcasts_clean$model_name <- sapply(strsplit(all_hindcasts_clean$model_id, "_"), function(x) x[1])
}

# Add group classification
cat("Adding group classification...\n")
all_hindcasts_clean$group <- ifelse(grepl("bac|16S", all_hindcasts_clean$rank.name), "16S", "ITS")
all_hindcasts_clean$pretty_group <- ifelse(all_hindcasts_clean$group == "16S", "Bacteria", "Fungi")

# Add rank classification
cat("Adding rank classification...\n")
all_hindcasts_clean$rank_only <- sapply(strsplit(all_hindcasts_clean$rank.name, "_"), function(x) x[1])

# Add forecast type classification
cat("Adding forecast type classification...\n")
all_hindcasts_clean$fcast_type <- "CLR"  # All CLR models

# Add pretty names for ranks
cat("Adding pretty rank names...\n")
rank_names_map <- c(
  "phylum" = "Phylum",
  "class" = "Class", 
  "order" = "Order",
  "family" = "Family",
  "genus" = "Genus",
  "functional" = "Functional group"
)
all_hindcasts_clean$pretty_name <- rank_names_map[all_hindcasts_clean$rank_only]

# Add year and month columns for temporal analysis
cat("Adding temporal columns...\n")
all_hindcasts_clean$year <- lubridate::year(all_hindcasts_clean$dates)
all_hindcasts_clean$month <- lubridate::month(all_hindcasts_clean$dates)
all_hindcasts_clean$season <- case_when(
  all_hindcasts_clean$month %in% c(12, 1, 2) ~ "Winter",
  all_hindcasts_clean$month %in% c(3, 4, 5) ~ "Spring", 
  all_hindcasts_clean$month %in% c(6, 7, 8) ~ "Summer",
  all_hindcasts_clean$month %in% c(9, 10, 11) ~ "Fall"
)

# Add uncertainty measures
cat("Adding uncertainty measures...\n")
all_hindcasts_clean$uncertainty <- all_hindcasts_clean$hi - all_hindcasts_clean$lo
all_hindcasts_clean$uncertainty_50 <- all_hindcasts_clean$hi_75 - all_hindcasts_clean$lo_25

# Add forecast horizon (days from start of forecast period)
cat("Adding forecast horizon...\n")
all_hindcasts_clean <- all_hindcasts_clean %>%
  group_by(model_id, plotID) %>%
  mutate(
    min_date = min(dates, na.rm=TRUE),
    forecast_horizon = as.numeric(difftime(dates, min_date, units="days"))
  ) %>%
  ungroup()

# Quality checks
cat("\nQuality Checks:\n")
cat("Forecasts with valid uncertainty:", sum(all_hindcasts_clean$uncertainty > 0, na.rm=TRUE), "\n")
cat("Forecasts with valid 50% uncertainty:", sum(all_hindcasts_clean$uncertainty_50 > 0, na.rm=TRUE), "\n")
cat("Forecasts with reasonable bounds (lo < med < hi):", 
    sum(all_hindcasts_clean$lo < all_hindcasts_clean$med & 
         all_hindcasts_clean$med < all_hindcasts_clean$hi, na.rm=TRUE), "\n")

# Summary statistics by model
cat("\nSummary by Model:\n")
model_summary <- all_hindcasts_clean %>%
  group_by(model_id, model_name, rank.name, taxon, group) %>%
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

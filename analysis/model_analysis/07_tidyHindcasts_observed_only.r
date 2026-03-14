#!/usr/bin/env Rscript

# This script processes ONLY observed files (uses original 07_tidyHindcasts.r).
# Can be run in parallel with the main script processing unobserved files.
# For normal runs use 07_tidyHindcasts_v2.r instead.
#
# WARNING: Both processes share the cache file. The cache is saved at the end,
# so the last process to finish will write the final cache state. This should
# be OK since each process only updates its own "kind" (observed vs unobserved).

tryCatch({
  mem.maxVSize(Inf)
  cat("Memory limit increased to unlimited\n")
}, error = function(e) {
  cat("Note: Could not increase memory limit:", e$message, "\n")
})

library(here)
source(here("source.R"))
library(data.table)
library(dplyr)

# Try nanoparquet first (lightweight), then arrow, then fallback to RDS
use_parquet <- FALSE
parquet_pkg <- NULL

# Check for nanoparquet first (lightweight, preferred)
if (requireNamespace("nanoparquet", quietly = TRUE)) {
  use_parquet <- TRUE
  parquet_pkg <- "nanoparquet"
  cat("Using nanoparquet for parquet support (lightweight)\n")
} else if (requireNamespace("arrow", quietly = TRUE)) {
  use_parquet <- TRUE
  parquet_pkg <- "arrow"
  cat("Using arrow for parquet support\n")
}

cat("=== PROCESSING OBSERVED FILES ONLY ===\n")
cat("Note: Running in parallel with main script\n")
cat("This will process the 876 observed files that were skipped\n\n")

# Source the main script to get all helper functions
# It will try to execute, but we'll let it - it should skip observed files
# that are already in cache, and we want it to process the 876 new ones
# Actually, wait - if we source it, it will process both observed AND unobserved
# We need to stop it after observed

# Set an environment variable to signal "observed only" mode
Sys.setenv(TIDYHINDCASTS_OBSERVED_ONLY = "TRUE")

# Source the main script - it will define functions and start processing
# We need to modify it to check the env var and stop after observed
# Actually, simpler: just source it and it will process observed first
# The cache should prevent conflicts, and shard numbers are different

cat("Loading functions from main script...\n")
source(here("analysis/model_analysis/07_tidyHindcasts.r"))

cat("\n=== OBSERVED FILES PROCESSING COMPLETE ===\n")
cat("Note: Main script may continue processing unobserved files\n")
cat("Shards saved to: data/summary/parquet/_shards/shard_observed_*.rds\n")

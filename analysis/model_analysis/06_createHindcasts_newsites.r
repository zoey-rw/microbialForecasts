#!/usr/bin/env Rscript

# Script for generating REAL hindcasts for unobserved sites
# Key fixes: Remove aggressive environmental filtering, use full time series, simplify timepoint logic
# CRITICAL FIX: Added fallback mechanism to read individual site files if site_output_list is incomplete
# This prevents truth data loss when some sites fail during processing but individual files are generated

cat("=== GENERATING HINDCASTS FOR UNOBSERVED SITES ===\n")

# Setup
source("../../source.R")
library(dplyr)
library(here)

## --- Helper Functions ---

# Deduplicate by keys with optional prioritization
dedup_key <- function(df, keys, prefer = NULL, arrange_by = NULL) {
  keys <- keys[keys %in% names(df)]
  if (!length(keys) || !nrow(df)) return(df)
  
  if (!is.null(prefer) && prefer %in% names(df)) {
    df <- df %>%
      dplyr::mutate(.prefer = !is.na(.data[[prefer]]) & is.finite(.data[[prefer]])) %>%
      dplyr::arrange(dplyr::desc(.prefer))
  }
  if (!is.null(arrange_by) && all(arrange_by %in% names(df))) {
    df <- df %>% dplyr::arrange(dplyr::across(dplyr::all_of(arrange_by)))
  }
  df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(keys)), .keep_all = TRUE) %>%
    dplyr::select(-dplyr::any_of(".prefer"))
}

# Ensure dateID and dates columns are properly formatted
ensure_dates <- function(df) {
  if ("dateID" %in% names(df)) {
    df$dateID <- suppressWarnings(as.numeric(as.character(df$dateID)))
  }
  if ("dates" %in% names(df) && !inherits(df$dates, "Date")) {
    df$dates <- tryCatch(
      as.Date(df$dates),
      error = function(e) {
        if ("dateID" %in% names(df)) {
          as.Date(paste0(df$dateID, "01"), "%Y%m%d")
        } else {
          as.Date(NA)
        }
      }
    )
  }
  if (!"dates" %in% names(df) && "dateID" %in% names(df)) {
    df$dates <- as.Date(paste0(df$dateID, "01"), "%Y%m%d")
  }
  df
}

# Ensure columns exist with default values
ensure_cols <- function(df, cols, default = NA_real_) {
  for (nm in cols) {
    if (!nm %in% names(df)) df[[nm]] <- default
  }
  df
}

# Debug printing helper
dbg <- function(..., .on = FALSE) {
  if (.on) {
    cat(...)
    cat("\n")
    try(flush.console(), silent = TRUE)
  }
}

# Safe readRDS with retry logic
safe_readRDS <- function(path, retry = 1, sleep = 0.5) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      if (retry > 0 && grepl("connection|cannot open", e$message, ignore.case = TRUE)) {
        Sys.sleep(sleep)
        safe_readRDS(path, retry = retry - 1, sleep = sleep)
      } else {
        stop(e)
      }
    }
  )
}

# Source the fixed prepBetaRegData function directly
source(here::here("microbialForecast/R/prepBetaRegData.r"))
source(here::here("microbialForecast/R/run_hindcast.r"))

# Load required packages for visualization and optimization
library(ggplot2)
library(gridExtra)
library(data.table)
library(stringr)
library(tidyr)
library(lubridate)
library(padr)
library(parallel)
library(doParallel)
library(foreach)

# CRITICAL: Limit thread usage to prevent system overload
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1") 
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv("OMP_THREAD_LIMIT" = "1")  # avoids OpenMP device contention in parallel workers

## Define stable absolute project root for workers
project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
cat("PROJECT ROOT:", project_root, "\n")

# Test mode controls
is_local_test <- identical(tolower(Sys.getenv("LOCAL_TEST", "false")), "true")
if (is_local_test) {
  cat("🧪 LOCAL TEST MODE ENABLED\n")
}

# Diagnostic figures flag
make_figs <- identical(tolower(Sys.getenv("FIGS", "true")), "true")
if (make_figs) {
  cat("📊 DIAGNOSTIC FIGURES ENABLED\n")
}

all_results <- list()
successful_taxa <- 0
failed_taxa <- 0

# Discover available models directly from file system (ignore CSV)
cat("Discovering available models from file system...\n")

# Search in multiple model output directories
# Override with MODEL_DIRS env var to read from a different output directory (e.g., rerun)
env_model_dirs <- Sys.getenv("MODEL_DIRS", "")
if (nzchar(env_model_dirs)) {
  model_dirs <- strsplit(env_model_dirs, ",")[[1]]
  cat("MODEL_DIRS override:", paste(model_dirs, collapse=", "), "\n")
} else {
  model_dirs <- c("cloglog_beta_driver_uncertainty")
}

all_files <- character(0)
for (model_dir in model_dirs) {
  model_output_dir <- here("data/model_outputs", model_dir)
  if (dir.exists(model_output_dir)) {
    cat("Searching for summaries in:", model_dir, "\n")
    files <- list.files(
      path = model_output_dir,
      pattern = "summary_.*_beta_regression.rds",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(files) > 0) {
      cat("Found", length(files), "files in", model_dir, "\n")
      all_files <- c(all_files, files)
    }
  }
}
cat("Found", length(all_files), "summary files total\n")

# Filter out individual chain files - only keep combined files
combined_files <- all_files[!grepl("chain[0-9]", all_files)]
cat("Found", length(all_files), "total files (including chains)\n")
cat("Found", length(combined_files), "combined model files (excluding individual chains)\n")

# Throttle in local test
if (is_local_test && length(combined_files) > 0) {
  # Keep a small representative subset
  combined_files <- head(sort(combined_files), 6)
  cat("LOCAL_TEST: Limiting combined files to", length(combined_files), "\n")
}

# Extract model information from file paths
model_info <- data.frame(
  file_path = combined_files,
  stringsAsFactors = FALSE
)

# Extract model subdirectory from file path
extract_model_subdir <- function(file_path) {
  # Extract the subdirectory name from path like: .../model_outputs/SUBDIR/model_name/...
  path_parts <- strsplit(file_path, "/")[[1]]
  model_outputs_idx <- which(path_parts == "model_outputs")
  if (length(model_outputs_idx) > 0 && model_outputs_idx < length(path_parts)) {
    return(path_parts[model_outputs_idx + 1])
  }
  return("cloglog_beta_driver_uncertainty")  # Fallback
}

model_info$model_subdir <- sapply(model_info$file_path, extract_model_subdir)

# Get file modification time for deduplication
file_info <- file.info(model_info$file_path)
model_info$mtime <- file_info$mtime

# Parse model ID from filename
model_info$model_id <- gsub("^.*summary_", "", basename(model_info$file_path))
model_info$model_id <- gsub("_beta_regression\\.rds$", "", model_info$model_id)
model_info$model_id <- gsub("_with_legacy_covariate_beta_regression\\.rds$", "", model_info$model_id)

# Parse model name and taxon from model_id
# CRITICAL FIX: Handle double-prefix models (e.g., cycl_only_cycl_only_...)
# Model ID format: env_cycl_herbicide_stress_20130601_20180101_with_legacy_covariate
# OR: cycl_only_cycl_only_dissim_nitrate_reduction_20130601_20180101_with_legacy_covariate
model_info$model_name <- ifelse(
  grepl("^cycl_only_cycl_only_", model_info$model_id), "cycl_only",
  ifelse(
    grepl("^env_cycl_env_cycl_", model_info$model_id), "env_cycl",
    ifelse(
      grepl("^env_cov_env_cov_", model_info$model_id), "env_cov",
      gsub("^([^_]+_[^_]+)_.*$", "\\1", model_info$model_id)
    )
  )
)

model_info$taxon <- ifelse(
  grepl("^cycl_only_cycl_only_", model_info$model_id), 
  gsub("^cycl_only_cycl_only_", "", model_info$model_id),
  ifelse(
    grepl("^env_cycl_env_cycl_", model_info$model_id),
    gsub("^env_cycl_env_cycl_", "", model_info$model_id),
    ifelse(
      grepl("^env_cov_env_cov_", model_info$model_id),
      gsub("^env_cov_env_cov_", "", model_info$model_id),
      gsub("^[^_]+_[^_]+_", "", model_info$model_id)
    )
  )
)

# Extract taxon (remove date suffix)
model_info$taxon <- gsub("_[0-9]{8}_.*$", "", model_info$taxon)

# Parse date range from model_id
# Extract min_date (first 8-digit date)
model_info$min_date <- gsub("^.*_([0-9]{8})_[0-9]{8}_.*$", "\\1", model_info$model_id)
# Extract max_date (second 8-digit date)
model_info$max_date <- gsub("^.*_[0-9]{8}_([0-9]{8})_.*$", "\\1", model_info$model_id)

# Model path is the subdirectory extracted from the file path
model_info$model_path <- file.path("data/model_outputs", model_info$model_subdir)

# Filter for models we want to process
available_models <- model_info %>%
  filter(
    model_name %in% c("env_cycl", "env_cov", "cycl_only"),
    min_date %in% c("20130601"),
    max_date %in% c("20180101")  # CRITICAL: Only 2018 is the accepted max (matching observed script)
  )

# CRITICAL: Exclude reconstructed_from_checkpoints models
n_before <- nrow(available_models)
available_models <- available_models %>%
  filter(!grepl("reconstructed_from_checkpoints", file_path, ignore.case = TRUE) &
         !grepl("reconstructed_from_checkpoints", model_id, ignore.case = TRUE))
if (n_before > nrow(available_models)) {
  cat("Excluded", n_before - nrow(available_models), "reconstructed_from_checkpoints models\n")
}

# CRITICAL: Ensure all models are from driver_uncertainty directories only
n_before <- nrow(available_models)
available_models <- available_models %>%
  filter(grepl("driver_uncertainty", file_path, ignore.case = TRUE))
if (n_before > nrow(available_models)) {
  cat("Excluded", n_before - nrow(available_models), "non-driver_uncertainty models\n")
}

# CRITICAL: Deduplicate by model_id (keep newest file based on modification time)
# There may be multiple summary files for the same model_id (e.g., from different chains or runs)
# Newer files are typically in parent directories, older ones in taxon-specific subdirectories
n_before_dedup <- nrow(available_models)
if (n_before_dedup > 0) {
  # Sort by model_id and mtime (newest first), then keep first occurrence of each model_id
  available_models <- available_models %>%
    arrange(model_id, desc(mtime)) %>%
    distinct(model_id, .keep_all = TRUE)
  # Remove mtime column as it's no longer needed
  available_models$mtime <- NULL
}
if (n_before_dedup > nrow(available_models)) {
  cat("Deduplicated", n_before_dedup - nrow(available_models), "duplicate model_ids (keeping newest file by modification time)\n")
  cat("Final unique models:", nrow(available_models), "\n")
}

# Prioritize env_cycl models first, then env_cov, then cycl_only
available_models <- available_models %>%
  arrange(
    factor(model_name, levels = c("env_cycl", "env_cov", "cycl_only"), ordered = TRUE),
    taxon, min_date, max_date
  )

# Optional taxon filter via env var: TAXON_FILTER="cellulolytic,saprotroph,acetate_simple"
taxon_filter_str <- Sys.getenv("TAXON_FILTER", "")
if (nchar(taxon_filter_str) > 0) {
  taxon_filter_vec <- trimws(strsplit(taxon_filter_str, ",")[[1]])
  available_models <- available_models %>% filter(taxon %in% taxon_filter_vec)
  cat("TAXON_FILTER: Limiting to", paste(taxon_filter_vec, collapse=", "),
      "->", nrow(available_models), "models\n")
}

if (is_local_test && nrow(available_models) > 0) {
  # Ensure coverage of model types where possible
  available_models <- available_models %>% group_by(model_name) %>% slice_head(n = 1) %>% ungroup()
  cat("LOCAL_TEST: Limiting available_models to", nrow(available_models), "(one per model type)\n")
}

# Add required columns for compatibility
available_models$species <- available_models$taxon
available_models$min.date <- available_models$min_date
available_models$max.date <- available_models$max_date
available_models$fcast_type <- "unobserved"

# CRITICAL FIX: Determine correct rank.name for each taxon
cat("Determining correct rank names for each taxon...\n")

# Load rank data to create taxon -> rank mapping
rank_files <- c(
  "data/clean/groupAbundances_16S_2023.rds",
  "data/clean/groupAbundances_ITS_2023.rds"
)

taxon_to_rank <- list()
for (file in rank_files) {
  if (file.exists(file)) {
    data <- safe_readRDS(file)
    for (rank_name in names(data)) {
      rank_data <- data[[rank_name]]
      taxon_cols <- colnames(rank_data)[!colnames(rank_data) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", "other")]
      for (taxon in taxon_cols) {
        taxon_to_rank[[taxon]] <- rank_name
      }
    }
  }
}

# Function to get correct rank name
get_rank_name <- function(taxon) {
  if (taxon %in% names(taxon_to_rank)) {
    return(taxon_to_rank[[taxon]])
  } else {
    return("functional_group")  # Fallback
  }
}

# Apply correct rank names
available_models$rank.name <- sapply(available_models$taxon, get_rank_name)
cat("Rank name mapping complete for", nrow(available_models), "models\n")

# Path utilities for bullet-proof filesystem operations
nz_path <- function(...) normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
exists_file <- function(...) file.exists(nz_path(...))

# Helper function to determine if a model is a driver uncertainty model
is_driver_model <- function(source_path, model_id = "") {
  # All cloglog models are driver uncertainty models
  any(grepl("/cloglog_beta_driver_uncertainty/", source_path)) || 
    any(grepl("/logit_beta_driver_uncertainty/", source_path)) || 
    grepl("driver_uncertainty", model_id)
}

# Centralized model path resolution - uses model_path parameter
resolve_model_paths <- function(project_root, model_name, taxon, model_id, model_path) {
  # Use the provided model_path (extracted from discovered file paths)
  # If model_path is a relative path, make it absolute
  if (!is.null(model_path) && !is.na(model_path) && model_path != "") {
    if (startsWith(model_path, "data/model_outputs/")) {
      base_dir <- nz_path(project_root, model_path)
    } else {
      base_dir <- nz_path(project_root, "data/model_outputs", model_path)
    }
  } else {
    # Fallback to cloglog driver uncertainty if model_path not provided
    base_dir <- nz_path(project_root, "data/model_outputs/cloglog_beta_driver_uncertainty")
  }

  # Search recursively in subdirectories for summary files
  all_summaries <- list.files(
    file.path(base_dir, model_name),
    pattern = paste0("summary_", gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", basename(model_id)), ".*\\.rds$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Search recursively in subdirectories for sample files
  all_samples <- list.files(
    file.path(base_dir, model_name),
    pattern = paste0("samples_", gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", basename(model_id)), ".*\\.rds$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Filter out chain files
  all_samples <- all_samples[!grepl("_chain[0-9]", all_samples)]

  list(
    summary = Filter(file.exists, all_summaries),
    samples = Filter(file.exists, all_samples),
    searched = list(summary = all_summaries, samples = all_samples)
  )
}

# Cache existing files to avoid repeated list.files calls (set globally before function)
existing_hindcast_files_cache <- NULL

# Get reference time from site effects files (when they were last regenerated)
get_site_effects_reference_time <- function() {
  site_effects_files <- c(
    here("data/summary/site_effects_unobserved_env_cycl.rds"),
    here("data/summary/site_effects_unobserved_env_cov.rds"),
    here("data/summary/site_effects_unobserved_cycl_only.rds")
  )
  existing <- site_effects_files[file.exists(site_effects_files)]
  if (length(existing) > 0) {
    return(max(file.info(existing)$mtime))
  }
  return(NULL)
}

# Function to check if all required hindcast files exist for a taxon AND are newer than site effects
check_existing_hindcasts <- function(model_id, required_sites) {
  # Use cached files (should be set before first call)
  if (is.null(existing_hindcast_files_cache) || length(existing_hindcast_files_cache) == 0) {
    return(FALSE)
  }
  
  # Get reference time from site effects files
  reference_time <- get_site_effects_reference_time()
  
  # Escape model_id for regex
  model_id_escaped <- gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_id)
  
  # Check for individual site files with both possible naming patterns
  # Pattern 1: original model_id
  pattern_1 <- paste0("^hindcasts_", model_id_escaped, "_")
  existing_files_1 <- existing_hindcast_files_cache[grepl(pattern_1, existing_hindcast_files_cache)]
  
  # Pattern 2: model_id with _with_legacy_covariate suffix  
  pattern_2 <- paste0("^hindcasts_", model_id_escaped, "_with_legacy_covariate_")
  existing_files_2 <- existing_hindcast_files_cache[grepl(pattern_2, existing_hindcast_files_cache)]
  
  # Check if all required sites have files in either pattern
  expected_files_1 <- paste0("hindcasts_", model_id, "_", required_sites, "_unobserved.rds")
  expected_files_2 <- paste0("hindcasts_", model_id, "_with_legacy_covariate_", required_sites, "_unobserved.rds")
  
  missing_sites_1 <- required_sites[!expected_files_1 %in% existing_files_1]
  missing_sites_2 <- required_sites[!expected_files_2 %in% existing_files_2]
  
  # Use the pattern that has fewer missing sites (or the first one if equal)
  if (length(missing_sites_1) <= length(missing_sites_2)) {
    missing_sites <- missing_sites_1
    pattern_files <- existing_files_1
    expected_files <- expected_files_1
  } else {
    missing_sites <- missing_sites_2
    pattern_files <- existing_files_2
    expected_files <- expected_files_2
  }
  
  # If files are missing, need to regenerate
  if (length(missing_sites) > 0) {
    return(FALSE)
  }
  
  # If we have a reference time, check that all files are newer than site effects
  if (!is.null(reference_time)) {
    hindcast_dir <- file.path(project_root, "data", "hindcasts", "driver_uncertainty")
    file_paths <- file.path(hindcast_dir, expected_files)
    existing_paths <- file_paths[file.exists(file_paths)]
    
    if (length(existing_paths) > 0) {
      file_times <- file.info(existing_paths)$mtime
      # All files must be newer than or equal to reference time
      if (any(file_times < reference_time)) {
        return(FALSE)  # Need to regenerate - some files are older than site effects
      }
    }
  }
  
  return(TRUE)  # All files exist and are up to date
}

# Define required sites for unobserved sites
required_unobserved_sites <- c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", "WREF", "YELL")

if (is_local_test) {
  required_unobserved_sites <- c("BARR", "BONA", "HEAL")
  cat("LOCAL_TEST: Limiting required unobserved sites to:", paste(required_unobserved_sites, collapse=", "), "\n")
}

# Pre-create output directory and cache existing files once (expensive operation)
hindcast_dir <- file.path(project_root, "data", "hindcasts", "driver_uncertainty")
dir.create(hindcast_dir, recursive = TRUE, showWarnings = FALSE)

# Cache existing hindcast files once at the start to avoid repeated list.files calls
existing_hindcast_files_cache <<- if (dir.exists(hindcast_dir)) {
  list.files(hindcast_dir, pattern = "_unobserved\\.rds$", full.names = FALSE)
} else {
  character(0)
}

# Filter out taxa that already have complete hindcast files
filtered_models <- list()
skipped_count <- 0

for (i in seq_len(nrow(available_models))) {
  taxon_config <- available_models[i, ]
  model_id <- taxon_config$model_id
  
  # Check if hindcasts already exist for all required sites
  if (check_existing_hindcasts(model_id, required_unobserved_sites)) {
    skipped_count <- skipped_count + 1
  } else {
    filtered_models[[length(filtered_models) + 1]] <- taxon_config
  }
}

cat("Total available models:", nrow(available_models), "\n")
cat("Models with complete files (skipped):", skipped_count, "\n")
cat("Models needing processing:", length(filtered_models), "\n")

# Show priority breakdown
if (length(filtered_models) > 0) {
  priority_breakdown <- table(sapply(filtered_models, function(x) x$model_name))
  cat("Priority breakdown (env_cycl first):\n")
  for (model_type in names(priority_breakdown)) {
    cat("  ", model_type, ":", priority_breakdown[[model_type]], "taxa\n")
  }
}

# Load data once before parallel processing
cat("Loading data for parallel processing...\n")
bacteria_file <- here("data/clean/groupAbundances_16S_2023.rds")
fungi_file <- here("data/clean/groupAbundances_ITS_2023.rds")

if (!file.exists(bacteria_file)) {
  stop("Bacteria data file not found: ", bacteria_file)
}
if (!file.exists(fungi_file)) {
  stop("Fungi data file not found: ", fungi_file)
}

bacteria <- safe_readRDS(bacteria_file)
fungi <- safe_readRDS(fungi_file)
all_ranks <- c(bacteria, fungi)

# Load environmental data
env_data <- safe_readRDS(here("data/clean/all_predictor_data.rds"))

# Load predicted site effects for unobserved sites - will load model-specific files per taxon
# Note: pred_effects are loaded per-model-type in process_single_taxon to ensure correct filtering

cat("✅ All data loaded successfully\n")

# Now convert to list format for processing with correct rank names
cat("Preparing taxa list with correct rank names...\n")
models_to_process_df <- dplyr::bind_rows(filtered_models)
taxa_to_process <- vector("list", nrow(models_to_process_df))
for (i in seq_len(nrow(models_to_process_df))) {
  row <- models_to_process_df[i, ]
  
  actual_model_id <- row$model_id
  
    # Determine the correct rank name based on the species
    species_name <- row$species
    
    # First, try to find the species in the available ranks
    rank_name <- NULL
    for (rank in names(all_ranks)) {
      if (species_name %in% colnames(all_ranks[[rank]])) {
        rank_name <- rank
        break
      }
    }
    
    # If not found, use fallback logic
    if (is.null(rank_name)) {
      if (grepl("_bac$", species_name)) {
        rank_name <- "genus_bac"
      } else if (grepl("_fun$", species_name)) {
        rank_name <- "genus_fun"
      } else if (species_name %in% names(all_ranks)) {
        rank_name <- species_name
      } else {
        # Final fallback to original rank.name
        rank_name <- row$rank.name
        cat("WARNING: Could not determine rank name for", species_name, ". Using fallback:", rank_name, "\n")
      }
    }
  
  taxa_to_process[[i]] <- list(
    model_id = actual_model_id,
    rank.name = rank_name,
    taxon = species_name,
    model_name = row$model_name,
    model_path = row$model_path,
    min.date = row$min.date,
    max.date = row$max.date,
    fcast_type = row$fcast_type
  )
}

cat("Prepared", length(taxa_to_process), "taxa for processing\n")

# Function to process a single taxon (for unobserved sites)
process_single_taxon <- function(taxon_config, all_ranks, env_data) {
  model_id <- taxon_config$model_id
  rank.name <- taxon_config$rank.name
  taxon <- taxon_config$taxon
  model_name <- taxon_config$model_name
  model_path <- taxon_config$model_path
  min.date <- taxon_config$min.date
  max.date <- taxon_config$max.date
  fcast_type <- taxon_config$fcast_type
  
  # Get the rank data from all_ranks
  rank.df <- all_ranks[[as.character(rank.name)]]
  if (is.null(rank.df)) {
    stop("Rank data not found for rank.name: ", rank.name)
  }
  
  # Create rank.df_spec for this taxon
  rank.df_spec <- rank.df %>%
    dplyr::select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", dplyr::all_of(taxon))
  rank.df_spec$other <- 1-rank.df_spec[[taxon]]
  
    dbg("\n=== PROCESSING TAXON:", taxon, "===", .on = TRUE)
  
  tryCatch({
    # Load model summary and samples using bullet-proof path resolution
    paths <- resolve_model_paths(project_root, model_name, taxon, model_id, model_path)
    
    if (length(paths$summary) == 0L) {
      dbg("❌ No summary file for", model_id, "\nSearched:\n  ", paste(paths$searched$summary, collapse="\n  "), .on = TRUE)
      return(list(success = FALSE, error = "Model summary not found"))
    }
    if (length(paths$samples) == 0L) {
      dbg("❌ No samples file for", model_id, "\nSearched:\n  ", paste(paths$searched$samples, collapse="\n  "), .on = TRUE)
      return(list(success = FALSE, error = "Model samples not found"))
    }

    summary_path <- normalizePath(paths$summary[[1]], winslash = "/", mustWork = TRUE)
    samples_path <- normalizePath(paths$samples[[1]], winslash = "/", mustWork = TRUE)

    model_summary <- safe_readRDS(summary_path)
    model.dat <- safe_readRDS(samples_path)
    
    # 🔐 Robust param_samples extraction
    param_samples <- NULL
    candidates <- list(
      model.dat$samples,
      model.dat$param_samples,
      if (!is.null(model.dat$fit) && !is.null(model.dat$fit$samples)) model.dat$fit$samples else NULL,
      if (is.list(model.dat) && length(model.dat) > 0 && is.list(model.dat[[1]]) && !is.null(model.dat[[1]]$samples)) model.dat[[1]]$samples else NULL
    )
    for (cand in candidates) if (!is.null(cand)) { param_samples <- cand; break }
    if (is.null(param_samples)) {
      stop("Could not locate parameter samples in model.dat (checked $samples, $param_samples, $fit$samples)")
    }
    # Convert mcmc.list to matrix once (avoids repeated conversion in fcast_logit_beta)
    if (inherits(param_samples, "mcmc.list")) {
      param_samples <- do.call(rbind, lapply(param_samples, as.matrix))
    } else if (!is.matrix(param_samples)) {
      param_samples <- as.matrix(param_samples)
    }
    
    # Compute driver flag once and reuse for all output paths
    driver_flag <- is_driver_model(samples_path, model_id)
    
      # Element map (from summarizeBetaRegModels.r):
      # [[1]] parameter summaries (summary_df)
      # [[2]] plot-level means (pred.means) - Mean, SD columns, plotID, siteID, timepoint, plot_num, site_num
      # [[3]] plot-level quantiles (pred.quantiles) - 2.5%, 25%, 50%, 75%, 97.5% columns, plotID, siteID, timepoint, plot_num, site_num
      # [[4]] gelman diagnostics (gd)
      # NOTE: This script uses [[2]] (pred.means) which is the existing approach for unobserved sites.
      # The plot alignment fix in run_hindcast.r uses plotID (not plot_num) for filtering, which works
      # with both [[2]] and [[3]] as long as plotID exists. For new sites, plotID won't exist in cal_quants
      # (which only contains observed sites), so calibration extraction correctly returns empty.

    stopifnot(length(model_summary) >= 2L, is.data.frame(model_summary[[2]]))

    cal_quants <- model_summary[[2]]
    
    # 🔐 Validate/normalize calibration frame
    if (!is.data.frame(cal_quants) || nrow(cal_quants) == 0) {
      cat("⚠️  model_summary[[2]] empty or not a data.frame for", taxon, ". Skipping taxon.\n")
      return(list(success = FALSE, error = "Empty model summary", taxon = taxon))
    }

    req_cols <- c("plotID","siteID","timepoint")
    missing  <- setdiff(req_cols, names(cal_quants))
    if (length(missing)) {
      cat("⚠️  cal_quants missing columns for", taxon, ":", paste(missing, collapse=", "), ". Skipping taxon.\n")
      return(list(success = FALSE, error = paste("Missing columns:", paste(missing, collapse=", ")), taxon = taxon))
    }

    # dateID sometimes absent — we can reconstruct later, but keep type safe if present
    cal_quants <- ensure_dates(cal_quants)

    cal_quants$plotID   <- as.character(cal_quants$plotID)
    cal_quants$siteID   <- as.character(cal_quants$siteID)
    cal_quants$timepoint<- suppressWarnings(as.integer(cal_quants$timepoint))
    
    # Validate timepoint - if NA, this indicates a problem in model summary generation
    if (any(is.na(cal_quants$timepoint))) {
      na_count <- sum(is.na(cal_quants$timepoint))
      dbg("⚠️  ERROR: Found", na_count, "NA timepoints in cal_quants for", taxon, .on = TRUE)
      dbg("  This indicates a problem in model summary generation. Cannot proceed.", .on = TRUE)
      return(list(success = FALSE, error = paste("NA timepoints in model summary:", na_count), taxon = taxon))
    }

    # IMPORTANT: use [[2]] for fcast as well — it has plot_num/site_num/timepoint you need
    plot_summary <- cal_quants
    
    # Get the rank data from pre-loaded data
    rank.df <- all_ranks[[as.character(rank.name)]]
    if (is.null(rank.df)) {
      dbg("⚠️  ERROR: Rank data not found for rank.name:", rank.name, "taxon:", taxon, .on = TRUE)
      return(list(success = FALSE, error = paste("Rank data not found:", rank.name), taxon = taxon))
    }
    
  rank.df_spec <- rank.df %>%
      dplyr::select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", dplyr::all_of(taxon))
    
    # Validate and fix dateID upstream: reconstruct from dates if missing
    rank.df_spec <- ensure_dates(rank.df_spec)
    
    # If dateID is missing, try to reconstruct from dates column
    if (any(is.na(rank.df_spec$dateID)) && "dates" %in% names(rank.df_spec)) {
      na_dateid_mask <- is.na(rank.df_spec$dateID)
      dates_vals <- rank.df_spec$dates[na_dateid_mask]
      
      # Try to extract dateID from dates (handles Date objects and character strings)
      if (inherits(dates_vals, "Date") || (is.character(dates_vals) && any(!is.na(dates_vals)))) {
        tryCatch({
          date_parsed <- as.Date(dates_vals)
          reconstructed_dateID <- as.numeric(format(date_parsed, "%Y%m"))
          rank.df_spec$dateID[na_dateid_mask] <- reconstructed_dateID
          cat("  Reconstructed", sum(na_dateid_mask), "missing dateID values from dates column\n")
        }, error = function(e) {
          # If reconstruction fails, keep NA and filter later
        })
      }
    }
    
    # Filter out rows with missing dateID - this is invalid data that shouldn't be in the raw files
    if (any(is.na(rank.df_spec$dateID))) {
      na_count <- sum(is.na(rank.df_spec$dateID))
      cat("⚠️  WARNING: Found", na_count, "rows with NA dateID in raw data for", taxon, "\n")
      cat("  This indicates invalid data in groupAbundances files - these rows should not exist.\n")
      cat("  Filtering out invalid rows.\n")
      
      rows_before <- nrow(rank.df_spec)
      rank.df_spec <- rank.df_spec[!is.na(rank.df_spec$dateID), ]
      
      if (nrow(rank.df_spec) == 0) {
        dbg("⚠️  ERROR: No valid dateIDs remaining after filtering for", taxon, ". Skipping taxon.", .on = TRUE)
        return(list(success = FALSE, error = "No valid dateIDs after filtering", taxon = taxon))
      }
      
      dbg("  Removed", rows_before - nrow(rank.df_spec), "invalid rows. Remaining:", 
          length(unique(rank.df_spec$siteID)), "sites,",
          length(unique(rank.df_spec$plotID)), "plots,",
          length(unique(rank.df_spec$dateID)), "unique dateIDs", .on = TRUE)
    }
    
    # Prepare full time series model inputs (align with observed script)
  keep_vec <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", taxon)
  
    # Hindcast end (keep your full window, but we'll still respect trained horizon)
    hindcast_end <- "20200101"

    # Build list of UNOBSERVED sites that actually exist in covariates
    correct_unobserved_sites <- c("ABBY","BARR","BONA","DEJU","HEAL","KONA","LAJA","LENO","MLBS","RMNP","SOAP","TOOL","WREF","YELL")
    available_sites <- rownames(env_data$temp)
    unobserved_sites <- intersect(correct_unobserved_sites, available_sites)
    if (length(unobserved_sites) == 0) {
      dbg("⚠️ No unobserved sites found in covariates for", taxon, .on = TRUE)
    }

    # For unobserved sites, we need to construct the data structure manually
    # because prepBetaRegData is designed for observed sites
    full.ts.model.inputs <- list()
    
    # Use the environmental data directly
    full.ts.model.inputs$temp <- env_data$temp
    full.ts.model.inputs$mois <- env_data$mois
    full.ts.model.inputs$relEM <- env_data$relEM_plot
    full.ts.model.inputs$LAI <- as.matrix(env_data$LAI)
    full.ts.model.inputs$pH <- env_data$pH
    full.ts.model.inputs$pC <- env_data$pC
    full.ts.model.inputs$pH_sd <- env_data$pH_sd
    full.ts.model.inputs$pC_sd <- env_data$pC_sd
    full.ts.model.inputs$temp_sd <- env_data$temp_sd
    full.ts.model.inputs$mois_sd <- env_data$mois_sd
    # sin_mo and cos_mo computed after trained_time_map is built (below)
    
    # Create truth.plot.long structure for unobserved sites
    # Get all plots for unobserved sites from the rank data
    unobserved_plots <- rank.df_spec %>%
      dplyr::filter(siteID %in% unobserved_sites) %>%
      dplyr::select(siteID, plotID) %>%
      dplyr::distinct()
    
    full.ts.model.inputs$truth.plot.long <- unobserved_plots %>%
      dplyr::mutate(
        plot_num = row_number(),
        site_num = as.numeric(as.factor(siteID))
      )

    # Build trained time map from summary [[2]]
    stopifnot(length(model_summary) >= 2L, is.data.frame(model_summary[[2]]))
    cal_quants <- model_summary[[2]]

    cal_quants <- ensure_dates(cal_quants)
    if ("dateID" %in% names(cal_quants)) {
      time_train <- cal_quants %>%
        dplyr::select(timepoint, dateID) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(timepoint), !is.na(dateID)) %>%
        dplyr::arrange(timepoint) %>%
        dplyr::mutate(
          trained_date_num = timepoint,
          dates = as.Date(paste0(dateID, "01"), "%Y%m%d")
        )
    } else {
      # 🔐 reconstruct from min.date + length of unique timepoints
      tp <- sort(unique(cal_quants$timepoint))
      base <- as.Date(min.date, "%Y%m%d")
      seq_dates <- base + lubridate::months(0:(length(tp) - 1))
      time_train <- data.frame(
        timepoint = tp,
        dateID    = as.numeric(format(seq_dates, "%Y%m")),
        trained_date_num = tp,
        dates     = as.Date(paste0(as.numeric(format(seq_dates, "%Y%m")), "01"), "%Y%m%d"),
        stringsAsFactors = FALSE
      )
    }
    if (nrow(time_train) == 0) stop("time_train is empty")

    cal_end_dateID   <- max(time_train$dateID, na.rm = TRUE)
    max_trained_num  <- max(time_train$trained_date_num, na.rm = TRUE)


    # For unobserved sites, use ALL their data since they weren't in the calibration period
    # All their data is considered "hindcast" data for evaluation purposes
    hindcast_truth_data <- rank.df_spec %>%
      dplyr::filter(siteID %in% unobserved_sites) %>%
      dplyr::select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", dplyr::all_of(taxon)) %>%
      dplyr::rename(truth = !!rlang::sym(taxon)) %>%
      dplyr::mutate(species = !!taxon) %>%
      ensure_dates()
    
    if (!nrow(hindcast_truth_data)) {
      cat("⚠️  No truth data for unobserved sites for", taxon, ". Skipping taxon.\n")
      return(list(success = FALSE, error = "No truth data for unobserved sites", taxon = taxon))
    }

    # For unobserved sites, create trained_time_map from their actual truth data
    # Since these sites weren't in calibration, all their dates are hindcast dates
    trained_time_map <- hindcast_truth_data %>%
      dplyr::distinct(dateID) %>%
      dplyr::arrange(dateID) %>%
      ensure_dates() %>%
      dplyr::mutate(
        trained_date_num = seq.int(from = 1L, length.out = dplyr::n()),
        dates = as.Date(paste0(dateID, "01"), "%Y%m%d")
      )

    # Integrity checks
    if (anyDuplicated(trained_time_map$dateID)) stop("trained_time_map has duplicate dateID")
    if (any(is.na(trained_time_map$trained_date_num))) stop("trained_time_map has NA trained_date_num")

    # Compute seasonal covariates from actual dates
    sin_cos <- get_sin_cos(as.character(trained_time_map$dateID))
    full.ts.model.inputs$sin_mo <- sin_cos$sin
    full.ts.model.inputs$cos_mo <- sin_cos$cos

    # Subset covariate matrices to columns matching trained_time_map dateIDs
    # This ensures positional indexing in create_covariate_samples_fixed reads correct time periods
    date_cols <- as.character(trained_time_map$dateID)
    subset_to_dates <- function(mat, cols) {
      if (is.null(mat)) return(NULL)
      if (!is.matrix(mat)) mat <- as.matrix(mat)
      available_cols <- intersect(cols, colnames(mat))
      if (length(available_cols) == 0) return(mat)  # fallback to original if no column names match
      mat[, available_cols, drop = FALSE]
    }
    full.ts.model.inputs$temp    <- subset_to_dates(full.ts.model.inputs$temp,    date_cols)
    full.ts.model.inputs$mois    <- subset_to_dates(full.ts.model.inputs$mois,    date_cols)
    full.ts.model.inputs$temp_sd <- subset_to_dates(full.ts.model.inputs$temp_sd, date_cols)
    full.ts.model.inputs$mois_sd <- subset_to_dates(full.ts.model.inputs$mois_sd, date_cols)
    full.ts.model.inputs$pH      <- subset_to_dates(full.ts.model.inputs$pH,      date_cols)
    full.ts.model.inputs$pH_sd   <- subset_to_dates(full.ts.model.inputs$pH_sd,   date_cols)
    full.ts.model.inputs$pC      <- subset_to_dates(full.ts.model.inputs$pC,      date_cols)
    full.ts.model.inputs$pC_sd   <- subset_to_dates(full.ts.model.inputs$pC_sd,   date_cols)
    full.ts.model.inputs$relEM   <- subset_to_dates(full.ts.model.inputs$relEM,   date_cols)
    full.ts.model.inputs$LAI     <- subset_to_dates(full.ts.model.inputs$LAI,     date_cols)

    # Expand covariates to match horizon if any matrices are still too short
    total_horizon <- max(trained_time_map$trained_date_num)

    ensure_cols <- function(mat, target_cols) {
      if (is.null(mat)) return(NULL)
      if (!is.matrix(mat)) mat <- as.matrix(mat)
      if (ncol(mat) == 0L) {
        return(matrix(NA_real_, nrow = nrow(mat), ncol = target_cols))
      }
      if (ncol(mat) >= target_cols) return(mat)
      add <- matrix(mat[, 1, drop = FALSE], nrow = nrow(mat), ncol = target_cols - ncol(mat))
      cbind(mat, add)
    }

    full.ts.model.inputs$N.date  <- total_horizon
    full.ts.model.inputs$temp    <- ensure_cols(full.ts.model.inputs$temp,    total_horizon)
    full.ts.model.inputs$mois    <- ensure_cols(full.ts.model.inputs$mois,    total_horizon)
    full.ts.model.inputs$relEM   <- ensure_cols(full.ts.model.inputs$relEM,   total_horizon)
    full.ts.model.inputs$LAI     <- ensure_cols(full.ts.model.inputs$LAI,     total_horizon)
    full.ts.model.inputs$pH      <- ensure_cols(full.ts.model.inputs$pH,      total_horizon)
    full.ts.model.inputs$pC      <- ensure_cols(full.ts.model.inputs$pC,      total_horizon)
    full.ts.model.inputs$pH_sd   <- ensure_cols(full.ts.model.inputs$pH_sd,   total_horizon)
    full.ts.model.inputs$pC_sd   <- ensure_cols(full.ts.model.inputs$pC_sd,   total_horizon)
    full.ts.model.inputs$temp_sd <- ensure_cols(full.ts.model.inputs$temp_sd, total_horizon)
    full.ts.model.inputs$mois_sd <- ensure_cols(full.ts.model.inputs$mois_sd, total_horizon)

    # Plot list for those sites
    plot_site_key <- full.ts.model.inputs$truth.plot.long %>%
      dplyr::select(siteID, plotID) %>%
      dplyr::distinct()
    plot_site_key <- plot_site_key %>% dplyr::filter(siteID %in% unobserved_sites)

    # We need non-colliding plot/site numbers for the new plots
    trained_plot_map <- cal_quants %>%
      dplyr::select(plotID, siteID, plot_num, site_num) %>%
      dplyr::distinct()
    max_plot_num <- suppressWarnings(max(trained_plot_map$plot_num, na.rm = TRUE))
    max_site_num <- suppressWarnings(max(trained_plot_map$site_num, na.rm = TRUE))
    if (!is.finite(max_plot_num)) max_plot_num <- 0L
    if (!is.finite(max_site_num)) max_site_num <- 0L

    # helper to assemble the truth for one plot (hindcast only)
    assemble_unobserved_truth <- function(plotID, siteID) {
      # Get the actual truth data for this plot/site
      site_data <- hindcast_truth_data %>%
        dplyr::filter(plotID == !!plotID, siteID == !!siteID)
      
      if (nrow(site_data) == 0) {
        return(data.frame())
      }
      
      # Use the actual truth data directly - no complex mapping needed
      truth_rows <- site_data %>%
        dplyr::mutate(
          date_num = row_number(),  # Simple sequential numbering
          timepoint = NA_real_
        ) %>%
        dplyr::select(dplyr::any_of(c("siteID", "plotID", "dateID", "dates", "species", "truth", "date_num", "timepoint")))
      
      truth_rows
    }

    site_output_list <- list()

    # CRITICAL FIX: Load and validate pred_effects ONCE per taxon (outside site/plot loops)
    # This must exist - if it doesn't, that's a data integrity problem that must be fixed upstream
    pred_effects_file <- here("data/summary", paste0("site_effects_unobserved_", model_name, ".rds"))
    
    if (!file.exists(pred_effects_file)) {
      stop("CRITICAL ERROR: pred_effects file not found: ", basename(pred_effects_file), 
           "\n  This file must exist for model type '", model_name, "'.",
           "\n  Please generate it using the site effects prediction script (05_predictSiteEffects.r).",
           "\n  Expected path: ", pred_effects_file)
    }
    
    pred_effects_all <- safe_readRDS(pred_effects_file)
    
    # Validate structure
    required_cols <- c("siteID", "fit", "model_id", "taxon")
    missing_cols <- setdiff(required_cols, names(pred_effects_all))
    if (length(missing_cols) > 0) {
      stop("CRITICAL ERROR: pred_effects file missing required columns: ", paste(missing_cols, collapse = ", "),
           "\n  File: ", basename(pred_effects_file),
           "\n  Available columns: ", paste(names(pred_effects_all), collapse = ", "))
    }
    
    # Check for se_fit (required for uncertainty propagation)
    if (!"se_fit" %in% names(pred_effects_all)) {
      warning("WARNING: pred_effects file missing 'se_fit' column for uncertainty propagation. ",
              "This file may need to be regenerated using 05_predictSiteEffects.r. ",
              "Hindcasts will use default uncertainty (0.1) for site effects.")
    }
    
    # CRITICAL FIX: pred_effects files store model_id with _beta_regression suffix
    # The script's model_id doesn't have this suffix, so we need to add it when searching
    model_id_with_suffix <- paste0(model_id, "_beta_regression")
    
    # Filter by model_id (with suffix) and taxon to get only relevant site effects for this specific model
    pred_effects_filt <- pred_effects_all %>%
      dplyr::filter(model_id == !!model_id_with_suffix, taxon == !!taxon)
    
    if (nrow(pred_effects_filt) == 0) {
      # Check if model_id exists but taxon doesn't match
      model_exists <- any(pred_effects_all$model_id == model_id_with_suffix)
      taxon_exists <- any(pred_effects_all$taxon == taxon)
      
      error_msg <- paste0(
        "CRITICAL ERROR: No pred_effects found for model_id '", model_id_with_suffix, "' and taxon '", taxon, "'\n",
        "  File: ", basename(pred_effects_file), "\n",
        "  File contains ", nrow(pred_effects_all), " rows\n",
        "  Unique model_ids in file: ", length(unique(pred_effects_all$model_id)), "\n",
        "  Unique taxa in file: ", length(unique(pred_effects_all$taxon)), "\n"
      )
      
      if (model_exists && !taxon_exists) {
        matching_model_taxa <- unique(pred_effects_all$taxon[pred_effects_all$model_id == model_id_with_suffix])
        error_msg <- paste0(error_msg,
          "  Model_id '", model_id_with_suffix, "' exists but taxon '", taxon, "' does not match.\n",
          "  Available taxa for this model_id: ", paste(matching_model_taxa, collapse = ", "))
      } else if (!model_exists && taxon_exists) {
        matching_taxon_models <- unique(pred_effects_all$model_id[pred_effects_all$taxon == taxon])
        error_msg <- paste0(error_msg,
          "  Taxon '", taxon, "' exists but model_id '", model_id_with_suffix, "' does not match.\n",
          "  Available model_ids for this taxon: ", paste(head(matching_taxon_models, 5), collapse = ", "))
      } else if (!model_exists && !taxon_exists) {
        error_msg <- paste0(error_msg,
          "  Neither model_id nor taxon found in file.\n",
          "  Sample model_ids: ", paste(head(unique(pred_effects_all$model_id), 5), collapse = ", "), "\n",
          "  Sample taxa: ", paste(head(unique(pred_effects_all$taxon), 5), collapse = ", "))
      }
      
      stop(error_msg,
           "\n  This indicates a mismatch between model outputs and site effects predictions.",
           "\n  Script searched for model_id with suffix: '", model_id_with_suffix, "'",
           "\n  Please regenerate pred_effects using 05_predictSiteEffects.r with matching model_ids.")
    }
    
    # Validate we have site effects for all required unobserved sites
    required_sites <- intersect(c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", "WREF", "YELL"), 
                                pred_effects_filt$siteID)
    missing_sites <- setdiff(c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", "WREF", "YELL"), 
                             pred_effects_filt$siteID)
    
    if (length(missing_sites) > 0) {
      warning("WARNING: pred_effects missing site effects for some unobserved sites: ", 
              paste(missing_sites, collapse = ", "),
              "\n  Only ", length(required_sites), " of 14 required sites present.",
              "\n  Missing sites will use random site effects, but this should be fixed upstream.")
    }
    
    dbg("  ✅ Loaded ", nrow(pred_effects_filt), " pred_effects for ", model_id, " (", length(required_sites), " sites)", .on = TRUE)

    # Pre-compute column index once per taxon (avoids repeated grep in extract_all_parameters)
    col_idx <- build_col_index(param_samples, model_name)

    for (siteID in unobserved_sites) {
      plot_list <- unique(plot_site_key$plotID[plot_site_key$siteID == siteID])

      if (length(plot_list) == 0) {
        dbg("  Skipping site", siteID, "- no plots present in truth", .on = TRUE)
        next
      }

      plot_output_list <- list()
    
      for (plotID in plot_list) {
        tryCatch({
          # Assign NEW (non-colliding) plot/site numbers for this new site/plot
          max_plot_num <- max_plot_num + 1L
          max_site_num <- max_site_num + 1L
          assigned_plot_num <- max_plot_num
          assigned_site_num <- max_site_num

          truth_df <- assemble_unobserved_truth(plotID, siteID)

          if (nrow(truth_df) == 0) {
            dbg("    ⚠️  Skipping plot", plotID, "- no hindcast truth rows", .on = TRUE)
            next
          }

          # Validate & coerce truth
          truth_df$truth <- suppressWarnings(as.numeric(truth_df$truth))
          bad <- is.na(truth_df$truth) | truth_df$truth < 0 | truth_df$truth > 1
          if (any(bad)) {
            looks_like_date <- mean(truth_df$truth > 200000, na.rm = TRUE) > 0.2
            stop(sprintf("truth values invalid for %s/%s (min=%.3f, max=%.3f). dateID-like? %s",
                         siteID, plotID,
                         min(truth_df$truth, na.rm = TRUE),
                         max(truth_df$truth, na.rm = TRUE),
                         if (looks_like_date) "YES" else "NO"))
          }

          truth_data <- truth_df %>%
            dplyr::mutate(
            plot_num = assigned_plot_num,
              site_num = assigned_site_num
            ) %>%
            dplyr::select(siteID, plotID, dateID, dates, species, truth,
                          plot_num, date_num, site_num, timepoint)

        # Validate truth_data before hindcast
        truth_data$truth <- suppressWarnings(as.numeric(truth_data$truth))
        req_cols <- c("siteID","plotID","dateID","dates","species","truth","plot_num","date_num","site_num","timepoint")
        missing <- setdiff(req_cols, names(truth_data))
        if (length(missing)) stop("truth_data missing columns: ", paste(missing, collapse=", "))
        if (any(is.na(truth_data$date_num))) stop("truth_data has NA date_num")
        bad <- is.na(truth_data$truth) | truth_data$truth < 0 | truth_data$truth > 1
        if (any(bad)) {
          stop(sprintf("truth_data has invalid truth (min=%.3f, max=%.3f)",
            min(truth_data$truth, na.rm=TRUE), max(truth_data$truth, na.rm=TRUE)))
        }
        stopifnot(is.finite(assigned_plot_num), is.finite(assigned_site_num))
    
    # pred_effects already loaded and validated outside the plot loop (above)
    # 🟢 CALL fcast_logit_beta TWICE - once for modeled effects, once for random effects
    # First call: modeled site effects (predict_site_effects = filtered pred_effects)
    # CRITICAL FIX: Filter pred_effects_filt to only the current siteID to avoid extraction issues
    pred_effects_site <- pred_effects_filt %>%
      dplyr::filter(siteID == !!siteID)
    
    # Validate we have site effects for this specific site
    if (nrow(pred_effects_site) == 0) {
      warning("WARNING: No site effects found for siteID ", siteID, " in pred_effects_filt. ",
              "This site will use random site effects. ",
              "Available siteIDs in pred_effects_filt: ", paste(unique(pred_effects_filt$siteID), collapse=", "))
      pred_effects_site <- NULL  # Will trigger random site effects
    } else if (nrow(pred_effects_site) > 1) {
      warning("WARNING: Multiple site effect rows found for siteID ", siteID, 
              " (", nrow(pred_effects_site), " rows). Using first row.")
      pred_effects_site <- pred_effects_site[1, , drop = FALSE]
    }
    
    tryCatch({
      hindcast.plot.modeled <- fcast_logit_beta(
        plotID,
        full.ts.model.inputs,
        param_samples,
        truth_data,
        plot_summary = cal_quants,
        Nmc = 250,
        predict_site_effects = pred_effects_site,
            rank.name = rank.name,
            model_id = model_id,
            col_idx = col_idx,
            metadata = list(
              species = taxon,
              cal_end_dateID = cal_end_dateID,
              trained_time_map = trained_time_map,
              newsite = TRUE
            )
          )
        }, error = function(e) {
          dbg("  ❌ fcast_logit_beta (modeled) failed for plot", plotID, ":", conditionMessage(e), .on = TRUE)
          hindcast.plot.modeled <- NULL
        })

        # Second call: random site effects (predict_site_effects = FALSE)
        tryCatch({
          hindcast.plot.random <- fcast_logit_beta(
            plotID,
            full.ts.model.inputs,
            param_samples,
            truth_data,
            plot_summary = cal_quants,
            Nmc = 250,
            predict_site_effects = NULL,
            rank.name = rank.name,
            model_id = model_id,
            col_idx = col_idx,
            metadata = list(
              species = taxon,
              cal_end_dateID = cal_end_dateID,
              trained_time_map = trained_time_map,
              newsite = TRUE
            )
          )
        }, error = function(e) {
          dbg("  ❌ fcast_logit_beta (random) failed for plot", plotID, ":", conditionMessage(e), .on = TRUE)
          hindcast.plot.random <- NULL
        })
        
        # Combine both results
        hindcast_results <- list()
        if (!is.null(hindcast.plot.modeled) && nrow(hindcast.plot.modeled) > 0) {
          hindcast.plot.modeled$predicted_site_effect <- TRUE
          hindcast_results[["modeled"]] <- hindcast.plot.modeled
        }
        if (!is.null(hindcast.plot.random) && nrow(hindcast.plot.random) > 0) {
          hindcast.plot.random$predicted_site_effect <- FALSE
          hindcast_results[["random"]] <- hindcast.plot.random
        }
        
        # Process each result type (modeled and random)
        plot_results <- list()
        
        for (result_type in names(hindcast_results)) {
          hindcast.plot <- hindcast_results[[result_type]]
          
          # Add metadata columns
          hindcast.plot <- hindcast.plot %>%
            dplyr::mutate(
              model_name = !!model_name,
              time_period = "20130601_20180101",
              species = !!taxon,
              rank_name = !!rank.name,
              newsite = TRUE,  # keep logical
              newsite_label = "New site",
              model_id = !!model_id,
              site_prediction = ifelse(predicted_site_effect, 
                                       "New time x site (modeled effect)",
                                       "New time x site (random effect)")
            )
          
          # Store result with type identifier
          plot_results[[paste0(plotID, "_", result_type)]] <- hindcast.plot
        }
        
        # Add all results to plot_output_list
        for (result_name in names(plot_results)) {
          plot_output_list[[result_name]] <- plot_results[[result_name]]
        }
        
        }, error = function(e) {
          # Silent error handling
        })
    }
    
    if (length(plot_output_list) > 0) {
      site_data <- rbindlist(plot_output_list, fill = T)
      
      # CRITICAL FIX: Clean up duplicate columns created by rbindlist (optimized)
      duplicate_cols <- grep("\\.(x|y)(\\.|$)", names(site_data), value = TRUE)
      if (length(duplicate_cols) > 0) {
        site_data <- site_data[, !names(site_data) %in% duplicate_cols, with = FALSE]
      }
      
      # CRITICAL FIX: Preserve correct truth values (only check if needed)
      if ("truth" %in% names(site_data)) {
        truth_values <- site_data$truth
        corrupted_mask <- !is.na(truth_values) & truth_values > 10000
        
        if (any(corrupted_mask)) {
          # Try to recover truth values from the taxon column
          if (taxon %in% names(site_data)) {
            site_data$truth[corrupted_mask] <- site_data[[taxon]][corrupted_mask]
          } else {
            site_data$truth[corrupted_mask] <- NA_real_
          }
        }
      }
      
      site_output_list[[siteID]] <- site_data
    }
  }
  
  # Create output directory once (outside loop) and save all site files in batch
  if (length(site_output_list) > 0) {
    # All models go to driver uncertainty directory
    output_dir <- file.path(project_root, "data", "hindcasts", "driver_uncertainty")
    
    # Create directory if it doesn't exist (only once per taxon)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save individual site files in batch
    for (siteID in names(site_output_list)) {
      site_output_file <- file.path(output_dir, paste0("hindcasts_", model_id, "_", siteID, "_unobserved.rds"))
      saveRDS(site_output_list[[siteID]], site_output_file)
    }
  }
  
  # Combine all results for this taxon
  if (length(site_output_list) > 0) {
    tax_output <- rbindlist(site_output_list, fill = T)
    
    # CRITICAL FIX: Clean up duplicate columns created by rbindlist (optimized)
    duplicate_cols <- grep("\\.(x|y)(\\.|$)", names(tax_output), value = TRUE)
    if (length(duplicate_cols) > 0) {
      tax_output <- tax_output[, !names(tax_output) %in% duplicate_cols, with = FALSE]
    }
    
    # CRITICAL FIX: Preserve correct truth values (only check if needed)
    if ("truth" %in% names(tax_output)) {
      truth_values <- tax_output$truth
      corrupted_mask <- !is.na(truth_values) & truth_values > 10000
      
      if (any(corrupted_mask)) {
        # Try to recover truth values from the taxon column
        if (taxon %in% names(tax_output)) {
          tax_output$truth[corrupted_mask] <- tax_output[[taxon]][corrupted_mask]
        } else {
          tax_output$truth[corrupted_mask] <- NA_real_
        }
      }
    }
  } else {
    tax_output <- data.frame()
  }

    # Ensure pretty_group is populated (canonical source: fill_pretty_group from package)
    if (nrow(tax_output) > 0) {
      tax_output <- fill_pretty_group(as.data.table(tax_output))
    }

    if (nrow(tax_output) > 0) {
      # All models go to driver uncertainty summary directory
      summary_dir <- file.path(project_root, "data", "summary", "driver_uncertainty")
      
      # Create directory if it doesn't exist
      dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  
      # Save results for this taxon
      out.path <- file.path(summary_dir, paste0("hindcasts_", taxon, "_unobserved.rds"))
      saveRDS(tax_output, out.path)
      
      # Generate diagnostic figures if enabled
      if (make_figs || is_local_test) {
        tryCatch({
          # Check for required columns
          need_cols <- c("plotID", "siteID", "dates", "med", "lo", "hi", "fcast_period")
          missing_cols <- setdiff(need_cols, names(tax_output))
          
          if (length(missing_cols) == 0) {
            # Create figure directory
            fig_root <- file.path(project_root, "figures", "hindcast_diagnostics", "unobserved_sites", "driver_uncertainty", model_id)
            dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)
            fig_root <- normalizePath(fig_root, winslash = "/", mustWork = TRUE)
            
            # Check for valid numeric data
            valid_data <- tax_output[!is.na(tax_output$med) & !is.na(tax_output$lo) & !is.na(tax_output$hi) & 
                                    is.finite(tax_output$med) & is.finite(tax_output$lo) & is.finite(tax_output$hi), ]
            
            if (nrow(valid_data) >= 5) {
              # Call diagnostic function
              if (exists("generate_hindcast_diagnostics")) {
                generate_hindcast_diagnostics(tax_output, model_id, taxon, out_dir = fig_root)
                cat("  ✅ Diagnostic figures generated for", taxon, "\n")
              } else {
                cat("  ⚠️  generate_hindcast_diagnostics function not found\n")
              }
            } else {
              cat("  ⚠️  Skipping diagnostic figures - insufficient valid data (", nrow(valid_data), " valid rows)\n")
            }
          } else {
            cat("  ⚠️  Skipping diagnostic figures - missing required columns:", paste(missing_cols, collapse = ", "), "\n")
          }
        }, error = function(e) {
          cat("  ⚠️  Error generating diagnostic figures:", conditionMessage(e), "\n")
        })
      }
      
      return(list(success = TRUE, data = tax_output, taxon = taxon, model_id = model_id))
  } else {
      return(list(success = FALSE, error = "No hindcasts generated", taxon = taxon, model_id = model_id))
  }
  
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    dbg("❌ ERROR processing", taxon, model_id, ":", error_msg, .on = TRUE)
    # Capture full error message
    return(list(success = FALSE, error = error_msg, taxon = taxon, model_id = model_id))
  })
}

# Process taxa in parallel with configurable cores (prioritizing env_cycl models)
# Allow override via environment variable, default to 2 cores
n_cores <- as.integer(Sys.getenv("HINDCAST_CORES", "2"))
test_taxa_count <- length(taxa_to_process)
cat("Processing", test_taxa_count, "taxa with", n_cores, "cores\n")

# Set up parallel processing
cat("Setting up parallel processing with", n_cores, "cores\n")
cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

# Export data to workers to avoid redundant loading
cat("Exporting data to parallel workers...\n")
clusterExport(cl, c("all_ranks", "env_data", "taxa_to_process", "project_root", 
                    "process_single_taxon", "resolve_model_paths", "nz_path", "exists_file", 
                    "is_driver_model", "required_unobserved_sites"))
clusterEvalQ(cl, { setwd(project_root); getwd() })   # pin worker wd to project root

# Load required packages in workers
clusterEvalQ(cl, {
  # CRITICAL: Set thread limits in workers to prevent nested parallelization
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1") 
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
  Sys.setenv(NUMEXPR_NUM_THREADS = "1")
  Sys.setenv("OMP_THREAD_LIMIT" = "1")
  
  library(here)
  source(here::here("source.R"))
  library(dplyr)
  # Source functions directly instead of loading package
  if(file.exists(here::here("microbialForecast/R/statsFunctions.r"))) {
    source(here::here("microbialForecast/R/statsFunctions.r"))
  }
  if(file.exists(here::here("microbialForecast/R/helperFunctions.r"))) {
    source(here::here("microbialForecast/R/helperFunctions.r"))
  }
  if(file.exists(here::here("microbialForecast/R/prepBetaRegData.r"))) {
    source(here::here("microbialForecast/R/prepBetaRegData.r"))
  }
  if(file.exists(here::here("microbialForecast/R/run_hindcast.r"))) {
    source(here::here("microbialForecast/R/run_hindcast.r"))
  }
  library(gridExtra)
  library(data.table)
  library(ggplot2)
  library(tidyr)
})

# Process taxa in parallel
cat("Starting parallel processing at", format(Sys.time(), "%H:%M:%S"), "\n")
start_time <- Sys.time()

# Process taxa with better error handling and progress tracking
results <- tryCatch({
  foreach(i = seq_along(taxa_to_process), 
          .packages = c("dplyr", "doParallel", "parallel", "gridExtra", "ggplot2", "data.table", "tidyr", "lubridate", "here"),
          .errorhandling = "pass",
          .export = c("process_single_taxon", "resolve_model_paths", "nz_path", "exists_file", "is_driver_model")) %dopar% {
    
    # Enhanced progress tracking
    if (i %% 10 == 0 || i == 1 || i == length(taxa_to_process)) {
      dbg("=== PROCESSING TAXON", i, "of", length(taxa_to_process), "at", format(Sys.time(), "%H:%M:%S"), "===", .on = TRUE)
    }
    
    # Get the taxon configuration
    taxon_config <- taxa_to_process[[i]]
    
    # Process the taxon with enhanced error handling and memory management
    result <- tryCatch({
      # Force garbage collection before processing each taxon
      gc()
      process_single_taxon(taxon_config, all_ranks, env_data)
    }, error = function(e) {
      cat("❌ ERROR processing", taxon_config$taxon, ":", conditionMessage(e), "\n")
      # Force garbage collection after error
      gc()
      return(list(success = FALSE, error = conditionMessage(e), taxon = taxon_config$taxon, model_id = taxon_config$model_id))
    })
    
    # Log result
    if (result$success) {
      cat("✅ SUCCESS:", result$taxon, "\n")
    } else {
      cat("❌ FAILED:", result$taxon, "-", result$error, "\n")
    }
    
    result
  }
}, error = function(e) {
  cat("ERROR in parallel processing:", e$message, "\n")
  cat("Error traceback:\n")
  print(traceback())
  # Stop the cluster and return empty results
  try(stopCluster(cl), silent = TRUE)
  return(list())
})

# Stop the cluster with error handling
tryCatch({
  stopCluster(cl)
  cat("Cluster stopped successfully\n")
}, error = function(e) {
  cat("Warning: Error stopping cluster:", e$message, "\n")
})

end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Parallel processing completed in", round(elapsed_time, 2), "seconds\n")

# Process results
cat("\n=== PROCESSING RESULTS ===\n")
successful_results <- list()
failed_results <- list()

for (i in seq_along(results)) {
  result <- results[[i]]
  if (is.null(result)) {
    cat("❌ FAILED: NULL RESULT for index", i, "\n")
    failed_results[[length(failed_results) + 1]] <- list(success = FALSE, error = "NULL result", taxon = "unknown")
    next
  }
  if (is.list(result) && "success" %in% names(result) && result$success) {
    successful_results[[length(successful_results) + 1]] <- result
    cat("✅ SUCCESS:", ifelse("taxon" %in% names(result), result$taxon, "unknown"), "\n")
  } else {
    failed_results[[length(failed_results) + 1]] <- result
    error_msg <- if ("error" %in% names(result)) result$error else 
                 if ("message" %in% names(result)) result$message else
                 if (is.character(result)) result else
                 if (inherits(result, "condition")) conditionMessage(result) else
                 "Unknown error"
    taxon_name <- if ("taxon" %in% names(result)) result$taxon else "unknown"
    model_id_name <- if ("model_id" %in% names(result)) result$model_id else "unknown"
    dbg("❌ FAILED:", taxon_name, "(", model_id_name, ") -", substr(error_msg, 1, 200), .on = TRUE)
    # Print full error message if it's longer than 200 chars
    if (nchar(error_msg) > 200) {
      error_lines <- strsplit(error_msg, "\n")[[1]]
      for (line in error_lines) {
        if (nchar(trimws(line)) > 0) {
          dbg("     ", line, .on = TRUE)
        }
      }
    }
  }
}

# Update counters
successful_taxa <- length(successful_results)
failed_taxa <- length(failed_results)

# Populate all_results from successful results
# Use model_id as key (not taxon) to avoid overwriting across model types
all_results <- list()
for (result in successful_results) {
  key <- if (!is.null(result$model_id)) result$model_id else result$taxon
  all_results[[key]] <- result$data
}

# Combine all taxa results
if (length(all_results) > 0) {
  combined_results <- rbindlist(all_results, fill = TRUE)
  
  # All models are driver uncertainty - save to driver uncertainty directory
  combined_path <- here("data", "summary", "driver_uncertainty", "combined_unobserved_hindcasts.rds")
  saveRDS(combined_results, combined_path)
  cat("Combined results saved to:", combined_path, "\n")
  cat("Total rows in combined results:", nrow(combined_results), "\n")
  cat("Number of unique taxa:", length(unique(combined_results$species)), "\n")
  cat("Number of unique sites:", length(unique(combined_results$siteID)), "\n")
  cat("Number of unique plots:", length(unique(combined_results$plotID)), "\n")
}

# Final summary
cat("\n=== FINAL SUMMARY ===\n")
cat("Total taxa available:", nrow(available_models) + skipped_count, "\n")
cat("Taxa with complete files (skipped):", skipped_count, "\n")
cat("Taxa processed:", length(taxa_to_process), "\n")
cat("Successful taxa:", length(successful_results), "\n")
cat("Failed taxa:", length(failed_results), "\n")

if (length(successful_results) > 0) {
  success_rate <- round(100 * length(successful_results) / length(taxa_to_process), 1)
  cat("Success rate:", success_rate, "%\n")
  
  # Combine all successful results
  cat("Combining successful results...\n")
  all_data <- list()
  for (result in successful_results) {
    key <- if (!is.null(result$model_id)) result$model_id else result$taxon
    all_data[[key]] <- result$data
  }
  
  combined_results <- rbindlist(all_data, fill = TRUE)
  
  # CRITICAL FIX: Clean up duplicate columns created by rbindlist
  cat("Cleaning up duplicate columns from rbindlist...\n")
  
  # Remove duplicate columns (.x, .y suffixes) that rbindlist creates
  duplicate_patterns <- c("\\.x$", "\\.y$", "\\.x\\.1$", "\\.y\\.1$")
  duplicate_cols <- character(0)
  
  for (pattern in duplicate_patterns) {
    duplicate_cols <- c(duplicate_cols, grep(pattern, names(combined_results), value = TRUE))
  }
  
  if (length(duplicate_cols) > 0) {
    dbg("Found duplicate columns:", paste(duplicate_cols, collapse=", "), .on = TRUE)
    # Remove duplicate columns
    combined_results <- combined_results[, !names(combined_results) %in% duplicate_cols, with = FALSE]
    dbg("Removed", length(duplicate_cols), "duplicate columns", .on = TRUE)
  }
  
  # CRITICAL FIX: Validate and fix truth column corruption
  if ("truth" %in% names(combined_results)) {
    # Check if truth values are corrupted with dateID values
    truth_values <- combined_results$truth
    corrupted_mask <- !is.na(truth_values) & truth_values > 10000
    
    if (any(corrupted_mask)) {
      dbg("WARNING: Found", sum(corrupted_mask), "corrupted truth values (dateID-like values)", .on = TRUE)
      dbg("Sample corrupted values:", head(truth_values[corrupted_mask]), .on = TRUE)
      
      # Set corrupted truth values to NA
      combined_results$truth[corrupted_mask] <- NA_real_
      dbg("Set corrupted truth values to NA", .on = TRUE)
    }
    
    # Validate remaining truth values are in [0,1] range
    valid_truth_mask <- !is.na(combined_results$truth) & 
                       combined_results$truth >= 0 & 
                       combined_results$truth <= 1
    
    if (any(!valid_truth_mask & !is.na(combined_results$truth))) {
      dbg("WARNING: Found truth values outside [0,1] range", .on = TRUE)
      combined_results$truth[!valid_truth_mask] <- NA_real_
    }
    
    dbg("Final truth data: ", sum(!is.na(combined_results$truth)), "valid values", .on = TRUE)
  }
  
  # Add site_prediction assignment (only if not already set, or fix based on logical newsite)
  if (!"site_prediction" %in% names(combined_results)) {
    combined_results$site_prediction <- dplyr::case_when(
      combined_results$predicted_site_effect == TRUE & combined_results$newsite == TRUE ~ "New time x site (modeled effect)",
      combined_results$predicted_site_effect == FALSE & combined_results$newsite == TRUE ~ "New time x site (random effect)",
      is.na(combined_results$newsite) | combined_results$newsite == FALSE ~ "New time (observed site)",
      TRUE ~ "New time (observed site)"
    )
  } else {
    # Fix any inconsistent site_prediction based on logical newsite and predicted_site_effect
    needs_fix <- !is.na(combined_results$newsite) & combined_results$newsite == TRUE &
                 ((combined_results$predicted_site_effect == TRUE & combined_results$site_prediction != "New time x site (modeled effect)") |
                  (combined_results$predicted_site_effect == FALSE & combined_results$site_prediction != "New time x site (random effect)"))
    if (any(needs_fix, na.rm = TRUE)) {
      combined_results$site_prediction[needs_fix & combined_results$predicted_site_effect == TRUE] <- "New time x site (modeled effect)"
      combined_results$site_prediction[needs_fix & combined_results$predicted_site_effect == FALSE] <- "New time x site (random effect)"
    }
  }
  
  # All models are driver uncertainty - save to driver uncertainty directory
  combined_path <- here("data", "summary", "driver_uncertainty", "combined_unobserved_hindcasts.rds")
  saveRDS(combined_results, combined_path)
  cat("Combined results saved to:", combined_path, "\n")
  cat("Total rows in combined results:", nrow(combined_results), "\n")
  cat("Number of unique taxa:", length(unique(combined_results$species)), "\n")
  cat("Number of unique sites:", length(unique(combined_results$siteID)), "\n")
  cat("Number of unique plots:", length(unique(combined_results$plotID)), "\n")
  
} else {
  cat("❌ No taxa processed successfully\n")
}

cat("\n=== PARALLEL HINDCAST GENERATION COMPLETE ===\n")

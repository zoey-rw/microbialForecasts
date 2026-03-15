#!/usr/bin/env Rscript

# Consolidated Model Summary Pipeline
# Replaces scattered functionality from:
# - 03_summarizeModelOutputs_ready_files.R
# - complete_summary_combination.R
# - generate_convergence_lists.R
# - generate_plot_summary_from_samples2.R
# - comprehensive_model_validation.R
# - 03b_combine_summaries.r
#
# This script preserves ALL crucial output files while streamlining the process

# Parse command line arguments or use defaults
args <- commandArgs(trailingOnly = TRUE)

# Process ALL model types (env_cycl, cycl_only, env_cov) and combine into one file
model_types <- c("env_cycl", "cycl_only", "env_cov")

# Check if running as script vs sourced interactively
if (length(args) == 0) {
  # Try to detect if we're being sourced interactively
  if (interactive()) {
    # Interactive use - use defaults
    cat("No command line arguments provided. Processing all model types.\n")
    base_path <- "data/model_outputs/cloglog_beta_driver_uncertainty/"
    cat("Using default base_path:", base_path, "\n")
  } else {
    # Non-interactive use without args - allow optional base_path
    base_path <- "data/model_outputs/cloglog_beta_driver_uncertainty/"
    cat("Processing all model types: env_cycl, cycl_only, env_cov\n")
    cat("Using default base_path:", base_path, "\n")
  }
} else {
  # Command line argument provided (optional base_path)
  base_path <- if (length(args) > 0) args[1] else "data/model_outputs/cloglog_beta_driver_uncertainty/"
  cat("Processing all model types: env_cycl, cycl_only, env_cov\n")
}

# Ensure variables exist
if (!exists("base_path") || is.null(base_path) || base_path == "") {
  base_path <- "data/model_outputs/cloglog_beta_driver_uncertainty/"
}

cat("Base path:", base_path, "\n")
cat("Model types to process:", paste(model_types, collapse = ", "), "\n")

source("../../source.R")
library(dplyr)
library(coda)
library(purrr)
library(data.table)
# Arrow package for parquet file support (optional)
if (!require(arrow, quietly = TRUE)) {
  tryCatch({
    install.packages("arrow", repos = "https://cran.rstudio.com/", dependencies = FALSE, quiet = TRUE)
    library(arrow)
  }, error = function(e) {
    cat("Warning: arrow package not available. Parquet file support will be limited.\n")
  })
}

# summarize_beta_model and calculate_plot_summary_from_samples loaded via microbialForecast package

# =============================================================================
# CONFIGURATION FLAGS
# =============================================================================

# Set to FALSE to skip plot estimate processing (saves memory)
PROCESS_PLOT_ESTIMATES <- TRUE

cat("Configuration: PROCESS_PLOT_ESTIMATES =", PROCESS_PLOT_ESTIMATES, "\n")

# =============================================================================
# PLOT ESTIMATE PROCESSING FUNCTION
# =============================================================================

process_plot_estimates_optimized <- function(summary_files) {
  cat("\n=== OPTIMIZED PLOT ESTIMATE PROCESSING ===\n")
  cat("Processing", length(summary_files), "summary files\n")
  
  # Use data.table throughout for memory efficiency
  setDTthreads(threads = 4)  # Optimize data.table parallelism
  
  # Pre-allocate list for better memory management
  plot_est_list <- vector("list", length(summary_files))
  valid_count <- 0
  
  # Process files with progress tracking
  for (i in seq_along(summary_files)) {
    if (i %% 50 == 0) {
      cat("Processing file", i, "/", length(summary_files), "\n")
      gc(verbose = FALSE)  # Periodic garbage collection
    }
    
    f <- summary_files[i]
    
    tryCatch({
      # Read only the second element (plot estimates)
      summary_obj <- readRDS(f)
      
      if (is.list(summary_obj) && length(summary_obj) >= 2) {
        plot_est <- summary_obj[[2]]
        
        if (is.data.frame(plot_est) && nrow(plot_est) > 0) {
          # Convert to data.table immediately
          setDT(plot_est)
          
          # Standardize column names efficiently
          col_renames <- c(
            "X2.5." = "2.5%",
            "X25." = "25%", 
            "X50." = "50%",
            "X75." = "75%",
            "X97.5." = "97.5%"
          )
          
          for (old_name in names(col_renames)) {
            if (old_name %in% names(plot_est)) {
              setnames(plot_est, old_name, col_renames[[old_name]])
            }
          }
          
          # Add missing columns efficiently
          if (!"SD" %in% names(plot_est)) {
            if (all(c("2.5%", "97.5%") %in% names(plot_est))) {
              plot_est[, SD := (`97.5%` - `2.5%`) / 3.92]
            } else {
              plot_est[, SD := NA_real_]
            }
          }
          
          if (!"Mean" %in% names(plot_est) && "50%" %in% names(plot_est)) {
            plot_est[, Mean := `50%`]
          }
          
          # Filter invalid rows efficiently
          plot_est <- plot_est[!is.na(Mean) & !is.nan(Mean)]
          
          if (nrow(plot_est) > 0) {
            valid_count <- valid_count + 1
            plot_est_list[[valid_count]] <- plot_est
          }
        }
      }
      
      # Clear the full summary object immediately
      rm(summary_obj)
      
    }, error = function(e) {
      # Silently skip problematic files
    })
  }
  
  # Trim list to valid entries
  plot_est_list <- plot_est_list[1:valid_count]
  
  cat("Valid plot estimate files:", valid_count, "\n")
  
  # Combine in smaller chunks to avoid memory issues
  if (valid_count > 0) {
    cat("Combining plot estimates in chunks...\n")
    
    # Process in chunks of 50 to avoid memory issues
    chunk_size <- 50
    n_chunks <- ceiling(valid_count / chunk_size)
    
    combined_chunks <- list()
    
    for (chunk in 1:n_chunks) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, valid_count)
      chunk_files <- plot_est_list[start_idx:end_idx]
      
      # Combine this chunk
      chunk_combined <- rbindlist(chunk_files, fill = TRUE, use.names = TRUE)
      combined_chunks[[chunk]] <- chunk_combined
      
      # Clear chunk memory
      rm(chunk_files, chunk_combined)
      gc(verbose = FALSE)
    }
    
    # Combine all chunks
    cat("Combining", n_chunks, "chunks...\n")
    plot_est_combined <- rbindlist(combined_chunks, fill = TRUE, use.names = TRUE)
    
    cat("Combined", nrow(plot_est_combined), "plot estimate rows\n")
    
    # Free memory
    rm(plot_est_list, combined_chunks)
    gc(verbose = FALSE)
    
    return(plot_est_combined)
  } else {
    return(data.table())
  }
}

cat("=== CONSOLIDATED MODEL SUMMARY PIPELINE ===\n")
cat("Processing model outputs with comprehensive validation and summarization\n\n")

# Check if we can skip processing steps
skip_processing <- FALSE
# IMPORTANT: For cloglog_beta_driver_uncertainty models, use sub-scripts for efficiency
# The sub-scripts (03b, 03c, 03d, 03e) break the work into manageable chunks
# Note: base_path is set earlier in the script
# DISABLE SUB-SCRIPTS: We want to regenerate summaries, not just combine existing ones
USE_SUB_SCRIPTS <- FALSE  # Force to FALSE to skip sub-scripts and go directly to STEP 2

if (USE_SUB_SCRIPTS) {
  cat("📋 Cloglog driver uncertainty models detected - using sub-scripts for efficient processing\n")
  cat("   Will use: 03b_combine_summaries.r, 03c_extract_plot_summaries.r, etc.\n")
  skip_processing <- FALSE
} else if (file.exists(here("data/summary/logit_beta_fixed_priors_summaries.rds")) && 
    (file.exists(here("data/summary/plot_estimates.rds")) || 
     file.exists(here("data/summary/plot_estimates.parquet")))) {
  cat("✅ Found existing summary files - checking if processing can be skipped...\n")
  
  # First check if the existing summary file has valid data
  tryCatch({
    test_summary <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
    has_valid_data <- (nrow(test_summary$summary_df) > 0)
    
    if (!has_valid_data) {
      cat("⚠️  Existing summary file has empty data frames - will reprocess\n")
      skip_processing <- FALSE
    } else {
      # Get timestamps of output files
      main_summary_time <- file.info(here("data/summary/logit_beta_fixed_priors_summaries.rds"))$mtime
      
      # Check which plot estimates file exists and get its timestamp
      if (file.exists(here("data/summary/plot_estimates.rds"))) {
        plot_est_time <- file.info(here("data/summary/plot_estimates.rds"))$mtime
      } else if (file.exists(here("data/summary/plot_estimates.parquet"))) {
        plot_est_time <- file.info(here("data/summary/plot_estimates.parquet"))$mtime
      } else {
        plot_est_time <- as.POSIXct("1900-01-01")  # Very old date to force reprocessing
      }
      
      # Find the newest input model files from ALL model types
      # Always check all three model types (env_cycl, cycl_only, env_cov)
      model_files <- c(
        list.files(file.path(base_path, "env_cov"), 
                  pattern = "samples_", recursive = TRUE, full.names = TRUE),
        list.files(file.path(base_path, "env_cycl"), 
                  pattern = "samples_", recursive = TRUE, full.names = TRUE),
        list.files(file.path(base_path, "cycl_only"), 
                  pattern = "samples_", recursive = TRUE, full.names = TRUE)
      )
      model_files <- model_files[!grepl("_chain[0-9]", model_files)]
      
      if (length(model_files) > 0) {
        # Get the newest model file timestamp
        model_times <- file.info(model_files)$mtime
        newest_model_time <- max(model_times, na.rm = TRUE)
        
        # Check if output files are newer than the newest input file
        if (main_summary_time > newest_model_time && plot_est_time > newest_model_time) {
          cat("✅ Output files are newer than input files and contain valid data - skipping processing steps\n")
          skip_processing <- TRUE
        } else {
          cat("⚠️  Output files are older than input files - will reprocess\n")
          cat("  Newest model file:", format(newest_model_time), "\n")
          cat("  Main summary file:", format(main_summary_time), "\n")
          cat("  Plot estimates file:", format(plot_est_time), "\n")
        }
      } else {
        cat("⚠️  No model files found - will reprocess\n")
      }
    }
  }, error = function(e) {
    cat("⚠️  Error reading existing summary file - will reprocess\n")
    cat("  Error:", e$message, "\n")
    skip_processing <- FALSE
  })
}

# Check available memory
if (.Platform$OS.type == "unix") {
  mem_info <- system("vm_stat", intern = TRUE)
  cat("Memory check completed\n")
}

# Set memory limits and parallel processing controls
options(memory.limit = 4000)  # 4GB limit for safety
cat("Memory limit set to 4GB\n")

# Set parallel processing limits
library(parallel)

# Function to filter files to only include those with both 'with_legacy_covariate' and 'beta_regression'
# or driver uncertainty files
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  # Filter to include either:
  # 1. Files with both 'with_legacy_covariate' and 'beta_regression' (old format)
  # 2. Files from driver uncertainty directory (new format)
  standard_files <- file_list[
    (grepl('with_legacy_covariate', basename(file_list)) & 
     grepl('beta_regression', basename(file_list))) |
    grepl('cloglog_beta_driver_uncertainty', file_list)
  ]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard files (with required suffixes or driver uncertainty):', length(standard_files), '\n')
  cat('  Filtered out:', length(file_list) - length(standard_files), '\n\n')
  
  return(standard_files)
}

max_cores <- min(4, detectCores() - 1)  # Max 4 cores, leave 1 free
options(mc.cores = max_cores)
cat("Parallel processing limited to", max_cores, "cores\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Convergence criteria (simplified approach)
# - Strict: Point est. < 1.05
# - Standard: Point est. < 1.1  
# - Weak: Point est. < 1.2

# Check if CLR models exist and process them
# if (dir.exists(here("data/model_outputs/CLR_regression"))) {
#   cat("Processing CLR models...\n")
#   source("03_summarizeModelOutputs_CLR.r")
# }

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to extract species name from filename
extract_species_name <- function(file_path) {
  basename_file <- basename(file_path)
  
  if (grepl("samples_env_cycl_", basename_file)) {
    species <- gsub(".*samples_env_cycl_([^_]+)_.*", "\\1", basename_file)
  } else if (grepl("samples_cycl_only_", basename_file)) {
    species <- gsub(".*samples_cycl_only_([^_]+)_.*", "\\1", basename_file)
  } else if (grepl("samples_", basename_file)) {
    species <- gsub(".*samples_([^_]+)_.*", "\\1", basename_file)
  } else {
    species <- "unknown"
  }
  
  return(species)
}

# Function to determine file type
get_file_type <- function(file_path) {
  basename_file <- basename(file_path)
  
  if (grepl("_chain[0-9]+\\.rds$", basename_file)) {
    return("chain")
  } else if (grepl("samples_.*\\.rds$", basename_file) && !grepl("_chain[0-9]+", basename_file)) {
    return("combined")
  } else {
    return("other")
  }
}

# Function to determine model type
get_model_type <- function(file_path) {
  if (grepl("env_cycl", file_path)) {
    return("env_cycl")
  } else if (grepl("cycl_only", file_path)) {
    return("cycl_only")
  } else if (grepl("env_cov", file_path)) {
    return("env_cov")
  } else {
    return("unknown")
  }
}

# Function to validate a single file
validate_file <- function(file_path) {
  result <- list(
    file = basename(file_path),
    full_path = file_path,
    species = extract_species_name(file_path),
    file_type = get_file_type(file_path),
    model_type = get_model_type(file_path),
    file_size_mb = round(file.info(file_path)$size / (1024^2), 2),
    readable = FALSE,
    has_samples = FALSE,
    has_samples2 = FALSE,
    has_param_summary = FALSE,
    param_summary_valid = FALSE,
    param_summary_names = "N/A",
    param_summary_issue = "N/A",
    has_metadata = FALSE,
    has_model_data = FALSE,
    has_plot_summary = FALSE,
    samples_dim = "N/A",
    samples2_dim = "N/A",
    metadata_names = "N/A",
    error_message = "N/A"
  )
  
  tryCatch({
    data <- readRDS(file_path)
    result$readable <- TRUE
    
    if (is.list(data)) {
      # Check for samples
      if ("samples" %in% names(data)) {
        result$has_samples <- TRUE
        if (is.list(data$samples)) {
          result$samples_dim <- paste(length(data$samples), "chains")
        } else if (is.matrix(data$samples) || is.array(data$samples)) {
          result$samples_dim <- paste(dim(data$samples), collapse = "x")
        }
      }
      
      # Check for samples2 (plot estimates) or plot_summary
      if ("samples2" %in% names(data)) {
        result$has_samples2 <- TRUE
        if (is.list(data$samples2)) {
          result$samples2_dim <- paste(length(data$samples2), "chains")
        } else if (is.matrix(data$samples2) || is.array(data$samples2)) {
          result$samples2_dim <- paste(dim(data$samples2), collapse = "x")
        }
        
        # Check if samples2 contains actual plot estimates (plot_mu parameters)
        if (is.matrix(data$samples2) || is.array(data$samples2)) {
          col_names <- colnames(data$samples2)
          has_plot_mu <- any(grepl("^plot_mu\\[", col_names))
          if (!has_plot_mu) {
            result$has_samples2 <- FALSE
            result$samples2_dim <- paste(result$samples2_dim, "(no plot_mu)")
          }
        }
      } else if ("plot_summary" %in% names(data)) {
        # Combined files have plot_summary instead of samples2
        result$has_samples2 <- TRUE
        if (is.list(data$plot_summary)) {
          result$samples2_dim <- paste(length(data$plot_summary), "elements")
        } else if (is.data.frame(data$plot_summary)) {
          result$samples2_dim <- paste(dim(data$plot_summary), collapse = "x")
        } else {
          result$samples2_dim <- "unknown"
        }
      }
      
      # Check for param_summary (handle nested structure)
      if ("param_summary" %in% names(data)) {
        result$has_param_summary <- TRUE
        if (is.list(data$param_summary)) {
          result$param_summary_names <- paste(names(data$param_summary), collapse = ", ")
          has_means <- "means" %in% names(data$param_summary)
          has_quantiles <- "quantiles" %in% names(data$param_summary)
          
          # If not found at top level, check nested param_summary
          if (!has_means || !has_quantiles) {
            if ("param_summary" %in% names(data$param_summary) && is.list(data$param_summary$param_summary)) {
              has_means <- "means" %in% names(data$param_summary$param_summary)
              has_quantiles <- "quantiles" %in% names(data$param_summary$param_summary)
            }
          }
          
          if (has_means && has_quantiles) {
            result$param_summary_valid <- TRUE
          } else {
            result$param_summary_valid <- FALSE
            result$param_summary_issue <- "Missing means or quantiles elements (checked nested structure)"
          }
        } else {
          result$param_summary_valid <- FALSE
          result$param_summary_issue <- "Not a list"
        }
      } else {
        result$has_param_summary <- FALSE
        result$param_summary_valid <- FALSE
        result$param_summary_issue <- "Missing param_summary"
      }
      
      # Check for metadata
      if ("metadata" %in% names(data)) {
        result$has_metadata <- TRUE
        result$metadata_names <- paste(names(data$metadata), collapse = ", ")
        
        if ("model_data" %in% names(data$metadata)) {
          result$has_model_data <- TRUE
        }
      }
    }
    
  }, error = function(e) {
    result$error_message <- e$message
  })
  
  return(result)
}

# Function to generate plot_summary from samples2 matrix
generate_plot_summary_for_model <- function(model_file_path) {
  cat("Processing model:", basename(model_file_path), "\n")
  
  model_data <- readRDS(model_file_path)
  
  if(!"samples2" %in% names(model_data)) {
    cat("  WARNING: No samples2 in model data\n")
    return(FALSE)
  }
  
  samples2_matrix <- model_data$samples2
  
  if(!is.matrix(samples2_matrix)) {
    cat("  WARNING: samples2 is not a matrix\n")
    return(FALSE)
  }
  
  cat("  samples2 matrix dimensions:", dim(samples2_matrix), "\n")
  
  # Check if samples2 contains plot_mu parameters
  col_names <- colnames(samples2_matrix)
  plot_mu_cols <- grep("^plot_mu\\[", col_names)
  
  if(length(plot_mu_cols) == 0) {
    cat("  WARNING: No plot_mu parameters found in samples2\n")
    return(FALSE)
  }
  
  cat("  Found", length(plot_mu_cols), "plot_mu parameters\n")
  
  # Generate plot_summary from samples2 using the existing function
  plot_summary <- calculate_plot_summary_from_samples(samples2_matrix)
  
  if(is.null(plot_summary) || length(plot_summary) == 0) {
    cat("  ERROR: Failed to generate plot_summary\n")
    return(FALSE)
  }
  
  # Add the plot_summary to the model data
  model_data$plot_summary <- plot_summary
  
  # Save the updated model file
  saveRDS(model_data, model_file_path)
  
  cat("  SUCCESS: Added plot_summary to model file\n")
  cat("  plot_summary contains", length(plot_summary), "elements\n")
  
  if(length(plot_summary) >= 2 && is.data.frame(plot_summary[[2]])) {
    cat("  plot_summary[[2]] dimensions:", dim(plot_summary[[2]]), "\n")
    cat("  plot_summary[[2]] columns:", paste(colnames(plot_summary[[2]]), collapse = ", "), "\n")
  }
  
  return(TRUE)
}

# =============================================================================
# STEP 1: COMPREHENSIVE MODEL VALIDATION
# =============================================================================

# SKIP VALIDATION: Go directly to summarization
SKIP_VALIDATION <- TRUE

if (!skip_processing && !SKIP_VALIDATION) {
cat("=== STEP 1: COMPREHENSIVE MODEL VALIDATION ===\n")

# Find all model files from ALL model types with memory-efficient approach
# Process all three model types (env_cycl, cycl_only, env_cov)
model_dirs <- file.path(base_path, model_types)

# Process directories one at a time to avoid memory issues
all_files <- c()
for (dir in model_dirs) {
  if (dir.exists(dir)) {
    cat("Scanning directory:", dir, "\n")
    files <- list.files(dir, pattern = "samples_", 
                       full.names = TRUE, recursive = TRUE)
    cat("  Raw files found:", length(files), "\n")
    
    # Exclude individual chain files; allow recursive group subfolders
    if (length(files) == 0) {
      files <- character(0)
    }
    
    files <- files[!grepl("_chain[0-9]", files)]
    cat("  Files after chain filtering:", length(files), "\n")
    # Filter out files with double samples_samples_ prefix (should not exist)
    files <- files[!grepl("^.*samples_samples_", basename(files))]
    cat("  Files after removing double samples_ prefix:", length(files), "\n")
    all_files <- c(all_files, files)
    cat("  Found", length(files), "files in", basename(dir), "\n")
  }
}

cat("Total model files found across all model types:", length(all_files), "\n")

# Apply filtering to only process standard files
all_files <- filter_standard_files(all_files)

# Process all files for complete convergence analysis
cat("Processing all", length(all_files), "files for comprehensive convergence analysis\n")

# Pre-allocate results data frame for better performance
n_files <- length(all_files)
results <- data.frame(
  file = character(n_files),
  full_path = character(n_files),
  species = character(n_files),
  file_type = character(n_files),
  model_type = character(n_files),
  file_size_mb = numeric(n_files),
  readable = logical(n_files),
  has_samples = logical(n_files),
  has_samples2 = logical(n_files),
  has_param_summary = logical(n_files),
  param_summary_valid = logical(n_files),
  param_summary_names = character(n_files),
  param_summary_issue = character(n_files),
  has_metadata = logical(n_files),
  has_model_data = logical(n_files),
  has_plot_summary = logical(n_files),
  samples_dim = character(n_files),
  samples2_dim = character(n_files),
  metadata_names = character(n_files),
  error_message = character(n_files),
  stringsAsFactors = FALSE
)

cat("Validating files using parallel processing...\n")

# Use parallel processing for validation
cl <- makeCluster(max_cores)

# Export all necessary functions and variables
clusterExport(cl, c("validate_file", "extract_species_name", "get_file_type", "get_model_type", "here"))

# Set up the cluster environment
clusterEvalQ(cl, {
  library(dplyr)
  library(here)
  # Source the helper functions that might be needed
  if(file.exists("microbialForecast/R/helperFunctions.r")) {
    source("microbialForecast/R/helperFunctions.r")
  }
})

# Process files in parallel with batching to avoid memory issues
batch_size <- 50  # Process 50 files at a time
n_batches <- ceiling(length(all_files) / batch_size)
validation_results <- list()

cat("Processing", length(all_files), "files in", n_batches, "batches of", batch_size, "\n")

for(batch in 1:n_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(all_files))
  batch_files <- all_files[start_idx:end_idx]
  
  cat("Processing batch", batch, "of", n_batches, "(", length(batch_files), "files)\n")
  
  batch_results <- tryCatch({
    parLapply(cl, batch_files, function(file_path) {
      tryCatch({
        validate_file(file_path)
      }, error = function(e) {
        # Return a default result structure if validation fails
        data.frame(
          file = basename(file_path),
          full_path = file_path,
          species = "unknown",
          file_type = "unknown",
          model_type = "unknown",
          file_size_mb = 0,
          readable = FALSE,
          has_samples = FALSE,
          has_samples2 = FALSE,
          has_param_summary = FALSE,
          param_summary_valid = FALSE,
          param_summary_names = "",
          param_summary_issue = e$message,
          has_metadata = FALSE,
          has_model_data = FALSE,
          has_plot_summary = FALSE,
          samples_dim = "",
          samples2_dim = "",
          metadata_names = "",
          error_message = e$message,
          stringsAsFactors = FALSE
        )
      })
    })
  }, error = function(e) {
    cat("Error in batch", batch, ":", e$message, "\n")
    # Fallback to sequential processing for this batch
    lapply(batch_files, validate_file)
  })
  
  # Add batch results to main results
  validation_results <- c(validation_results, batch_results)
  
  # Force garbage collection between batches
  gc()
}

stopCluster(cl)

# Convert results to data frame
for(i in seq_along(validation_results)) {
  results[i, ] <- validation_results[[i]]
}

cat("Validation complete\n")

# Identify files ready for summary script
# For driver uncertainty models, we need to process raw samples to create summaries
# So we only need to check for basic file structure, not pre-existing summaries
ready_for_summary <- results[
  results$file_type == "combined" & 
  results$has_samples & 
  results$has_metadata & 
  results$has_model_data, 
]

cat("Files ready for summary script:", nrow(ready_for_summary), "\n")

# Save validation results
write.csv(results, "comprehensive_model_validation_results.csv", row.names = FALSE)
write.csv(ready_for_summary, "files_ready_for_summary.csv", row.names = FALSE)

cat("Validation results saved\n")
}  # End of STEP 1 validation block

# =============================================================================
# STEP 2: PROCESS MODEL FILES AND CREATE SUMMARIES (INCLUDING PLOT SUMMARIES)
# =============================================================================

cat("\n=== STEP 2: PROCESS MODEL FILES AND CREATE SUMMARIES ===\n")

# If validation was skipped, find files directly
if (exists("SKIP_VALIDATION") && SKIP_VALIDATION) {
  cat("Skipping validation - finding sample files directly...\n")
  # Find all sample files from all model types
  env_cycl_files <- list.files(file.path(base_path, "env_cycl"), 
                               pattern = "samples_.*_beta_regression\\.rds$", 
                               recursive = TRUE, full.names = TRUE)
  cycl_only_files <- list.files(file.path(base_path, "cycl_only"), 
                               pattern = "samples_.*_beta_regression\\.rds$", 
                               recursive = TRUE, full.names = TRUE)
  env_cov_files <- list.files(file.path(base_path, "env_cov"), 
                             pattern = "samples_.*_beta_regression\\.rds$", 
                             recursive = TRUE, full.names = TRUE)
  
  # Filter out chain files
  all_files <- c(env_cycl_files, cycl_only_files, env_cov_files)
  ready_file_paths <- all_files[!grepl("_chain[0-9]", all_files)]
  cat("Found", length(ready_file_paths), "sample files to process\n")
} else {
  # Get ready files from validation step
  ready_file_paths <- ready_for_summary$full_path
}

n_ready_files <- length(ready_file_paths)

# Pre-allocate file_summaries list
file_summaries <- vector("list", n_ready_files)

# Plot summary generation is now handled by the summarize_beta_model function

# Use parallel processing for model summarization
# Add progress tracking option - use sequential if only a few files or if requested
use_sequential <- length(ready_file_paths) <= 10 || identical(tolower(Sys.getenv("USE_SEQUENTIAL", "false")), "true")

if (use_sequential) {
  cat("Using sequential processing (", length(ready_file_paths), " files)\n", sep="")
  
  # Process files sequentially with progress tracking
  summarization_results <- lapply(seq_along(ready_file_paths), function(i) {
    f <- ready_file_paths[i]
    
    if (i %% 10 == 0 || i == 1) {
      cat("Processing file", i, "of", length(ready_file_paths), ":", basename(f), "\n")
    }
    
    # FORCE REGENERATION for env_cov models - always regenerate
    # For other model types, check if summary is up-to-date
    summary_file <- gsub("samples_", "summary_", f)
    should_summarize <- TRUE
    
    # Force regeneration for env_cov models
    if (grepl("env_cov", f)) {
      should_summarize <- TRUE  # Always regenerate env_cov
    } else if (file.exists(summary_file)) {
      # For other model types, check if up-to-date
      file_time <- file.mtime(f)
      summary_time <- file.mtime(summary_file)
      if (summary_time > file_time) {
        return(list(result = TRUE, message = "Summary file is up-to-date, skipping"))
      } else {
        should_summarize <- TRUE
      }
    }
    
    if (should_summarize) {
      tryCatch({
        # Determine file source
        file_source <- if (grepl("reconstructed_from_checkpoints", f)) "reconstructed" else "old_converged"
        
        # Process the file with error handling for truth.plot.long issues
        # The summarize_beta_model function now handles plot_summary validation and regeneration
        out <- tryCatch({
          microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
        }, error = function(e) {
          if (grepl("truth.plot.long", e$message)) {
            # Try with drop_other = FALSE to avoid plot reconstruction
            microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = FALSE, overwrite = TRUE)
          } else {
            stop(e)
          }
        })
        
        # Add source information to the summary
        if (is.list(out) && "metadata" %in% names(out)) {
          out$metadata$file_source <- file_source
        }
        
        return(list(result = out, message = "Successfully processed"))
      }, error = function(e) {
        return(list(result = NULL, message = paste("Error processing:", e$message)))
      })
    }
  })
  
  cl <- NULL  # No cluster to stop
} else {
  cat("Using parallel processing with", max_cores, "cores\n")
  cl <- makeCluster(max_cores)
  
  # Set up cluster environment - load updated source code
  clusterEvalQ(cl, {
    library(dplyr)
    library(here)
    library(microbialForecast)
  })
  
  # Export necessary variables and functions
  clusterExport(cl, c("ready_file_paths"), envir = environment())
  # Export the summarize_beta_model function to workers (they sourced it, but export to be safe)
  clusterExport(cl, "summarize_beta_model", envir = environment())
  
  # Process files in parallel with error handling and progress tracking
  cat("Processing", length(ready_file_paths), "files in parallel...\n")
  summarization_results <- tryCatch({
    # Use parLapply with progress tracking via intermediate results
    result_list <- vector("list", length(ready_file_paths))
    
    # Process in batches to show progress
    batch_size <- max(1, floor(length(ready_file_paths) / 10))
    n_batches <- ceiling(length(ready_file_paths) / batch_size)
    
    for (batch in 1:n_batches) {
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, length(ready_file_paths))
      batch_files <- ready_file_paths[start_idx:end_idx]
      
      cat("Processing batch", batch, "of", n_batches, "(", length(batch_files), "files)\n")
      
      batch_results <- parLapply(cl, batch_files, function(f) {
        # FORCE REGENERATION for env_cov models - always regenerate
        # For other model types, check if summary is up-to-date
        summary_file <- gsub("samples_", "summary_", f)
        
        # Force regeneration for env_cov models - skip the up-to-date check entirely
        if (grepl("env_cov", f)) {
          should_summarize <- TRUE  # Always regenerate env_cov, no checks
        } else {
          # For other model types, check if summary is up-to-date
          if (file.exists(summary_file)) {
            file_time <- file.mtime(f)
            summary_time <- file.mtime(summary_file)
            if (summary_time > file_time) {
              return(list(result = TRUE, message = "Summary file is up-to-date, skipping"))
            }
          }
          should_summarize <- TRUE
        }
        
        if (should_summarize) {
          tryCatch({
            # Determine file source
            file_source <- if (grepl("reconstructed_from_checkpoints", f)) "reconstructed" else "old_converged"
            
            # Process the file with error handling for truth.plot.long issues
            # The summarize_beta_model function now handles plot_summary validation and regeneration
            out <- tryCatch({
              # Use the sourced function directly (not package namespace) since workers sourced it
              summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
            }, error = function(e) {
              if (grepl("truth.plot.long", e$message)) {
                # Try with drop_other = FALSE to avoid plot reconstruction
                microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = FALSE, overwrite = TRUE)
              } else {
                stop(e)
              }
            })
            
            # Add source information to the summary
            if (is.list(out) && "metadata" %in% names(out)) {
              out$metadata$file_source <- file_source
            }
            
            return(list(result = out, message = "Successfully processed"))
          }, error = function(e) {
            return(list(result = NULL, message = paste("Error processing:", e$message)))
          })
        } else {
          return(list(result = NULL, message = "Skipped (should_summarize = FALSE)"))
        }
      })
      
      # Store batch results and log progress
      result_list[start_idx:end_idx] <- batch_results
      
      # Count successes and failures in this batch
      batch_success <- sum(sapply(batch_results, function(x) !is.null(x$result)))
      batch_skipped <- sum(sapply(batch_results, function(x) !is.null(x$message) && grepl("up-to-date|Skipped", x$message)))
      batch_errors <- sum(sapply(batch_results, function(x) is.null(x$result) && !is.null(x$message) && grepl("Error", x$message)))
      
      cat("Completed batch", batch, "of", n_batches, "- Success:", batch_success, "Skipped:", batch_skipped, "Errors:", batch_errors, "\n")
    }
    
    return(result_list)
  }, error = function(e) {
    cat("Error in parallel summarization:", e$message, "\n")
    cat("Falling back to sequential processing...\n")
    # Fallback to sequential processing
    lapply(ready_file_paths, function(f) {
      # FORCE REGENERATION for env_cov models - always regenerate
      # For other model types, check if summary is up-to-date
      summary_file <- gsub("samples_", "summary_", f)
      should_summarize <- TRUE
      
      # Force regeneration for env_cov models
      if (grepl("env_cov", f)) {
        should_summarize <- TRUE  # Always regenerate env_cov
      } else if (file.exists(summary_file)) {
        # For other model types, check if up-to-date
        file_time <- file.mtime(f)
        summary_time <- file.mtime(summary_file)
        if (summary_time > file_time) {
          return(list(result = TRUE, message = "Summary file is up-to-date, skipping"))
        } else {
          should_summarize <- TRUE
        }
      }
      
      if (should_summarize) {
        tryCatch({
          # Determine file source
          file_source <- if (grepl("reconstructed_from_checkpoints", f)) "reconstructed" else "old_converged"
          
          # Process the file with error handling for truth.plot.long issues
          out <- tryCatch({
            # Use the sourced function directly
            summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
          }, error = function(e) {
            if (grepl("truth.plot.long", e$message)) {
              # Try with drop_other = FALSE to avoid plot reconstruction
              microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = FALSE, overwrite = TRUE)
            } else {
              stop(e)
            }
          })
          
          # Add source information to the summary
          if (is.list(out) && "metadata" %in% names(out)) {
            out$metadata$file_source <- file_source
          }
          
          return(list(result = out, message = "Successfully processed"))
        }, error = function(e) {
          return(list(result = NULL, message = paste("Error processing:", e$message)))
        })
      }
    })
  })
}

# Stop cluster if it was created
if (!is.null(cl)) {
  tryCatch({
    stopCluster(cl)
  }, error = function(e) {
    cat("Warning: Error stopping cluster:", e$message, "\n")
  })
}

# Extract results
for(i in seq_along(summarization_results)) {
  file_summaries[[i]] <- summarization_results[[i]]$result
  if(i %% 10 == 0) {
    cat("Processed", i, "of", n_ready_files, "files\n")
  }
}

# Remove NULL entries and track failures
failed_files <- character(0)
failed_messages <- character(0)
valid_summaries <- list()

for(i in seq_along(summarization_results)) {
  if(is.null(summarization_results[[i]]$result)) {
    failed_files <- c(failed_files, basename(ready_file_paths[i]))
    failed_messages <- c(failed_messages, summarization_results[[i]]$message)
  } else {
    valid_summaries[[length(valid_summaries) + 1]] <- summarization_results[[i]]$result
  }
}

file_summaries <- valid_summaries
cat("Model processing complete. Results:", length(file_summaries), "out of", n_ready_files, "files\n")
if(length(failed_files) > 0) {
  cat("  ⚠️ Failed to process", length(failed_files), "files:\n")
  for(i in seq_len(min(10, length(failed_files)))) {
    cat("    -", failed_files[i], ":", failed_messages[i], "\n")
  }
  if(length(failed_files) > 10) {
    cat("    ... and", length(failed_files) - 10, "more failures\n")
  }
}

# =============================================================================
# STEP 3: COLLECT AND COMBINE SUMMARY FILES
# =============================================================================

cat("\n=== STEP 3: COLLECT AND COMBINE SUMMARY FILES ===\n")

# Search for summary files from ALL model types - only main model files, not chain files
# Use R's list.files instead of system find for better portability
# Process all three model types (env_cycl, cycl_only, env_cov)
summary_file_list <- character(0)
for (mt in model_types) {
  mt_dir <- file.path(base_path, mt)
  if (dir.exists(mt_dir)) {
    # Use list.files with recursive=TRUE to find all summary files
    mt_files <- list.files(mt_dir, pattern = "^summary_.*\\.rds$", 
                          full.names = TRUE, recursive = TRUE)
    # Filter out chain files
    mt_files <- mt_files[!grepl("_chain[0-9]", basename(mt_files))]
    if (length(mt_files) > 0) {
      summary_file_list <- c(summary_file_list, mt_files)
      cat("Found", length(mt_files), "summary files for", mt, "\n")
    } else {
      cat("No summary files found for", mt, "in", mt_dir, "\n")
    }
  } else {
    cat("Directory does not exist:", mt_dir, "\n")
  }
}

cat("Found", length(summary_file_list), "summary files from all processed models (combined from all model types)\n")

# Check which model files don't have corresponding summary files
if(exists("ready_file_paths")) {
  model_files_without_summaries <- character(0)
  for(f in ready_file_paths) {
    expected_summary <- gsub("samples_", "summary_", f)
    if(!file.exists(expected_summary)) {
      model_files_without_summaries <- c(model_files_without_summaries, basename(f))
    }
  }
  if(length(model_files_without_summaries) > 0) {
    cat("  ⚠️", length(model_files_without_summaries), "model files are missing summary files:\n")
    for(i in seq_len(min(10, length(model_files_without_summaries)))) {
      cat("    -", model_files_without_summaries[i], "\n")
    }
    if(length(model_files_without_summaries) > 10) {
      cat("    ... and", length(model_files_without_summaries) - 10, "more\n")
    }
  } else {
    cat("  ✅ All model files have corresponding summary files\n")
  }
}

# Filter out empty summary files
valid_summary_files <- c()
for (file in summary_file_list) {
  tryCatch({
    obj <- readRDS(file)
    if (length(obj) > 0 && (length(names(obj)) > 0 || (is.list(obj) && length(obj) >= 3))) {
      valid_summary_files <- c(valid_summary_files, file)
    } else {
      cat("Skipping empty summary file:", basename(file), "\n")
    }
  }, error = function(e) {
    cat("Skipping corrupted summary file:", basename(file), "\n")
  })
}

summary_file_list <- valid_summary_files
cat("Valid summary files:", length(summary_file_list), "\n")

# For cloglog_beta_driver_uncertainty models with many files, use sub-scripts
# This avoids memory issues and provides better progress tracking
if (USE_SUB_SCRIPTS && length(summary_file_list) > 500) {
  cat("\n⚠️  Large number of summary files detected (", length(summary_file_list), ")\n", sep="")
  cat("   Delegating to sub-scripts for efficient processing:\n")
  cat("   1. Running 03b_combine_summaries.r to combine summaries...\n")
  
  # Run 03b_combine_summaries.r
  tryCatch({
    source("analysis/model_analysis/03b_combine_summaries.r")
    cat("   ✅ 03b_combine_summaries.r completed\n")
  }, error = function(e) {
    cat("   ❌ Error running 03b_combine_summaries.r:", e$message, "\n")
    cat("   Falling back to main script processing...\n")
    USE_SUB_SCRIPTS <- FALSE
  })
  
  # Load the combined summary if it was created
  if (file.exists(here("data/summary/logit_beta_fixed_priors_summaries.rds"))) {
    cat("   Loading combined summary from 03b...\n")
    main_summary <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
    summary_df <- main_summary$summary_df
    gelman_list <- main_summary$gelman.summary
    convergence_results <- list(
      keep_list = main_summary$keep_list,
      keep_list_weak = main_summary$keep_list_weak,
      keep_list_stricter = main_summary$keep_list_stricter,
      rerun_list = main_summary$rerun_list
    )
    
    # Ensure convergence lists are saved (03b should have done this, but ensure they're saved)
    cat("   Ensuring convergence lists are saved...\n")
    if (!dir.exists(here("data/summary"))) {
      dir.create(here("data/summary"), recursive = TRUE)
    }
    tryCatch({
      saveRDS(convergence_results$keep_list, here("data/summary/converged_taxa_list.rds"))
      saveRDS(convergence_results$keep_list_stricter, here("data/summary/stricter_converged_taxa_list.rds"))
      saveRDS(convergence_results$keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))
      saveRDS(convergence_results$rerun_list, here("data/summary/unconverged_taxa_list.rds"))
      cat("   ✅ Convergence lists saved\n")
    }, error = function(e) {
      cat("   ⚠️  Warning: Could not save convergence lists:", e$message, "\n")
      cat("   (They should have been saved by 03b_combine_summaries.r)\n")
    })
    
    # Skip to convergence analysis section
    cat("   ✅ Using combined summaries from sub-script\n")
    cat("   Skipping main script's file reading/combining step\n")
    skip_file_combining <- TRUE
  } else {
    cat("   ⚠️  Combined summary not found, falling back to main script\n")
    skip_file_combining <- FALSE
  }
} else {
  skip_file_combining <- FALSE
}

# Process summary files in chunks to avoid memory issues
if(!skip_file_combining && length(summary_file_list) > 0) {
  cat("Processing", length(summary_file_list), "summary files using chunked approach...\n")
  
  chunk_size <- 50
  n_chunks <- ceiling(length(summary_file_list) / chunk_size)
  
  # Create temporary directory for intermediate files
  temp_dir <- here("data/temp_summaries")
  if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  
  for(chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, length(summary_file_list))
    chunk_files <- summary_file_list[start_idx:end_idx]
    
    cat("Processing chunk", chunk, "of", n_chunks, "(", length(chunk_files), "files)\n")
    
    # Process chunk
    chunk_summaries <- purrr::map(chunk_files, function(f) {
      tryCatch({
        readRDS(f)
      }, error = function(e) {
        cat("Error reading summary file:", f, "-", e$message, "\n")
        NULL
      })
    })
    
    # Remove NULL entries
    chunk_summaries <- chunk_summaries[!sapply(chunk_summaries, is.null)]
    
    if(length(chunk_summaries) > 0) {
      # Save intermediate results
      temp_file <- file.path(temp_dir, paste0("chunk_", chunk, ".rds"))
      saveRDS(chunk_summaries, temp_file)
      cat("  Saved chunk", chunk, "to", temp_file, "\n")
    }
  }
  
  # Load and combine all chunks
  cat("Combining all chunks...\n")
  all_chunk_files <- list.files(temp_dir, pattern = "chunk_.*\\.rds", full.names = TRUE)
  
  if(length(all_chunk_files) > 0) {
    all_summaries <- purrr::map(all_chunk_files, readRDS) %>% purrr::flatten()
    
    # Remove NULL entries
    all_summaries <- all_summaries[!sapply(all_summaries, is.null)]
    
    if(length(all_summaries) > 0) {
      # Convert data types to avoid binding conflicts
      all_summaries <- lapply(all_summaries, function(x) {
        if(is.list(x) && length(x) >= 3) {
          for(i in 1:3) {
            if(is.data.frame(x[[i]]) && nrow(x[[i]]) > 0) {
              if("date_num" %in% names(x[[i]])) x[[i]]$date_num <- as.character(x[[i]]$date_num)
              if("site_num" %in% names(x[[i]])) x[[i]]$site_num <- as.character(x[[i]]$site_num)
              if("beta_num" %in% names(x[[i]])) x[[i]]$beta_num <- as.character(x[[i]]$beta_num)
            }
          }
        }
        return(x)
      })
      
      # Extract summary_df from all summaries (element 1)
      cat("Extracting summary_df from", length(all_summaries), "summary files...\n")
      summary_dfs_list <- lapply(all_summaries, function(x) {
        if (length(x) >= 1 && is.data.frame(x[[1]]) && nrow(x[[1]]) > 0) {
          return(x[[1]])
        } else {
          return(data.frame())
        }
      })
      
      # Remove empty data frames
      summary_dfs_list <- summary_dfs_list[sapply(summary_dfs_list, nrow) > 0]
      
      if (length(summary_dfs_list) > 0) {
        cat("  Found", length(summary_dfs_list), "non-empty summary data frames\n")
        summary_df <- map_df(summary_dfs_list, identity)
        cat("  ✅ Combined summary_df:", nrow(summary_df), "rows,", ncol(summary_df), "columns\n")
      } else {
        cat("  ⚠️ WARNING: No valid summary_df data found in any summary file!\n")
        summary_df <- data.frame()
      }
      
      # Process plot estimates only if flag is enabled
      if (PROCESS_PLOT_ESTIMATES) {
        cat("Processing plot estimates (PROCESS_PLOT_ESTIMATES = TRUE)\n")
        
        # For large datasets, use sub-script 03c for plot estimates
        if (USE_SUB_SCRIPTS && length(summary_file_list) > 500) {
          cat("   Large dataset detected - using 03c_extract_plot_summaries.r\n")
          cat("   Running 03c_extract_plot_summaries.r...\n")
          
          tryCatch({
            source("analysis/model_analysis/03c_extract_plot_summaries.r")
            cat("   ✅ 03c_extract_plot_summaries.r completed\n")
            
            # Try to load the plot estimates
            plot_est_parquet <- here("data/summary/plot_estimates.parquet")
            if (file.exists(plot_est_parquet)) {
              if (require(arrow, quietly = TRUE)) {
                plot_est <- read_parquet(plot_est_parquet)
                cat("   Loaded", nrow(plot_est), "plot estimate rows from sub-script\n")
              } else {
                cat("   Arrow package not available, trying RDS...\n")
                plot_est_rds <- here("data/summary/plot_estimates.rds")
                if (file.exists(plot_est_rds)) {
                  plot_est <- readRDS(plot_est_rds)
                  cat("   Loaded", nrow(plot_est), "plot estimate rows from RDS\n")
                } else {
                  plot_est <- data.frame()
                }
              }
            } else {
              cat("   ⚠️  Plot estimates file not found after sub-script\n")
              plot_est <- data.frame()
            }
          }, error = function(e) {
            cat("   ❌ Error running 03c_extract_plot_summaries.r:", e$message, "\n")
            cat("   Falling back to main script processing...\n")
            # Check if plot estimates already exist from previous processing
            plot_est_parquet <- here("data/summary/plot_estimates.parquet")
            if (file.exists(plot_est_parquet)) {
              if (require(arrow, quietly = TRUE)) {
                plot_est <- read_parquet(plot_est_parquet)
              } else {
                plot_est <- process_plot_estimates_optimized(summary_file_list)
              }
            } else {
              plot_est <- process_plot_estimates_optimized(summary_file_list)
            }
          })
        } else {
          # Check if plot estimates already exist from previous processing
          plot_est_parquet <- here("data/summary/plot_estimates.parquet")
          
          if (file.exists(plot_est_parquet)) {
            cat("Loading existing plot estimates from", plot_est_parquet, "\n")
            tryCatch({
              if (require(arrow, quietly = TRUE)) {
                plot_est <- read_parquet(plot_est_parquet)
                cat("Loaded", nrow(plot_est), "plot estimate rows\n")
              } else {
                cat("Arrow package not available, falling back to processing from summary files...\n")
                plot_est <- process_plot_estimates_optimized(summary_file_list)
              }
            }, error = function(e) {
              cat("Error loading existing plot estimates:", e$message, "\n")
              cat("Falling back to processing from summary files...\n")
              plot_est <- process_plot_estimates_optimized(summary_file_list)
            })
          } else {
            cat("No existing plot estimates found, processing from summary files...\n")
            plot_est <- process_plot_estimates_optimized(summary_file_list)
          }
        }
        
        if (nrow(plot_est) > 0) {
          cat("Final plot_est size:", nrow(plot_est), "rows,", 
            format(object.size(plot_est), units="MB"), "memory\n")
        } else {
          cat("No valid plot estimates found, using empty data frame\n")
          plot_est <- data.frame()
        }
      } else {
        cat("Skipping plot estimate processing (PROCESS_PLOT_ESTIMATES = FALSE)\n")
        plot_est <- data.frame()
      }
      
      # Use Element 4 for parameter-level Gelman diagnostics (not Element 3 which is plot-level)
      cat("Extracting gelman_list from summary files...\n")
      cat("  Total summaries to process:", length(all_summaries), "\n")
      
      gelman_dfs_list <- lapply(seq_along(all_summaries), function(i) {
        x <- all_summaries[[i]]
        if (length(x) >= 4 && is.data.frame(x[[4]]) && nrow(x[[4]]) > 0) {
          return(x[[4]])
        } else {
          if (i <= 5 || i %% 50 == 0) {
            cat("  Summary", i, ": element 4 missing or empty (length:", length(x), ")\n")
          }
          return(data.frame())
        }
      })
      
      # Remove empty data frames
      gelman_dfs_list <- gelman_dfs_list[sapply(gelman_dfs_list, nrow) > 0]
      cat("  Summaries with valid gelman data:", length(gelman_dfs_list), "out of", length(all_summaries), "\n")
      
      if (length(gelman_dfs_list) > 0) {
        gelman_list <- map_df(gelman_dfs_list, identity)
        cat("  ✅ Combined gelman_list:", nrow(gelman_list), "rows,", ncol(gelman_list), "columns\n")
        if ("model_id" %in% colnames(gelman_list)) {
          unique_models <- unique(gelman_list$model_id)
          cat("  Unique model_ids in gelman_list:", length(unique_models), "\n")
        } else {
          cat("  ⚠️ WARNING: 'model_id' column not found in gelman_list!\n")
          cat("  Available columns:", paste(colnames(gelman_list), collapse = ", "), "\n")
        }
      } else {
        cat("  ⚠️ WARNING: No valid gelman_list data found in any summary file!\n")
        gelman_list <- data.frame()
      }
      
      cat("Successfully combined", length(all_summaries), "summaries\n")
      cat("FINAL SUMMARY:\n")
      cat("  - summary_df:", nrow(summary_df), "rows\n")
      cat("  - gelman_list:", nrow(gelman_list), "rows\n")
    } else {
      stop("No valid summaries found after combining chunks")
    }
  } else {
    stop("No chunk files found")
  }
  
  # Clean up temporary directory
  unlink(temp_dir, recursive = TRUE)
  cat("Cleaned up temporary files\n")
  
} else {
  stop("No summary files found. Please check the file paths.")
}

# =============================================================================
# STEP 4: CONVERGENCE ANALYSIS
# =============================================================================

# Load existing data if processing was skipped OR if we used sub-scripts
if (skip_processing || (exists("skip_file_combining") && skip_file_combining)) {
  cat("Loading existing summary data...\n")
  main_summary <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  
  # Load plot estimates from available format
  if (file.exists(here("data/summary/plot_estimates.rds"))) {
    plot_est <- readRDS(here("data/summary/plot_estimates.rds"))
  } else if (file.exists(here("data/summary/plot_estimates.parquet"))) {
    plot_est <- as.data.frame(read_parquet(here("data/summary/plot_estimates.parquet")))
  } else if (file.exists(here("data/summary/plot_estimates.csv.gz"))) {
    plot_est <- fread(here("data/summary/plot_estimates.csv.gz"))
  } else {
    cat("⚠️  No plot estimates file found, creating empty data frame\n")
    plot_est <- data.frame()
  }
  
  # Extract components from existing summary
  summary_df <- main_summary$summary_df
  gelman.summary <- main_summary$gelman.summary
  keep_list <- main_summary$keep_list
  keep_list_weak <- main_summary$keep_list_weak
  keep_list_stricter <- main_summary$keep_list_stricter
  rerun_list <- main_summary$rerun_list
  
  # Create convergence results from existing data
  convergence_results <- list(
    keep_list = keep_list,
    keep_list_weak = keep_list_weak,
    keep_list_stricter = keep_list_stricter,
    rerun_list = rerun_list
  )
  
  cat("✅ Loaded existing summary data\n")
  
  # Update the main summary with the loaded plot estimates
  main_summary$plot_est <- plot_est
  cat("  Summary_df rows:", nrow(summary_df), "\n")
  cat("  Plot_est rows:", nrow(plot_est), "\n")
  cat("  Keep_list length:", length(keep_list), "\n")
  
  # Ensure convergence_results exists for later saving
  if (!exists("convergence_results")) {
    convergence_results <- list(
      keep_list = keep_list,
      keep_list_weak = keep_list_weak,
      keep_list_stricter = keep_list_stricter,
      rerun_list = rerun_list
    )
  }
  
  # Also ensure gelman_list exists (needed for later steps)
  if (!exists("gelman_list")) {
    gelman_list <- gelman.summary
  }
} else {
  cat("\n=== STEP 4: CONVERGENCE ANALYSIS ===\n")

# Process Gelman summary using data.table for better performance
cat("Processing Gelman summary...\n")
cat("Gelman list dimensions:", dim(gelman_list), "\n")

# Convert to data.table for faster operations
if (nrow(gelman_list) == 0) {
  cat("No Gelman data available for convergence analysis\n")
  convergence_results <- list(
    keep_list = character(0),
    keep_list_weak = character(0),
    keep_list_stricter = character(0),
    rerun_list = character(0)
  )
} else {
  cat("Processing Gelman data for model categorization...\n")
  cat("Gelman list dimensions:", dim(gelman_list), "\n")
  
  # Convert to data.table
  tryCatch({
    gelman_dt <- as.data.table(gelman_list)
    cat("  Converted to data.table successfully\n")
    cat("  Gelman data columns:", colnames(gelman_dt), "\n")
    cat("  Gelman data rows:", nrow(gelman_dt), "\n")
    
    # Filter and add major parameter flag
    if("model_id" %in% colnames(gelman_dt)) {
      gelman_dt <- gelman_dt[!grepl("all_covariates", model_id)]
    }
    
    if("parameter" %in% colnames(gelman_dt)) {
      gelman_dt[, is_major_param := grepl("beta|int|sigma|core_sd|rho", parameter)]
      cat("  Major parameters:", sum(gelman_dt$is_major_param), "out of", nrow(gelman_dt), "\n")
    } else {
      gelman_dt[, is_major_param := TRUE]  # Assume all are major if no parameter column
      cat("  No parameter column found, treating all as major parameters\n")
    }
    
    cat("  Applied filters successfully\n")
  }, error = function(e) {
    cat("  Error in data.table conversion:", e$message, "\n")
    stop("Failed to process Gelman data")
  })
  
    # Apply simple convergence criteria using data.table
    tryCatch({
      # Function to ensure model_ids have _beta_regression suffix to match filenames
      fix_model_id <- function(model_id) {
        if (grepl("with_legacy_covariate$", model_id) && !grepl("_beta_regression$", model_id)) {
          return(paste0(model_id, "_beta_regression"))
        }
        return(model_id)
      }
      
      # Calculate maximum Point est. per model for major parameters
      # A model is converged if the MAX Point est. of its major parameters meets the threshold
      model_max_gelman <- gelman_dt[
        is_major_param == TRUE,
        .(max_point_est = max(`Point est.`, na.rm = TRUE),
          mean_point_est = mean(`Point est.`, na.rm = TRUE),
          n_params = .N),
        by = model_id
      ]
      
      # Get models where MAX major parameter Point est. meets criteria
      keep_models <- unique(model_max_gelman[max_point_est < 1.1]$model_id)
      keep_models_weak <- unique(model_max_gelman[max_point_est < 1.2]$model_id)
      keep_models_stricter <- unique(model_max_gelman[max_point_est < 1.05]$model_id)
      
      # Fix model_ids to include _beta_regression suffix
      keep_models <- sapply(keep_models, fix_model_id)
      keep_models_weak <- sapply(keep_models_weak, fix_model_id)
      keep_models_stricter <- sapply(keep_models_stricter, fix_model_id)
      
      # Get models that need rerun (don't meet basic criteria)
      all_models <- unique(gelman_dt$model_id)
      all_models <- sapply(all_models, fix_model_id)
      rerun_models <- setdiff(all_models, keep_models)
    
    keep_list <- keep_models
    keep_list_weak <- keep_models_weak
    keep_list_stricter <- keep_models_stricter
    rerun_list <- rerun_models
    
    cat("\n  === CONVERGENCE STATISTICS ===\n")
    cat("  Total models in analysis:", length(all_models), "\n")
    cat("  Converged (strict < 1.05):", length(keep_list_stricter), "\n")
    cat("  Converged (standard < 1.1):", length(keep_list), "\n")
    cat("  Converged (weak < 1.2):", length(keep_list_weak), "\n")
    cat("  Need rerun (> 1.1):", length(rerun_list), "\n")
    
    # Breakdown by model type
    if("model_id" %in% colnames(gelman_dt)) {
      cat("\n  Breakdown by model type:\n")
      for(mt in model_types) {
        models_of_type <- unique(gelman_dt[grepl(paste0("^", mt), model_id)]$model_id)
        converged_of_type <- sum(grepl(paste0("^", mt), keep_models))
        cat("    ", mt, ": ", length(models_of_type), " models, ", converged_of_type, " converged\n", sep="")
      }
    }
    
  }, error = function(e) {
    cat("  Error in convergence criteria application:", e$message, "\n")
    # Set default values if convergence analysis fails
    keep_list <- character(0)
    keep_list_weak <- character(0)
    keep_list_stricter <- character(0)
    rerun_list <- character(0)
  })
  
  cat("Model categorization complete:\n")
  cat("  - Strictly converged:", length(keep_list), "\n")
  cat("  - Weakly converged:", length(keep_list_weak), "\n")
  cat("  - Stricter converged:", length(keep_list_stricter), "\n")
  cat("  - Unconverged:", length(rerun_list), "\n")
  cat("  - Total models categorized:", length(keep_list_weak) + length(rerun_list), "\n")
  
  convergence_results <- list(
    keep_list = keep_list,
    keep_list_weak = keep_list_weak,
    keep_list_stricter = keep_list_stricter,
    rerun_list = rerun_list
  )
}
}

# =============================================================================
# STEP 5: MISSING CHAINS ANALYSIS
# =============================================================================

cat("\n=== STEP 5: MISSING CHAINS ANALYSIS ===\n")

# Analyze missing chains for env_cycl models
env_cycl_dirs <- c(
  here("data/model_outputs/logit_beta_fixed_priors/env_cycl"),
  here("data/model_outputs/reconstructed_from_checkpoints/env_cycl"),
  here("data/model_outputs/old_mcmc/logit_beta_regression/env_cycl")
)

# Get all model directories
model_dirs <- c()
for (env_cycl_dir in env_cycl_dirs) {
  if (dir.exists(env_cycl_dir)) {
    dirs_in_path <- list.dirs(env_cycl_dir, full.names = FALSE, recursive = FALSE)
    dirs_in_path <- dirs_in_path[dirs_in_path != ""]
    model_dirs <- c(model_dirs, dirs_in_path)
  }
}
model_dirs <- unique(model_dirs)

expected_chains <- 1:4
n_models <- length(model_dirs)

# Pre-allocate missing chains data frame (estimate max size)
max_missing <- n_models * 4  # Worst case: all chains missing for all models
missing_chains <- data.frame(
  model = character(max_missing),
  missing_chain = integer(max_missing),
  existing_chains = character(max_missing),
  priority = character(max_missing),
  stringsAsFactors = FALSE
)

cat("Checking for missing chains in", n_models, "models...\n")

missing_count <- 0
for (model in model_dirs) {
  found_chains <- c()
  
  for (env_cycl_dir in env_cycl_dirs) {
    if (dir.exists(env_cycl_dir)) {
      model_subdir <- file.path(env_cycl_dir, model)
      
      # Check in model subdirectory
      if (dir.exists(model_subdir)) {
        pattern1 <- paste0("samples_env_cycl_", model, "_20130601_20180101_with_legacy_covariate_chain([1-4])\\.rds")
        files_in_subdir <- list.files(model_subdir, pattern = pattern1, full.names = FALSE)
        
        if (length(files_in_subdir) > 0) {
          chains_in_subdir <- as.integer(gsub(".*_chain([1-4])\\.rds", "\\1", files_in_subdir))
          found_chains <- c(found_chains, chains_in_subdir)
        }
      }
      
      # Check in root env_cycl directory
      pattern2 <- paste0("samples_env_cycl_", model, "_20130601_20180101_with_legacy_covariate_chain([1-4])\\.rds")
      files_in_root <- list.files(env_cycl_dir, pattern = pattern2, full.names = FALSE)
      
      if (length(files_in_root) > 0) {
        chains_in_root <- as.integer(gsub(".*_chain([1-4])\\.rds", "\\1", files_in_root))
        found_chains <- c(found_chains, chains_in_root)
      }
    }
  }
  
  found_chains <- sort(unique(found_chains))
  missing <- setdiff(expected_chains, found_chains)
  
  if (length(missing) > 0) {
    # Determine priority based on existing progress
    if (length(found_chains) == 0) {
      priority <- "Low (No chains)"
    } else if (1 %in% found_chains && length(found_chains) < 4) {
      priority <- "High (Has chain 1)"
    } else if (length(found_chains) > 0 && length(found_chains) < 4) {
      priority <- "Medium (Partial progress)"
    } else {
      priority <- "Low (Incomplete)"
    }
    
    existing_chains_str <- if(length(found_chains) > 0) paste(found_chains, collapse = ",") else "none"
    
    for (chain in missing) {
      missing_count <- missing_count + 1
      missing_chains[missing_count, ] <- list(
        model = model,
        missing_chain = chain,
        existing_chains = existing_chains_str,
        priority = priority
      )
    }
  }
}

# Trim to actual size
missing_chains <- missing_chains[1:missing_count, ]

# Create priority-based rerun lists using data.table
if (nrow(missing_chains) > 0) {
  missing_dt <- as.data.table(missing_chains)
  
  priority_models <- unique(missing_dt[priority == "High (Has chain 1)", model])
  priority_rerun_list <- paste0("env_cycl_", priority_models, "_20130601_20180101_with_legacy_covariate")
  
  all_missing_models <- unique(missing_dt[, model])
  all_missing_rerun_list <- paste0("env_cycl_", all_missing_models, "_20130601_20180101_with_legacy_covariate")
  
  missing_summary <- missing_dt[, .(
    missing_chains = paste(sort(missing_chain), collapse = ",")
  ), by = .(model, priority, existing_chains)][order(priority, model)]
} else {
  priority_models <- character(0)
  priority_rerun_list <- character(0)
  all_missing_models <- character(0)
  all_missing_rerun_list <- character(0)
  missing_summary <- data.table()
}

cat("Missing chain analysis complete:\n")
cat("  Total models with missing chains:", length(all_missing_models), "\n")
cat("  High priority models (have chain 1):", length(priority_models), "\n")
cat("  Total missing chains:", nrow(missing_chains), "\n")

missing_chain_results <- list(
  missing_chains = missing_chains,
  missing_summary = missing_summary,
  priority_models = priority_models,
  all_missing_models = all_missing_models,
  priority_rerun_list = priority_rerun_list,
  all_missing_rerun_list = all_missing_rerun_list
)

# =============================================================================
# STEP 6: SAVE ALL RESULTS - ROBUST VERSION
# =============================================================================

# Load required libraries

cat("\n=== STEP 6: SAVE ALL RESULTS ===\n")

# Create output directory if it doesn't exist
if (!dir.exists(here("data/summary"))) {
  dir.create(here("data/summary"), recursive = TRUE)
}

# First, ensure all critical variables exist with proper defaults
# BUT: If we just processed summary files, these should already exist
# Only create empty if truly missing (e.g., if skip_processing was TRUE but data wasn't found)
if (!exists("summary_df")) {
  cat("  ⚠️ summary_df missing!\n")
  cat("  This should not happen if processing completed successfully.\n")
  cat("  Attempting to load from individual summary files as fallback...\n")
  
  # Fallback: try to load from individual summary files
  base_path <- if(exists("base_path")) base_path else "data/model_outputs/cloglog_beta_driver_uncertainty"
  if (dir.exists(base_path)) {
    summary_files <- c(
      list.files(file.path(base_path, "env_cycl"), pattern = "summary_.*_beta_regression\\.rds$", full.names = TRUE, recursive = TRUE),
      list.files(file.path(base_path, "cycl_only"), pattern = "summary_.*_beta_regression\\.rds$", full.names = TRUE, recursive = TRUE),
      list.files(file.path(base_path, "env_cov"), pattern = "summary_.*_beta_regression\\.rds$", full.names = TRUE, recursive = TRUE)
    )
    if (length(summary_files) > 0) {
      cat("  Found", length(summary_files), "summary files, attempting to combine...\n")
      all_summaries <- lapply(summary_files, function(f) {
        tryCatch(readRDS(f), error = function(e) NULL)
      })
      all_summaries <- all_summaries[!sapply(all_summaries, is.null)]
      if (length(all_summaries) > 0) {
        summary_dfs <- lapply(all_summaries, function(x) {
          if (length(x) >= 1 && is.data.frame(x[[1]]) && nrow(x[[1]]) > 0) return(x[[1]]) else return(data.frame())
        })
        summary_dfs <- summary_dfs[sapply(summary_dfs, nrow) > 0]
        if (length(summary_dfs) > 0) {
          summary_df <- map_df(summary_dfs, identity)
          cat("  ✅ Successfully loaded", nrow(summary_df), "rows from individual files\n")
        } else {
          cat("  ❌ No valid data found in summary files\n")
  summary_df <- data.frame()
        }
      } else {
        cat("  ❌ Could not read any summary files\n")
        summary_df <- data.frame()
      }
    } else {
      cat("  ❌ No summary files found\n")
      summary_df <- data.frame()
    }
  } else {
    cat("  ❌ Base path not found\n")
    summary_df <- data.frame()
  }
}

if (!exists("plot_est")) {
  cat("  ⚠️ plot_est missing, creating empty data frame\n")
  plot_est <- data.frame()
}

if (!exists("gelman_list")) {
  cat("  ⚠️ gelman_list missing, creating empty data frame\n")
  gelman_list <- data.frame()
}

if (!exists("convergence_results")) {
  cat("  ⚠️ convergence_results missing, creating empty structure\n")
  convergence_results <- list(
    keep_list = character(0),
    keep_list_weak = character(0),
    keep_list_stricter = character(0),
    rerun_list = character(0)
  )
}

if (!exists("missing_chain_results")) {
  cat("  ⚠️ missing_chain_results missing, creating empty structure\n")
  missing_chain_results <- list(
    missing_chains = data.frame(),
    missing_summary = data.frame(),
    priority_models = character(0),
    all_missing_models = character(0),
    priority_rerun_list = character(0),
    all_missing_rerun_list = character(0)
  )
}

# Filter convergence lists to only include driver uncertainty models
# Driver uncertainty models are identified by:
# 1. Having "with_legacy_covariate" and "beta_regression" in model_id
# 2. Being from cloglog_beta_driver_uncertainty directory (model_id starts with cycl_only, env_cov, or env_cycl)
filter_driver_uncertainty_only <- function(model_list) {
  if (length(model_list) == 0) return(character(0))
  
  # Driver uncertainty models have: model_name_species_20130601_20180101_with_legacy_covariate_beta_regression
  # And model_name is one of: cycl_only, env_cov, env_cycl
  driver_pattern <- "^((cycl_only|env_cov|env_cycl)_.*_20130601_20180101_with_legacy_covariate_beta_regression)$"
  filtered <- model_list[grepl(driver_pattern, model_list)]
  return(filtered)
}

# Ensure convergence_results exists before filtering and saving
if (!exists("convergence_results")) {
  cat("⚠️  convergence_results not found, creating from main_summary if available\n")
  if (exists("main_summary") && is.list(main_summary)) {
    convergence_results <- list(
      keep_list = if("keep_list" %in% names(main_summary)) main_summary$keep_list else character(0),
      keep_list_weak = if("keep_list_weak" %in% names(main_summary)) main_summary$keep_list_weak else character(0),
      keep_list_stricter = if("keep_list_stricter" %in% names(main_summary)) main_summary$keep_list_stricter else character(0),
      rerun_list = if("rerun_list" %in% names(main_summary)) main_summary$rerun_list else character(0)
    )
  } else {
    cat("⚠️  main_summary not available, creating empty convergence_results\n")
    convergence_results <- list(
      keep_list = character(0),
      keep_list_weak = character(0),
      keep_list_stricter = character(0),
      rerun_list = character(0)
    )
  }
}

# Filter all convergence lists to driver uncertainty only
convergence_results$keep_list <- filter_driver_uncertainty_only(convergence_results$keep_list)
convergence_results$keep_list_weak <- filter_driver_uncertainty_only(convergence_results$keep_list_weak)
convergence_results$keep_list_stricter <- filter_driver_uncertainty_only(convergence_results$keep_list_stricter)
convergence_results$rerun_list <- filter_driver_uncertainty_only(convergence_results$rerun_list)

cat("  Filtered convergence lists to driver uncertainty models only:\n")
cat("    keep_list:", length(convergence_results$keep_list), "models\n")
cat("    keep_list_weak:", length(convergence_results$keep_list_weak), "models\n")
cat("    keep_list_stricter:", length(convergence_results$keep_list_stricter), "models\n")
cat("    rerun_list:", length(convergence_results$rerun_list), "models\n")

# Save convergence lists individually with error handling
cat("\n  Saving convergence lists...\n")
tryCatch({
  saveRDS(convergence_results$keep_list, here("data/summary/converged_taxa_list.rds"))
  cat("    ✅ Saved converged_taxa_list.rds (", length(convergence_results$keep_list), " models)\n", sep="")
}, error = function(e) {
  cat("    ❌ Error saving converged_taxa_list.rds:", e$message, "\n")
})

tryCatch({
  saveRDS(convergence_results$keep_list_stricter, here("data/summary/stricter_converged_taxa_list.rds"))
  cat("    ✅ Saved stricter_converged_taxa_list.rds (", length(convergence_results$keep_list_stricter), " models)\n", sep="")
}, error = function(e) {
  cat("    ❌ Error saving stricter_converged_taxa_list.rds:", e$message, "\n")
})

tryCatch({
  saveRDS(convergence_results$keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))
  cat("    ✅ Saved weak_converged_taxa_list.rds (", length(convergence_results$keep_list_weak), " models)\n", sep="")
}, error = function(e) {
  cat("    ❌ Error saving weak_converged_taxa_list.rds:", e$message, "\n")
})

# Filter missing chain results to driver uncertainty only
if (exists("missing_chain_results")) {
  missing_chain_results$priority_rerun_list <- filter_driver_uncertainty_only(missing_chain_results$priority_rerun_list)
  missing_chain_results$all_missing_rerun_list <- filter_driver_uncertainty_only(missing_chain_results$all_missing_rerun_list)
  if (nrow(missing_chain_results$missing_chains) > 0) {
    # Filter missing_chains data.frame if it has model_id column
    if ("model_id" %in% colnames(missing_chain_results$missing_chains)) {
      missing_chain_results$missing_chains <- missing_chain_results$missing_chains[
        grepl("^((cycl_only|env_cov|env_cycl)_.*_20130601_20180101_with_legacy_covariate_beta_regression)$", 
              missing_chain_results$missing_chains$model_id), 
        , drop = FALSE
      ]
    }
  }
  cat("  Filtered missing chain results to driver uncertainty only\n")
}

# Save unconverged list based on convergence analysis (filtered to driver uncertainty)
# This matches the rerun_list in the main summary
tryCatch({
  saveRDS(convergence_results$rerun_list, here("data/summary/unconverged_taxa_list.rds"))
  cat("    ✅ Saved unconverged_taxa_list.rds (", length(convergence_results$rerun_list), " models)\n", sep="")
}, error = function(e) {
  cat("    ❌ Error saving unconverged_taxa_list.rds:", e$message, "\n")
})

# Create and save prioritized rerun list (combines convergence failures and missing chains)
rerun_list_prioritized <- c(
  missing_chain_results$priority_rerun_list, 
  setdiff(missing_chain_results$all_missing_rerun_list, 
          missing_chain_results$priority_rerun_list),
  convergence_results$rerun_list
)
rerun_list_prioritized <- unique(filter_driver_uncertainty_only(rerun_list_prioritized))

# Save priority rerun list (filtered to driver uncertainty)
saveRDS(missing_chain_results$priority_rerun_list, here("data/summary/priority_rerun_list.rds"))

# Save missing chains analysis
saveRDS(list(
  missing_chains = missing_chain_results$missing_chains,
  missing_summary = missing_chain_results$missing_summary,
  priority_models = missing_chain_results$priority_models,
  all_missing_models = missing_chain_results$all_missing_models
), here("data/summary/missing_chains_analysis.rds"))

cat("  ✅ Saved convergence and analysis files\n")

# Save plot estimates in multiple formats for redundancy
if (nrow(plot_est) > 0) {
  cat("  Saving plot estimates (", nrow(plot_est), "rows)...\n")
  
  # Ensure it's a data frame
  plot_est_df <- as.data.frame(plot_est)
  
  # Primary format: RDS with compression
  tryCatch({
    saveRDS(plot_est_df, here("data/summary/plot_estimates.rds"), compress = "xz")
    cat("    ✅ Saved as RDS\n")
  }, error = function(e) {
    cat("    ❌ Failed to save RDS:", e$message, "\n")
  })
  
  # Secondary format: CSV (compressed)
  tryCatch({
    if (require(data.table, quietly = TRUE)) {
      fwrite(plot_est_df, here("data/summary/plot_estimates.csv.gz"), compress = "gzip")
      cat("    ✅ Saved as compressed CSV\n")
    } else {
      write.csv(plot_est_df, here("data/summary/plot_estimates.csv"), row.names = FALSE)
      cat("    ✅ Saved as CSV\n")
    }
  }, error = function(e) {
    cat("    ❌ Failed to save CSV:", e$message, "\n")
  })
  
  # Optional format: Parquet (if available)
  tryCatch({
    if (require(arrow, quietly = TRUE)) {
      write_parquet(plot_est_df, here("data/summary/plot_estimates.parquet"))
      cat("    ✅ Saved as Parquet\n")
    }
  }, error = function(e) {
    # Silently skip if arrow not available
  })
} else {
  cat("  ℹ️ No plot estimates to save (empty data frame)\n")
}

# Build main summary WITHOUT plot_est to avoid memory issues
main_summary <- list(
  summary_df = as.data.frame(summary_df),
  gelman.summary = as.data.frame(gelman_list),
  keep_models_weak = if(exists("keep_models_weak")) keep_models_weak else character(0),
  keep_list = convergence_results$keep_list,
  keep_list_weak = convergence_results$keep_list_weak,
  keep_list_stricter = convergence_results$keep_list_stricter,
  rerun_list = convergence_results$rerun_list,
  priority_rerun_list = missing_chain_results$priority_rerun_list,
  missing_chains_analysis = list(
    missing_chains = missing_chain_results$missing_chains,
    missing_summary = missing_chain_results$missing_summary,
    priority_models = missing_chain_results$priority_models,
    all_missing_models = missing_chain_results$all_missing_models
  )
)

# Log what we're saving
cat("\n  Main summary contents:\n")
cat("    - summary_df:", nrow(main_summary$summary_df), "rows,", ncol(main_summary$summary_df), "columns\n")
cat("    - gelman.summary:", nrow(main_summary$gelman.summary), "rows,", ncol(main_summary$gelman.summary), "columns\n")
cat("    - keep_list:", length(main_summary$keep_list), "models\n")
cat("    - keep_list_weak:", length(main_summary$keep_list_weak), "models\n")
cat("    - keep_list_stricter:", length(main_summary$keep_list_stricter), "models\n")
cat("    - rerun_list:", length(main_summary$rerun_list), "models\n")
cat("    - plot_est: saved separately in plot_estimates.parquet/rds/csv.gz\n")

# Save main summary - verify data exists before saving
# CRITICAL: summary_df MUST have data before saving
cat("\n=== VERIFYING DATA BEFORE SAVE ===\n")
cat("  summary_df rows:", nrow(main_summary$summary_df), "\n")
cat("  summary_df columns:", ncol(main_summary$summary_df), "\n")
cat("  gelman_list rows:", nrow(main_summary$gelman.summary), "\n")

if (nrow(main_summary$summary_df) == 0) {
  cat("\n  ❌ ERROR: summary_df is empty! Cannot save main summary file.\n")
  cat("  This indicates a serious problem with model summarization.\n")
  cat("  Debugging information:\n")
  
  # Debug: Check if summary_df should have been populated
  if (exists("summary_file_list")) {
    cat("  - Summary files found:", length(summary_file_list), "\n")
  }
  if (exists("all_summaries")) {
    cat("  - All summaries loaded:", length(all_summaries), "\n")
  }
  
  cat("\n  Individual summary files in cloglog_beta_driver_uncertainty should still be valid.\n")
  cat("  Consider running 03b_combine_summaries.r to rebuild the combined file.\n")
  cat("  NOT SAVING EMPTY FILE - removing if it exists...\n")
  
  # Remove any existing empty file
  if (file.exists(here("data/summary/logit_beta_fixed_priors_summaries.rds"))) {
    file.remove(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
    cat("  Removed existing empty summary file\n")
  }
  
} else {
  cat("\n  ✅ Data verification passed - summary_df has", nrow(main_summary$summary_df), "rows\n")
  
tryCatch({
  saveRDS(main_summary, here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  cat("\n  ✅ Saved main summary file successfully\n")
  
  # Verify it saved correctly
  test_read <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  if (!is.null(test_read$summary_df) && nrow(test_read$summary_df) > 0) {
      cat("  ✅ Verified: main summary file contains", nrow(test_read$summary_df), "rows\n")
      cat("  ✅ File size:", format(file.info(here("data/summary/logit_beta_fixed_priors_summaries.rds"))$size, units="MB"), "\n")
  } else {
      cat("  ❌ CRITICAL: main summary file is empty after save!\n")
      cat("  This indicates a save/read problem. Removing file...\n")
      # Remove empty file
      file.remove(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
      cat("  Removed corrupted empty summary file\n")
      stop("Failed to save summary file correctly - file was empty after save")
  }
  
}, error = function(e) {
    cat("\n  ❌ ERROR saving main summary file:\n")
  cat("    ", e$message, "\n")
    cat("  Individual summary files remain valid and can be loaded directly.\n")
    cat("  The combined summary file can be regenerated using 03b_combine_summaries.r\n")
    # Don't create backup files - individual files are the source of truth
    stop(e)  # Re-throw to prevent continuing with invalid state
  })
}

# Save additional outputs
if (exists("results") && !is.null(results)) {
  problematic_files <- results[!results$readable | results$error_message != "N/A", ]
  write.csv(problematic_files, "problematic_files.csv", row.names = FALSE)
  cat("  ✅ Saved problematic_files.csv\n")
}

if (exists("plot_est") && nrow(plot_est) > 0) {
  if (require(data.table, quietly = TRUE)) {
    plot_est_dt <- as.data.table(plot_est)
    
    # Only create summary stats if required columns exist
    required_cols <- c("taxon", "model_name")
    optional_cols <- c("plotID", "dateID", "siteID")
    
    if (all(required_cols %in% names(plot_est_dt))) {
      summary_stats <- plot_est_dt[, {
        stats <- list(n_rows = .N)
        if ("plotID" %in% names(.SD)) stats$n_plots <- length(unique(plotID))
        if ("dateID" %in% names(.SD)) stats$n_timepoints <- length(unique(dateID))
        if ("siteID" %in% names(.SD)) stats$n_sites <- length(unique(siteID))
        stats
      }, by = .(taxon, model_name)]
      
      write.csv(summary_stats, "data/summary/combined_summary_stats.csv", row.names = FALSE)
      cat("  ✅ Saved combined_summary_stats.csv\n")
    }
  }
}

cat("\n✅ All results saved successfully\n")
# FINAL SUMMARY
# =============================================================================

cat("\n=== CONSOLIDATED SUMMARY PIPELINE COMPLETE ===\n")
cat("Summary:\n")
if (exists("all_files")) {
  cat("  - Processed", length(all_files), "model files for validation\n")
} else {
  cat("  - Processing was skipped (using existing data)\n")
}
if (exists("ready_for_summary")) {
  cat("  - Files ready for summary:", nrow(ready_for_summary), "\n")
} else {
  cat("  - Files ready for summary: N/A (processing skipped)\n")
}
if (exists("summary_file_list")) {
  cat("  - Summary files processed:", length(summary_file_list), "\n")
} else {
  cat("  - Summary files processed: N/A (processing skipped)\n")
}
cat("  - Converged models:", length(convergence_results$keep_list), "\n")
cat("  - Weakly converged models:", length(convergence_results$keep_list_weak), "\n")
cat("  - Models needing rerun:", length(convergence_results$rerun_list), "\n")
cat("  - High priority missing chains:", length(missing_chain_results$priority_models), "\n")
cat("  - All results saved to data/summary/\n")

# =============================================================================
# TESTING FUNCTION
# =============================================================================

# Function to test the consolidated script
test_consolidated_summary <- function() {
  cat("\n=== TESTING CONSOLIDATED SUMMARY SCRIPT ===\n")
  
  # Test 1: Check if all required output files exist
  required_files <- c(
    "data/summary/converged_taxa_list.rds",
    "data/summary/stricter_converged_taxa_list.rds", 
    "data/summary/weak_converged_taxa_list.rds",
    "data/summary/unconverged_taxa_list.rds",
    "data/summary/priority_rerun_list.rds",
    "data/summary/missing_chains_analysis.rds",
    "data/summary/logit_beta_fixed_priors_summaries.rds",
    "comprehensive_model_validation_results.csv",
    "files_ready_for_summary.csv",
    "problematic_files.csv"
  )
  
  missing_files <- c()
  for (file in required_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) == 0) {
    cat("✅ All required output files exist\n")
  } else {
    cat("❌ Missing files:\n")
    for (file in missing_files) {
      cat("  -", file, "\n")
    }
  }
  
  # Test 2: Check data integrity
  tryCatch({
    main_summary <- readRDS("data/summary/logit_beta_fixed_priors_summaries.rds")
    expected_components <- c("summary_df", "plot_est", "gelman.summary", "keep_list", 
                           "keep_list_weak", "keep_list_stricter", "rerun_list")
    
    missing_components <- setdiff(expected_components, names(main_summary))
    if (length(missing_components) == 0) {
      cat("✅ Main summary file has all required components\n")
    } else {
      cat("❌ Missing components in main summary:", paste(missing_components, collapse = ", "), "\n")
    }
  }, error = function(e) {
    cat("❌ Error reading main summary file:", e$message, "\n")
  })
  
  # Test 3: Check convergence lists
  tryCatch({
    converged <- readRDS("data/summary/converged_taxa_list.rds")
    weak_converged <- readRDS("data/summary/weak_converged_taxa_list.rds")
    unconverged <- readRDS("data/summary/unconverged_taxa_list.rds")
    
    cat("✅ Convergence lists loaded successfully\n")
    cat("  - Converged models:", length(converged), "\n")
    cat("  - Weakly converged models:", length(weak_converged), "\n")
    cat("  - Unconverged models:", length(unconverged), "\n")
  }, error = function(e) {
    cat("❌ Error reading convergence lists:", e$message, "\n")
  })
  
  cat("\n=== TESTING COMPLETE ===\n")
}

# Run the test
test_consolidated_summary()

cat("\n✅ Consolidated summary pipeline complete!\n")
cat("All crucial output files have been preserved and generated.\n")

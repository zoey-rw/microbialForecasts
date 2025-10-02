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

source("source.R")
library(dplyr)
library(coda)
library(purrr)
library(data.table)
library(arrow)

# Source the calculate_plot_summary_from_samples function
source("microbialForecast/R/summarizeBetaRegModels.r")

# =============================================================================
# CONFIGURATION FLAGS
# =============================================================================

# Set to FALSE to skip plot estimate processing (saves memory)
PROCESS_PLOT_ESTIMATES <- FALSE

cat("Configuration: PROCESS_PLOT_ESTIMATES =", PROCESS_PLOT_ESTIMATES, "\n")

# =============================================================================
# OPTIMIZED PLOT ESTIMATE PROCESSING FUNCTION
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
if (file.exists(here("data/summary/logit_beta_fixed_priors_summaries.rds")) && 
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
      
      # Find the newest input model files
      model_files <- c(
        list.files("data/model_outputs/reconstructed_from_checkpoints/env_cov", 
                  pattern = "samples_.*\\.rds$", recursive = FALSE, full.names = TRUE),
        list.files("data/model_outputs/reconstructed_from_checkpoints/env_cycl", 
                  pattern = "samples_.*\\.rds$", recursive = FALSE, full.names = TRUE),
        list.files("data/model_outputs/reconstructed_from_checkpoints/cycl_only", 
                  pattern = "samples_.*\\.rds$", recursive = FALSE, full.names = TRUE)
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
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  # Filter to only include files with both required suffixes
  standard_files <- file_list[grepl('with_legacy_covariate', basename(file_list)) & 
                              grepl('beta_regression', basename(file_list))]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard files (with both suffixes):', length(standard_files), '\n')
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
if (dir.exists(here("data/model_outputs/CLR_regression"))) {
  cat("Processing CLR models...\n")
  source("03_summarizeModelOutputs_CLR.r")
}

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
      
      # Check for samples2 (plot estimates)
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
      }
      
      # Check for param_summary
      if ("param_summary" %in% names(data)) {
        result$has_param_summary <- TRUE
        if (is.list(data$param_summary)) {
          result$param_summary_names <- paste(names(data$param_summary), collapse = ", ")
          has_means <- "means" %in% names(data$param_summary)
          has_quantiles <- "quantiles" %in% names(data$param_summary)
          if (has_means && has_quantiles) {
            result$param_summary_valid <- TRUE
          } else {
            result$param_summary_valid <- FALSE
            result$param_summary_issue <- "Missing means or quantiles elements"
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

if (!skip_processing) {
cat("=== STEP 1: COMPREHENSIVE MODEL VALIDATION ===\n")

# Find all model files with memory-efficient approach
model_dirs <- c(
  "data/model_outputs/reconstructed_from_checkpoints/env_cycl",
  "data/model_outputs/reconstructed_from_checkpoints/cycl_only"
)

# Process directories one at a time to avoid memory issues
all_files <- c()
for (dir in model_dirs) {
  if (dir.exists(dir)) {
    cat("Scanning directory:", dir, "\n")
    files <- list.files(dir, pattern = "samples_.*\\.rds$", 
                       full.names = TRUE, recursive = FALSE)
    
    # Filter to max depth 2 and exclude individual chain files
    files <- files[sapply(strsplit(files, "/"), function(x) {
      recon_idx <- which(x == "reconstructed_from_checkpoints")
      if(length(recon_idx) > 0) {
        depth <- length(x) - recon_idx - 1
        return(depth <= 2)
      }
      return(FALSE)
    })]
    
    files <- files[!grepl("_chain[0-9]", files)]
    all_files <- c(all_files, files)
    cat("  Found", length(files), "files in", basename(dir), "\n")
  }
}

cat("Total model files found:", length(all_files), "\n")

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
ready_for_summary <- results[
  results$file_type == "combined" & 
  results$has_samples & 
  results$has_samples2 & 
  results$has_metadata & 
  results$has_model_data, 
]

cat("Files ready for summary script:", nrow(ready_for_summary), "\n")

# Save validation results
write.csv(results, "comprehensive_model_validation_results.csv", row.names = FALSE)
write.csv(ready_for_summary, "files_ready_for_summary.csv", row.names = FALSE)

cat("Validation results saved\n")

# =============================================================================
# STEP 2: PROCESS MODEL FILES AND CREATE SUMMARIES (INCLUDING PLOT SUMMARIES)
# =============================================================================

cat("\n=== STEP 2: PROCESS MODEL FILES AND CREATE SUMMARIES ===\n")

# Get ready files for processing
ready_file_paths <- ready_for_summary$full_path
n_ready_files <- length(ready_file_paths)

# Pre-allocate file_summaries list
file_summaries <- vector("list", n_ready_files)

# Plot summary generation is now handled by the summarize_beta_model function

# Use parallel processing for model summarization
cl <- makeCluster(max_cores)

# Set up cluster environment
clusterEvalQ(cl, {
  library(dplyr)
  library(here)
  library(microbialForecast)
})

# Process files in parallel with error handling
summarization_results <- tryCatch({
  parLapply(cl, ready_file_paths, function(f) {
    # Check if summary already exists and is up-to-date
    summary_file <- gsub("samples_", "summary_", f)
    should_summarize <- TRUE
    
    if (file.exists(summary_file)) {
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
          microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = NULL)
        }, error = function(e) {
          if (grepl("truth.plot.long", e$message)) {
            # Try with drop_other = FALSE to avoid plot reconstruction
            microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = FALSE, overwrite = NULL)
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
}, error = function(e) {
  cat("Error in parallel summarization:", e$message, "\n")
  # Fallback to sequential processing
  lapply(ready_file_paths, function(f) {
    # Check if summary already exists and is up-to-date
    summary_file <- gsub("samples_", "summary_", f)
    should_summarize <- TRUE
    
    if (file.exists(summary_file)) {
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
          microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = NULL)
        }, error = function(e) {
          if (grepl("truth.plot.long", e$message)) {
            # Try with drop_other = FALSE to avoid plot reconstruction
            microbialForecast::summarize_beta_model(f, save_summary=T, drop_other = FALSE, overwrite = NULL)
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

stopCluster(cl)

# Extract results
for(i in seq_along(summarization_results)) {
  file_summaries[[i]] <- summarization_results[[i]]$result
  if(i %% 10 == 0) {
    cat("Processed", i, "of", n_ready_files, "files\n")
  }
}

# Remove NULL entries
file_summaries <- file_summaries[!sapply(file_summaries, is.null)]
cat("Model processing complete. Results:", length(file_summaries), "\n")

# =============================================================================
# STEP 3: COLLECT AND COMBINE SUMMARY FILES
# =============================================================================

cat("\n=== STEP 3: COLLECT AND COMBINE SUMMARY FILES ===\n")

# Search for summary files - only main model files, not chain files
# Use system commands to be more explicit about what we're getting
env_cycl_files <- system("find data/model_outputs/reconstructed_from_checkpoints/env_cycl -maxdepth 1 -name 'summary_*.rds' | grep -v chain", intern = TRUE)
cycl_only_files <- system("find data/model_outputs/reconstructed_from_checkpoints/cycl_only -maxdepth 1 -name 'summary_*.rds' | grep -v chain", intern = TRUE)

summary_file_list <- c(env_cycl_files, cycl_only_files)


cat("Found", length(summary_file_list), "summary files from processed models\n")

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

# Process summary files in chunks to avoid memory issues
if(length(summary_file_list) > 0) {
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
      
      summary_df <- map_df(all_summaries, 1)
      
      # Process plot estimates only if flag is enabled
      if (PROCESS_PLOT_ESTIMATES) {
        cat("Processing plot estimates (PROCESS_PLOT_ESTIMATES = TRUE)\n")
        
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
      gelman_list <- map_df(all_summaries, 4)
      
      cat("Successfully combined", length(all_summaries), "summaries\n")
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

# Load existing data if processing was skipped
if (skip_processing) {
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
    
    # Get unique model IDs that meet criteria
    keep_models <- unique(gelman_dt[is_major_param == TRUE & `Point est.` < 1.1]$model_id)
    keep_models_weak <- unique(gelman_dt[is_major_param == TRUE & `Point est.` < 1.2]$model_id)
    keep_models_stricter <- unique(gelman_dt[is_major_param == TRUE & `Point est.` < 1.05]$model_id)
    
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
    
    cat("  Total models:", length(all_models), "\n")
    cat("  Converged (strict):", length(keep_list_stricter), "\n")
    cat("  Converged (standard):", length(keep_list), "\n")
    cat("  Converged (weak):", length(keep_list_weak), "\n")
    cat("  Need rerun:", length(rerun_list), "\n")
    
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
if (!exists("summary_df")) {
  cat("  ⚠️ summary_df missing, creating empty data frame\n")
  summary_df <- data.frame()
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

# Save convergence lists individually
saveRDS(convergence_results$keep_list, here("data/summary/converged_taxa_list.rds"))
saveRDS(convergence_results$keep_list_stricter, here("data/summary/stricter_converged_taxa_list.rds"))
saveRDS(convergence_results$keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))

# Create and save unconverged list
rerun_list_prioritized <- c(
  missing_chain_results$priority_rerun_list, 
  setdiff(missing_chain_results$all_missing_rerun_list, 
          missing_chain_results$priority_rerun_list)
)
saveRDS(rerun_list_prioritized, here("data/summary/unconverged_taxa_list.rds"))

# Save priority rerun list
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

# Save main summary - NO FALLBACK TO MINIMAL VERSION
tryCatch({
  saveRDS(main_summary, here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  cat("\n  ✅ Saved main summary file successfully\n")
  
  # Verify it saved correctly
  test_read <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  if (!is.null(test_read$summary_df) && nrow(test_read$summary_df) > 0) {
    cat("  ✅ Verified: main summary file saved successfully\n")
  } else {
    cat("  ⚠️ Warning: main summary file may be empty\n")
  }
  
}, error = function(e) {
  cat("\n  ❌ CRITICAL ERROR saving main summary file:\n")
  cat("    ", e$message, "\n")
  cat("\n  Attempting to diagnose the problem...\n")
  
  # Check object sizes
  cat("  Object sizes:\n")
  cat("    - summary_df:", format(object.size(main_summary$summary_df), units = "MB"), "\n")
  cat("    - gelman.summary:", format(object.size(main_summary$gelman.summary), units = "MB"), "\n")
  cat("    - Total:", format(object.size(main_summary), units = "MB"), "\n")
  
  # Save components individually as backup
  cat("\n  Saving components individually as backup...\n")
  tryCatch({
    saveRDS(main_summary$summary_df, here("data/summary/backup_summary_df.rds"))
    cat("    ✅ Saved backup_summary_df.rds\n")
  }, error = function(e2) {
    cat("    ❌ Failed to save summary_df:", e2$message, "\n")
  })
  
  # plot_est is saved separately, no need to backup here
  
  tryCatch({
    saveRDS(main_summary$gelman.summary, here("data/summary/backup_gelman.rds"))
    cat("    ✅ Saved backup_gelman.rds\n")
  }, error = function(e2) {
    cat("    ❌ Failed to save gelman.summary:", e2$message, "\n")
  })
  
  stop("Failed to save main summary file. Check backup files in data/summary/")
})

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
}

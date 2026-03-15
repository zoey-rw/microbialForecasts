#!/usr/bin/env Rscript

# 03c_extract_plot_summaries.r
# Standalone script to extract and combine plot summaries from summary files
# This script focuses specifically on the plot summary extraction that was causing crashes
# in the main 03_summarizeModelOutputs.r script
# Updated to handle driver uncertainty model types

source("../../source.R")
library(data.table)
if (require(arrow, quietly = TRUE)) {
  library(arrow)
  use_parquet <- TRUE
} else {
  use_parquet <- FALSE
  cat("Note: arrow package not available, using RDS format instead of Parquet\n")
}
library(here)

cat("=== EXTRACTING AND COMBINING PLOT SUMMARIES ===\n")
cat("Standalone script for plot summary extraction\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Memory management settings
options(memory.limit = 2000)  # 2GB limit for safety (more conservative)
setDTthreads(threads = 2)     # Reduce data.table parallelism to save memory

# Processing settings
CHUNK_SIZE <- 10              # Smaller chunks to avoid memory issues
MIN_ROWS_PER_FILE <- 1        # Minimum rows required to include a file
TEMP_DIR <- here("data/temp_plot_summaries")

cat("Configuration:\n")
cat("  - Chunk size:", CHUNK_SIZE, "files per chunk\n")
cat("  - Memory limit: 2GB\n")
cat("  - Data.table threads:", getDTthreads(), "\n")
cat("  - Temp directory:", TEMP_DIR, "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to standardize column names in plot estimates
standardize_plot_columns <- function(dt) {
  if (!is.data.table(dt)) {
    setDT(dt)
  }
  
  # Standardize quantile column names
  col_renames <- c(
    "X2.5." = "2.5%",
    "X25." = "25%", 
    "X50." = "50%",
    "X75." = "75%",
    "X97.5." = "97.5%"
  )
  
  for (old_name in names(col_renames)) {
    if (old_name %in% names(dt)) {
      setnames(dt, old_name, col_renames[[old_name]])
    }
  }
  
  # Add missing columns if needed
  if (!"SD" %in% names(dt)) {
    if (all(c("2.5%", "97.5%") %in% names(dt))) {
      dt[, SD := (`97.5%` - `2.5%`) / 3.92]
    } else {
      dt[, SD := NA_real_]
    }
  }
  
  if (!"Mean" %in% names(dt) && "50%" %in% names(dt)) {
    dt[, Mean := `50%`]
  }
  
  return(dt)
}

# Function to validate plot estimate data - simplified for sparse matrices
validate_plot_estimate <- function(plot_est) {
  if (!is.data.frame(plot_est) || nrow(plot_est) == 0) {
    return(FALSE)
  }
  
  # Check for some basic plot-related columns
  plot_cols <- c("plot_num", "timepoint", "siteID", "plotID")
  has_plot_info <- any(plot_cols %in% names(plot_est))
  
  if (!has_plot_info) {
    return(FALSE)
  }
  
  # Check for either quantile columns OR Mean column (sparse matrices are OK)
  has_quantiles <- any(c("2.5%", "25%", "50%", "75%", "97.5%", "X2.5.", "X25.", "X50.", "X75.", "X97.5.") %in% names(plot_est))
  has_mean <- "Mean" %in% names(plot_est)
  
  # Accept if we have either quantiles OR mean column
  return(has_quantiles || has_mean)
}

# Function to process a single summary file
process_single_summary_file <- function(file_path) {
  tryCatch({
    # Read the summary file
    summary_obj <- readRDS(file_path)
    
    # Check if it's a valid summary object
    if (!is.list(summary_obj) || length(summary_obj) < 2) {
      return(NULL)
    }
    
    # Extract plot estimates from element 2 (element 3 is just placeholder structure)
    plot_est <- summary_obj[[2]]
    
    # Clear the full summary object immediately to save memory
    rm(summary_obj)
    gc(verbose = FALSE)
    
    # Validate and process plot estimates
    if (!validate_plot_estimate(plot_est)) {
      return(NULL)
    }
    
    # Convert to data.table and standardize columns
    plot_est_dt <- as.data.table(plot_est)
    plot_est_dt <- standardize_plot_columns(plot_est_dt)
    
  # For sparse matrices, filter out completely empty rows (all NA)
  # This will dramatically reduce memory usage
  if (nrow(plot_est_dt) > 0) {
    # Check which rows have at least one non-NA value in key columns
    key_cols <- c("Mean", "2.5%", "25%", "50%", "75%", "97.5%", "X2.5.", "X25.", "X50.", "X75.", "X97.5.")
    available_cols <- key_cols[key_cols %in% names(plot_est_dt)]
    
    if (length(available_cols) > 0) {
      # Create a logical vector for rows with at least one non-NA value
      has_data <- rowSums(!is.na(plot_est_dt[, ..available_cols])) > 0
      
      # BUT: Keep rows that have metadata even if they don't have plot estimates
      # This preserves the structure needed for downstream analysis
      has_metadata <- !is.na(plot_est_dt$siteID) & !is.na(plot_est_dt$plotID) & !is.na(plot_est_dt$model_id)
      
      # Keep rows that either have data OR have metadata
      keep_rows <- has_data | has_metadata
      plot_est_dt <- plot_est_dt[keep_rows]
    }
  }
  
  # Check if we still have enough data after filtering
  if (nrow(plot_est_dt) < MIN_ROWS_PER_FILE) {
    return(NULL)
  }

  return(plot_est_dt)
    
  }, error = function(e) {
    cat("  Error processing file:", basename(file_path), "-", e$message, "\n")
    return(NULL)
  })
}

# Function to process files in chunks and save intermediate results
process_chunk_and_save_model_type <- function(chunk_files, chunk_idx, model_type, model_temp_dir) {
  cat("  Processing chunk", chunk_idx, "with", length(chunk_files), "files\n")
  
  # Process all files in this chunk
  chunk_results <- list()
  valid_count <- 0
  
  for (i in seq_along(chunk_files)) {
    file_path <- chunk_files[i]
    
    if (i %% 10 == 0) {
      cat("    File", i, "/", length(chunk_files), "\n")
    }
    
    result <- process_single_summary_file(file_path)
    
    if (!is.null(result)) {
      valid_count <- valid_count + 1
      chunk_results[[valid_count]] <- result
    }
  }
  
  # Combine results for this chunk
  if (valid_count > 0) {
    cat("    Combining", valid_count, "valid files in chunk", chunk_idx, "\n")
    
    # Use rbindlist for efficient combining
    chunk_combined <- rbindlist(chunk_results, fill = TRUE, use.names = TRUE)
    
    # Save intermediate result with model type prefix
    if (use_parquet) {
      temp_file <- file.path(model_temp_dir, paste0("plot_chunk_", model_type, "_", chunk_idx, ".parquet"))
      write_parquet(chunk_combined, temp_file)
    } else {
      temp_file <- file.path(model_temp_dir, paste0("plot_chunk_", model_type, "_", chunk_idx, ".rds"))
      saveRDS(chunk_combined, temp_file)
    }
    
    cat("    Saved chunk", chunk_idx, "with", nrow(chunk_combined), "rows\n")
    
    # Clear memory
    rm(chunk_results, chunk_combined)
    gc(verbose = FALSE)
    
    return(temp_file)
  } else {
    cat("    No valid data in chunk", chunk_idx, "\n")
    return(NULL)
  }
}

# =============================================================================
# MAIN PROCESSING
# =============================================================================

# Process by model type to avoid memory issues and allow resumability
model_types <- c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only")

# Check if specific model type is requested via command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  requested_type <- args[1]
  if (requested_type %in% model_types) {
    model_types <- requested_type
    cat("Processing only model type:", requested_type, "\n")
  } else {
    cat("Unknown model type:", requested_type, ". Processing all types.\n")
  }
}

# Find summary files by model type
cat("Searching for summary files by model type...\n")
all_summary_files <- list()

for (model_type in model_types) {
  # Handle different directory structures for different model types
  if (grepl("driver_uncertainty_", model_type)) {
    model_path <- file.path("data/model_outputs/cloglog_beta_driver_uncertainty", model_type)
  } else {
    model_path <- file.path("data/model_outputs/cloglog_beta_driver_uncertainty", model_type)
  }
  
  if (dir.exists(model_path)) {
    files <- list.files(
      path = model_path,
      pattern = "^summary_.*\\.rds$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    # Filter out chain files and other non-standard files
    files <- files[!grepl("_chain[0-9]+\\.rds$", files)]
    files <- files[!grepl("_chain[0-9]+_beta_regression", files)]
    
    all_summary_files[[model_type]] <- files
    cat("Found", length(files), "files for", model_type, "\n")
  } else {
    cat("Directory not found:", model_path, "\n")
    all_summary_files[[model_type]] <- character(0)
  }
}

# Collect all chunks from all model types for final combination
all_model_chunks <- list()

# Process each model type separately
for (model_type in model_types) {
  summary_files <- all_summary_files[[model_type]]
  
  if (length(summary_files) == 0) {
    cat("No files found for", model_type, "- skipping\n")
    next
  }
  
  cat("\n=== PROCESSING MODEL TYPE:", toupper(model_type), "===\n")
  cat("Files to process:", length(summary_files), "\n")

  # Create temporary directory for this model type
  model_temp_dir <- file.path(TEMP_DIR, model_type)
  if (!dir.exists(model_temp_dir)) {
    dir.create(model_temp_dir, recursive = TRUE)
  }
  
  # Process files in chunks
  cat("Processing files in chunks of", CHUNK_SIZE, "...\n")
  n_chunks <- ceiling(length(summary_files) / CHUNK_SIZE)
  chunk_files <- list()
  
  for (chunk_idx in 1:n_chunks) {
    start_idx <- (chunk_idx - 1) * CHUNK_SIZE + 1
    end_idx <- min(chunk_idx * CHUNK_SIZE, length(summary_files))
    chunk_file_list <- summary_files[start_idx:end_idx]
    
    temp_file <- process_chunk_and_save_model_type(chunk_file_list, chunk_idx, model_type, model_temp_dir)
    
    if (!is.null(temp_file)) {
      chunk_files[[length(chunk_files) + 1]] <- temp_file
    }
    
    # Force garbage collection between chunks
    gc(verbose = FALSE)
  }
  
  cat("Processed", length(chunk_files), "valid chunks for", model_type, "\n")

  # Store chunks for later combination across all model types
  if (length(chunk_files) > 0) {
    all_model_chunks[[model_type]] <- chunk_files
  }
  
  # Clean up temporary directory for this model type
  unlink(model_temp_dir, recursive = TRUE)
  cat("  ✅ Cleaned up temporary files for", model_type, "\n")
}  # End of model type processing loop

# =============================================================================
# COMBINE ALL CHUNKS FROM ALL MODEL TYPES (matching main script output)
# =============================================================================

cat("\n=== COMBINING ALL PLOT ESTIMATES FROM ALL MODEL TYPES ===\n")

# Collect all chunk files from all model types
all_chunk_files <- unlist(all_model_chunks)

if (length(all_chunk_files) == 0) {
  cat("No plot estimate chunks found to combine\n")
} else {
  cat("Found", length(all_chunk_files), "chunks across all model types\n")
  cat("Combining in batches to avoid memory issues...\n")
  
  # Process chunks in batches for final combination
  COMBINE_BATCH_SIZE <- 10  # Smaller batches for final combination
  n_combine_batches <- ceiling(length(all_chunk_files) / COMBINE_BATCH_SIZE)
  combined_batches <- list()
  
  for (batch_idx in 1:n_combine_batches) {
    start_idx <- (batch_idx - 1) * COMBINE_BATCH_SIZE + 1
    end_idx <- min(batch_idx * COMBINE_BATCH_SIZE, length(all_chunk_files))
    batch_files <- all_chunk_files[start_idx:end_idx]
    
    cat("  Combining batch", batch_idx, "of", n_combine_batches, "(", length(batch_files), "chunks)\n")
    
    # Load and combine this batch
    batch_data <- list()
    for (chunk_file in batch_files) {
      if (grepl("\\.parquet$", chunk_file) && use_parquet) {
        chunk_data <- read_parquet(chunk_file)
      } else {
        chunk_data <- readRDS(chunk_file)
      }
      batch_data[[length(batch_data) + 1]] <- chunk_data
    }
    
    if (length(batch_data) > 0) {
      batch_combined <- rbindlist(batch_data, fill = TRUE, use.names = TRUE)
      combined_batches[[batch_idx]] <- batch_combined
      
      # Clear memory
      rm(batch_data, batch_combined)
      gc(verbose = FALSE)
    }
  }
  
  # Combine all batches
  if (length(combined_batches) > 0) {
    cat("  Combining", length(combined_batches), "batches into final plot estimates...\n")
    plot_est_combined <- rbindlist(combined_batches, fill = TRUE, use.names = TRUE)
    
    cat("  ✅ Combined", nrow(plot_est_combined), "plot estimate rows from all model types\n")
    
    # Clear memory
    rm(combined_batches)
    gc(verbose = FALSE)
    
    # Save in same formats as main script
    if (nrow(plot_est_combined) > 0) {
      cat("\n  Saving plot estimates in multiple formats (matching main script)...\n")
      
      # Ensure it's a data frame
      plot_est_df <- as.data.frame(plot_est_combined)
      
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
      
      # Clear memory
      rm(plot_est_combined, plot_est_df)
      gc(verbose = FALSE)
      
      cat("\n  ✅ Plot estimates saved to data/summary/plot_estimates.{rds,parquet,csv.gz}\n")
      cat("  ✅ Output matches main script format\n")
    } else {
      cat("  ⚠️ No plot estimates to save (empty data frame)\n")
    }
  }
}

# =============================================================================
# CLEANUP
# =============================================================================

cat("\nCleaning up temporary files...\n")
if (dir.exists(TEMP_DIR)) {
  unlink(TEMP_DIR, recursive = TRUE)
  cat("  ✅ Cleaned up temporary directory\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=== PLOT SUMMARY EXTRACTION COMPLETE ===\n")
cat("Summary:\n")

total_files_processed <- 0
for (model_type in model_types) {
  files_count <- length(all_summary_files[[model_type]])
  total_files_processed <- total_files_processed + files_count
  cat("  -", model_type, "files processed:", files_count, "\n")
}

cat("  - Total input files processed:", total_files_processed, "\n")
cat("  - Output files saved to:", here("data/summary"), "\n")
cat("  - Model types processed:", paste(model_types, collapse = ", "), "\n")

cat("\n✅ SUCCESS: Plot summaries extracted and combined successfully!\n")
cat("Output files (matching main script):\n")
cat("  - data/summary/plot_estimates.rds\n")
cat("  - data/summary/plot_estimates.parquet\n")
cat("  - data/summary/plot_estimates.csv.gz\n")

cat("\nThis script can be run independently of the main 03_summarizeModelOutputs.r\n")
cat("The plot estimates are now available in the same format as the main script.\n")
cat("\nTo process a specific model type only, run:\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r env_cycl\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r cycl_only\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r env_cov\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r driver_uncertainty_env_cycl\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r driver_uncertainty_env_cov\n")
cat("  Rscript analysis/model_analysis/03c_extract_plot_summaries.r driver_uncertainty_cycl_only\n")

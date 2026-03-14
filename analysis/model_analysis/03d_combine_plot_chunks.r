#!/usr/bin/env Rscript

# 03d_combine_plot_chunks.r
# Script to combine plot estimate chunks on-demand
# This allows other scripts to load the combined data when needed

source("source.R")
library(data.table)
library(arrow)
library(here)

cat("=== COMBINING PLOT ESTIMATE CHUNKS ===\n")
cat("On-demand combination of plot estimate chunks\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Memory management settings
options(memory.limit = 1000)  # 1GB limit for safety
setDTthreads(threads = 1)     # Reduce data.table parallelism to save memory

# Processing settings
BATCH_SIZE <- 2  # Process chunks in very small batches
TEMP_DIR <- here("data/temp_plot_combination")

cat("Configuration:\n")
cat("  - Batch size:", BATCH_SIZE, "chunks per batch\n")
cat("  - Memory limit: 2GB\n")
cat("  - Data.table threads:", getDTthreads(), "\n")
cat("  - Temp directory:", TEMP_DIR, "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to combine chunks for a specific model type
combine_chunks_for_model <- function(model_type) {
  cat("=== COMBINING CHUNKS FOR", toupper(model_type), "===\n")
  
  # Find chunk directory
  chunk_dir <- file.path(here("data/summary"), paste0("plot_estimates_", model_type))
  
  if (!dir.exists(chunk_dir)) {
    cat("  ❌ Chunk directory not found:", chunk_dir, "\n")
    return(NULL)
  }
  
  # Find all chunk files
  chunk_files <- list.files(chunk_dir, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(chunk_files) == 0) {
    cat("  ❌ No chunk files found in:", chunk_dir, "\n")
    return(NULL)
  }
  
  cat("  Found", length(chunk_files), "chunk files\n")
  
  # Create temporary directory
  if (!dir.exists(TEMP_DIR)) {
    dir.create(TEMP_DIR, recursive = TRUE)
  }
  
  # Process chunks in batches
  n_batches <- ceiling(length(chunk_files) / BATCH_SIZE)
  cat("  Processing", length(chunk_files), "chunks in", n_batches, "batches of", BATCH_SIZE, "\n")
  
  # Create a temporary file for the final combined result
  final_temp_file <- file.path(TEMP_DIR, paste0("final_combined_", model_type, ".parquet"))
  
  # Process first batch to initialize the final file
  first_batch <- chunk_files[1:min(BATCH_SIZE, length(chunk_files))]
  cat("  Processing first batch...\n")
  
  all_chunks <- list()
  for (i in seq_along(first_batch)) {
    chunk_data <- read_parquet(first_batch[[i]])
    all_chunks[[i]] <- chunk_data
  }
  
  combined_data <- rbindlist(all_chunks, fill = TRUE, use.names = TRUE)
  write_parquet(combined_data, final_temp_file)
  
  # Clear memory
  rm(all_chunks, combined_data)
  gc(verbose = FALSE)
  
  # Process remaining batches and append to the final file
  if (n_batches > 1) {
    for (batch_idx in 2:n_batches) {
      start_idx <- (batch_idx - 1) * BATCH_SIZE + 1
      end_idx <- min(batch_idx * BATCH_SIZE, length(chunk_files))
      batch_files <- chunk_files[start_idx:end_idx]
      
      cat("  Processing batch", batch_idx, "/", n_batches, "(", length(batch_files), "files)\n")
    
    # Load batch
    all_chunks <- list()
    for (i in seq_along(batch_files)) {
      chunk_data <- read_parquet(batch_files[[i]])
      all_chunks[[i]] <- chunk_data
    }
    
    batch_combined <- rbindlist(all_chunks, fill = TRUE, use.names = TRUE)
    
    # Read existing data and combine
    existing_data <- read_parquet(final_temp_file)
    combined_data <- rbindlist(list(existing_data, batch_combined), fill = TRUE, use.names = TRUE)
    
    # Save updated combined data
    write_parquet(combined_data, final_temp_file)
    
    # Clear memory
    rm(all_chunks, batch_combined, existing_data, combined_data)
    gc(verbose = FALSE)
    }
  }
  
  # Get row count without loading full data
  cat("  Getting final row count for", model_type, "...\n")
  temp_data <- read_parquet(final_temp_file)
  total_rows <- nrow(temp_data)
  rm(temp_data)
  gc(verbose = FALSE)
  
  # Move the temporary file to the final location
  final_file <- file.path(here("data/summary"), paste0("plot_estimates_", model_type, ".parquet"))
  file.rename(final_temp_file, final_file)
  
  cat("  ✅ Combined successfully:", total_rows, "rows for", model_type, "\n")
  cat("  ✅ Saved to:", final_file, "\n")
  
  return(NULL)  # Don't return the data to save memory
}

# Function to save combined data
save_combined_data <- function(plot_estimates, model_type) {
  cat("\nSaving combined plot estimates for", model_type, "...\n")
  
  # Create output directory
  output_dir <- here("data/summary")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (nrow(plot_estimates) > 0) {
    cat("Saving", nrow(plot_estimates), "plot estimate rows for", model_type, "...\n")
    
    # Save in multiple formats for redundancy
    success_count <- 0
    
    # Format 1: Parquet (fastest, smallest)
    tryCatch({
      parquet_file <- file.path(output_dir, paste0("plot_estimates_", model_type, ".parquet"))
      write_parquet(plot_estimates, parquet_file)
      file_size_mb <- round(file.info(parquet_file)$size / 1048576, 2)
      cat("  ✅ Saved as Parquet:", file_size_mb, "MB\n")
      success_count <- success_count + 1
    }, error = function(e) {
      cat("  ❌ Failed to save Parquet:", e$message, "\n")
    })
    
    # Format 2: Compressed CSV (good compatibility)
    tryCatch({
      csv_file <- file.path(output_dir, paste0("plot_estimates_", model_type, ".csv.gz"))
      fwrite(plot_estimates, csv_file, compress = "gzip")
      file_size_mb <- round(file.info(csv_file)$size / 1048576, 2)
      cat("  ✅ Saved as compressed CSV:", file_size_mb, "MB\n")
      success_count <- success_count + 1
    }, error = function(e) {
      cat("  ❌ Failed to save CSV:", e$message, "\n")
    })
    
    # Format 3: RDS (R compatibility)
    tryCatch({
      rds_file <- file.path(output_dir, paste0("plot_estimates_", model_type, ".rds"))
      saveRDS(as.data.frame(plot_estimates), rds_file, compress = "xz")
      file_size_mb <- round(file.info(rds_file)$size / 1048576, 2)
      cat("  ✅ Saved as RDS:", file_size_mb, "MB\n")
      success_count <- success_count + 1
    }, error = function(e) {
      cat("  ❌ Failed to save RDS:", e$message, "\n")
    })
    
    if (success_count == 0) {
      cat("  ⚠️ Failed to save plot estimates for", model_type, "in any format\n")
    }
    
    # Print summary statistics
    cat("\n=== SUMMARY STATISTICS FOR", toupper(model_type), "===\n")
    cat("Total rows:", nrow(plot_estimates), "\n")
    cat("Memory usage:", format(object.size(plot_estimates), units = "MB"), "\n")
    
    if ("taxon" %in% names(plot_estimates)) {
      cat("Unique taxa:", length(unique(plot_estimates$taxon)), "\n")
    }
    if ("model_name" %in% names(plot_estimates)) {
      cat("Unique models:", length(unique(plot_estimates$model_name)), "\n")
    }
    if ("siteID" %in% names(plot_estimates)) {
      cat("Unique sites:", length(unique(plot_estimates$siteID)), "\n")
    }
    if ("plotID" %in% names(plot_estimates)) {
      cat("Unique plots:", length(unique(plot_estimates$plotID)), "\n")
    }
    
  } else {
    cat("No plot estimates to save for", model_type, "(empty dataset)\n")
  }
}

# =============================================================================
# MAIN PROCESSING
# =============================================================================

# Check if specific model type is requested via command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  requested_type <- args[1]
  if (requested_type %in% c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only")) {
    model_types <- requested_type
    cat("Combining only model type:", requested_type, "\n")
  } else {
    cat("Unknown model type:", requested_type, ". Combining all types.\n")
    model_types <- c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only")
  }
} else {
  model_types <- c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only")
  cat("Combining all model types\n")
}

# Process each model type
all_model_type_files <- list()
for (model_type in model_types) {
  cat("\n")
  result <- combine_chunks_for_model(model_type)
  
  # Collect combined files for each model type
  combined_file <- file.path(here("data/summary"), paste0("plot_estimates_", model_type, ".parquet"))
  if (file.exists(combined_file)) {
    all_model_type_files[[model_type]] <- combined_file
  }
  
  # Clear memory between model types
  gc(verbose = FALSE)
}

# =============================================================================
# OPTIONALLY COMBINE ALL MODEL TYPES INTO ONE FILE (matching main script output)
# =============================================================================

cat("\n=== OPTIONAL: COMBINING ALL MODEL TYPES INTO ONE FILE ===\n")
cat("This will create plot_estimates.rds/parquet/csv.gz matching main script output\n")
cat("Combining", length(all_model_type_files), "model type files...\n")

if (length(all_model_type_files) > 0) {
  # Load and combine all model types
  all_model_data <- list()
  for (model_type in names(all_model_type_files)) {
    file_path <- all_model_type_files[[model_type]]
    cat("  Loading", model_type, "...\n")
    model_data <- read_parquet(file_path)
    all_model_data[[model_type]] <- model_data
    rm(model_data)
    gc(verbose = FALSE)
  }
  
  # Combine all model types
  cat("  Combining all model types...\n")
  plot_est_all <- rbindlist(all_model_data, fill = TRUE, use.names = TRUE)
  
  cat("  ✅ Combined", nrow(plot_est_all), "plot estimate rows from all model types\n")
  
  # Clear memory
  rm(all_model_data)
  gc(verbose = FALSE)
  
  # Save in same formats as main script
  if (nrow(plot_est_all) > 0) {
    cat("\n  Saving plot estimates in multiple formats (matching main script)...\n")
    
    # Ensure it's a data frame
    plot_est_df <- as.data.frame(plot_est_all)
    
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
    rm(plot_est_all, plot_est_df)
    gc(verbose = FALSE)
    
    cat("\n  ✅ Combined plot estimates saved to data/summary/plot_estimates.{rds,parquet,csv.gz}\n")
    cat("  ✅ Output matches main script format\n")
  } else {
    cat("  ⚠️ No plot estimates to save (empty data frame)\n")
  }
} else {
  cat("  ⚠️ No model type files found to combine\n")
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

cat("\n=== PLOT CHUNK COMBINATION COMPLETE ===\n")
cat("Combined plot estimates are now available in data/summary/\n")
cat("Files created (matching main script):\n")
cat("  - plot_estimates.rds (combined from all model types)\n")
cat("  - plot_estimates.parquet (combined from all model types)\n")
cat("  - plot_estimates.csv.gz (combined from all model types)\n")
cat("\nPer-model-type files (if needed):\n")
for (model_type in model_types) {
  cat("  - plot_estimates_", model_type, ".parquet\n", sep = "")
}

cat("\nThese files can now be used by:\n")
cat("  - Script 6: Hindcast initialization\n")
cat("  - Script 7: Scoring metrics\n")
cat("  - Script 9: Phenology functionality\n")
cat("  - All scripts expecting plot_estimates.rds/parquet/csv.gz\n")

cat("\nTo combine specific model types only, run:\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r env_cycl\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r cycl_only\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r env_cov\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r driver_uncertainty_env_cycl\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r driver_uncertainty_env_cov\n")
cat("  Rscript analysis/model_analysis/03d_combine_plot_chunks.r driver_uncertainty_cycl_only\n")


#!/usr/bin/env Rscript

# Script to combine existing summary files from logit_beta_regression directories
# This script efficiently combines existing summary files with proper memory management

source("../../source.R")
library(dplyr)
library(tidyr)

cat("=== COMBINING EXISTING SUMMARY FILES ===\n")
cat("Collecting summary files from cloglog_beta_driver_uncertainty directories...\n")

# Base path
base_path <- here("data/model_outputs/cloglog_beta_driver_uncertainty")

# Find all summary files in subdirectories (recursive search to handle subdirectories)
# Use recursive = FALSE to avoid picking up stale duplicates in per-taxon subdirectories
env_cycl_files <- list.files(file.path(base_path, "env_cycl"),
                            pattern = "summary_env_cycl_.*_beta_regression\\.rds$",
                            recursive = FALSE,
                            full.names = TRUE)
cycl_only_files <- list.files(file.path(base_path, "cycl_only"),
                             pattern = "summary_cycl_only_.*_beta_regression\\.rds$",
                             recursive = FALSE,
                             full.names = TRUE)
env_cov_files <- list.files(file.path(base_path, "env_cov"),
                           pattern = "summary_env_cov_.*_beta_regression\\.rds$",
                           recursive = FALSE,
                           full.names = TRUE)

all_summary_files <- c(env_cycl_files, cycl_only_files, env_cov_files)
cat("Found", length(all_summary_files), "summary files from processed models\n")

# Quick read and combine function with chunking for memory efficiency
combine_summaries <- function(files) {
  cat("Reading and combining", length(files), "files...\n")
  
  # Process in chunks to avoid memory issues with large file sets
  CHUNK_SIZE <- 50  # Process 50 files at a time
  n_chunks <- ceiling(length(files) / CHUNK_SIZE)
  
  if (n_chunks > 1) {
    cat("Processing in", n_chunks, "chunks of", CHUNK_SIZE, "files each...\n")
  }
  
  all_summaries <- list()
  
  for (chunk_idx in 1:n_chunks) {
    start_idx <- (chunk_idx - 1) * CHUNK_SIZE + 1
    end_idx <- min(chunk_idx * CHUNK_SIZE, length(files))
    chunk_files <- files[start_idx:end_idx]
    
    if (n_chunks > 1) {
      cat("  Processing chunk", chunk_idx, "of", n_chunks, "(", length(chunk_files), "files)\n")
    }
    
    # Read files in this chunk
    chunk_summaries <- lapply(chunk_files, function(f) {
    tryCatch({
      readRDS(f)
    }, error = function(e) {
      cat("Error reading", basename(f), "\n")
      NULL
    })
  })
  
    # Remove NULL entries and add to master list
    chunk_summaries <- chunk_summaries[!sapply(chunk_summaries, is.null)]
    all_summaries <- c(all_summaries, chunk_summaries)
    
    # Clean up chunk memory
    rm(chunk_summaries)
    gc(verbose = FALSE)
  }
  
  summaries <- all_summaries
  cat("Successfully read", length(summaries), "files\n")
  
  # Extract and combine summary_df (element 1) in chunks
  cat("Extracting and combining summary_df...\n")
  summary_dfs <- list()
  
  chunk_size_df <- 50
  n_chunks_df <- ceiling(length(summaries) / chunk_size_df)
  
  for (chunk_idx in 1:n_chunks_df) {
    start_idx <- (chunk_idx - 1) * chunk_size_df + 1
    end_idx <- min(chunk_idx * chunk_size_df, length(summaries))
    chunk_summaries <- summaries[start_idx:end_idx]
    
    chunk_summary_dfs <- lapply(chunk_summaries, function(x) {
    if (length(x) >= 1 && is.data.frame(x[[1]]) && nrow(x[[1]]) > 0) {
      return(x[[1]])
    } else {
      return(data.frame())
    }
  })
  
    chunk_combined <- bind_rows(chunk_summary_dfs)
    if (nrow(chunk_combined) > 0) {
      summary_dfs[[length(summary_dfs) + 1]] <- chunk_combined
    }
    
    rm(chunk_summaries, chunk_summary_dfs, chunk_combined)
    gc(verbose = FALSE)
  }
  
  # Combine all chunks
  combined_summary_df <- bind_rows(summary_dfs)
  cat("Combined summary_df:", nrow(combined_summary_df), "rows\n")
  
  # Extract and combine gelman summaries (element 4) in chunks
  cat("Extracting and combining gelman.summary...\n")
  gelman_dfs <- list()
  
  for (chunk_idx in 1:n_chunks_df) {
    start_idx <- (chunk_idx - 1) * chunk_size_df + 1
    end_idx <- min(chunk_idx * chunk_size_df, length(summaries))
    chunk_summaries <- summaries[start_idx:end_idx]
    
    chunk_gelman_dfs <- lapply(chunk_summaries, function(x) {
    if (length(x) >= 4 && is.data.frame(x[[4]]) && nrow(x[[4]]) > 0) {
      return(x[[4]])
    } else {
      return(data.frame())
    }
  })
  
    chunk_combined <- bind_rows(chunk_gelman_dfs)
    if (nrow(chunk_combined) > 0) {
      gelman_dfs[[length(gelman_dfs) + 1]] <- chunk_combined
    }
    
    rm(chunk_summaries, chunk_gelman_dfs, chunk_combined)
    gc(verbose = FALSE)
  }
  
  # Combine all gelman chunks
  combined_gelman <- bind_rows(gelman_dfs)
  cat("Combined gelman.summary:", nrow(combined_gelman), "rows\n")
  
  # Clean up intermediate lists
  rm(summaries, summary_dfs, gelman_dfs)
  gc(verbose = FALSE)
  
  # Process gelman data for model categorization (matching main script criteria)
  cat("Processing gelman data for model categorization...\n")
  
  # Convert to data.table for faster operations (matching main script)
  library(data.table)
  gelman_dt <- as.data.table(combined_gelman)
  
  # Filter and add major parameter flag (matching main script)
  if("model_id" %in% colnames(gelman_dt)) {
    gelman_dt <- gelman_dt[!grepl("all_covariates", model_id)]
  }
  
  if("parameter" %in% colnames(gelman_dt)) {
    gelman_dt[, is_major_param := grepl("beta|int|sigma|core_sd|rho", parameter)]
  } else {
    gelman_dt[, is_major_param := TRUE]
  }
  
  # Function to ensure model_ids have _beta_regression suffix to match filenames (matching main script)
  fix_model_id <- function(model_id) {
    if (grepl("with_legacy_covariate$", model_id) && !grepl("_beta_regression$", model_id)) {
      return(paste0(model_id, "_beta_regression"))
    }
    return(model_id)
  }
  
  # Apply simple convergence criteria (matching main script)
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
  
  # Convert back to data.frame for gelman.summary (for compatibility)
  gelman.summary <- as.data.frame(gelman_dt)
  
  # Remove is_major_param column if it was added (to match main script output)
  if("is_major_param" %in% colnames(gelman.summary)) {
    gelman.summary <- gelman.summary[, colnames(gelman.summary) != "is_major_param"]
  }

  cat("Model categorization complete:\n")
  cat("keep_list length:", length(keep_list), "\n")
  cat("keep_list_weak length:", length(keep_list_weak), "\n")
  cat("rerun_list length:", length(rerun_list), "\n")
  cat("Total models categorized:", length(keep_list_weak) + length(rerun_list), "\n")
  
  # Skip plot estimates for now - will be processed in scripts 3c and 3d
  cat("Skipping plot estimates - will be processed separately\n")
  
  # Create final object with the correct structure (matching main script)
  result <- list(
    summary_df = combined_summary_df,
    plot_est = data.frame(),  # Empty - will be processed in 3c
    gelman.summary = gelman.summary,
    keep_models_weak = if(exists("keep_models_weak")) keep_models_weak else character(0),
    keep_list = keep_list,
    keep_list_weak = keep_list_weak,
    keep_list_stricter = keep_list_stricter,
    rerun_list = rerun_list,
    priority_rerun_list = character(0),
    missing_chains_analysis = list()
  )
  
  result
}

# Combine summary files for all taxa calibration models
if(length(all_summary_files) > 0) {
  main_summary <- combine_summaries(all_summary_files)
  
  # Save the combined summaries
  cat("Saving combined summaries...\n")
  saveRDS(main_summary, here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  cat("Main summaries saved as logit_beta_fixed_priors_summaries.rds\n")
  
  # Save convergence lists individually (matching main script)
  cat("Saving convergence lists...\n")
  if (!dir.exists(here("data/summary"))) {
    dir.create(here("data/summary"), recursive = TRUE)
  }
  
  saveRDS(main_summary$keep_list, here("data/summary/converged_taxa_list.rds"))
  saveRDS(main_summary$keep_list_stricter, here("data/summary/stricter_converged_taxa_list.rds"))
  saveRDS(main_summary$keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))
  saveRDS(main_summary$rerun_list, here("data/summary/unconverged_taxa_list.rds"))
  
  cat("  ✅ Saved converged_taxa_list.rds\n")
  cat("  ✅ Saved stricter_converged_taxa_list.rds\n")
  cat("  ✅ Saved weak_converged_taxa_list.rds\n")
  cat("  ✅ Saved unconverged_taxa_list.rds\n")
} else {
  cat("No valid summary data found to combine\n")
}

cat("Summary combination complete!\n")
cat("Summary:\n")
cat("  - Total summary files processed:", length(all_summary_files), "\n")
cat("  - Models categorized:", length(main_summary$keep_list_weak) + length(main_summary$rerun_list), "\n")
cat("  - All required outputs saved (matching main script format)\n")
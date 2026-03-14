#!/usr/bin/env Rscript

# 03e_load_plot_estimates.r
# Functions to load plot estimates on-demand from chunks
# This avoids memory issues by loading only what's needed

source("source.R")
library(data.table)
library(arrow)
library(here)

# =============================================================================
# LOAD PLOT ESTIMATES FUNCTIONS
# =============================================================================

#' Load plot estimates for a specific model type
#' @param model_type Character string: "env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", or "driver_uncertainty_cycl_only"
#' @param subset_cols Character vector: specific columns to load (optional)
#' @param filter_expr Expression: data.table filter expression (optional)
#' @param max_chunks Integer: maximum number of chunks to load (optional)
#' @return data.table with plot estimates
load_plot_estimates <- function(model_type, subset_cols = NULL, filter_expr = NULL, max_chunks = NULL) {
  
  # Find chunk directory
  chunk_dir <- file.path(here("data/summary"), paste0("plot_estimates_", model_type))
  
  if (!dir.exists(chunk_dir)) {
    stop("Chunk directory not found: ", chunk_dir)
  }
  
  # Find all chunk files
  chunk_files <- list.files(chunk_dir, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(chunk_files) == 0) {
    stop("No chunk files found in: ", chunk_dir)
  }
  
  # Limit chunks if requested
  if (!is.null(max_chunks) && max_chunks < length(chunk_files)) {
    chunk_files <- chunk_files[1:max_chunks]
    cat("Loading only first", max_chunks, "chunks out of", length(chunk_files), "total\n")
  }
  
  cat("Loading plot estimates for", model_type, "from", length(chunk_files), "chunks...\n")
  
  # Load chunks one by one and combine
  all_data <- list()
  
  for (i in seq_along(chunk_files)) {
    if (i %% 5 == 0) {
      cat("  Loading chunk", i, "/", length(chunk_files), "\n")
    }
    
    # Load chunk
    chunk_data <- read_parquet(chunk_files[[i]])
    
    # Subset columns if requested
    if (!is.null(subset_cols)) {
      available_cols <- subset_cols[subset_cols %in% names(chunk_data)]
      if (length(available_cols) > 0) {
        chunk_data <- chunk_data[, ..available_cols]
      }
    }
    
    # Apply filter if requested
    if (!is.null(filter_expr)) {
      chunk_data <- chunk_data[eval(filter_expr)]
    }
    
    # Only keep if there's data
    if (nrow(chunk_data) > 0) {
      all_data[[length(all_data) + 1]] <- chunk_data
    }
    
    # Clear memory
    rm(chunk_data)
    gc(verbose = FALSE)
  }
  
  # Combine all chunks
  if (length(all_data) > 0) {
    cat("  Combining", length(all_data), "chunks...\n")
    combined_data <- rbindlist(all_data, fill = TRUE, use.names = TRUE)
    rm(all_data)
    gc(verbose = FALSE)
    
    cat("  ✅ Loaded", nrow(combined_data), "rows for", model_type, "\n")
    return(combined_data)
  } else {
    cat("  ⚠️ No data found for", model_type, "\n")
    return(data.table())
  }
}

#' Load plot estimates for multiple model types
#' @param model_types Character vector: model types to load
#' @param subset_cols Character vector: specific columns to load (optional)
#' @param filter_expr Expression: data.table filter expression (optional)
#' @param max_chunks Integer: maximum number of chunks to load per model type (optional)
#' @return data.table with plot estimates from all model types
load_all_plot_estimates <- function(model_types = c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only"), 
                                   subset_cols = NULL, 
                                   filter_expr = NULL, 
                                   max_chunks = NULL) {
  
  all_data <- list()
  
  for (model_type in model_types) {
    cat("\n=== LOADING", toupper(model_type), "===\n")
    
    model_data <- load_plot_estimates(model_type, subset_cols, filter_expr, max_chunks)
    
    if (nrow(model_data) > 0) {
      # Add model type column
      model_data[, model_type := model_type]
      all_data[[length(all_data) + 1]] <- model_data
    }
    
    # Clear memory
    rm(model_data)
    gc(verbose = FALSE)
  }
  
  if (length(all_data) > 0) {
    cat("\n=== COMBINING ALL MODEL TYPES ===\n")
    combined_data <- rbindlist(all_data, fill = TRUE, use.names = TRUE)
    rm(all_data)
    gc(verbose = FALSE)
    
    cat("✅ Total combined rows:", nrow(combined_data), "\n")
    return(combined_data)
  } else {
    cat("⚠️ No data found for any model type\n")
    return(data.table())
  }
}

#' Get summary statistics for plot estimates without loading all data
#' @param model_type Character string: "env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", or "driver_uncertainty_cycl_only"
#' @return list with summary statistics
get_plot_estimates_summary <- function(model_type) {
  
  # Find chunk directory
  chunk_dir <- file.path(here("data/summary"), paste0("plot_estimates_", model_type))
  
  if (!dir.exists(chunk_dir)) {
    return(list(error = "Chunk directory not found"))
  }
  
  # Find all chunk files
  chunk_files <- list.files(chunk_dir, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(chunk_files) == 0) {
    return(list(error = "No chunk files found"))
  }
  
  # Load first chunk to get structure
  first_chunk <- read_parquet(chunk_files[1])
  
  # Get file sizes
  file_sizes <- file.info(chunk_files)$size
  total_size_mb <- sum(file_sizes) / 1048576
  
  # Estimate total rows
  rows_per_chunk <- nrow(first_chunk)
  estimated_total_rows <- rows_per_chunk * length(chunk_files)
  
  return(list(
    model_type = model_type,
    n_chunks = length(chunk_files),
    estimated_total_rows = estimated_total_rows,
    rows_per_chunk = rows_per_chunk,
    total_size_mb = round(total_size_mb, 2),
    columns = names(first_chunk),
    n_columns = length(names(first_chunk))
  ))
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if (interactive() || length(commandArgs(trailingOnly = TRUE)) >= 0) {
  
  cat("=== PLOT ESTIMATES LOADING FUNCTIONS ===\n")
  cat("Available functions:\n")
  cat("  - load_plot_estimates(model_type, subset_cols, filter_expr, max_chunks)\n")
  cat("  - load_all_plot_estimates(model_types, subset_cols, filter_expr, max_chunks)\n")
  cat("  - get_plot_estimates_summary(model_type)\n\n")
  
  # Show summary for both model types
  for (model_type in c("env_cycl", "cycl_only", "env_cov", "driver_uncertainty_env_cycl", "driver_uncertainty_env_cov", "driver_uncertainty_cycl_only")) {
    cat("=== SUMMARY FOR", toupper(model_type), "===\n")
    summary_info <- get_plot_estimates_summary(model_type)
    
    if ("error" %in% names(summary_info)) {
      cat("  ❌", summary_info$error, "\n")
    } else {
      cat("  Chunks:", summary_info$n_chunks, "\n")
      cat("  Estimated rows:", summary_info$estimated_total_rows, "\n")
      cat("  Rows per chunk:", summary_info$rows_per_chunk, "\n")
      cat("  Total size:", summary_info$total_size_mb, "MB\n")
      cat("  Columns:", summary_info$n_columns, "\n")
      cat("  Sample columns:", paste(head(summary_info$columns, 5), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  cat("Example usage:\n")
  cat("  # Load all data for env_cycl\n")
  cat("  env_data <- load_plot_estimates('env_cycl')\n\n")
  
  cat("  # Load only specific columns\n")
  cat("  subset_data <- load_plot_estimates('env_cycl', subset_cols = c('siteID', 'plotID', 'Mean'))\n\n")
  
  cat("  # Load with filter\n")
  cat("  filtered_data <- load_plot_estimates('env_cycl', filter_expr = quote(siteID == 'BART'))\n\n")
  
  cat("  # Load only first 5 chunks\n")
  cat("  sample_data <- load_plot_estimates('env_cycl', max_chunks = 5)\n\n")
  
  cat("  # Load all model types\n")
  cat("  all_data <- load_all_plot_estimates()\n\n")
}

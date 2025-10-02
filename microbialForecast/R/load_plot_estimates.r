#' Load plot estimates from chunked Parquet files
#' 
#' This function loads plot estimates from the chunked Parquet files created by
#' the plot estimates extraction pipeline. It handles memory efficiently by
#' loading chunks incrementally and combining them.
#' 
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @param subset_cols Character vector: specific columns to load (optional)
#' @param filter_expr Expression: data.table filter expression (optional)
#' @param max_chunks Integer: maximum number of chunks to load (optional)
#' @return data.table with plot estimates
#' @export
load_plot_estimates <- function(model_type, subset_cols = NULL, filter_expr = NULL, max_chunks = NULL) {
  
  # Find chunk directory
  chunk_dir <- file.path(here::here("data/summary"), paste0("plot_estimates_", model_type))
  
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
    chunk_data <- arrow::read_parquet(chunk_files[[i]])
    
    # Subset columns if requested
    if (!is.null(subset_cols)) {
      available_cols <- subset_cols[subset_cols %in% names(chunk_data)]
      if (length(available_cols) > 0) {
        # Convert to data.frame for simple column selection, then back to data.table
        chunk_data <- as.data.frame(chunk_data)
        chunk_data <- chunk_data[, available_cols, drop = FALSE]
        chunk_data <- data.table::as.data.table(chunk_data)
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
    combined_data <- data.table::rbindlist(all_data, fill = TRUE, use.names = TRUE)
    rm(all_data)
    gc(verbose = FALSE)
    
    cat("  ✅ Loaded", nrow(combined_data), "rows for", model_type, "\n")
    return(combined_data)
  } else {
    cat("  ⚠️ No data found for", model_type, "\n")
    return(data.table::data.table())
  }
}

#' Load plot estimates for multiple model types
#' 
#' @param model_types Character vector: model types to load
#' @param subset_cols Character vector: specific columns to load (optional)
#' @param filter_expr Expression: data.table filter expression (optional)
#' @param max_chunks Integer: maximum number of chunks to load per model type (optional)
#' @return data.table with plot estimates from all model types
#' @export
load_all_plot_estimates <- function(model_types = c("env_cycl", "cycl_only"), 
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
    combined_data <- data.table::rbindlist(all_data, fill = TRUE, use.names = TRUE)
    rm(all_data)
    gc(verbose = FALSE)
    
    cat("✅ Total combined rows:", nrow(combined_data), "\n")
    return(combined_data)
  } else {
    cat("⚠️ No data found for any model type\n")
    return(data.table::data.table())
  }
}

#' Get summary statistics for plot estimates without loading all data
#' 
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return list with summary statistics
#' @export
get_plot_estimates_summary <- function(model_type) {
  
  # Find chunk directory
  chunk_dir <- file.path(here::here("data/summary"), paste0("plot_estimates_", model_type))
  
  if (!dir.exists(chunk_dir)) {
    return(list(error = "Chunk directory not found"))
  }
  
  # Find all chunk files
  chunk_files <- list.files(chunk_dir, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(chunk_files) == 0) {
    return(list(error = "No chunk files found"))
  }
  
  # Load first chunk to get structure
  first_chunk <- arrow::read_parquet(chunk_files[1])
  
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

#' Load plot estimates for hindcast initialization
#' 
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return data.table with plot estimates for hindcast initialization
#' @export
load_plot_estimates_for_hindcast <- function(model_type = "env_cycl") {
  cat("Loading plot estimates for hindcast initialization...\n")
  
  # Load only essential columns for hindcast initialization
  essential_cols <- c("siteID", "plotID", "timepoint", "Mean", "date_num", "species")
  
  plot_data <- load_plot_estimates(
    model_type = model_type,
    subset_cols = essential_cols,
    max_chunks = NULL  # Load all chunks
  )
  
  cat("✅ Loaded", nrow(plot_data), "rows for hindcast initialization\n")
  return(plot_data)
}

#' Load plot estimates for scoring metrics
#' 
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return data.table with plot estimates for scoring
#' @export
load_plot_estimates_for_scoring <- function(model_type = "env_cycl") {
  cat("Loading plot estimates for scoring metrics...\n")
  
  # Load columns needed for scoring
  scoring_cols <- c("siteID", "plotID", "timepoint", "Mean", "truth", "species", "model_name")
  
  plot_data <- load_plot_estimates(
    model_type = model_type,
    subset_cols = scoring_cols,
    max_chunks = NULL  # Load all chunks
  )
  
  cat("✅ Loaded", nrow(plot_data), "rows for scoring metrics\n")
  return(plot_data)
}

#' Load plot estimates for phenology analysis
#' 
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return data.table with plot estimates for phenology
#' @export
load_plot_estimates_for_phenology <- function(model_type = "env_cycl") {
  cat("Loading plot estimates for phenology analysis...\n")
  
  # Load columns needed for phenology
  phenology_cols <- c("siteID", "plotID", "timepoint", "Mean", "date_num", "species", "taxon", "model_name", "model_id", "dates", "fcast_type", "time_period", "pretty_group", "rank_only")
  
  plot_data <- load_plot_estimates(
    model_type = model_type,
    subset_cols = phenology_cols,
    max_chunks = NULL  # Load all chunks
  )
  
  cat("✅ Loaded", nrow(plot_data), "rows for phenology analysis\n")
  return(plot_data)
}

#' Load plot estimates for a specific site
#' 
#' @param site_id Character string: site ID to filter by
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return data.table with plot estimates for the specified site
#' @export
load_plot_estimates_for_site <- function(site_id, model_type = "env_cycl") {
  cat("Loading plot estimates for site:", site_id, "\n")
  
  plot_data <- load_plot_estimates(
    model_type = model_type,
    filter_expr = quote(siteID == site_id),
    max_chunks = NULL  # Load all chunks
  )
  
  cat("✅ Loaded", nrow(plot_data), "rows for site", site_id, "\n")
  return(plot_data)
}

#' Load plot estimates for a specific species
#' 
#' @param species_name Character string: species name to filter by
#' @param model_type Character string: "env_cycl" or "cycl_only"
#' @return data.table with plot estimates for the specified species
#' @export
load_plot_estimates_for_species <- function(species_name, model_type = "env_cycl") {
  cat("Loading plot estimates for species:", species_name, "\n")
  
  plot_data <- load_plot_estimates(
    model_type = model_type,
    filter_expr = quote(species == species_name),
    max_chunks = NULL  # Load all chunks
  )
  
  cat("✅ Loaded", nrow(plot_data), "rows for species", species_name, "\n")
  return(plot_data)
}




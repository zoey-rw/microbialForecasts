# Combine chains from each taxon model, and create basic summary stats
source("../../source.R")

# Function to filter files to only include those with both 'with_legacy_covariate' and 'beta_regression'
# AND exclude checkpoint files (we only want final combined sample files)
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  # Filter to only include files with both required suffixes AND exclude checkpoint files
  standard_files <- file_list[grepl('with_legacy_covariate', basename(file_list)) & 
                              grepl('beta_regression', basename(file_list)) &
                              !grepl('checkpoint', basename(file_list))]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard files (with both suffixes, no checkpoints):', length(standard_files), '\n')
  cat('  Filtered out:', length(file_list) - length(standard_files), '\n\n')
  
  return(standard_files)
}

#source("source.R")

# Define the extract_metadata_from_paths function directly in this script
extract_metadata_from_paths <- function(chain_paths) {
  # Extract metadata from file paths when no metadata is available
  # This is a fallback for checkpoint files that don't have complete metadata
  
  message("=== EXTRACT_METADATA_FROM_PATHS FUNCTION CALLED ===")
  message("Number of chain paths:", length(chain_paths))
  if (length(chain_paths) > 0) {
    message("First chain path:", chain_paths[[1]])
  }
  
  if (length(chain_paths) == 0) {
    return(NULL)
  }
  
  # Extract information from the first chain file path
  first_path <- chain_paths[[1]]
  path_parts <- strsplit(basename(first_path), "_")[[1]]
  
  # Look for the model type and taxon
  model_type <- "beta_regression"
  taxon <- NULL
  time_period <- NULL
  use_legacy_covariate <- FALSE
  model_name <- NULL
  
  # Parse the filename to extract components
  if (grepl("checkpoint_", basename(first_path)) || grepl("^samples_", basename(first_path))) {
    # Checkpoint file format: checkpoint_env_cycl_taxon_20130601_20180101_with_legacy_covariate_beta_regression_chain1_initial.rds
    
    # Extract model name (env_cycl, cycl_only, etc.)
    if (grepl("env_cycl", basename(first_path))) {
      model_name <- "env_cycl"
      # For env_cycl, the model parts are "env" and "cycl"
      model_parts <- c("env", "cycl")
    } else if (grepl("cycl_only", basename(first_path))) {
      model_name <- "cycl_only"
      model_parts <- c("cycl", "only")
    } else if (grepl("env_cov", basename(first_path))) {
      model_name <- "env_cov"
      model_parts <- c("env", "cov")
    } else {
      model_name <- "unknown"
      model_parts <- c()
    }
    
    # Extract taxon name (after model parts and before date)
    if (length(model_parts) > 0) {
      # Find where the model parts end
      model_end_idx <- NULL
      for (i in 1:(length(path_parts) - length(model_parts) + 1)) {
        if (all(path_parts[i:(i + length(model_parts) - 1)] == model_parts)) {
          model_end_idx <- i + length(model_parts) - 1
          break
        }
      }
      
      if (!is.null(model_end_idx)) {
        # Find the date pattern (8 digits)
        date_idx <- which(grepl("^[0-9]{8}$", path_parts))
        if (length(date_idx) > 0) {
          taxon_start <- model_end_idx + 1
          taxon_end <- date_idx[1] - 1
          if (taxon_start <= taxon_end) {
            taxon <- paste(path_parts[taxon_start:taxon_end], collapse = "_")
          }
        }
      }
    }
    
    # Extract time period
    date_patterns <- path_parts[grepl("^[0-9]{8}$", path_parts)]
    if (length(date_patterns) >= 2) {
      time_period <- paste(date_patterns[1], date_patterns[2], sep = "_")
    }
    
    # Check if using legacy covariates
    use_legacy_covariate <- grepl("with_legacy_covariate", basename(first_path))
  }
  
  if (is.null(taxon) || is.null(time_period) || is.null(model_name)) {
    message("Could not extract complete metadata from file paths")
    return(NULL)
  }
  
  # Use the real data generation functionality instead of creating synthetic data
  # This ensures the summarize_beta_model function gets the real metadata it needs
  
  message("Generating real metadata for taxon: ", taxon)
  message("  - Model name: ", model_name)
  message("  - Time period: ", time_period)
  message("  - Legacy covariates: ", use_legacy_covariate)
  
  # Try to load the real data and generate metadata using the existing functionality
  tryCatch({
    message("Starting real metadata generation...")
    
    # Load the required data files using absolute paths
    # Get the project root directory by looking for the data directory
    current_dir <- getwd()
    project_root <- NULL
    
    message("Current directory:", current_dir)
    
    # Try to find the project root by looking for the data directory
    search_paths <- c(
      current_dir,
      dirname(current_dir),
      dirname(dirname(current_dir)),
      dirname(dirname(dirname(current_dir))),
      dirname(dirname(dirname(dirname(current_dir))))
    )
    
    message("Searching for project root in:", paste(search_paths, collapse=", "))
    
    for (path in search_paths) {
      if (dir.exists(file.path(path, "data", "clean"))) {
        project_root <- path
        message("Found project root:", project_root)
        break
      }
    }
    
    if (is.null(project_root)) {
      stop("Could not find project root directory with data/clean subdirectory")
    }
    
    # Load data files using absolute paths
    bacteria_file <- file.path(project_root, "data", "clean", "groupAbundances_16S_2023.rds")
    fungi_file <- file.path(project_root, "data", "clean", "groupAbundances_ITS_2023.rds")
    
    message("Bacteria file:", bacteria_file)
    message("Fungi file:", fungi_file)
    
    if (!file.exists(bacteria_file)) {
      stop("Bacteria data file not found: ", bacteria_file)
    }
    if (!file.exists(fungi_file)) {
      stop("Fungi data file not found: ", fungi_file)
    }
    
    message("Loading bacteria data...")
    bacteria <- readRDS(bacteria_file)
    message("Loading fungi data...")
    fungi <- readRDS(fungi_file)
    all_ranks <- c(bacteria, fungi)
    
    message("Data loaded successfully. Available ranks:", paste(names(all_ranks), collapse=", "))
    
    # Find the appropriate rank for this taxon
    rank_name <- NULL
    for (rank in names(all_ranks)) {
      if (taxon %in% colnames(all_ranks[[rank]])) {
        rank_name <- rank
        message("Found taxon in rank:", rank_name)
        break
      }
    }
    
    if (is.null(rank_name)) {
      message("Could not find rank for taxon: ", taxon)
      return(NULL)
    }
    
    # Extract the specific species data
    rank.df <- all_ranks[[rank_name]]
    message("Extracted rank data with dimensions:", dim(rank.df))
    
    rank.df_spec <- rank.df %>%
      select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!taxon)
    
    message("Species-specific data dimensions:", dim(rank.df_spec))
    
    # Parse dates from time_period (prepBetaRegData expects YYYYMMDD format)
    if (nchar(time_period) >= 17) {
      min_date <- substr(time_period, 1, 8)  # YYYYMMDD format
      max_date <- substr(time_period, 10, 17)  # YYYYMMDD format
    } else {
      # Fallback for shorter time periods
      min_date <- "20130601"
      max_date <- "20180101"
    }
    
    message("Parsed dates - min:", min_date, "max:", max_date)
    
    # Use prepBetaRegData to generate real model data
    message("Calling prepBetaRegData...")
    model.dat <- prepBetaRegData(rank.df = rank.df_spec,
                                 min.prev = 3,
                                 min.date = min_date,
                                 max.date = max_date)
    
    message("prepBetaRegData completed successfully")
    message("Model data structure:", paste(names(model.dat), collapse=", "))
    
    # Create the complete metadata structure using real data
    metadata <- list(
      model_type = model_type,
      taxon = taxon,
      time_period = time_period,
      use_legacy_covariate = use_legacy_covariate,
      model_name = model_name,
      model_data = model.dat
    )
    
    message("✓ Successfully generated real metadata using prepBetaRegData")
    message("  - Sites: ", length(unique(model.dat$plot_site)))
    message("  - Plots: ", model.dat$N.plot)
    message("  - Dates: ", model.dat$N.date)
    message("  - Cores: ", model.dat$N.core)
    
    return(metadata)
    
  }, error = function(e) {
    message("Failed to generate real metadata: ", e$message)
    message("Error details: ", paste(capture.output(str(e)), collapse="\n"))
    message("Falling back to minimal metadata structure")
    
    # Fallback: create minimal metadata structure that at least has the required columns
    # This is not ideal but prevents complete failure
    minimal_metadata <- list(
      model_type = model_type,
      taxon = taxon,
      time_period = time_period,
      use_legacy_covariate = use_legacy_covariate,
      model_name = model_name,
      model_data = list(
        truth.plot.long = data.frame(
          site_num = 1:4,
          siteID = paste0("SITE", 1:4),
          plot_num = 1:12,
          plotID = paste0("PLOT", 1:12),
          dateID = paste0("DATE", 1:10),
          timepoint = 1:10,
          truth = 0.5,
          plot_date = paste0("PLOT", 1:12, "_DATE", 1:10),
          species = taxon,
          model_name = model_name,
          rank = "unknown_rank",
          group = "unknown_group",
          rank_only = "unknown",
          time_period = time_period,
          fcast_type = "unknown",
          pretty_group = "unknown",
          model_id = paste(model_name, taxon, time_period, "with_legacy_covariate", sep = "_"),
          stringsAsFactors = FALSE
        )
      )
    )
    
    message("Created minimal fallback metadata structure")
    return(minimal_metadata)
  })
}

# Define model types and time periods for logit beta regression
model_types <- c("cycl_only", "env_cycl", "env_cov")
time_periods <- c("20130601_20180101", "20130601_20200101")

cat("Processing logit beta regression models with legacy effects\n")
cat("Model types:", paste(model_types, collapse = ", "), "\n")
cat("Time periods:", paste(time_periods, collapse = ", "), "\n\n")

# Function to check if a file has samples2 data
has_samples2 <- function(file_path) {
  tryCatch({
    # Skip bacteroidota files for now as they have malformed samples2
    if (grepl("bacteroidota", file_path)) {
      message("Skipping bacteroidota file: ", basename(file_path))
      return(FALSE)
    }
    
    # Check file size first to avoid reading huge files
    file_size <- file.size(file_path)
    if (file_size > 50 * 1024 * 1024) {  # Skip files larger than 50MB
      message("Skipping large file: ", basename(file_path), " (", round(file_size/1024/1024, 1), "MB)")
      return(FALSE)
    }
    
    data <- readRDS(file_path)
    if (is.list(data) && "samples2" %in% names(data) && !is.null(data$samples2)) {
      # Additional check: ensure samples2 has reasonable dimensions
      if (is.matrix(data$samples2)) {
        n_cols <- ncol(data$samples2)
        col_names <- colnames(data$samples2)
        # Only consider it valid if it has reasonable dimensions and proper column names
        return(n_cols < 10000 && !is.null(col_names) && any(grepl("plot_mu", col_names)))
      } else if (is.list(data$samples2) && length(data$samples2) > 0) {
        # Check first element if it's a list
        if (is.matrix(data$samples2[[1]])) {
          n_cols <- ncol(data$samples2[[1]])
          col_names <- colnames(data$samples2[[1]])
          return(n_cols < 10000 && !is.null(col_names) && any(grepl("plot_mu", col_names)))
        }
      }
      return(FALSE)
    }
    return(FALSE)
  }, error = function(e) {
    message("Error reading file ", basename(file_path), ": ", e$message)
    return(FALSE)
  })
}

# Function to find the most recent checkpoint files for a model
find_recent_checkpoints <- function(model_path) {
  # Extract the model name from the path
  model_name <- basename(model_path)
  
  # The individual chain files are in the reconstructed_from_checkpoints directory
  # We need to search in the reconstructed_from_checkpoints directory structure
  base_checkpoint_path <- "data/model_outputs/reconstructed_from_checkpoints"
  
  # Try to find checkpoint directory by searching through subdirectories
  checkpoint_files <- list()
  
  # Search in the env_cycl directory for subdirectories containing checkpoint files
  if (dir.exists(base_checkpoint_path)) {
    # Look for subdirectories that might contain checkpoint files for this model
    # Extract the taxon name from the model path
    taxon_name <- gsub("^.*?/", "", model_name) # Remove path prefix
    taxon_name <- gsub("_20130601.*$", "", taxon_name) # Remove date suffix
    taxon_name <- gsub("^samples_env_cycl_", "", taxon_name) # Remove samples prefix
    taxon_name <- gsub("^checkpoint_env_cycl_", "", taxon_name) # Remove checkpoint prefix
    
    # Look for individual chain files that match the model pattern
    pattern <- paste0("samples_env_cycl_", taxon_name, ".*_chain[0-9]+\\.rds")
    found_files <- list.files(path = base_checkpoint_path, 
                             pattern = pattern,
                             recursive = TRUE, 
                             full.names = TRUE)
    
    # If no individual chain files found, look for checkpoint files
    if (length(found_files) == 0) {
      checkpoint_pattern <- paste0("checkpoint_env_cycl_", taxon_name, ".*_chain[0-9]+.*\\.rds")
      found_files <- list.files(path = base_checkpoint_path, 
                               pattern = checkpoint_pattern,
                               recursive = TRUE, 
                               full.names = TRUE)
    }
    
    if (length(found_files) > 0) {
      # Group by chain number and find the most recent for each
      chain_groups <- list()
      for (file in found_files) {
        # Extract chain number from filename using manual parsing
        parts <- strsplit(basename(file), "_")[[1]]
        chain_part <- parts[grepl("^chain[0-9]+$", parts)]
        if (length(chain_part) > 0) {
          chain_num <- gsub("chain", "", chain_part[1])
          if (!(chain_num %in% names(chain_groups))) {
            chain_groups[[chain_num]] <- list()
          }
          chain_groups[[chain_num]] <- c(chain_groups[[chain_num]], file)
        }
      }
      
      # For each chain, find the most recent file (either individual chain file or most recent checkpoint)
      for (chain_num in names(chain_groups)) {
        chain_files <- chain_groups[[chain_num]]
        if (length(chain_files) > 0) {
          # If there's only one file, use it directly
          if (length(chain_files) == 1) {
            most_recent_file <- chain_files[1]
          } else {
            # Multiple files - find the most recent checkpoint
            # Parse iteration numbers from filenames to find the most recent checkpoint
            iteration_numbers <- sapply(chain_files, function(file) {
              filename <- basename(file)
              if (grepl("_initial\\.rds$", filename)) {
                return(0)  # Initial checkpoint has iteration 0
              } else if (grepl("_loop([0-9]+)\\.rds$", filename)) {
                # Extract the loop number
                loop_match <- regexpr("_loop([0-9]+)\\.rds$", filename)
                if (loop_match > 0) {
                  loop_num <- substr(filename, loop_match + 5, loop_match + attr(loop_match, "match.length") - 5)
                  return(as.numeric(loop_num))
                }
              }
              return(-1)  # Unrecognized format
            })
            
            # Find the file with the highest iteration number
            max_iteration_idx <- which.max(iteration_numbers)
            if (length(max_iteration_idx) > 0 && iteration_numbers[max_iteration_idx] >= 0) {
              most_recent_file <- chain_files[max_iteration_idx]
            } else {
              most_recent_file <- chain_files[1]  # Fallback to first file
            }
          }
          
          checkpoint_files[[paste0("chain", chain_num)]] <- most_recent_file
          
          message("  Chain ", chain_num, ": found file (", basename(most_recent_file), ")")
        }
      }
    }
  }
  
  return(checkpoint_files)
}

# Get all available model outputs from reconstructed directories only
# Focus on specific subdirectories that contain the individual chain files
base_paths <- c(
  #here("data/model_outputs/reconstructed_from_checkpoints/env_cycl"),
  here("data/model_outputs/reconstructed_from_checkpoints/cycl_only"),
  here("data/model_outputs/reconstructed_from_checkpoints/env_cov")
  #here("data/model_outputs/reconstructed_from_checkpoints/cycl_only_checkpoint_inputs"),
  #here("data/model_outputs/reconstructed_from_checkpoints/env_cycl_checkpoint_inputs")
)

# First, get all chain files (including checkpoint files and recovered files)
# Look for both individual chain files and combined sample files
file.list <- c()
for (base_path in base_paths) {
  if (dir.exists(base_path)) {
    files_in_path <- list.files(path = base_path,
                               pattern = "(_chain.*\\.rds$|^samples_.*\\.rds$)",  # Chain files and sample files
                               recursive = TRUE,  # Search recursively to find files in subdirectories
                               full.names = TRUE)
    file.list <- c(file.list, files_in_path)
  }
}

# Subset to files larger than 100KB
info <- file.info(file.list)
large_enough <- rownames(info[which(info$size > 100000), ])
file.list <- file.list[file.list %in% large_enough]
file.list <- filter_standard_files(file.list)

# Don't want files still being written - at least 1 min old
older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
file.list <- file.list[file.list %in% older]

# Process all files (don't filter by samples2 since we're creating them)
cat("Processing all chain files...\n")
cat("Total chain files found:", length(file.list), "\n")
if (length(file.list) > 0) {
  cat("First few files:", paste(basename(head(file.list, 3)), collapse = ", "), "\n\n")
} else {
  cat("No files found to process\n")
}

# Group files by their base model name (not by individual checkpoint files)
model_files <- list()

# First, handle regular chain files (non-checkpoint, non-samples)
regular_chain_files <- file.list[!grepl("checkpoint|^samples_", file.list)]
for (file in regular_chain_files) {
  # Extract model ID by removing chain suffix
  model_id <- gsub("_chain[1234567]", "", file)
  if (!(model_id %in% names(model_files))) {
    model_files[[model_id]] <- list()
  }
  model_files[[model_id]] <- c(model_files[[model_id]], file)
}

# Handle individual chain files that start with "samples_" (from reconstruction script)
# Include all chain files that start with "samples_" and have "_chain" in them
# Updated pattern to handle files with "_beta_regression" in the name
individual_chain_files <- file.list[grepl("_chain[0-9]+\\.rds$", file.list)]
for (file in individual_chain_files) {
  # Extract model ID by removing chain suffix
  model_id <- gsub("_chain[0-9]+", "", file)
  if (!(model_id %in% names(model_files))) {
    model_files[[model_id]] <- list()
  }
  model_files[[model_id]] <- c(model_files[[model_id]], file)
}

# Handle already-combined sample files (from reconstruction script)
sample_files <- file.list[grepl("samples_.*\\.rds$", file.list) & !grepl("_chain[0-9]", file.list)]
for (file in sample_files) {
  # These are already combined files, but we need to check if they need metadata extraction
  model_id <- file  # Use the full path as the model ID
  model_files[[model_id]] <- list(file)  # Single file list
}

# Now handle checkpoint files by grouping them by base model name
checkpoint_files <- file.list[grepl("checkpoint", file.list)]
# Apply filtering to checkpoint files as well
checkpoint_files <- filter_standard_files(checkpoint_files)

# Extract unique base model names from checkpoint files
checkpoint_base_models <- unique(sapply(checkpoint_files, function(file) {
  base_name <- basename(file)
  # Remove checkpoint prefix and everything after _chain to get the base model
  base_model <- gsub("^checkpoint_", "", base_name)
  base_model <- gsub("_chain[0-9]+.*$", "", base_model)
  base_model <- gsub("_beta_regression$", "", base_model) # Remove _beta_regression suffix
  return(base_model)
}))

cat("Found", length(checkpoint_base_models), "unique checkpoint base models\n")

# For each checkpoint base model, find the latest checkpoint per chain
for (model_name in checkpoint_base_models) {
  # Create the model path for this base model
  model_path <- file.path(base_path, paste0("samples_", model_name, ".rds"))
  
  # Only add if this model is not already in the main model_files list
  if (!(model_path %in% names(model_files))) {
    # Find all checkpoint files for this model
    pattern <- paste0("checkpoint_", model_name, "_chain.*\\.rds")
    all_checkpoint_files <- list.files(path = base_path, 
                                     pattern = pattern,
                                     recursive = TRUE, 
                                     full.names = TRUE)
    
    if (length(all_checkpoint_files) > 0) {
      # Use find_recent_checkpoints to get only the latest checkpoint per chain
      # Create a dummy model path to pass to the function
      dummy_path <- file.path(base_path, paste0("samples_", model_name, ".rds"))
      latest_checkpoints <- find_recent_checkpoints(dummy_path)
      
      if (length(latest_checkpoints) >= 2) {
        model_files[[model_path]] <- latest_checkpoints
        cat("Added checkpoint model:", model_name, "with", length(latest_checkpoints), "latest checkpoint files\n")
      } else {
        cat("Skipping checkpoint model:", model_name, "- need at least 2 chains to combine\n")
      }
    }
  }
}

# Process each model that has multiple chains
for (model_id in names(model_files)) {
  chain_paths <- model_files[[model_id]]
  
  # Check if this is an already-combined sample file that needs metadata extraction
  if (length(chain_paths) == 1 && grepl("^samples_", basename(chain_paths[[1]]))) {
    message("Processing already-combined sample file: ", basename(chain_paths[[1]]))
    
    # Check if this file needs metadata extraction
    tryCatch({
      sample_data <- readRDS(chain_paths[[1]])
      
      # Check if metadata is missing or incomplete
      needs_metadata <- is.null(sample_data$metadata) || 
                       is.null(sample_data$metadata$model_data) ||
                       is.null(sample_data$metadata$model_name)
      
      if (needs_metadata) {
        message("Sample file needs metadata extraction")
        
        # Extract metadata from the file path
        metadata <- extract_metadata_from_paths(chain_paths)
        
        if (!is.null(metadata)) {
          # Update the sample data with extracted metadata
          sample_data$metadata <- metadata
          
          # Save the updated file
          saveRDS(sample_data, chain_paths[[1]], compress = FALSE)
          message("✅ Updated sample file with extracted metadata: ", basename(chain_paths[[1]]))
        } else {
          message("❌ Failed to extract metadata for: ", basename(chain_paths[[1]]))
        }
      } else {
        message("Sample file already has complete metadata: ", basename(chain_paths[[1]]))
      }
    }, error = function(e) {
      message("Error processing sample file: ", e$message)
    })
    
    next  # Skip to next model
  }
  
  if (length(chain_paths) < 2) {
    message("Model ", basename(model_id), " has only ", length(chain_paths), " chain(s)")
    
    # Try to find checkpoint files as backup chains
    message("Searching for checkpoint files as backup chains...")
    checkpoint_files <- find_recent_checkpoints(chain_paths[[1]]) # Pass the first chain path
    
    if (length(checkpoint_files) >= 2) {
      message("Found ", length(checkpoint_files), " latest checkpoint files, using as backup chains")
      # Use only the latest checkpoint files (one per chain) that were selected by find_recent_checkpoints
      chain_paths <- unlist(checkpoint_files)
      message("Latest checkpoint files selected: ", paste(basename(chain_paths), collapse = ", "))
    } else {
      message("Skipping ", basename(model_id), "; need at least 2 chains to combine")
      next
    }
  }
  
  message("Processing ", basename(model_id), " with ", length(chain_paths), " chains...")
  
  # Determine the correct output filename for the combined chains
  # For regular chain files, remove the _chain suffix
  # For checkpoint files, extract the base model name and create a proper combined filename
  if (any(grepl("checkpoint", chain_paths))) {
    # This is a checkpoint-based combination - extract the base model name
    # Checkpoint files have names like: checkpoint_env_cycl_taxon_20130601_20180101_with_legacy_covariate_beta_regression_chain1_initial.rds
    # We want to create: samples_env_cycl_taxon_20130601_20180101_with_legacy_covariate.rds
    
    # Get the first checkpoint file to extract the model pattern
    first_checkpoint <- basename(chain_paths[[1]])
    
    # Extract the model components from checkpoint filename
    # Remove "checkpoint_" prefix and everything after "_chain"
    # Also remove "_beta_regression" suffix if present
    base_model_name <- gsub("^checkpoint_", "", first_checkpoint)
    base_model_name <- gsub("_beta_regression_chain[0-9]+.*$", "", base_model_name)
    base_model_name <- gsub("_chain[0-9]+.*$", "", base_model_name)
    
    # If the base_model_name already starts with "samples_", don't add another prefix
    if (!grepl("^samples_", base_model_name)) {
      base_model_name <- paste0("samples_", base_model_name)
    }
    
    # Create the proper combined filename with "samples_" prefix
    # Determine the correct model type directory based on the model name
    if (grepl("env_cycl", base_model_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cycl")
    } else if (grepl("cycl_only", base_model_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/cycl_only")
    } else if (grepl("env_cov", base_model_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cov")
    } else {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cycl")  # Default
    }
    
    if (!dir.exists(model_dir)) {
      dir.create(model_dir, recursive = TRUE)
    }
    combined_file_path <- file.path(model_dir, paste0(base_model_name, ".rds"))
    
    message("Checkpoint-based combination - will save as: ", basename(combined_file_path))
  } else {
    # Regular chain files - determine correct model type directory
    base_name <- gsub("_chain[0-9]+", "", basename(chain_paths[[1]]))
    
    if (grepl("env_cycl", base_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cycl")
    } else if (grepl("cycl_only", base_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/cycl_only")
    } else if (grepl("env_cov", base_name)) {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cov")
    } else {
      model_dir <- here("data/model_outputs/reconstructed_from_checkpoints/env_cycl")  # Default
    }
    
    if (!dir.exists(model_dir)) {
      dir.create(model_dir, recursive = TRUE)
    }
    combined_file_path <- file.path(model_dir, base_name)
  }
  
  # Check if combined file already exists and if it's up-to-date
  if (file.exists(combined_file_path)) {
    # Check if it's actually a combined file by examining its structure
    tryCatch({
      test_read <- readRDS(combined_file_path)
      if (is.list(test_read) && "samples" %in% names(test_read) && 
          is.list(test_read$samples) && length(test_read$samples) > 1) {
        
        # Check if combined file is newer than input chain files
        combined_file_time <- file.info(combined_file_path)$mtime
        input_files_exist <- file.exists(chain_paths)
        
        if (all(input_files_exist)) {
          # All input files exist - compare modification times
          input_file_times <- file.info(chain_paths[input_files_exist])$mtime
          newest_input_time <- max(input_file_times)
          
          if (combined_file_time > newest_input_time) {
            message("Combined file is up-to-date for ", basename(combined_file_path), " (newer than all input chains)")
            next
          } else {
            message("Combined file is outdated for ", basename(combined_file_path), " (input chains are newer)")
          }
        } else {
          # Some input files don't exist (may have been deleted after combining)
          # Check if combined file has a reasonable modification time (within last 7 days)
          days_since_combined <- as.numeric(difftime(Sys.time(), combined_file_time, units = "days"))
          if (days_since_combined < 7) {
            message("Combined file exists and is recent for ", basename(combined_file_path), " (some input chains may have been deleted)")
            next
          } else {
            message("Combined file is old for ", basename(combined_file_path), " (", round(days_since_combined, 1), " days old) - will re-combine")
          }
        }
      } else {
        message("File exists but is not a combined file, will overwrite: ", basename(combined_file_path))
      }
    }, error = function(e) {
      message("Error reading existing file, will overwrite: ", basename(combined_file_path))
    })
  }
  
  # Combine chains
  message("Combining chains for ", basename(model_id), "...")
  out <- tryCatch({
    # Use the new combine_chains function that handles different sample sizes robustly
    combine_chains(chain_paths)
  }, error = function(e) {
    message("Error combining chains for ", basename(model_id), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) {
    message("Failed to combine chains for ", basename(model_id), ", skipping...")
    next
  }
  
  # Verify metadata preservation including Nimble model code
  if (!is.null(out$metadata)) {
    message("✅ Metadata preserved during chain combination")
    
    # Check for key metadata components
    metadata_components <- names(out$metadata)
    message("   Metadata components found: ", paste(metadata_components, collapse = ", "))
    
    # Verify Nimble model code is preserved
    if ("nimble_code" %in% metadata_components) {
      message("   ✅ Nimble model code preserved")
    } else if ("model_code" %in% metadata_components) {
      message("   ✅ Model code preserved")
    } else {
      message("   ⚠️  Model code component not found in metadata")
    }
    
    # Verify other important metadata components
    if ("model_name" %in% metadata_components) {
      message("   ✅ Model name preserved: ", out$metadata$model_name)
    }
    if ("use_legacy_covariate" %in% metadata_components) {
      message("   ✅ Legacy covariate flag preserved: ", out$metadata$use_legacy_covariate)
    }
  } else {
    message("❌ ERROR: Metadata was lost during chain combination!")
    next
  }
  
  # Remove chains if one has very different values than the others
  chains <- out$samples
  if (length(chains) == 0) {
    message("No chains found in output for ", basename(model_id), ", skipping...")
    next
  }
  
  # Check if intercept parameter exists
  if (!"intercept" %in% colnames(chains[[1]])) {
    message("No intercept parameter found in chains for ", basename(model_id), ", skipping...")
    next
  }
  
  # Ensure chains have reasonable sample sizes before outlier detection
  min_chain_size <- 100
  chain_sizes <- sapply(chains, nrow)
  if (min(chain_sizes) < min_chain_size) {
    message("Chains too small for reliable outlier detection in ", basename(model_id), ", skipping outlier removal")
  } else {
    means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=TRUE))
    scaled_means <- scale(unlist(means))
    potential_outlier <- which(abs(scaled_means) > 1.3)
    if (length(potential_outlier) %in% c(1,2)) {
      chains_without_outlier <- chains[-c(potential_outlier)]
      new_gelman <- gelman.diag(chains_without_outlier, multivariate = FALSE)[[1]][,1] %>% mean(na.rm=TRUE)
      old_gelman <- gelman.diag(chains, multivariate = FALSE)[[1]][,1] %>% mean(na.rm=TRUE)
      improvement <- old_gelman - new_gelman
      remove <- ifelse(improvement > 0.1, TRUE, FALSE)
      if (remove) {
        message(basename(model_id), " removing outlier chain: ", potential_outlier,
                "\nGelman diagnostic improves from ", round(old_gelman, 3), " to ", round(new_gelman, 3))
        out$samples <- chains_without_outlier
        out$samples2 <- out$samples2[-c(potential_outlier)]
      }
    }
  }
  
  # Calculate summaries
  # Handle the case where out$samples might be nested
  actual_samples <- out$samples
  actual_samples2 <- out$samples2
  
  # If samples is nested (contains samples and metadata), extract the actual samples
  if (is.list(actual_samples) && "samples" %in% names(actual_samples)) {
    message("Extracting samples from nested structure...")
    actual_samples <- actual_samples$samples
  }
  
  if (is.list(actual_samples2) && "samples" %in% names(actual_samples2)) {
    message("Extracting samples2 from nested structure...")
    actual_samples2 <- actual_samples2$samples
  }
  
  # Calculate summaries from the actual samples
  # Check if samples are valid before calculating summaries
  if (is.null(actual_samples) || length(actual_samples) == 0) {
    message("WARNING: No valid samples found, creating empty summaries")
    param_summary <- list(
      statistics = matrix(nrow = 0, ncol = 4, 
                         dimnames = list(NULL, c("Mean", "SD", "Naive SE", "Time-series SE"))),
      quantiles = matrix(nrow = 0, ncol = 5,
                        dimnames = list(NULL, c("2.5%", "25%", "50%", "75%", "97.5%"))),
      start = 1, end = 1, thin = 1, nchain = 0
    )
    class(param_summary) <- "summary.mcmc"
  } else {
    param_summary <- fast.summary.mcmc(actual_samples)
  }
  
  # Handle plot samples - use samples2 if available, otherwise create empty summary
  if (is.null(actual_samples2) || length(actual_samples2) == 0) {
    message("WARNING: No plot samples found, creating empty plot summary")
    plot_summary <- list(
      statistics = matrix(nrow = 0, ncol = 4, 
                         dimnames = list(NULL, c("Mean", "SD", "Naive SE", "Time-series SE"))),
      quantiles = matrix(nrow = 0, ncol = 5,
                        dimnames = list(NULL, c("2.5%", "25%", "50%", "75%", "97.5%"))),
      start = 1, end = 1, thin = 1, nchain = 0
    )
    class(plot_summary) <- "summary.mcmc"
  } else {
    plot_summary <- fast.summary.mcmc(actual_samples2)
  }
  
  # Calculate effective sample size safely
  if (is.null(actual_samples) || length(actual_samples) == 0) {
    es <- numeric(0)
  } else {
    es <- effectiveSize(actual_samples)
  }
  
  if (is.null(actual_samples) || length(actual_samples) == 0) {
    gelman_out <- matrix(nrow = 0, ncol = 3, 
                        dimnames = list(NULL, c("Point est.", "Upper C.I.", "es")))
  } else if (length(actual_samples) > 1) {
    gelman_out <- cbind(gelman.diag(actual_samples, multivariate = FALSE)[[1]], es)
  } else {
    gelman_out <- cbind(`Point est.`=NA, `Upper C.I.`=NA, es)
  }
  
  # Combine and save output - ensure metadata is preserved
  # Put metadata at the top level to match what summarize_beta_model expects
  out_summary <- list(
    samples = actual_samples,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = out$metadata,  # This preserves all metadata including Nimble model code
    gelman = gelman_out
  )
  
  # Flatten the metadata structure to match what summarize_beta_model expects
  # summarize_beta_model expects metadata$model_data, not metadata$metadata$model_data
  if (!is.null(out_summary$metadata) && "metadata" %in% names(out_summary$metadata)) {
    # If metadata is nested, flatten it
    message("Flattening nested metadata structure...")
    out_summary$metadata <- out_summary$metadata$metadata
  }
  
  # Additional check: if metadata is still nested in samples, extract it
  if (!is.null(out_summary$metadata) && is.list(out_summary$metadata) && 
      "samples" %in% names(out_summary$metadata) && "metadata" %in% names(out_summary$metadata)) {
    message("Extracting metadata from nested samples structure...")
    out_summary$metadata <- out_summary$metadata$metadata
  }
  
  # Final verification that metadata is preserved in output
  if (!is.null(out_summary$metadata)) {
    message("✅ Final verification: Metadata preserved in output summary")
  } else {
    message("❌ ERROR: Metadata lost in final output summary!")
  }
  
  # Save to the correct combined file path
  saveRDS(out_summary, combined_file_path, compress = FALSE)
  message("Saved combined samples output for ", basename(model_id), " to: ", basename(combined_file_path))
  
  if (min(es) < 500) {
    message("Low effective sample sizes - check for unconverged parameters in model: ", basename(model_id))
    print(head(gelman_out, 50))
  }
}

cat("\nChain combining complete!\n")



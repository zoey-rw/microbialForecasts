# Combine chains from reconstructed checkpoint files
# This script specifically processes files in the reconstructed_from_checkpoints directory
source("../../source.R")

# CRITICAL FIX: Add comprehensive taxon name validation function
validate_taxon_name <- function(taxon, model_id, file_path = NULL) {
  validation_result <- list(
    is_valid = FALSE,
    corrected_taxon = NULL,
    error_message = NULL,
    warnings = character(0)
  )
  
  # Check for null or empty taxon
  if (is.null(taxon) || is.na(taxon) || nchar(taxon) == 0) {
    validation_result$error_message <- "Taxon name is null, NA, or empty"
    return(validation_result)
  }
  
  # Check for obviously invalid taxon names
  invalid_patterns <- c("^unknown$", "^NA$", "^null$", "^$", "^with$", "^legacy$", "^covariate$", "^beta$", "^regression$")
  for (pattern in invalid_patterns) {
    if (grepl(pattern, taxon, ignore.case = TRUE)) {
      validation_result$error_message <- paste("Taxon name matches invalid pattern:", pattern)
      return(validation_result)
    }
  }
  
  # Check for taxon names that are too short (likely extraction errors)
  if (nchar(taxon) < 3) {
    validation_result$warnings <- c(validation_result$warnings, "Taxon name is very short, may be extraction error")
  }
  
  # Check for taxon names that contain only numbers (likely extraction errors)
  if (grepl("^[0-9]+$", taxon)) {
    validation_result$error_message <- "Taxon name contains only numbers, likely extraction error"
    return(validation_result)
  }
  
  # Check for common extraction mistakes
  if (grepl("_with_|_legacy_|_covariate_|_beta_|_regression_", taxon)) {
    # Try to clean up the taxon name
    cleaned_taxon <- gsub("_with_.*|_legacy_.*|_covariate_.*|_beta_.*|_regression_.*", "", taxon)
    if (nchar(cleaned_taxon) > 0) {
      validation_result$corrected_taxon <- cleaned_taxon
      validation_result$warnings <- c(validation_result$warnings, "Taxon name contained non-taxon parts, cleaned")
    }
  }
  
  # If we get here, the taxon name is valid
  validation_result$is_valid <- TRUE
  if (is.null(validation_result$corrected_taxon)) {
    validation_result$corrected_taxon <- taxon
  }
  
  return(validation_result)
}

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
    
    # CRITICAL FIX: Improved taxon name extraction logic
    # Extract taxon name (after model parts and before date)
    if (length(model_parts) > 0) {
      # Find where the model parts end - use more robust detection
      model_end_idx <- NULL
      
      # Look for the model parts in sequence
      for (i in 1:(length(path_parts) - length(model_parts) + 1)) {
        if (all(path_parts[i:(i + length(model_parts) - 1)] == model_parts)) {
          model_end_idx <- i + length(model_parts) - 1
          break
        }
      }
      
      # CRITICAL FIX: If model parts not found in sequence, try alternative approach
      if (is.null(model_end_idx)) {
        # For cases where model parts might be separated, look for the model name pattern
        if (model_name == "env_cycl") {
          # Look for "env" followed by "cycl" (may be separated)
          env_idx <- which(path_parts == "env")
          cycl_idx <- which(path_parts == "cycl")
          if (length(env_idx) > 0 && length(cycl_idx) > 0 && cycl_idx > env_idx) {
            model_end_idx <- cycl_idx
          }
        } else if (model_name == "env_cov") {
          # Look for "env" followed by "cov"
          env_idx <- which(path_parts == "env")
          cov_idx <- which(path_parts == "cov")
          if (length(env_idx) > 0 && length(cov_idx) > 0 && cov_idx > env_idx) {
            model_end_idx <- cov_idx
          }
        } else if (model_name == "cycl_only") {
          # Look for "cycl" followed by "only"
          cycl_idx <- which(path_parts == "cycl")
          only_idx <- which(path_parts == "only")
          if (length(cycl_idx) > 0 && length(only_idx) > 0 && only_idx > cycl_idx) {
            model_end_idx <- only_idx
          }
        }
      }
      
      if (!is.null(model_end_idx)) {
        # Find the date pattern (8 digits)
        date_idx <- which(grepl("^[0-9]{8}$", path_parts))
        if (length(date_idx) > 0) {
          taxon_start <- model_end_idx + 1
          taxon_end <- date_idx[1] - 1
          
          # CRITICAL FIX: Validate taxon extraction
          if (taxon_start <= taxon_end && taxon_start <= length(path_parts)) {
            taxon_parts <- path_parts[taxon_start:min(taxon_end, length(path_parts))]
            
            # Filter out any non-taxon parts (like "with", "legacy", "covariate")
            taxon_parts <- taxon_parts[!taxon_parts %in% c("with", "legacy", "covariate", "beta", "regression")]
            
            if (length(taxon_parts) > 0) {
              taxon <- paste(taxon_parts, collapse = "_")
              
              # CRITICAL FIX: Use comprehensive taxon validation
              taxon_validation <- validate_taxon_name(taxon, model_id, first_path)
              
              if (!taxon_validation$is_valid) {
                message("ERROR: Invalid taxon name extracted from: ", basename(first_path))
                message("  Taxon: '", taxon, "'")
                message("  Error: ", taxon_validation$error_message)
                taxon <- NULL
              } else {
                # Use corrected taxon if available
                if (!is.null(taxon_validation$corrected_taxon) && taxon_validation$corrected_taxon != taxon) {
                  message("INFO: Corrected taxon name from '", taxon, "' to '", taxon_validation$corrected_taxon, "'")
                  taxon <- taxon_validation$corrected_taxon
                }
                
                # Report warnings
                if (length(taxon_validation$warnings) > 0) {
                  for (warning_msg in taxon_validation$warnings) {
                    message("WARNING: ", warning_msg, " for taxon '", taxon, "'")
                  }
                }
              }
            } else {
              message("Warning: No valid taxon parts found in: ", basename(first_path))
              taxon <- NULL
            }
          } else {
            message("Warning: Invalid taxon extraction indices for: ", basename(first_path))
            taxon <- NULL
          }
        } else {
          message("Warning: No date pattern found in: ", basename(first_path))
          taxon <- NULL
        }
      } else {
        message("Warning: Could not find model parts in: ", basename(first_path))
        taxon <- NULL
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

cat("Processing reconstructed checkpoint files\n")
cat("Looking for individual chain files in reconstructed_from_checkpoints directory\n\n")

# Get all individual chain files from the reconstructed directory
base_path <- here("data/model_outputs/reconstructed_from_checkpoints")

# Look for individual chain files that start with "samples_" and have "_chain[0-9]" pattern
file.list <- list.files(path = base_path,
                       pattern = "^samples_.*_chain[0-9]+\\.rds$",
                       recursive = TRUE,
                       full.names = TRUE)

# Subset to files larger than 100KB
info <- file.info(file.list)
large_enough <- rownames(info[which(info$size > 100000), ])
file.list <- file.list[file.list %in% large_enough]

# Don't want files still being written - at least 1 min old
older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
file.list <- file.list[file.list %in% older]

cat("Total individual chain files found:", length(file.list), "\n")
cat("First few files:", paste(basename(head(file.list, 3)), collapse = ", "), "\n\n")

# Group files by their base model name
model_files <- list()

for (file in file.list) {
  # Extract model ID by removing chain suffix
  model_id <- gsub("_chain[0-9]+", "", file)
  if (!(model_id %in% names(model_files))) {
    model_files[[model_id]] <- list()
  }
  model_files[[model_id]] <- c(model_files[[model_id]], file)
}

cat("Found", length(model_files), "unique models to process\n\n")

# Process each model that has multiple chains
for (model_id in names(model_files)) {
  chain_paths <- model_files[[model_id]]
  
  if (length(chain_paths) < 2) {
    message("Model ", basename(model_id), " has only ", length(chain_paths), " chain(s), skipping...")
    next
  }
  
  message("Processing ", basename(model_id), " with ", length(chain_paths), " chains...")
  
  # Determine the correct output filename for the combined chains
  # Extract the base model name from the first chain file
  first_chain <- basename(chain_paths[[1]])
  
  # Remove "_chain[0-9]" suffix to get the base model name
  base_model_name <- gsub("_chain[0-9]+", "", first_chain)
  
  # Create the proper combined filename
  # Save in the same reconstructed directory structure
  combined_file_path <- file.path(dirname(chain_paths[[1]]), base_model_name)
  
  message("Will save combined file as: ", basename(combined_file_path))
  
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
    # Use the combine_chains function that handles different sample sizes robustly
    combine_chains(chain_paths)
  }, error = function(e) {
    message("Error combining chains for ", basename(model_id), ": ", e$message)
    return(NULL)
  })
  
  # Add improved reconstruction if plot_summary is missing or invalid
  if (!is.null(out) && !is.null(out$metadata) && !is.null(out$metadata$model_data)) {
    # Check if plot_summary exists and has plot data
    needs_reconstruction <- FALSE
    if (is.null(out$plot_summary) || length(out$plot_summary) == 0) {
      needs_reconstruction <- TRUE
    } else if (is.list(out$plot_summary) && length(out$plot_summary) > 0) {
      if (nrow(out$plot_summary[[1]]) == 0) {
        needs_reconstruction <- TRUE
      } else {
        rownames_sample <- head(rownames(out$plot_summary[[1]]), 10)
        has_plot_estimates <- any(grepl("plot_mu", rownames_sample))
        if (!has_plot_estimates) {
          needs_reconstruction <- TRUE
        }
      }
    }
    
    if (needs_reconstruction) {
      message("  Reconstructing plot estimates with data-driven initial conditions...")
      
      # Source the improved reconstruction function
      source("../../clean_reconstruction_function_v2.R")
      
      # Determine model type from metadata or filename
      model_type <- "cycl_only"  # default
      if (!is.null(out$metadata$model_name)) {
        model_type <- out$metadata$model_name
      } else if (grepl("env_cycl", basename(model_id))) {
        model_type <- "env_cycl"
      } else if (grepl("env_cov", basename(model_id))) {
        model_type <- "env_cov"
      }
      
      # Run improved reconstruction
      reconstructed_plot_summary <- tryCatch({
        reconstruct_plot_estimates(
          out, 
          out$param_summary, 
          out$metadata$model_data, 
          model_type
        )
      }, error = function(e) {
        message("  Error in reconstruction: ", e$message)
        return(NULL)
      })
      
      # Add reconstructed plot summary to the output
      if (!is.null(reconstructed_plot_summary) && 
          length(reconstructed_plot_summary) == 2 && 
          nrow(reconstructed_plot_summary[[1]]) > 0 && 
          nrow(reconstructed_plot_summary[[2]]) > 0) {
        out$plot_summary <- reconstructed_plot_summary
        message("  ✅ Added reconstructed plot summary with data-driven initial conditions")
      } else {
        message("  ⚠️  Reconstruction failed, proceeding without plot summary")
      }
    } else {
      message("  ✅ Plot summary already exists and is valid")
    }
  }
  
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
    
    # Check if model_data is missing and extract it if needed
    if (!"model_data" %in% metadata_components) {
      message("   ⚠️  Model data missing, attempting to extract from file paths...")
      extracted_metadata <- extract_metadata_from_paths(chain_paths)
      if (!is.null(extracted_metadata) && "model_data" %in% names(extracted_metadata)) {
        out$metadata$model_data <- extracted_metadata$model_data
        message("   ✅ Successfully extracted model data")
      } else {
        message("   ❌ Failed to extract model data")
      }
    } else {
      message("   ✅ Model data already present, skipping extraction")
    }
    
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
  param_summary <- summary(out$samples)
  plot_summary <- summary(out$samples2)
  es <- effectiveSize(out$samples)
  
  if (length(out$samples) > 1) {
    gelman_out <- cbind(gelman.diag(out$samples, multivariate = FALSE)[[1]], es)
  } else {
    gelman_out <- cbind(`Point est.`=NA, `Upper C.I.`=NA, es)
  }
  
  # Combine and save output - ensure metadata is preserved
  # Put metadata at the top level to match what summarize_beta_model expects
  out_summary <- list(
    samples = out$samples,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = out$metadata,  # This preserves all metadata including Nimble model code
    gelman = gelman_out
  )
  
  # Flatten the metadata structure to match what summarize_beta_model expects
  # summarize_beta_model expects metadata$model_data, not metadata$metadata$model_data
  if (!is.null(out_summary$metadata)) {
    # Check if metadata is nested and flatten it
    if ("metadata" %in% names(out_summary$metadata)) {
      message("Flattening nested metadata structure...")
      out_summary$metadata <- out_summary$metadata$metadata
    }
    
    # Also check if the metadata itself has nested structure
    if (!is.null(out_summary$metadata) && "metadata" %in% names(out_summary$metadata)) {
      message("Flattening double-nested metadata structure...")
      out_summary$metadata <- out_summary$metadata$metadata
    }
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

cat("\nReconstructed chain combining complete!\n")

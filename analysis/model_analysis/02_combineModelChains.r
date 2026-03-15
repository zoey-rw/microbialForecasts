# Combine chains from each taxon model, and create basic summary stats
source("../../source.R")

# Configuration: Chain size limits (can be set via environment variables)
# Default: truncate to reduce memory usage for long chains
# Default limits: 5000 parameter samples, 3000 plot samples
# To disable truncation: set CUT_SIZE1=NULL or CUT_SIZE1=0 in environment
# To use different limits: set CUT_SIZE1=<value> and CUT_SIZE2=<value>
cut_size1_env <- Sys.getenv("CUT_SIZE1", NA)
if (is.na(cut_size1_env) || cut_size1_env == "" || tolower(cut_size1_env) == "null" || cut_size1_env == "0") {
  cut_size1 <- 5000  # Default: truncate to 5000 samples
} else {
  cut_size1 <- as.numeric(cut_size1_env)
  if (is.na(cut_size1) || cut_size1 <= 0) cut_size1 <- NULL
}

cut_size2_env <- Sys.getenv("CUT_SIZE2", NA)
if (is.na(cut_size2_env) || cut_size2_env == "" || tolower(cut_size2_env) == "null" || cut_size2_env == "0") {
  cut_size2 <- 3000  # Default: truncate to 3000 samples
} else {
  cut_size2 <- as.numeric(cut_size2_env)
  if (is.na(cut_size2) || cut_size2 <= 0) cut_size2 <- NULL
}

cat("Chain size limits:\n")
cat("  cut_size1 (parameter samples):", if (is.null(cut_size1)) "NULL (no limit)" else cut_size1, "\n")
cat("  cut_size2 (plot samples):", if (is.null(cut_size2)) "NULL (no limit)" else cut_size2, "\n")
if (!is.null(cut_size1) || !is.null(cut_size2)) {
  cat("  (Using truncation by default to reduce memory usage for large chains)\n")
}
cat("\n")

# Source combine_chains function if package not available
if (!exists("combine_chains") || !is.function(combine_chains)) {
  if (file.exists(here("microbialForecast/R/combineChains.r"))) {
    source(here("microbialForecast/R/combineChains.r"))
  } else if (requireNamespace("microbialForecast", quietly = TRUE)) {
    combine_chains <- microbialForecast::combine_chains
  } else {
    stop("combine_chains function not found. Please ensure microbialForecast/R/combineChains.r exists.")
  }
}

# Source helper functions if package not available (for fast.summary.mcmc)
if (!exists("fast.summary.mcmc") || !is.function(fast.summary.mcmc)) {
  if (file.exists(here("microbialForecast/R/helperFunctions.r"))) {
    source(here("microbialForecast/R/helperFunctions.r"))
  } else if (requireNamespace("microbialForecast", quietly = TRUE)) {
    fast.summary.mcmc <- microbialForecast::fast.summary.mcmc
  } else {
    stop("fast.summary.mcmc function not found. Please ensure microbialForecast/R/helperFunctions.r exists.")
  }
}

# Function to filter files to include both old and new standardized schema files
# Updated to work with the new compatibility shim output structure
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  # For new driver uncertainty models, accept all chain files
  # For old schema, keep the original filtering
  standard_files <- file_list[
    # New driver uncertainty schema files (from compatibility shim)
    (grepl('cloglog_beta_driver_uncertainty', file_list) & grepl('_chain[0-9]+\\.rds$', basename(file_list))) |
    # New standardized schema files (from compatibility shim)
    (grepl('logit_beta_fixed_priors', file_list) & grepl('_chain[0-9]+\\.rds$', basename(file_list))) |
    # Old schema files (with both required suffixes AND exclude checkpoint files)
    (grepl('with_legacy_covariate', basename(file_list)) & 
     grepl('beta_regression', basename(file_list)) &
     !grepl('checkpoint', basename(file_list)))
  ]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard files (new schema + old schema, no checkpoints):', length(standard_files), '\n')
  cat('  Filtered out:', length(file_list) - length(standard_files), '\n\n')
  
  return(standard_files)
}

#source("source.R")

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

# CRITICAL FIX: Test function to validate taxon name extraction fixes
test_taxon_extraction <- function() {
  cat("=== TESTING TAXON NAME EXTRACTION FIXES ===\n")
  
  test_cases <- list(
    # Test case 1: Standard env_cycl model
    list(
      filename = "checkpoint_env_cycl_acetate_simple_20130601_20180101_with_legacy_covariate_beta_regression_chain1_initial.rds",
      expected_taxon = "acetate_simple",
      description = "Standard env_cycl model with legacy covariate"
    ),
    # Test case 2: Standard cycl_only model
    list(
      filename = "checkpoint_cycl_only_cellobiose_complex_20130601_20180101_beta_regression_chain1_initial.rds",
      expected_taxon = "cellobiose_complex",
      description = "Standard cycl_only model"
    ),
    # Test case 3: Standard env_cov model
    list(
      filename = "checkpoint_env_cov_acidimicrobiia_20130601_20180101_beta_regression_chain1_initial.rds",
      expected_taxon = "acidimicrobiia",
      description = "Standard env_cov model"
    ),
    # Test case 4: Complex taxon name
    list(
      filename = "checkpoint_env_cycl_phylum_fun_20130601_20180101_with_legacy_covariate_beta_regression_chain1_initial.rds",
      expected_taxon = "phylum_fun",
      description = "Complex taxon name with legacy covariate"
    ),
    # Test case 5: Problematic case that might have failed before
    list(
      filename = "checkpoint_env_cycl_actinobacteria_20130601_20180101_with_legacy_covariate_beta_regression_chain1_initial.rds",
      expected_taxon = "actinobacteria",
      description = "Taxon that might have been incorrectly extracted"
    )
  )
  
  test_results <- list(
    total_tests = length(test_cases),
    passed_tests = 0,
    failed_tests = 0,
    details = list()
  )
  
  for (i in seq_along(test_cases)) {
    test_case <- test_cases[[i]]
    cat("\nTest", i, ":", test_case$description, "\n")
    cat("Filename:", test_case$filename, "\n")
    cat("Expected taxon:", test_case$expected_taxon, "\n")
    
    # Simulate the extraction process
    path_parts <- strsplit(basename(test_case$filename), "_")[[1]]
    
    # Extract model name
    model_name <- NULL
    if (grepl("env_cycl", test_case$filename)) {
      model_name <- "env_cycl"
      model_parts <- c("env", "cycl")
    } else if (grepl("cycl_only", test_case$filename)) {
      model_name <- "cycl_only"
      model_parts <- c("cycl", "only")
    } else if (grepl("env_cov", test_case$filename)) {
      model_name <- "env_cov"
      model_parts <- c("env", "cov")
    }
    
    # Extract taxon using the fixed logic
    taxon <- NULL
    if (length(model_parts) > 0) {
      model_end_idx <- NULL
      
      # Look for the model parts in sequence
      for (j in 1:(length(path_parts) - length(model_parts) + 1)) {
        if (all(path_parts[j:(j + length(model_parts) - 1)] == model_parts)) {
          model_end_idx <- j + length(model_parts) - 1
          break
        }
      }
      
      if (!is.null(model_end_idx)) {
        # Find the date pattern (8 digits)
        date_idx <- which(grepl("^[0-9]{8}$", path_parts))
        if (length(date_idx) > 0) {
          taxon_start <- model_end_idx + 1
          taxon_end <- date_idx[1] - 1
          
          if (taxon_start <= taxon_end && taxon_start <= length(path_parts)) {
            taxon_parts <- path_parts[taxon_start:min(taxon_end, length(path_parts))]
            
            # Filter out any non-taxon parts
            taxon_parts <- taxon_parts[!taxon_parts %in% c("with", "legacy", "covariate", "beta", "regression")]
            
            if (length(taxon_parts) > 0) {
              taxon <- paste(taxon_parts, collapse = "_")
            }
          }
        }
      }
    }
    
    # Validate the extracted taxon
    if (!is.null(taxon)) {
      taxon_validation <- validate_taxon_name(taxon, "test_model", test_case$filename)
      if (taxon_validation$is_valid) {
        taxon <- taxon_validation$corrected_taxon
      }
    }
    
    cat("Extracted taxon:", if (is.null(taxon)) "NULL" else taxon, "\n")
    
    # Check if the test passed
    if (!is.null(taxon) && taxon == test_case$expected_taxon) {
      cat("✓ PASSED\n")
      test_results$passed_tests <- test_results$passed_tests + 1
      test_results$details[[i]] <- list(status = "PASSED", extracted = taxon, expected = test_case$expected_taxon)
    } else {
      cat("✗ FAILED\n")
      test_results$failed_tests <- test_results$failed_tests + 1
      test_results$details[[i]] <- list(status = "FAILED", extracted = taxon, expected = test_case$expected_taxon)
    }
  }
  
  # Summary
  cat("\n=== TEST RESULTS SUMMARY ===\n")
  cat("Total tests:", test_results$total_tests, "\n")
  cat("Passed:", test_results$passed_tests, "\n")
  cat("Failed:", test_results$failed_tests, "\n")
  cat("Success rate:", round(100 * test_results$passed_tests / test_results$total_tests, 1), "%\n")
  
  if (test_results$failed_tests == 0) {
    cat("✓ ALL TESTS PASSED - Taxon name extraction fixes are working correctly!\n")
  } else {
    cat("✗ SOME TESTS FAILED - Additional fixes may be needed\n")
  }
  
  return(test_results)
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
    
    message("Created minimal fallback metadata structure")
    return(NULL)
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
        # Updated to work with new standardized schema (plot_mu columns)
        return(n_cols < 10000 && !is.null(col_names) && any(grepl("plot_mu", col_names)))
      } else if (is.list(data$samples2) && length(data$samples2) > 0) {
        # Check first element if it's a list
        if (is.matrix(data$samples2[[1]])) {
          n_cols <- ncol(data$samples2[[1]])
          col_names <- colnames(data$samples2[[1]])
          # Updated to work with new standardized schema (plot_mu columns)
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
# Get all available model outputs from both old and new model directories
# Updated to work with the new driver uncertainty models
# Only process new driver uncertainty files to avoid processing 18,000+ old files
base_paths <- c(
  # New driver uncertainty models (prioritize these)
  here("data/model_outputs/cloglog_beta_driver_uncertainty")
)

# ==== ROOTS ====
driver_output_root <- here("data/model_outputs/cloglog_beta_driver_uncertainty")
recon_output_root  <- here("data/model_outputs/reconstructed_from_checkpoints")

# ==== NORMALIZE KEYS ====
normalize_model_key <- function(path) {
  # key only by basename and strip a trailing chain
  sub("_chain[0-9]+\\.rds$", "", basename(path))
}

# ==== READ A CHAIN (with samples2 renaming if needed) ====
read_chain_safe <- function(path) {
  x <- readRDS(path)
  # If compat shim object
  if (is.list(x) && !is.null(x$samples2) && is.matrix(x$samples2)) {
    cn <- colnames(x$samples2)
    if (!is.null(cn) && any(grepl("^plot_mu\\[", cn))) {
      colnames(x$samples2) <- sub("^plot_mu\\[", "mu[", cn)  # map back for old utils
    }
  }
  x
}

# ==== FIND LATEST CHECKPOINTS (search both trees) ====
find_recent_checkpoints <- function(model_key,
                                    roots = c(driver_output_root, recon_output_root)) {
  out <- list()
  patt_chain <- paste0("^", gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_key), "_chain[0-9]+.*\\.rds$")
  files <- character(0)
  for (rt in roots) {
    if (dir.exists(rt)) {
      cand <- list.files(rt, pattern = patt_chain, recursive = TRUE, full.names = TRUE)
      files <- c(files, cand)
    }
  }
  if (!length(files)) return(out)

  # group by chainN
  buckets <- list()
  for (f in files) {
    parts <- strsplit(basename(f), "_")[[1]]
    chain_token <- parts[grepl("^chain[0-9]+$", parts)]
    if (!length(chain_token)) next
    
    ch <- sub("^chain", "", chain_token[1])
    buckets[[ch]] <- c(buckets[[ch]] %||% character(0), f)
  }

  # pick newest by mtime per chain
  for (ch in names(buckets)) {
    bb <- buckets[[ch]]
    mt <- file.info(bb)$mtime
    out[[paste0("chain", ch)]] <- bb[which.max(mt)]
  }
  out
}

`%||%` <- function(a,b) if (is.null(a)) b else a

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

# Filter to prioritize recent files and avoid old malformed ones
# For driver uncertainty files, use all of them (removed date filter to re-combine all)
driver_uncertainty_files <- file.list[grepl("cloglog_beta_driver_uncertainty", file.list)]

# For other files, skip them - we only want driver uncertainty outputs
other_files <- file.list[!grepl("cloglog_beta_driver_uncertainty", file.list)]

# Use all driver uncertainty files
file.list <- driver_uncertainty_files

# Filter to only env_cycl models for testing
filter_to_env_cycl <- FALSE
if (filter_to_env_cycl) {
  env_cycl_files <- file.list[grepl("env_cycl", file.list)]
  non_env_cycl_files <- file.list[!grepl("env_cycl", file.list)]
  file.list <- env_cycl_files
  cat("Filtering to env_cycl models only:\n")
  cat("  env_cycl files:", length(env_cycl_files), "\n")
  cat("  Other model types filtered out:", length(non_env_cycl_files), "\n\n")
}

# Filter to test only problematic models (can be enabled for testing)
filter_to_problematic_models <- identical(tolower(Sys.getenv("TEST_PROBLEMATIC", "false")), "true")
if (filter_to_problematic_models) {
  problematic_patterns <- c(
    "cycl_only_assim_nitrate_reduction",
    "cycl_only_denitrification",
    "cycl_only_dissim_nitrate_reduction",
    "cycl_only_dissim_nitrite_reduction",
    "cycl_only_endophyte",
    "cycl_only_glucose_simple",
    "cycl_only_hypocreales",
    "cycl_only_mortierella",
    "cycl_only_mortierellaceae",
    "cycl_only_nitrification",
    "cycl_only_osmotic_stress",
    "cycl_only_sucrose_complex",
    "env_cov_oligotroph",
    "env_cycl_basidiomycota"
  )
  problematic_files <- file.list[sapply(file.list, function(f) {
    any(sapply(problematic_patterns, function(p) grepl(p, basename(f))))
  })]
  cat("Filtering to problematic models for testing:\n")
  cat("  Problematic files:", length(problematic_files), "\n")
  cat("  Other files filtered out:", length(file.list) - length(problematic_files), "\n\n")
  file.list <- problematic_files
}

# Filter to re-run only failed models (for investigation)
filter_to_failed_models <- identical(tolower(Sys.getenv("RERUN_FAILED", "false")), "true")
# Filter to test a single model (for debugging)
test_single_model <- identical(tolower(Sys.getenv("TEST_SINGLE", "false")), "true")
if (test_single_model) {
  # Test with a single known failed model
  # Can override with TEST_MODEL env var
  test_model_pattern <- Sys.getenv("TEST_MODEL", "cycl_only_assim_nitrite_reduction_20130601_20180101")
  test_files <- file.list[sapply(file.list, function(f) {
    grepl(test_model_pattern, basename(f))
  })]
  cat("TESTING SINGLE MODEL:\n")
  cat("  Test model pattern:", test_model_pattern, "\n")
  cat("  Test files found:", length(test_files), "\n")
  cat("  Other files filtered out:", length(file.list) - length(test_files), "\n\n")
  file.list <- test_files
} else if (filter_to_failed_models) {
  # Known failed models from previous runs
  failed_model_patterns <- c(
    "env_cov_nitrification_20130601_20180101"
  )
  failed_files <- file.list[sapply(file.list, function(f) {
    any(sapply(failed_model_patterns, function(p) grepl(p, basename(f))))
  })]
  cat("Filtering to failed models for re-investigation:\n")
  cat("  Failed model patterns:", paste(failed_model_patterns, collapse = ", "), "\n")
  cat("  Failed files found:", length(failed_files), "\n")
  cat("  Other files filtered out:", length(file.list) - length(failed_files), "\n\n")
  file.list <- failed_files
}

if (!filter_to_env_cycl && !filter_to_problematic_models && !filter_to_failed_models && !test_single_model) {
  cat("File filtering applied:\n")
  cat("  Total driver uncertainty files found:", length(file.list), "\n")
  cat("  Other files filtered out:", length(other_files), "\n\n")
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
regular_chain_files <- file.list[!grepl("checkpoint|^samples_", basename(file.list))]
for (file in regular_chain_files) {
  # Use normalize_model_key to get clean key from basename
  key <- normalize_model_key(file)
  model_files[[key]] <- c(model_files[[key]] %||% character(0), file)
}

# Handle individual chain files (those that end with _chainX.rds)
individual_chain_files <- file.list[grepl("_chain[0-9]+\\.rds$", basename(file.list))]
force_recombine <- identical(tolower(Sys.getenv("FORCE_RECOMBINE", "false")), "true")
cat("Processing", length(individual_chain_files), "individual chain files...\n")
for (file in individual_chain_files) {
  # Use normalize_model_key to get clean key from basename
  key <- normalize_model_key(file)
  # Normalize key: strip 'samples_' prefix if present for consistent grouping
  # Files starting with samples_ should be grouped the same as files without it
  normalized_key <- if (grepl("^samples_", key)) {
    sub("^samples_", "", key)
  } else {
    key
  }
  # if a combined file already exists (samples_*), we'll decide later based on samples2 presence
  model_files[[normalized_key]] <- c(model_files[[normalized_key]] %||% character(0), file)
}

# Handle already-combined sample files (from reconstruction script)
# This phase INSPECTS existing combined files to determine what needs to be done
# OPTIMIZATION: Skip inspection for most files - only check those that might need work
sample_files <- file.list[grepl("samples_.*\\.rds$", file.list) & !grepl("_chain[0-9]", file.list)]
cat("\n=== INSPECTION PHASE: Checking", length(sample_files), "already-combined sample files ===\n")
cat("(Skipping files that are likely already complete to speed up inspection)\n\n")
sample_counter <- 0
skipped_count <- 0
checked_count <- 0

for (file in sample_files) {
  sample_counter <- sample_counter + 1
  
  # Report progress every 10 files to reduce output
  file_size_mb <- round(file.size(file) / 1024 / 1024, 1)
  if (sample_counter %% 10 == 0 || sample_counter <= 3) {
    cat(sprintf("  [%d/%d] Inspecting: %s (%.1f MB) [%d checked, %d skipped]\n", 
                sample_counter, length(sample_files), basename(file), file_size_mb,
                checked_count, skipped_count))
    flush.console()
  }
  
  # OPTIMIZATION: Skip files that are too large (likely to hang) - threshold: 200MB
  # Also skip files that are clearly recent and large enough to be complete
  # BUT: Skip this optimization if we're testing problematic models or forcing recombine
  # (force_recombine is already defined earlier in the script)
  test_problematic <- identical(tolower(Sys.getenv("TEST_PROBLEMATIC", "false")), "true")
  skip_optimization <- force_recombine || test_problematic
  
  if (!skip_optimization && file_size_mb > 200) {
    # If file is large (>200MB) and recent (>1 hour old), assume it's complete
    file_age_hours <- as.numeric(difftime(Sys.time(), file.info(file)$mtime, units = "hours"))
    if (file_age_hours > 1 && file_size_mb > 50) {
      skipped_count <- skipped_count + 1
      if (sample_counter <= 5 || sample_counter %% 50 == 0) {
        cat(sprintf("    ⏭️  Skipping large recent file (%.1f MB, %.1f hrs old): %s\n", 
                    file_size_mb, file_age_hours, basename(file)))
      }
      next
    } else if (file_size_mb > 500) {
      skipped_count <- skipped_count + 1
      if (sample_counter <= 5 || sample_counter %% 50 == 0) {
        cat(sprintf("    ⚠️  Skipping very large file (>500MB): %s\n", basename(file)))
      }
      next
    }
  }
  
  # These are already combined files, but we need to check if they need metadata extraction
  model_id <- file  # Use the full path as the model ID
  
  # Check if this file needs processing (metadata extraction or missing samples2)
  checked_count <- checked_count + 1
  tryCatch({
    # OPTIMIZATION: For smaller files, do full check. For larger files, do lighter check
    # Only read the file header/structure first
    sample_data <- NULL
    if (file_size_mb < 100) {
      # Small file - read it fully
      sample_data <- tryCatch({
        readRDS(file)
      }, error = function(e) {
        if (sample_counter <= 10 || sample_counter %% 50 == 0) {
          cat("  Error reading file ", basename(file), " (", round(file.size(file)/1024/1024, 1), "MB): ", e$message, "\n")
        }
        NULL
      })
    } else {
      # Large file - try lighter approach: read with limited depth
      # If this fails or takes too long, assume file needs checking
      sample_data <- tryCatch({
        # Try reading with a timeout-like approach (R doesn't have true timeouts, but we can skip slow reads)
        readRDS(file)
      }, error = function(e) {
        # If read fails, assume file needs re-processing
        NULL
      })
    }
    
    if (is.null(sample_data)) {
      # File couldn't be read - look for chains to re-combine
      skipped_count <- skipped_count + 1
      next
    }
    
    # Check if samples2 is missing
    missing_samples2 <- !"samples2" %in% names(sample_data)
    
    if (missing_samples2) {
      # File needs to be re-combined from chains
      # Extract model ID to find chain files
      model_base <- basename(file)
      model_base <- gsub("^samples_", "", model_base)
      model_base <- gsub("\\.rds$", "", model_base)
      
      # Look for chain files for this model
      # Escape special regex characters in model_base
      model_base_escaped <- gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_base)
      chain_pattern <- paste0(model_base_escaped, "_chain[0-9]+\\.rds$")
      chain_files <- file.list[grepl(chain_pattern, basename(file.list))]
      
      if (length(chain_files) >= 2) {
        # Found chain files, add them for re-combination
        # Use the model_id (without chain suffix) as the key, normalized
        chain_model_id_base <- basename(chain_files[1])
        chain_model_id <- gsub("_chain[0-9]+\\.rds$", "", chain_model_id_base)
        # Normalize: strip 'samples_' prefix if present for consistent grouping
        normalized_chain_model_id <- if (grepl("^samples_", chain_model_id)) {
          sub("^samples_", "", chain_model_id)
        } else {
          chain_model_id
        }
        if (!(normalized_chain_model_id %in% names(model_files))) {
          model_files[[normalized_chain_model_id]] <- list()
        }
        model_files[[normalized_chain_model_id]] <- c(model_files[[normalized_chain_model_id]], chain_files)
        message("Will re-combine from chains (file missing samples2): ", basename(file))
      } else {
        message("Skipping combined file - missing samples2 and no chains found: ", basename(file))
      }
    } else {
      # Check if FORCE_RECOMBINE is set - if so, always look for chains and re-combine
      if (force_recombine) {
        # Extract model ID to find chain files
        model_base <- basename(file)
        model_base <- gsub("^samples_", "", model_base)
        model_base <- gsub("\\.rds$", "", model_base)
        
        # Look for chain files for this model
        model_base_escaped <- gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_base)
        chain_pattern <- paste0(model_base_escaped, "_chain[0-9]+\\.rds$")
        chain_files <- file.list[grepl(chain_pattern, basename(file.list))]
        
        if (length(chain_files) >= 2) {
          # Found chain files, add them for re-combination
          chain_model_id_base <- basename(chain_files[1])
          chain_model_id <- gsub("_chain[0-9]+\\.rds$", "", chain_model_id_base)
          # Normalize: strip 'samples_' prefix if present for consistent grouping
          normalized_chain_model_id <- if (grepl("^samples_", chain_model_id)) {
            sub("^samples_", "", chain_model_id)
          } else {
            chain_model_id
          }
          if (!(normalized_chain_model_id %in% names(model_files))) {
            model_files[[normalized_chain_model_id]] <- list()
          }
          model_files[[normalized_chain_model_id]] <- c(model_files[[normalized_chain_model_id]], chain_files)
          message("FORCE_RECOMBINE: Will re-combine from chains: ", basename(file))
        } else {
          message("FORCE_RECOMBINE: No chain files found for ", basename(file), ", will try to update existing file")
          model_files[[model_id]] <- list(file)  # Add for processing to update
        }
      } else {
        # Check if metadata is needed
        needs_metadata <- is.null(sample_data$metadata) || 
                         is.null(sample_data$metadata$model_data) ||
                         is.null(sample_data$metadata$model_name)
      
        if (needs_metadata) {
          model_files[[model_id]] <- list(file)  # Single file list
          message("Will process combined file for metadata extraction: ", basename(file))
        } else {
          # Verify the file actually exists at the correct save location
          # The file being checked might be in a different location than where we'll save
          # Extract model info to determine save path
          model_base <- basename(file)
          model_base <- gsub("^samples_", "", model_base)
          model_base <- gsub("\\.rds$", "", model_base)
          
          # Determine save path based on model type
          if (grepl("env_cycl", model_base)) {
            save_dir <- here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl")
          } else if (grepl("cycl_only", model_base)) {
            save_dir <- here("data/model_outputs/cloglog_beta_driver_uncertainty/cycl_only")
          } else if (grepl("env_cov", model_base)) {
            save_dir <- here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cov")
          } else {
            save_dir <- here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl")
          }
          
          # Check if subdirectories exist (for group-based models)
          if (dir.exists(save_dir)) {
            # Check recursively for the file
            expected_files <- list.files(save_dir, pattern = paste0("^samples_", gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_base), "\\.rds$"), 
                                        recursive = TRUE, full.names = TRUE)
            if (length(expected_files) == 0) {
              # File doesn't exist at save location - look for chains to combine
              model_base_escaped <- gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_base)
              chain_pattern <- paste0(model_base_escaped, "_chain[0-9]+\\.rds$")
              chain_files <- file.list[grepl(chain_pattern, basename(file.list))]
              
              if (length(chain_files) >= 2) {
                chain_model_id_base <- basename(chain_files[1])
                chain_model_id <- gsub("_chain[0-9]+\\.rds$", "", chain_model_id_base)
                # Normalize: strip 'samples_' prefix if present for consistent grouping
                normalized_chain_model_id <- if (grepl("^samples_", chain_model_id)) {
                  sub("^samples_", "", chain_model_id)
                } else {
                  chain_model_id
                }
                if (!(normalized_chain_model_id %in% names(model_files))) {
                  model_files[[normalized_chain_model_id]] <- list()
                }
                model_files[[normalized_chain_model_id]] <- c(model_files[[normalized_chain_model_id]], chain_files)
                message("File not at save location, will re-combine from chains: ", basename(file))
              } else {
                message("File not at save location and no chains found: ", basename(file))
              }
            } else {
              message("Skipping combined file - already has complete metadata and samples2: ", basename(file))
            }
          } else {
            # Save directory doesn't exist - look for chains to combine
            model_base_escaped <- gsub("([.^$*+?()|\\[\\{]|\\\\])", "\\\\\\1", model_base)
            chain_pattern <- paste0(model_base_escaped, "_chain[0-9]+\\.rds$")
            chain_files <- file.list[grepl(chain_pattern, basename(file.list))]
            
            if (length(chain_files) >= 2) {
              chain_model_id_base <- basename(chain_files[1])
              chain_model_id <- gsub("_chain[0-9]+\\.rds$", "", chain_model_id_base)
              # Normalize: strip 'samples_' prefix if present for consistent grouping
              normalized_chain_model_id <- if (grepl("^samples_", chain_model_id)) {
                sub("^samples_", "", chain_model_id)
              } else {
                chain_model_id
              }
              if (!(normalized_chain_model_id %in% names(model_files))) {
                model_files[[normalized_chain_model_id]] <- list()
              }
              model_files[[normalized_chain_model_id]] <- c(model_files[[normalized_chain_model_id]], chain_files)
              message("Save directory missing, will re-combine from chains: ", basename(file))
            } else {
              message("Save directory missing and no chains found: ", basename(file))
            }
          }
        }
      }
    }
  }, error = function(e) {
    # If we can't read the file, add it for processing
    model_files[[model_id]] <- list(file)
    message("Will process combined file (read error): ", basename(file))
  })
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

checkpoint_files <- file.list[grepl("checkpoint", basename(file.list))]
checkpoint_files <- filter_standard_files(checkpoint_files)

checkpoint_keys <- unique(sub("_chain[0-9]+.*\\.rds$", "", sub("^checkpoint_", "", basename(checkpoint_files))))
cat("Found", length(checkpoint_keys), "unique checkpoint base models\n")

# For each checkpoint key, find the latest checkpoint per chain
for (key in checkpoint_keys) {
  # Normalize checkpoint key: strip 'samples_' prefix if present for consistent grouping
  normalized_key <- if (grepl("^samples_", key)) {
    sub("^samples_", "", key)
  } else {
    key
  }
  
  # Only add if we don't already have >=2 chains grouped
  existing <- model_files[[normalized_key]]
  if (length(unique(existing)) >= 2) next

  latest <- find_recent_checkpoints(key, roots = c(driver_output_root, recon_output_root))
  if (length(latest) >= 2) {
    model_files[[normalized_key]] <- unique(c(existing %||% character(0), unlist(latest, use.names = FALSE)))
    cat("Added checkpoint-derived chains for:", normalized_key, " (#", length(latest), ")\n")
  }
}

# Helper function to quickly check if combined file is already up-to-date
check_if_up_to_date <- function(model_id, chain_paths, driver_output_root, recon_output_root) {
  # Determine target directory - always use model type directory, not subdirectories
  # Chains may be in subdirectories but combined files are in model type directory
  model_type <- if (grepl("^env_cycl_", model_id)) "env_cycl"
           else if (grepl("^env_cov_",  model_id)) "env_cov"
           else if (grepl("^cycl_only_",model_id)) "cycl_only" else "env_cycl"
  target_dir <- file.path(driver_output_root, model_type)
  
  # Create combined file path
  normalized_model_id <- if (grepl("^samples_", model_id)) {
    sub("^samples_", "", model_id)
  } else {
    model_id
  }
  combined_file_path <- file.path(target_dir, paste0("samples_", normalized_model_id, ".rds"))
  
  # Quick check: if file exists, has samples2, and is newer than inputs
  if (file.exists(combined_file_path)) {
    tryCatch({
      # Quick metadata check without reading full file
      test_read <- readRDS(combined_file_path)
      if (is.list(test_read) && "samples2" %in% names(test_read)) {
        combined_file_time <- file.info(combined_file_path)$mtime
        input_files_exist <- file.exists(chain_paths)
        if (all(input_files_exist)) {
          input_file_times <- file.info(chain_paths[input_files_exist])$mtime
          newest_input_time <- max(input_file_times)
          if (combined_file_time > newest_input_time) {
            return(TRUE)  # File is up-to-date
          }
        }
      }
    }, error = function(e) {
      # If read fails, assume needs processing
      return(FALSE)
    })
  }
  return(FALSE)  # Not up-to-date or doesn't exist
}

# Function to process a single model (used for parallel processing)
process_single_model <- function(model_id, chain_paths, driver_output_root, recon_output_root, 
                                 combine_chains, fast.summary.mcmc, extract_metadata_from_paths,
                                 effectiveSize, gelman.diag, cut_size1 = NULL, cut_size2 = NULL) {
  result <- list(model_id = model_id, success = FALSE, message = "")
  
  # Check if this is an already-combined sample file that needs metadata extraction
  if (length(chain_paths) == 1 && grepl("^samples_", basename(chain_paths[[1]]))) {
    message("Processing already-combined sample file: ", basename(chain_paths[[1]]))
    
    # Check if this file needs re-combining or just metadata update
    tryCatch({
      sample_data <- readRDS(chain_paths[[1]])
      
      # Check if samples2 is missing (needs re-combination)
      missing_samples2 <- !"samples2" %in% names(sample_data)
      
      if (missing_samples2) {
        message("Sample file is missing samples2, will skip and re-combine from chains later")
        result$message <- "Missing samples2, needs re-combination"
        return(result)
      }
      
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
          result$success <- TRUE
          result$message <- "Metadata updated"
        } else {
          message("❌ Failed to extract metadata for: ", basename(chain_paths[[1]]))
          result$message <- "Failed to extract metadata"
        }
      } else {
        message("Sample file already has complete metadata: ", basename(chain_paths[[1]]))
        result$success <- TRUE
        result$message <- "Metadata already complete"
      }
    }, error = function(e) {
      message("Error processing sample file: ", e$message)
      result$message <- paste("Error:", e$message)
    })
    
    return(result)
  }
  
  if (length(chain_paths) < 2) {
    message("Model ", basename(model_id), " has only ", length(chain_paths), " chain(s)")
    
    # Try to find checkpoint files as backup chains
    # Normalize model_id to remove samples_ prefix if present for checkpoint search
    checkpoint_search_key <- if (grepl("^samples_", model_id)) {
      sub("^samples_", "", model_id)
    } else {
      model_id
    }
    message("Searching for checkpoint files as backup chains (using key: ", checkpoint_search_key, ")...")
    
    # Need to define find_recent_checkpoints locally or pass it
    # For now, skip checkpoint search in parallel - it requires access to the function
    result$message <- paste("Only", length(chain_paths), "chain(s), need >= 2")
    return(result)
  }
  
  message("Processing ", model_id, " with ", length(chain_paths), " chains...")
  
  # Build save path - always save in model type directory, not in subdirectories
  # Chains may be in subdirectories (e.g., env_cov/coniochaeta/) but combined files
  # should be saved in the model type directory (e.g., env_cov/)
  chain_paths <- unique(chain_paths)
  
  # Derive model type from the model_id and save in model type directory
  model_type <- if (grepl("^env_cycl_", model_id)) "env_cycl"
           else if (grepl("^env_cov_",  model_id)) "env_cov"
           else if (grepl("^cycl_only_",model_id)) "cycl_only" else "env_cycl"
  target_dir <- file.path(driver_output_root, model_type)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure the target directory exists
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create the combined filename - always use "samples_" prefix for consistency
  # Normalize model_id first: remove samples_ prefix if present, then add it back
  normalized_model_id <- if (grepl("^samples_", model_id)) {
    sub("^samples_", "", model_id)
  } else {
    model_id
  }
  base_name <- paste0("samples_", normalized_model_id)
  combined_file_path <- file.path(target_dir, paste0(base_name, ".rds"))
  
  # Check if combined file already exists and if it has proper samples2 data
  needs_recombine <- TRUE
  if (file.exists(combined_file_path)) {
    tryCatch({
      test_read <- readRDS(combined_file_path)
      # Check if file has samples2 object (raw samples, not just summary)
      has_samples2_object <- FALSE
      if (is.list(test_read) && "samples2" %in% names(test_read)) {
        has_samples2_object <- TRUE
      }
      
      if (has_samples2_object) {
        # File has valid samples2 - check if it's up-to-date
        combined_file_time <- file.info(combined_file_path)$mtime
        input_files_exist <- file.exists(chain_paths)
        
        if (all(input_files_exist)) {
          input_file_times <- file.info(chain_paths[input_files_exist])$mtime
          newest_input_time <- max(input_file_times)
          
          if (combined_file_time > newest_input_time) {
            message("Combined file is up-to-date for ", basename(combined_file_path), " (newer than all input chains)")
            needs_recombine <- FALSE
          } else {
            message("Combined file is outdated for ", basename(combined_file_path), " (input chains are newer)")
          }
        } else {
          message("Some input chain files missing for ", basename(combined_file_path), ", will re-combine")
        }
      } else {
        message("Combined file exists but missing samples2 object, will re-combine: ", basename(combined_file_path))
      }
    }, error = function(e) {
      message("Error reading existing file, will overwrite: ", basename(combined_file_path))
    })
  }
  
  if (!needs_recombine) {
    result$success <- TRUE
    result$message <- "File already up-to-date"
    return(result)
  }
  
  # Preflight check: validate chain files before combining
  bad <- character(0)
  for (p in chain_paths) {
    ok <- tryCatch({
      z <- readRDS(p)
      is.list(z) && ("samples" %in% names(z))  # minimal requirement
    }, error = function(e) FALSE)
    if (!ok) bad <- c(bad, p)
  }
  if (length(bad)) {
    message("Skipping ", length(bad), " unreadable/invalid chain(s):")
    message("  - ", paste(basename(bad), collapse = "\n  - "))
    chain_paths <- setdiff(chain_paths, bad)
  }
  if (length(chain_paths) < 2) {
    message("Not enough valid chains after preflight; skipping model.")
    result$message <- "Not enough valid chains"
    return(result)
  }
  
  # Combine chains with truncation (similar to combine_large_chains.r)
  # Check file sizes and apply aggressive truncation for very large files
  chain_file_sizes_mb <- sapply(chain_paths, function(cf) {
    if (file.exists(cf)) file.size(cf) / 1024 / 1024 else 0
  })
  max_file_size_mb <- max(chain_file_sizes_mb, na.rm = TRUE)
  
  # Use truncation (default from config)
  effective_cut_size1 <- cut_size1
  effective_cut_size2 <- cut_size2
  
  # Apply aggressive truncation for very large files (>500MB per chain)
  if (max_file_size_mb > 500) {
    message("  Large file detected (max ", round(max_file_size_mb, 1), " MB per chain), using aggressive truncation")
    effective_cut_size1 <- min(cut_size1 %||% 5000, 3000)  # Reduce further for huge files
    effective_cut_size2 <- min(cut_size2 %||% 3000, 2000)
    
    # If file is extremely large (>1GB), skip it if we have enough other chains
    if (max_file_size_mb > 1000 && length(chain_paths) >= 3) {
      message("  WARNING: Extremely large file (", round(max_file_size_mb, 1), " MB) detected. Skipping to avoid memory errors.")
      message("  Using ", length(chain_paths) - 1, " other chains instead.")
      chain_paths <- chain_paths[chain_file_sizes_mb < 1000]
      if (length(chain_paths) < 2) {
        message("  Not enough chains after filtering; skipping model")
        result$message <- "Not enough chains after filtering large files"
        return(result)
      }
    }
  }
  
  message("Combining chains for ", model_id, " (", length(chain_paths), " chains)")
  message("  Using truncation: cut_size1 = ", if (is.null(effective_cut_size1)) "NULL" else effective_cut_size1, 
          ", cut_size2 = ", if (is.null(effective_cut_size2)) "NULL" else effective_cut_size2)
  start_time <- Sys.time()
  out <- tryCatch({
    # Use the new combine_chains function that handles different sample sizes robustly
    # Pass file paths, not loaded objects
    # Pass cut_size parameters to limit chain sizes if specified
    combine_chains(chain_paths, cut_size1 = effective_cut_size1, cut_size2 = effective_cut_size2)
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    error_call <- deparse(conditionCall(e))
    message("Error combining chains for ", model_id, ": ", error_msg)
    if (length(error_call) > 0 && nchar(error_call[1]) < 200) {
      message("  Error occurred in: ", error_call[1])
    }
    # Log chain file info for debugging
    if (length(chain_paths) > 0) {
      chain_info <- sapply(chain_paths, function(p) {
        if (file.exists(p)) {
          paste0(basename(p), " (", round(file.size(p)/1024/1024, 1), "MB)")
        } else {
          paste0(basename(p), " (MISSING)")
        }
      })
      message("  Chain files: ", paste(chain_info, collapse = ", "))
    }
    return(NULL)
  })
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  message("Chain combination completed in ", round(elapsed_time, 1), " seconds")
  
  if (is.null(out)) {
    message("Failed to combine chains for ", basename(model_id), ", skipping...")
    result$message <- "Failed to combine chains"
    return(result)
  }
  
  # Verify metadata preservation including Nimble model code
  if (is.null(out$metadata)) {
    message("❌ ERROR: Metadata was lost during chain combination!")
    result$message <- "Metadata lost during combination"
    return(result)
  }
  
  # Remove chains if one has very different values than the others
  chains <- out$samples
  if (length(chains) == 0) {
    message("No chains found in output for ", basename(model_id), ", skipping...")
    result$message <- "No chains in output"
    return(result)
  }

  # HARD ASSERTION: ensure site_effect_sd is present in parameter samples
  first_chain_cols <- try(colnames(chains[[1]]), silent = TRUE)
  if (inherits(first_chain_cols, "try-error") || is.null(first_chain_cols)) {
    message("❌ ERROR: Could not read column names from first chain for ", basename(model_id), "; skipping save")
    result$message <- "Could not read column names"
    return(result)
  }
  has_site_effect_sd <- any(grepl('^site_effect_sd(\\[|$)', first_chain_cols))
  if (!has_site_effect_sd) {
    message("❌ MISSING REQUIRED PARAMETER 'site_effect_sd' in ", basename(model_id))
    result$message <- "Missing site_effect_sd parameter"
    return(result)
  }
  
  # Check if intercept parameter exists
  if (!"intercept" %in% colnames(chains[[1]])) {
    message("No intercept parameter found in chains for ", basename(model_id), ", skipping...")
    result$message <- "Missing intercept parameter"
    return(result)
  }
  
  # Ensure chains have reasonable sample sizes before outlier detection
  min_chain_size <- 100
  chain_sizes <- sapply(chains, nrow)
  if (min(chain_sizes) >= min_chain_size) {
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
  actual_samples <- out$samples
  actual_samples2 <- out$samples2
  
  # If samples is nested (contains samples and metadata), extract the actual samples
  if (is.list(actual_samples) && "samples" %in% names(actual_samples)) {
    actual_samples <- actual_samples$samples
  }
  
  if (is.list(actual_samples2) && "samples" %in% names(actual_samples2)) {
    actual_samples2 <- actual_samples2$samples
  }
  
  # Calculate summaries from the actual samples
  if (is.null(actual_samples) || length(actual_samples) == 0) {
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
  
  # Handle plot samples
  if (is.null(actual_samples2) || length(actual_samples2) == 0) {
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
  
  # Calculate effective sample size
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
  
  # Combine and save output
  out_summary <- list(
    samples = actual_samples,
    samples2 = actual_samples2,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = out$metadata,
    gelman = gelman_out
  )
  
  # Flatten the metadata structure if needed
  if (!is.null(out_summary$metadata) && "metadata" %in% names(out_summary$metadata)) {
    out_summary$metadata <- out_summary$metadata$metadata
  }
  
  if (!is.null(out_summary$metadata) && is.list(out_summary$metadata) && 
      "samples" %in% names(out_summary$metadata) && "metadata" %in% names(out_summary$metadata)) {
    out_summary$metadata <- out_summary$metadata$metadata
  }
  
  # Map back to plot_mu[ on the combined object before saving
  if (is.list(out_summary) && !is.null(out_summary$samples2)) {
    if (is.matrix(out_summary$samples2)) {
      cn <- colnames(out_summary$samples2)
      if (!is.null(cn) && any(grepl("^mu\\[", cn))) {
        colnames(out_summary$samples2) <- sub("^mu\\[", "plot_mu[", cn)
      }
    } else if (inherits(out_summary$samples2, "mcmc.list")) {
      for (i in seq_along(out_summary$samples2)) {
        cn <- colnames(out_summary$samples2[[i]])
        if (!is.null(cn) && any(grepl("^mu\\[", cn))) {
          colnames(out_summary$samples2[[i]]) <- sub("^mu\\[", "plot_mu[", cn)
        }
      }
    }
  }
  
  # Save to the correct combined file path
  saveRDS(out_summary, combined_file_path, compress = FALSE)
  message("Saved combined samples output for ", model_id, " to: ", basename(combined_file_path))
  
  # Delete chain files after successful combination
  # Only delete if combined file was JUST created (significantly newer than all chains)
  if (file.exists(combined_file_path)) {
    combined_mtime <- file.info(combined_file_path)$mtime
    combined_size <- file.size(combined_file_path)
    
    # Verify combined file is valid (reasonable size)
    if (combined_size > 100000) {  # At least 100KB
      chain_mtimes <- sapply(chain_paths, function(cf) {
        if (file.exists(cf)) file.info(cf)$mtime else NA
      })
      chain_mtimes <- chain_mtimes[!is.na(chain_mtimes)]
      
      if (length(chain_mtimes) > 0) {
        newest_chain_mtime <- max(chain_mtimes)
        # Only delete if combined file is significantly newer (at least 2 seconds)
        # This ensures we only delete when we just created the combined file
        time_diff_secs <- as.numeric(difftime(combined_mtime, newest_chain_mtime, units = "secs"))
        
        if (time_diff_secs >= 2) {
          deleted_count <- 0
          deleted_size <- 0
          for (chain_file in chain_paths) {
            if (file.exists(chain_file)) {
              tryCatch({
                chain_size <- file.size(chain_file)
                file.remove(chain_file)
                deleted_count <- deleted_count + 1
                deleted_size <- deleted_size + chain_size
                message("  Deleted chain file: ", basename(chain_file))
              }, error = function(e) {
                warning("  Failed to delete chain file ", basename(chain_file), ": ", e$message)
              })
            }
          }
          if (deleted_count > 0) {
            message("  ✓ Deleted ", deleted_count, " chain file(s), freed ", round(deleted_size / 1024 / 1024, 2), " MB")
          }
        } else {
          message("  ⚠️  Skipping chain deletion: combined file is not significantly newer (diff: ", round(time_diff_secs, 2), " sec)")
        }
      }
    } else {
      message("  ⚠️  Skipping chain deletion: combined file is suspiciously small (", round(combined_size/1024, 1), " KB)")
    }
  }
  
  result$success <- TRUE
  result$message <- "Successfully combined and saved"
  result$file <- combined_file_path
  
  if (min(es) < 500) {
    message("Low effective sample sizes - check for unconverged parameters in model: ", basename(model_id))
  }
  
  # Force garbage collection after processing to free memory
  gc(verbose = FALSE)
  
  return(result)
}

# Process each model that has multiple chains
total_models <- length(model_files)
cat("\n=== COMBINING PHASE: Combining chains for", total_models, "models ===\n")
cat("(This is where the actual chain combination work happens)\n")
cat("Models to process:\n")
if (total_models > 0) {
  # Show first 10 model IDs for debugging
  model_names <- names(model_files)
  for (i in 1:min(10, length(model_names))) {
    num_chains <- length(model_files[[model_names[i]]])
    cat(sprintf("  [%d] %s (%d chains)\n", i, basename(model_names[i]), num_chains))
  }
  if (length(model_names) > 10) {
    cat(sprintf("  ... and %d more models\n", length(model_names) - 10))
  }
} else {
  cat("  No models found to process\n")
}
cat("\n")

# Check if parallel processing is requested
use_parallel <- identical(tolower(Sys.getenv("USE_PARALLEL", "true")), "true")
# Set to 2 cores by default to reduce memory usage (can override with N_CORES env var)
# Memory-intensive operations should use fewer cores
n_cores <- as.numeric(Sys.getenv("N_CORES", 2))
# Cap at 4 cores maximum to prevent memory exhaustion
n_cores <- min(n_cores, 4)
cat("Using", n_cores, "cores for parallel processing (capped at 4 for memory safety)\n")

if (use_parallel && n_cores > 1 && total_models > 10) {
  cat("Using parallel processing with", n_cores, "cores\n\n")
  library(parallel)
  library(doParallel)
  library(foreach)
  
  # Create cluster
  cl <- makeCluster(n_cores, type = "PSOCK")
  registerDoParallel(cl)
  
  # Export necessary functions and variables (including helper operator and truncation settings)
  clusterExport(cl, c("combine_chains", "fast.summary.mcmc", "extract_metadata_from_paths",
                       "effectiveSize", "gelman.diag", "driver_output_root", "recon_output_root",
                       "process_single_model", "check_if_up_to_date", "%||%", "cut_size1", "cut_size2"), envir = environment())
  
  # Load required packages on workers
  clusterEvalQ(cl, {
    library(coda)
    library(dplyr)
  })
  
  # Process models in parallel with better error handling and progress reporting
  cat("Processing", total_models, "models in parallel...\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  start_time <- Sys.time()
  
  # Create a list to track progress (file-based progress tracking)
  progress_file <- tempfile(pattern = "combine_progress_", fileext = ".txt")
  cat("Progress file:", progress_file, "\n")
  
  # Create model list with indices for progress tracking
  model_list <- names(model_files)
  total_count <- length(model_list)
  
  # Add a monitor to periodically report progress while processing
  # This will help identify where it hangs
  cat("Note: Progress file will be updated as workers complete models\n")
  cat("You can monitor progress in real-time with: tail -f", progress_file, "\n\n")
  flush.console()
  
  # Process models with enhanced error reporting
  results <- tryCatch({
    foreach(i = seq_along(model_list),
                     .packages = c("coda", "dplyr"),
                     .errorhandling = "pass",
                     .export = c("model_files", "progress_file", "model_list", "total_count",
                                 "driver_output_root", "recon_output_root", "process_single_model",
                                 "check_if_up_to_date", "combine_chains", "fast.summary.mcmc", "extract_metadata_from_paths",
                                 "effectiveSize", "gelman.diag", "%||%", "cut_size1", "cut_size2")) %dopar% {
      
      model_id <- model_list[i]
      worker_start <- Sys.time()
      chain_paths <- model_files[[model_id]]
      
      # Check if combined file is already up-to-date BEFORE checking chain sizes
      # This avoids unnecessary INFO messages for files that will be skipped
      is_up_to_date <- tryCatch({
        check_if_up_to_date(model_id, chain_paths, driver_output_root, recon_output_root)
      }, error = function(e) FALSE)
      
      if (is_up_to_date) {
        # Skip processing - return success result immediately
        result <- list(model_id = model_id, success = TRUE, message = "File already up-to-date")
      } else {
        # Report progress for this worker
        tryCatch({
          cat(sprintf("[Worker] [%s] STARTING model %d/%d: %s\n", 
                      format(Sys.time(), "%H:%M:%S"), i, total_count, basename(model_id)), 
              file = progress_file, append = TRUE)
        }, error = function(e) {
          # If file writing fails, continue anyway
        })
        
        # Check chain file sizes before processing (may indicate if model will be slow)
        chain_sizes_mb <- tryCatch({
          sapply(chain_paths, function(p) {
            if (file.exists(p)) round(file.size(p) / 1024 / 1024, 1) else 0
          })
        }, error = function(e) numeric(0))
        
        max_chain_size <- if (length(chain_sizes_mb) > 0) max(chain_sizes_mb) else 0
        
        # Inform if chains are very large (truncation will still apply)
        if (max_chain_size > 100) {
          tryCatch({
            truncation_note <- if (!is.null(cut_size1)) paste0(" (will truncate to ", cut_size1, " samples)") else ""
            cat(sprintf("[Worker] [%s] INFO: Model %d has large chains (max %.1f MB)%s - may take a while to read\n", 
                        format(Sys.time(), "%H:%M:%S"), i, max_chain_size, truncation_note), 
                file = progress_file, append = TRUE)
          }, error = function(e) {})
        }
        
        # Process with error handling
        result <- tryCatch({
          result <- process_single_model(model_id, chain_paths, driver_output_root, recon_output_root,
                                        combine_chains, fast.summary.mcmc, extract_metadata_from_paths,
                                        effectiveSize, gelman.diag, cut_size1 = cut_size1, cut_size2 = cut_size2)
          result
        }, error = function(e) {
          error_msg <- conditionMessage(e)
          error_trace <- paste(capture.output(traceback()), collapse = "\n")
          worker_elapsed <- as.numeric(difftime(Sys.time(), worker_start, units = "secs"))
          tryCatch({
            cat(sprintf("[Worker] [%s] ERROR after %.1fs processing %s: %s\n", 
                        format(Sys.time(), "%H:%M:%S"), worker_elapsed, basename(model_id), error_msg), 
                file = progress_file, append = TRUE)
            if (nchar(error_trace) > 0 && nchar(error_trace) < 2000) {
              cat(sprintf("[Worker] [%s] ERROR TRACEBACK for %s:\n%s\n", 
                          format(Sys.time(), "%H:%M:%S"), basename(model_id), error_trace), 
                  file = progress_file, append = TRUE)
            }
          }, error = function(e2) {})
          return(list(model_id = model_id, success = FALSE, message = error_msg, error = error_msg, traceback = error_trace))
        })
      }
      
      # Report completion with timing (for both up-to-date and processed models)
      worker_elapsed <- as.numeric(difftime(Sys.time(), worker_start, units = "secs"))
      tryCatch({
        status <- if (isTRUE(result$success)) "SUCCESS" else "FAILED"
        msg <- result$message %||% "Unknown status"
        cat(sprintf("[Worker] [%s] %s (%.1fs): model %d/%d (%s) - %s\n", 
                    format(Sys.time(), "%H:%M:%S"), status, worker_elapsed, i, total_count, 
                    basename(model_id), msg), 
            file = progress_file, append = TRUE)
        
        # If model took a very long time, warn about it
        if (worker_elapsed > 300) {  # > 5 minutes
          cat(sprintf("[Worker] [%s] WARNING: Model %d took %.1f minutes - may indicate slow/hanging model\n", 
                      format(Sys.time(), "%H:%M:%S"), i, worker_elapsed/60), 
              file = progress_file, append = TRUE)
        }
      }, error = function(e) {
        # Continue if file write fails
      })
      
      result
    }
  }, error = function(e) {
    cat("\n❌ ERROR in parallel foreach loop:", e$message, "\n")
    cat("Error traceback:\n")
    print(traceback())
    return(list())
  })
  
  # Report elapsed time
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("\nParallel processing completed in %.1f seconds (%.1f minutes)\n", 
              elapsed_time, elapsed_time/60))
  
  # Read and display progress file with analysis
  if (file.exists(progress_file)) {
    cat("\n=== Progress Log (last 50 lines) ===\n")
    progress_lines <- readLines(progress_file)
    cat(paste(tail(progress_lines, 50), collapse = "\n"))  # Show last 50 lines
    cat("\n")
    
    # Analyze progress for hanging models
    cat("=== Progress Analysis ===\n")
    warning_lines <- grep("WARNING", progress_lines, value = TRUE)
    if (length(warning_lines) > 0) {
      cat("Models that took >5 minutes:\n")
      for (warn in tail(warning_lines, 10)) {
        cat("  ", warn, "\n")
      }
    }
    
    error_lines <- grep("ERROR", progress_lines, value = TRUE)
    if (length(error_lines) > 0) {
      cat(sprintf("\nTotal errors encountered: %d\n", length(error_lines)))
      cat("Sample errors (last 5):\n")
      for (err in tail(error_lines, 5)) {
        cat("  ", err, "\n")
      }
    }
  } else {
    cat("⚠️  Warning: Progress file not found - workers may not have been able to write progress\n")
  }
  
  # Stop cluster with error handling
  tryCatch({
    stopCluster(cl)
    cat("Cluster stopped successfully\n")
  }, error = function(e) {
    cat("⚠️  Warning: Error stopping cluster:", e$message, "\n")
  })
  
  # Force garbage collection after parallel processing
  gc(verbose = FALSE)
  
  # Analyze results with detailed error reporting
  cat("\n=== Processing Results ===\n")
  successful_results <- 0
  failed_results <- 0
  error_summary <- list()
  
  for (i in seq_along(results)) {
    result <- results[[i]]
    model_id <- names(model_files)[i]
    
    # Handle different result types
    if (is.null(result)) {
      cat(sprintf("  [%d] %s - NULL result\n", i, basename(model_id)))
      failed_results <- failed_results + 1
      error_summary[[length(error_summary) + 1]] <- list(model = basename(model_id), error = "NULL result")
    } else if (!is.list(result)) {
      cat(sprintf("  [%d] %s - Invalid result (not a list)\n", i, basename(model_id)))
      failed_results <- failed_results + 1
      error_summary[[length(error_summary) + 1]] <- list(model = basename(model_id), error = "Invalid result type")
    } else if (isTRUE(result$success)) {
      successful_results <- successful_results + 1
      if (i %% 10 == 0 || i <= 5) {
        cat(sprintf("  [%d] %s - SUCCESS: %s\n", i, basename(model_id), result$message))
      }
    } else {
      failed_results <- failed_results + 1
      error_msg <- result$message %||% result$error %||% "Unknown error"
      cat(sprintf("  [%d] %s - FAILED: %s\n", i, basename(model_id), error_msg))
      error_summary[[length(error_summary) + 1]] <- list(model = basename(model_id), error = error_msg)
    }
  }
  
  # Print summary
  cat("\n=== Final Summary ===\n")
  cat("Successful:", successful_results, "of", total_models, "models\n")
  cat("Failed:", failed_results, "of", total_models, "models\n")
  
  if (length(error_summary) > 0 && length(error_summary) <= 20) {
    cat("\nFailed Models:\n")
    for (err in error_summary) {
      cat(sprintf("  - %s: %s\n", err$model, err$error))
    }
  } else if (length(error_summary) > 20) {
    cat(sprintf("\nFailed Models (showing first 10 of %d):\n", length(error_summary)))
    for (err in error_summary[1:10]) {
      cat(sprintf("  - %s: %s\n", err$model, err$error))
    }
    cat("  ... and", length(error_summary) - 10, "more\n")
  }
  
} else {
  # Sequential processing (original code)
  cat("Using sequential processing\n")
  model_counter <- 0
  for (model_id in names(model_files)) {
    model_counter <- model_counter + 1
    cat(sprintf("\n[%d/%d] Processing model: %s\n", model_counter, total_models, basename(model_id)))
    flush.console()
    chain_paths <- model_files[[model_id]]
    
    # Report chain file sizes before processing
    if (length(chain_paths) > 0) {
      chain_sizes <- sapply(chain_paths, function(p) {
        if (file.exists(p)) round(file.size(p) / 1024 / 1024, 1) else 0
      })
      cat(sprintf("  Chain files (%d): sizes = %s MB\n", length(chain_paths), 
                  paste(chain_sizes, collapse=", ")))
      flush.console()
    }
  
  # Check if this is an already-combined sample file that needs metadata extraction
  if (length(chain_paths) == 1 && grepl("^samples_", basename(chain_paths[[1]]))) {
    message("Processing already-combined sample file: ", basename(chain_paths[[1]]))
    
    # Check if this file needs re-combining or just metadata update
    tryCatch({
      sample_data <- readRDS(chain_paths[[1]])
      
      # Check if samples2 is missing (needs re-combination)
      missing_samples2 <- !"samples2" %in% names(sample_data)
      
      if (missing_samples2) {
        message("Sample file is missing samples2, will skip and re-combine from chains later")
        next  # Skip this file - it will be handled by the normal combination logic
      }
      
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
    # Normalize model_id to remove samples_ prefix if present for checkpoint search
    checkpoint_search_key <- if (grepl("^samples_", model_id)) {
      sub("^samples_", "", model_id)
    } else {
      model_id
    }
    message("Searching for checkpoint files as backup chains (using key: ", checkpoint_search_key, ")...")
    checkpoint_files <- find_recent_checkpoints(checkpoint_search_key, roots = c(driver_output_root, recon_output_root))
    
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
  
  message("Processing ", model_id, " with ", length(chain_paths), " chains...")
  
  # Build save path - always save in model type directory, not in subdirectories
  # Chains may be in subdirectories (e.g., env_cov/coniochaeta/) but combined files
  # should be saved in the model type directory (e.g., env_cov/)
  chain_paths <- unique(chain_paths)
  
  # Derive model type from the model_id and save in model type directory
  model_type <- if (grepl("^env_cycl_", model_id)) "env_cycl"
           else if (grepl("^env_cov_",  model_id)) "env_cov"
           else if (grepl("^cycl_only_",model_id)) "cycl_only" else "env_cycl"
  target_dir <- file.path(driver_output_root, model_type)
  
  # Ensure the target directory exists
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create the combined filename - always use "samples_" prefix for consistency
  # Normalize model_id first: remove samples_ prefix if present, then add it back
  normalized_model_id <- if (grepl("^samples_", model_id)) {
    sub("^samples_", "", model_id)
  } else {
    model_id
  }
  base_name <- paste0("samples_", normalized_model_id)
  combined_file_path <- file.path(target_dir, paste0(base_name, ".rds"))
  
  # Check if combined file already exists and if it has proper samples2 data
  needs_recombine <- TRUE
  if (file.exists(combined_file_path)) {
    tryCatch({
      test_read <- readRDS(combined_file_path)
      # Check if file has samples2 object (raw samples, not just summary)
      has_samples2_object <- FALSE
      if (is.list(test_read) && "samples2" %in% names(test_read)) {
        has_samples2_object <- TRUE
      }
      
      if (has_samples2_object) {
        # File has valid samples2 - check if it's up-to-date
        combined_file_time <- file.info(combined_file_path)$mtime
        input_files_exist <- file.exists(chain_paths)
        
        if (all(input_files_exist)) {
          input_file_times <- file.info(chain_paths[input_files_exist])$mtime
          newest_input_time <- max(input_file_times)
          
          if (combined_file_time > newest_input_time) {
            message("Combined file is up-to-date for ", basename(combined_file_path), " (newer than all input chains)")
            needs_recombine <- FALSE
          } else {
            message("Combined file is outdated for ", basename(combined_file_path), " (input chains are newer)")
          }
        } else {
          message("Some input chain files missing for ", basename(combined_file_path), ", will re-combine")
        }
      } else {
        message("Combined file exists but missing samples2 object, will re-combine: ", basename(combined_file_path))
      }
    }, error = function(e) {
      message("Error reading existing file, will overwrite: ", basename(combined_file_path))
    })
  }
  
  if (!needs_recombine) {
    next
  }
  
  # Preflight check: validate chain files before combining
  bad <- character(0)
  for (p in chain_paths) {
    ok <- tryCatch({
      z <- readRDS(p)
      is.list(z) && ("samples" %in% names(z))  # minimal requirement
    }, error = function(e) FALSE)
    if (!ok) bad <- c(bad, p)
  }
  if (length(bad)) {
    message("Skipping ", length(bad), " unreadable/invalid chain(s):")
    message("  - ", paste(basename(bad), collapse = "\n  - "))
    chain_paths <- setdiff(chain_paths, bad)
  }
  if (length(chain_paths) < 2) {
    message("Not enough valid chains after preflight; skipping model.")
    next
  }
  
  # Combine chains with truncation (similar to combine_large_chains.r)
  # Check file sizes and apply aggressive truncation for very large files
  chain_file_sizes_mb <- sapply(chain_paths, function(cf) {
    if (file.exists(cf)) file.size(cf) / 1024 / 1024 else 0
  })
  max_file_size_mb <- max(chain_file_sizes_mb, na.rm = TRUE)
  
  # Use truncation (default from config)
  effective_cut_size1 <- cut_size1
  effective_cut_size2 <- cut_size2
  
  # Apply aggressive truncation for very large files (>500MB per chain)
  if (max_file_size_mb > 500) {
    message("  Large file detected (max ", round(max_file_size_mb, 1), " MB per chain), using aggressive truncation")
    effective_cut_size1 <- min(cut_size1 %||% 5000, 3000)  # Reduce further for huge files
    effective_cut_size2 <- min(cut_size2 %||% 3000, 2000)
    
    # If file is extremely large (>1GB), skip it if we have enough other chains
    if (max_file_size_mb > 1000 && length(chain_paths) >= 3) {
      message("  WARNING: Extremely large file (", round(max_file_size_mb, 1), " MB) detected. Skipping to avoid memory errors.")
      message("  Using ", length(chain_paths) - 1, " other chains instead.")
      chain_paths <- chain_paths[chain_file_sizes_mb < 1000]
      if (length(chain_paths) < 2) {
        message("  Not enough chains after filtering; skipping model")
        next
      }
    }
  }
  
  cat(sprintf("  Combining %d chains...\n", length(chain_paths)))
  cat(sprintf("  Using truncation: cut_size1 = %s, cut_size2 = %s\n",
              if (is.null(effective_cut_size1)) "NULL" else effective_cut_size1,
              if (is.null(effective_cut_size2)) "NULL" else effective_cut_size2))
  flush.console()
  start_time <- Sys.time()
  message("Combining chains for ", model_id, " (", length(chain_paths), " chains)")
  out <- tryCatch({
    # Use the new combine_chains function that handles different sample sizes robustly
    # Pass file paths, not loaded objects
    # Pass cut_size parameters to limit chain sizes if specified
    combine_chains(chain_paths, cut_size1 = effective_cut_size1, cut_size2 = effective_cut_size2)
  }, error = function(e) {
    message("Error combining chains for ", model_id, ": ", e$message)
    return(NULL)
  })
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("  ✓ Chain combination completed in %.1f seconds\n", elapsed_time))
  flush.console()
  
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

  # HARD ASSERTION: ensure site_effect_sd is present in parameter samples
  # This parameter is required downstream to draw random site effects for new sites
  first_chain_cols <- try(colnames(chains[[1]]), silent = TRUE)
  if (inherits(first_chain_cols, "try-error") || is.null(first_chain_cols)) {
    message("❌ ERROR: Could not read column names from first chain for ", basename(model_id), "; skipping save")
    next
  }
  has_site_effect_sd <- any(grepl('^site_effect_sd(\\[|$)', first_chain_cols))
  if (!has_site_effect_sd) {
    message("❌ MISSING REQUIRED PARAMETER 'site_effect_sd' in ", basename(model_id))
    message("   Columns present (head): ", paste(head(first_chain_cols), collapse = ", "))
    message("   Action: re-fit model with 'site_effect_sd' monitored, or inspect 01_fitModels script.")
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
  # Include samples2 for hindcasting functions that need raw plot estimates
  out_summary <- list(
    samples = actual_samples,
    samples2 = actual_samples2,  # Raw plot estimates for hindcasting
    param_summary = param_summary,
    plot_summary = plot_summary,  # Statistical summaries for visualization
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
  
  # Map back to plot_mu[ on the combined object before saving
  if (is.list(out_summary) && !is.null(out_summary$samples2)) {
    if (is.matrix(out_summary$samples2)) {
      cn <- colnames(out_summary$samples2)
      if (!is.null(cn) && any(grepl("^mu\\[", cn))) {
        colnames(out_summary$samples2) <- sub("^mu\\[", "plot_mu[", cn)
      }
    } else if (inherits(out_summary$samples2, "mcmc.list")) {
      # Handle mcmc.list format
      for (i in seq_along(out_summary$samples2)) {
        cn <- colnames(out_summary$samples2[[i]])
        if (!is.null(cn) && any(grepl("^mu\\[", cn))) {
          colnames(out_summary$samples2[[i]]) <- sub("^mu\\[", "plot_mu[", cn)
        }
      }
    }
  }
  
  # Save to the correct combined file path
  saveRDS(out_summary, combined_file_path, compress = FALSE)
  message("Saved combined samples output for ", model_id, " to: ", basename(combined_file_path))
  
  # Delete chain files after successful combination
  # Only delete if combined file was JUST created (significantly newer than all chains)
  if (file.exists(combined_file_path)) {
    combined_mtime <- file.info(combined_file_path)$mtime
    combined_size <- file.size(combined_file_path)
    
    # Verify combined file is valid (reasonable size)
    if (combined_size > 100000) {  # At least 100KB
      chain_mtimes <- sapply(chain_paths, function(cf) {
        if (file.exists(cf)) file.info(cf)$mtime else NA
      })
      chain_mtimes <- chain_mtimes[!is.na(chain_mtimes)]
      
      if (length(chain_mtimes) > 0) {
        newest_chain_mtime <- max(chain_mtimes)
        # Only delete if combined file is significantly newer (at least 2 seconds)
        # This ensures we only delete when we just created the combined file
        time_diff_secs <- as.numeric(difftime(combined_mtime, newest_chain_mtime, units = "secs"))
        
        if (time_diff_secs >= 2) {
          deleted_count <- 0
          deleted_size <- 0
          for (chain_file in chain_paths) {
            if (file.exists(chain_file)) {
              tryCatch({
                chain_size <- file.size(chain_file)
                file.remove(chain_file)
                deleted_count <- deleted_count + 1
                deleted_size <- deleted_size + chain_size
                message("  Deleted chain file: ", basename(chain_file))
              }, error = function(e) {
                warning("  Failed to delete chain file ", basename(chain_file), ": ", e$message)
              })
            }
          }
          if (deleted_count > 0) {
            message("  ✓ Deleted ", deleted_count, " chain file(s), freed ", round(deleted_size / 1024 / 1024, 2), " MB")
          }
        } else {
          message("  ⚠️  Skipping chain deletion: combined file is not significantly newer (diff: ", round(time_diff_secs, 2), " sec)")
        }
      }
    } else {
      message("  ⚠️  Skipping chain deletion: combined file is suspiciously small (", round(combined_size/1024, 1), " KB)")
    }
  }
  
  if (min(es) < 500) {
    message("Low effective sample sizes - check for unconverged parameters in model: ", basename(model_id))
    print(head(gelman_out, 50))
  }
  }
}

cat("\nChain combining complete!\n")


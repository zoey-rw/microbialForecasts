#!/usr/bin/env Rscript

# RESTART CAPABILITY
# Can restart models using initial values from previous MCMC runs
# Automatically detects extreme values and applies fallback strategies

# Load required packages
library(microbialForecast)
    library(here)

# Set project root
here::i_am("analysis/model_analysis/01_restartFitModels.R")
project_root <- here()

cat("here() starts at", project_root, "\n")
cat("Project root set to:", getwd(), "\n")

# Load packages and create directories using package functions
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("logit_beta_fixed_priors", "logit_beta_fixed_priors/cycl_only")
)

cat("==================================================\n")
cat("Microbial Forecasts Environment Setup Complete!\n")
cat("Project root:", getwd(), "\n")
cat("Package status: microbialForecast loaded \n")
cat("Ready for analysis.\n")
cat("==================================================\n")

# Get arguments from the command line - handle both Rscript and R --slave -e
argv <- commandArgs(TRUE)
if (length(argv) == 0) {
    # Try alternative method for R --slave -e
    argv <- Sys.getenv("R_ARGS")
    if (argv != "") {
        argv <- strsplit(argv, " ")[[1]]
} else {
        argv <- character(0)
    }
}

if (length(argv) >= 2) {
    specific_model_idx <- as.numeric(argv[1])
    specific_chain_no <- as.numeric(argv[2])
    run_all_models <- FALSE
    cat("Running specific model index:", specific_model_idx, "chain:", specific_chain_no, "\n")
} else if (length(argv) == 1) {
    specific_model_idx <- as.numeric(argv[1])
    specific_chain_no <- 1  # Default to chain 1 if only model index provided
    run_all_models <- FALSE
    cat("Running specific model index:", specific_model_idx, "chain:", specific_chain_no, "\n")
} else {
    specific_model_idx <- 1
    specific_chain_no <- 1
    run_all_models <- TRUE
    cat("Running all models with all chains (default behavior)\n")
}

# Default restart settings
RESTART_ENABLED <- TRUE  # Re-enabled after confirming sequential functionality works
RESTART_MIN_ESS <- 10    # Minimum effective sample size for parameter extraction
RESTART_MAX_RHAT <- 1.1  # Maximum R-hat value allowed
RESTART_BURNIN_PROP <- 0.5  # Proportion of samples to discard as burnin
RESTART_FALLBACK_STRATEGY <- "median"  # "median", "random", or "zero"

# MCMC settings - configured for production convergence
burnin <- 200
thin <- 1
iter_per_chunk <- 200
init_iter <- 400
min_eff_size_perchain <- 10  # Target ESS per parameter
max_loops <- 10  # Production setting (10 * 200 = 2000 + 400 initial = 2400 iterations)
min_total_iterations <- 2000  # Production setting
max_save_size <- 20000  # Maximum samples to save

# Sequential execution settings
# NOTE: This should match batch_restart_models.R to avoid conflicts
# batch_restart_models.R uses nchains = 2, so we'll use the same here
nchains <- 4  # Number of chains per model (for HPC mode))

#### Run on all groups ----

source(here("source.R"))
# Restart functions are now defined inline for consistency with main script

# Load data early for filtering
cat("Loading data files for filtering...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully for", length(all_ranks), "ranks\n")

# MODEL IMPLEMENTATION
# This uses the beta regression framework with restart capability
#
# FILE PATH STRATEGY:
# - All file paths use here() function for consistency
# - here() resolves relative to the project root (microbialForecasts/)
# - No setwd() calls to avoid working directory confusion
# - All output files go to: here("data", "model_outputs", "logit_beta_regression", ...)

library(microbialForecast)
cat("✓ microbialForecast package loaded for consistent functions\n")

# Function to create restart initial values from previous model estimates
# This replaces the default initial values with actual parameter estimates (unless extreme)
create_restart_inits_from_estimates <- function(previous_samples, constants, model_name, use_fallback_for_extreme = TRUE) {
    cat("Creating restart initial values from previous model estimates...\n")
    
    # Get the final parameter estimates from previous samples
    if (is.list(previous_samples) && "samples" %in% names(previous_samples)) {
        samples_matrix <- previous_samples$samples
    } else {
        samples_matrix <- previous_samples
    }
    
    # Get the final row (latest estimates)
    final_estimates <- samples_matrix[nrow(samples_matrix), ]
    
    cat("  Previous model had", nrow(samples_matrix), "iterations\n")
    cat("  Using final parameter estimates for restart\n")
    
    # Check for extreme values and apply fallback if needed
    extreme_flags <- logical(length(final_estimates))
    names(extreme_flags) <- names(final_estimates)
    
    if (use_fallback_for_extreme) {
        # Define reasonable ranges for each parameter type
        for (param_name in names(final_estimates)) {
            param_value <- final_estimates[param_name]
            
            		if (grepl("precision", param_name)) {
			extreme_flags[param_name] <- param_value < 0.001 || param_value > 1000
		} else if (grepl("rho", param_name)) {
			extreme_flags[param_name] <- param_value < 0.001 || param_value > 0.999
		} else if (grepl("beta", param_name)) {
			extreme_flags[param_name] <- param_value < -50 || param_value > 50
		} else if (grepl("site_effect", param_name)) {
			extreme_flags[param_name] <- param_value < -10 || param_value > 10
		} else if (grepl("intercept", param_name)) {
			extreme_flags[param_name] <- param_value < -20 || param_value > 20
		} else if (grepl("legacy_effect", param_name)) {
			extreme_flags[param_name] <- param_value < -20 || param_value > 20
		} else if (grepl("site_effect_sd", param_name)) {
			extreme_flags[param_name] <- param_value < 0.001 || param_value > 10
		}
        }
        
        # Apply fallback for extreme values
        if (any(extreme_flags)) {
            cat("  ⚠️  Found", sum(extreme_flags), "extreme parameter values, applying fallback\n")
            
            			for (param_name in names(final_estimates)) {
				if (extreme_flags[param_name]) {
					old_value <- final_estimates[param_name]
					
					# Apply conservative fallback values
					if (grepl("precision", param_name)) {
						final_estimates[param_name] <- 2.0
					} else if (grepl("rho", param_name)) {
						final_estimates[param_name] <- 0.5
					} else if (grepl("beta", param_name)) {
						final_estimates[param_name] <- 0.0
					} else if (grepl("site_effect", param_name)) {
						final_estimates[param_name] <- 0.0
					} else if (grepl("intercept", param_name)) {
						final_estimates[param_name] <- 0.0
					} else if (grepl("legacy_effect", param_name)) {
						final_estimates[param_name] <- 0.0
					} else if (grepl("site_effect_sd", param_name)) {
						final_estimates[param_name] <- 0.1
					}
					
					cat("    ", param_name, ":", old_value, "→", final_estimates[param_name], "(fallback)\n")
				}
			}
        }
    }
    
    # Create the initial values list
    inits <- list()
    
    # Add all the estimated parameters
    for (param_name in names(final_estimates)) {
        if (!grepl("^Ex\\[|^mu\\[", param_name)) {  # Skip latent variables
            inits[[param_name]] <- final_estimates[param_name]
        }
    }
    
    # Add the latent variables using the stable approach (will be updated when constants are available)
    # For now, use placeholder dimensions that will be updated later
    inits$Ex <- matrix(0.3, nrow = 1, ncol = 1)  # Placeholder
    inits$mu <- matrix(0.3, nrow = 1, ncol = 1)  # Placeholder
    
    cat("  ✓ Created restart initial values from", length(final_estimates), "parameters\n")
    cat("  ✓ Fallback values applied:", sum(extreme_flags), "parameters\n")
    
    # Return the structure expected by the restart script
    return(list(
        initial_values = inits,
        extreme_flags = extreme_flags,
        fallback_used = extreme_flags  # All extreme values get fallback applied
    ))
}

# Function to update latent variable dimensions when constants become available
update_restart_inits_dimensions <- function(restart_inits, constants) {
    if (is.null(restart_inits) || is.null(constants)) {
        return(restart_inits)
    }
    
    cat("Updating restart initial values with proper dimensions...\n")
    
    # Update the latent variable matrices with correct dimensions
    restart_inits$initial_values$Ex <- matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date)
    restart_inits$initial_values$mu <- matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date)
    
    cat("  ✓ Updated Ex matrix dimensions:", dim(restart_inits$initial_values$Ex), "\n")
    cat("  ✓ Updated mu matrix dimensions:", dim(restart_inits$initial_values$mu), "\n")
    
    return(restart_inits)
}

# Wrapper function to handle parameter name mismatch
find_chain_files_wrapper <- function(model_name, species, min_date, max_date, use_legacy_covariate) {
  # Call the actual function with the correct parameter names
  return(find_chain_files(model_name, species, min_date, max_date, use_legacy_covariate))
}

# Enhanced check_continue function with restart capability
check_continue_with_restart <- function(samples, min_eff_size = 10, model_info = NULL) {
	# Validate input
	if (is.null(samples) || nrow(samples) == 0) {
		cat("  WARNING: Empty or NULL samples provided to check_continue, defaulting to continue\n")
		return(list(continue = TRUE, restart_recommended = FALSE, reason = "empty_samples"))
	}

	# Convert to mcmc object for effectiveSize calculation
	if (!inherits(samples, "mcmc")) {
		samples <- as.mcmc(samples)
	}

	# Calculate effective sample sizes for all parameters
	eff_sizes <- effectiveSize(samples)

	# Check if any parameter has ESS below threshold
	min_ess <- min(eff_sizes, na.rm = TRUE)

	# Continue if minimum ESS is below threshold
	continue <- min_ess < min_eff_size

	# Check for extreme values that might indicate restart is needed
	extreme_flags <- check_extreme_values(colMeans(samples))
	n_extreme <- sum(extreme_flags)

	# Recommend restart if we have extreme values or poor convergence
	restart_recommended <- n_extreme > 0 || min_ess < min_eff_size * 2

	cat("  ESS check - Min ESS:", round(min_ess, 1), "Target:", min_eff_size, "\n")
	cat("  Extreme values:", n_extreme, "/", length(extreme_flags), "\n")
	cat("  Continue sampling:", continue, "\n")
	cat("  Restart recommended:", restart_recommended, "\n")

	return(list(
		continue = continue,
		restart_recommended = restart_recommended,
		reason = if (n_extreme > 0) "extreme_values" else if (min_ess < min_eff_size * 2) "poor_convergence" else "normal"
	))
}

# Enhanced run_scenarios_fixed with restart capability
run_scenarios_with_restart <- function(j, chain_no) {
	# Initialize error tracking and logging
	start_time <- Sys.time()
	error_context <- list()

	tryCatch({
	cat("=== Starting model fitting (with restart capability) ===\n")
	cat("Model index:", j, "Chain:", chain_no, "\n")
	cat("Restart enabled:", RESTART_ENABLED, "\n")

		# Extract model parameters with validation
	rank.name <- params$rank.name[[j]]
	species <- params$species[[j]]
	model_id <- params$model_id[[j]]
	model_name <- params$model_name[[j]]
	min.date <- params$min.date[[j]]
	max.date <- params$max.date[[j]]
	scenario <- params$scenario[[j]]

	use_legacy_covariate <- grepl("Legacy with covariate", scenario)

	cat("Model configuration:\n")
	cat("  Species:", species, "(", rank.name, ")\n")
	cat("  Model:", model_name, "\n")
	cat("  Legacy covariate:", use_legacy_covariate, "\n")

	# Check if we should try to restart from previous chains
	restart_inits <- NULL
	if (RESTART_ENABLED) {
		cat("  Checking for existing chain files to restart from...\n")

		tryCatch({
			# Look for existing chain files for this model
			# Create wrapper to handle parameter name mismatch
			existing_chain_files <- find_chain_files_wrapper(model_name, species, min.date, max.date, use_legacy_covariate)

			if (length(existing_chain_files) > 0) {
				cat("  Found", length(existing_chain_files), "existing chain files\n")

				# Load the most recent chain file for restart
				latest_chain_file <- existing_chain_files[length(existing_chain_files)]
				cat("  Loading latest chain file:", basename(latest_chain_file), "\n")
				
				# Load the previous samples
				previous_samples <- readRDS(latest_chain_file)
				
				# Create restart initial values from previous estimates
				restart_inits <- create_restart_inits_from_estimates(
					previous_samples,
					constants = NULL,  # Will be set later when constants are available
					model_name = model_name,
					use_fallback_for_extreme = TRUE
				)

				cat("  ✓ Restart setup complete using previous parameter estimates\n")
				cat("    - Extreme values detected:", sum(restart_inits$extreme_flags), "/", length(restart_inits$extreme_flags), "\n")
				cat("    - Using fallback values:", sum(restart_inits$fallback_used), "\n")
				cat("    - Ready to restart with improved initial values\n")

			} else {
				cat("  No existing chain files found - starting fresh model\n")
			}
		}, error = function(e) {
			cat("  WARNING: Error setting up restart:", e$message, "\n")
			cat("  Falling back to fresh model run\n")
			restart_inits <- NULL
		})
	}

				# Validate data availability
		if (!exists("all_ranks") || is.null(all_ranks)) {
			stop("Data 'all_ranks' not available in environment")
		}

		# Get the specific group data with validation
		if (!(rank.name %in% names(all_ranks))) {
			stop("Rank name '", rank.name, "' not found in data. Available ranks: ",
				 paste(names(all_ranks), collapse=", "))
	}
	rank.df <- all_ranks[[rank.name]]

		# Check if species column exists
		if (!(species %in% colnames(rank.df))) {
			stop("Species '", species, "' not found in rank '", rank.name, "'. Available columns: ",
				 paste(colnames(rank.df), collapse=", "))
	}

	cat("Preparing model data for", rank.name, "\n")

	# Extract the specific species ONLY (no "other" column needed for beta regression)
	cat("DEBUG: Extracting species", species, "from", rank.name, "\n")

		# Validate required columns exist
		required_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", species)
		missing_cols <- required_cols[!required_cols %in% colnames(rank.df)]
		if (length(missing_cols) > 0) {
			stop("Missing required columns in rank data: ", paste(missing_cols, collapse=", "))
		}

	# For beta regression, we only need the species column (no "other" column)
	# The model will handle the constraint that proportions must sum to 1
	# Expected structure: 6 metadata columns + 1 species column = 7 total columns
	# This ensures prepBetaRegData creates a single-column y matrix
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!species)

	cat("DEBUG: rank.df_spec dimensions:", dim(rank.df_spec), "\n")
	cat("DEBUG: rank.df_spec columns:", colnames(rank.df_spec), "\n")
	cat("DEBUG: Expected columns: siteID, plotID, dateID, sampleID, dates, plot_date,", species, "\n")
	cat("DEBUG: Column count:", ncol(rank.df_spec), "(should be 7)\n")

	# Validate the data structure before calling prepBetaRegData
	if (ncol(rank.df_spec) != 7) {
		stop("Data structure error: expected 7 columns, got ", ncol(rank.df_spec), 
			 ". Columns: ", paste(colnames(rank.df_spec), collapse=", "))
	}

	# Extract and validate species data
		species_data <- rank.df_spec[[species]]
		if (all(is.na(species_data)) || all(species_data == 0) || all(species_data == 1)) {
			stop("Species '", species, "' has no valid variation (all NA, 0, or 1)")
		}

	# Note: For beta regression, we only need the species abundance column
	# The "other" column is not needed because:
	# 1. Beta regression models proportions directly (0-1 range)
	# 2. The model handles the constraint that proportions must sum to 1
	# 3. Including "other" would create a 2-column response, causing the error

	# Use prepBetaRegData with the species-specific data
		cat("Calling prepBetaRegData...\n")
	model.dat <- prepBetaRegData(rank.df = rank.df_spec,
								min.prev = 3,
								min.date = min.date,
								max.date = max.date)

	cat("Data prepared successfully\n")

	# Debug: Check data structure immediately after preparation
	cat("DEBUG: Data structure check:\n")
	cat("  model.dat$y dimensions:", dim(model.dat$y), "\n")
	cat("  model.dat$y class:", class(model.dat$y), "\n")
	cat("  model.dat$y column names:", if(is.matrix(model.dat$y)) paste(colnames(model.dat$y), collapse=", ") else "N/A", "\n")
	if (is.matrix(model.dat$y)) {
		cat("  model.dat$y first few rows:\n")
		print(head(model.dat$y, 3))
	}
	cat("  N.core calculated:", nrow(model.dat$y), "\n")
	cat("  N.spp calculated:", ncol(model.dat$y), "\n")
	cat("  Model type:", model_name, "\n")

	# CRITICAL: For beta regression, we expect exactly 1 column (the species abundance)
	if (ncol(model.dat$y) != 1) {
		cat("  ❌ ERROR: Beta regression requires exactly 1 column in response data\n")
		cat("  Current columns:", ncol(model.dat$y), "\n")
		cat("  Column names:", paste(colnames(model.dat$y), collapse=", "), "\n")
		cat("  This suggests the data preparation is including extra columns\n")
		stop("Response data must have exactly 1 column for beta regression")
	}

	# Validate model data structure
	if (!is.list(model.dat) || !("y" %in% names(model.dat))) {
		stop("Invalid model data structure from prepBetaRegData")
	}
	if (!is.matrix(model.dat$y) || nrow(model.dat$y) == 0) {
		stop("Model data 'y' is not a valid matrix or is empty")
	}

		# Prepare constants with validation
		cat("Preparing model constants...\n")
		required_constants <- c("plotID", "timepoint", "plot_site", "site_start", "plot_start",
							   "plot_index", "plot_num", "plot_site_num", "N.plot", "N.spp",
							   "N.core", "N.site", "N.date", "sin_mo", "cos_mo")

		missing_constants <- required_constants[!required_constants %in% names(model.dat)]
		if (length(missing_constants) > 0) {
			stop("Missing required constants in model data: ", paste(missing_constants, collapse=", "))
		}

		constants <- model.dat[required_constants]

	# Add environmental predictors with validation (ONLY for environmental models)
	if (model_name %in% c("env_cycl", "env_cov")) {
		env_predictors <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
		cat("Adding environmental predictors with validation...\n")

		for (pred in env_predictors) {
			if (pred %in% names(model.dat)) {
				constants[[pred]] <- model.dat[[pred]]

				# Validate predictor dimensions and structure
				pred_data <- model.dat[[pred]]
				if (is.matrix(pred_data)) {
					cat("  ✓ Added", pred, "predictor:", dim(pred_data), "matrix\n")
				} else if (is.vector(pred_data)) {
					cat("  ✓ Added", pred, "predictor:", length(pred_data), "vector\n")
				} else {
					cat("  ✓ Added", pred, "predictor:", class(pred_data), "object\n")
				}

				# Check for missing or extreme values
				if (is.numeric(pred_data)) {
					missing_pct <- mean(is.na(pred_data)) * 100
					if (missing_pct > 0) {
						cat("    WARNING:", pred, "has", round(missing_pct, 1), "% missing values\n")
					}

					if (is.matrix(pred_data)) {
						extreme_vals <- sum(abs(pred_data) > 10, na.rm = TRUE)
						if (extreme_vals > 0) {
							cat("    WARNING:", pred, "has", extreme_vals, "extreme values (>10)\n")
						}
					}
				}
			} else {
				cat("  ❌ ERROR:", pred, "predictor not found in model data\n")
				stop("Missing required environmental predictor: ", pred)
			}
		}
	} else {
		cat("Skipping environmental predictors for", model_name, "model\n")
		}

	# Add legacy covariate if needed
	if (use_legacy_covariate) {
		cat("Adding legacy covariate with enhanced validation...\n")

		# Create legacy covariate matrix properly for plot x time indexing (RESTORED from working model)
		# The legacy covariate should be 1 for legacy period (2013-2015), 0 for post-2015
		cat("Creating legacy covariate matrix for plot x time structure...\n")
		cat("Matrix dimensions needed:", constants$N.plot, "plots ×", constants$N.date, "time points\n")

		# Use time-based approach (simpler and more robust) - RESTORED from working model
		if ("timepoint" %in% names(model.dat)) {
			cat("Using timepoint for legacy calculation...\n")
			# Use timepoint-based legacy: assume early timepoints are legacy
			timepoints <- sort(unique(model.dat$timepoint))
			cat("Time points:", length(timepoints), "\n")
			cat("N.date constant:", constants$N.date, "\n")

			# Create legacy vector that matches N.date exactly
			if (length(timepoints) == constants$N.date) {
				legacy_by_time <- timepoints <= quantile(timepoints, 0.6)  # First 60% are legacy
			} else {
				# If dimensions don't match, create a vector of length N.date
				cat("WARNING: Timepoint length (", length(timepoints), ") != N.date (", constants$N.date, "), creating default pattern\n")
				legacy_by_time <- 1:constants$N.date <= quantile(1:constants$N.date, 0.6)
			}

			cat("Legacy time points:", sum(legacy_by_time), "\n")
			cat("Legacy vector length:", length(legacy_by_time), "\n")

			# Expand to full matrix (all plots have same legacy status per time point)
			constants$legacy <- matrix(as.numeric(legacy_by_time), nrow = constants$N.plot, ncol = constants$N.date, byrow = TRUE)

		} else {
			cat("WARNING: No time information available, creating default legacy pattern\n")
			# Default: first half of time period is legacy
			legacy_by_time <- rep(c(1, 0), length.out = constants$N.date)
			constants$legacy <- matrix(as.numeric(legacy_by_time), nrow = constants$N.plot, ncol = constants$N.date, byrow = TRUE)
		}

		cat("Legacy matrix created successfully\n")

		# Validate legacy covariate to prevent extreme values
		legacy_sum <- sum(constants$legacy)
		legacy_total <- length(constants$legacy)
		if (legacy_sum == 0 || legacy_sum == legacy_total) {
			cat("WARNING: Legacy covariate is all 0s or all 1s - this may cause numerical issues\n")
		}

		cat("Legacy covariate added:", legacy_sum, "legacy observations out of", legacy_total, "\n")
		cat("Legacy proportion:", round(legacy_sum/legacy_total, 3), "\n")
		cat("Legacy matrix dimensions:", nrow(constants$legacy), "x", ncol(constants$legacy), "\n")

		cat("  ✓ Legacy covariate added successfully\n")
	}
	
	# Update restart initial values with proper dimensions now that constants are available
	if (!is.null(restart_inits)) {
		restart_inits <- update_restart_inits_dimensions(restart_inits, constants)
	}

	# Model hyperparameters - adjust based on model type
		cat("Setting model hyperparameters...\n")
	if (model_name == "env_cycl") {
		constants$N.beta = 8
		constants$omega <- 0.1 * diag(8)
		constants$zeros <- rep(0, 8)
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
		constants$omega <- 0.1 * diag(6)
		constants$zeros <- rep(0, 6)
	} else {
		constants$N.beta = 2
		constants$omega <- 0.1 * diag(2)
		constants$zeros <- rep(0, 2)
	}



	cat("Constants prepared successfully\n")

	# Validate constants dimensions before creating initial values
	cat("Validating constants dimensions...\n")
	if (is.null(constants$N.plot) || constants$N.plot <= 0) {
		stop("Invalid N.plot constant:", constants$N.plot)
	}
	if (is.null(constants$N.date) || constants$N.date <= 0) {
		stop("Invalid N.date constant:", constants$N.date)
	}
	if (is.null(constants$N.site) || constants$N.site <= 0) {
		stop("Invalid N.site constant:", constants$N.site)
	}
	if (is.null(constants$N.core) || constants$N.core <= 0) {
		stop("Invalid N.core constant:", constants$N.core)
	}
	cat("  ✓ Constants validation passed\n")
	cat("    N.plot:", constants$N.plot, "N.date:", constants$N.date, "N.site:", constants$N.site, "N.core:", constants$N.core, "\n")



	# Create stable model using our function
	cat("Building Nimble model...\n")
	modelCode <- create_stable_model(model_name, use_legacy_covariate)

	# Create inits with restart capability
	cat("Creating initial values...\n")

	# Use restart initial values if available, otherwise create fresh inits
	if (!is.null(restart_inits)) {
		cat("  Using restart initial values from previous chains\n")

		# Get restart initial values
		inits <- restart_inits$initial_values
		
		# Debug: Check what we got from restart
		cat("  DEBUG: Restart inits structure:\n")
		cat("    Names:", paste(names(inits), collapse=", "), "\n")
		cat("    Lengths:", paste(sapply(inits, length), collapse=", "), "\n")
		
		# Convert individual parameter elements back to vectors for Nimble compatibility
		cat("  Converting restart inits to Nimble-compatible format...\n")
		
		# Convert beta parameters from individual elements to vector
		beta_names <- names(inits)[grepl("^beta\\[", names(inits))]
		if (length(beta_names) > 0) {
			cat("    Converting", length(beta_names), "beta parameters to vector\n")
			beta_values <- numeric(length(beta_names))
			for (i in 1:length(beta_names)) {
				param_name <- paste0("beta[", i, "]")
				if (param_name %in% names(inits)) {
					beta_values[i] <- inits[[param_name]]
				}
			}
			inits$beta <- beta_values
			# Remove individual beta elements
			inits <- inits[!grepl("^beta\\[", names(inits))]
			cat("    ✓ Beta vector created with", length(beta_values), "elements\n")
		}
		
		# Convert site_effect parameters from individual elements to vector
		site_effect_names <- names(inits)[grepl("^site_effect\\[", names(inits))]
		if (length(site_effect_names) > 0) {
			cat("    Converting", length(site_effect_names), "site_effect parameters to vector\n")
			site_effect_values <- numeric(length(site_effect_names))
			for (i in 1:length(site_effect_names)) {
				param_name <- paste0("site_effect[", i, "]")
				if (param_name %in% names(inits)) {
					site_effect_values[i] <- inits[[param_name]]
				}
			}
			inits$site_effect <- site_effect_values
			# Remove individual site_effect elements
			inits <- inits[!grepl("^site_effect\\[", names(inits))]
			cat("    ✓ Site_effect vector created with", length(site_effect_values), "elements\n")
		}
		
		cat("  ✓ Restart inits converted to Nimble format\n")
		
		# Debug: Show final restart inits structure
		cat("  FINAL restart inits structure:\n")
		cat("    Names:", paste(names(inits), collapse=", "), "\n")
		cat("    Lengths:", paste(sapply(inits, length), collapse=", "), "\n")
		if ("beta" %in% names(inits)) {
			cat("    Beta vector length:", length(inits$beta), "\n")
		}
		if ("site_effect" %in% names(inits)) {
			cat("    Site_effect vector length:", length(inits$site_effect), "\n")
		}

		# Ensure all required parameters are present
		set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model

		# Check for missing parameters, handling the case where restart inits have individual elements
		cat("  Checking for missing parameters in restart inits...\n")
		cat("  Available parameters:", paste(names(inits), collapse=", "), "\n")
		
		# Check if we have all required parameters (including individual beta and site_effect elements)
		has_precision <- "precision" %in% names(inits)
		has_rho <- "rho" %in% names(inits)
		has_beta <- any(grepl("^beta\\[", names(inits))) || "beta" %in% names(inits)
		has_intercept <- "intercept" %in% names(inits)
		has_site_effect <- any(grepl("^site_effect\\[", names(inits))) || "site_effect" %in% names(inits)
		has_site_effect_sd <- "site_effect_sd" %in% names(inits)
		has_legacy_effect <- if (use_legacy_covariate) "legacy_effect" %in% names(inits) else TRUE
		
		cat("  Parameter availability check:\n")
		cat("    precision:", has_precision, "\n")
		cat("    rho:", has_rho, "\n")
		cat("    beta:", has_beta, "\n")
		cat("    intercept:", has_intercept, "\n")
		cat("    site_effect:", has_site_effect, "\n")
		cat("    site_effect_sd:", has_site_effect_sd, "\n")
		cat("    legacy_effect:", has_legacy_effect, "\n")
		
		# Only add parameters that are truly missing
		if (!has_precision) {
			cat("  Adding missing parameter: precision\n")
			inits$precision <- rgamma(1, 2, 0.1)
		}
		if (!has_rho) {
			cat("  Adding missing parameter: rho\n")
			inits$rho <- runif(1, 0.1, 0.9)
		}
		if (!has_beta) {
			cat("  Adding missing parameter: beta\n")
			# Create beta vector with correct length for this model
			if (model_name == "env_cycl") {
				inits$beta <- rnorm(8, 0, 0.1)  # 8 beta parameters for env_cycl
			} else if (model_name == "env_cov") {
				inits$beta <- rnorm(6, 0, 0.1)  # 6 beta parameters for env_cov
			} else {
				inits$beta <- rnorm(2, 0, 0.1)  # 2 beta parameters for cycl_only
			}
		}
		if (!has_intercept) {
			cat("  Adding missing parameter: intercept\n")
			inits$intercept <- rnorm(1, 0, 0.5)
		}
		if (!has_site_effect) {
			cat("  Adding missing parameter: site_effect\n")
			inits$site_effect <- rnorm(constants$N.site, 0, 0.1)
		}
		if (!has_site_effect_sd) {
			cat("  Adding missing parameter: site_effect_sd\n")
			inits$site_effect_sd <- runif(1, 0.1, 1)
		}
		if (!has_legacy_effect && use_legacy_covariate) {
			cat("  Adding missing parameter: legacy_effect\n")
			inits$legacy_effect <- rnorm(1, 0, 0.1)
		}
		
		cat("  ✓ Parameter check complete\n")

		cat("  ✓ Restart initial values ready for chain", chain_no, "\n")
		cat("  RESTART MODE: Will use restart inits for MCMC\n")

	} else {
		cat("  Creating fresh initial values\n")
		cat("  FRESH MODE: Will create new inits for MCMC\n")

		# Create stable initial values using the working approach from sourced models
		cat("Creating initial values using working approach...\n")
		inits <- create_stable_inits(constants, model_name, model_data = model.dat)

		cat("  ✓ Fresh initial values created for chain", chain_no, "\n")
	}

	cat("Model built successfully\n")

	# Build model
		cat("Building Nimble model...\n")
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y = model.dat$y[,1,drop=FALSE]), inits = inits)

	# Debug: Check data dimensions
	cat("Data dimensions check:\n")
	cat("  model.dat$y dimensions:", dim(model.dat$y), "\n")
	cat("  constants$N.core:", constants$N.core, "\n")

	# Compile model
		cat("Compiling Nimble model...\n")
	cModel <- compileNimble(Rmodel)

	cat("Model compiled successfully\n")

	# Configure MCMC with proper sampler management
		cat("Configuring MCMC...\n")
	
	# Enhanced monitoring for comprehensive convergence analysis (matching working approach)
	monitored_params <- c(
		# Core parameters (primary focus) - MONITOR ALL ESSENTIAL PARAMETERS
		"precision", "rho", "beta", "site_effect_sd", "site_effect", "intercept", "legacy_effect"
	)

	# Use monitors2 for latent variables at different interval (matching working approach)
	monitored_latent_params <- c(
		# Monitor latent process variables at different interval for efficiency
		"Ex", "mu"
	)

	cat("Monitoring parameters for convergence analysis:\n")
	cat("  Core parameters:", paste(monitored_params, collapse=", "), "\n")
	cat("  Latent variables:", paste(monitored_latent_params, collapse=", "), "\n")
	cat("  Total beta parameters:", constants$N.beta, "\n")

	mcmcConf <- configureMCMC(
		model = cModel,
		monitors = monitored_params,
		monitors2 = monitored_latent_params,  # Use monitors2 for latent variables
		thin = 1,
		thin2 = 20,  # Sample latent variables every 20th iteration for efficiency
		enableWAIC = FALSE
	)

	# Add specialized samplers for better convergence of key parameters (matching working approach)
	cat("Adding specialized samplers for convergence improvement...\n")

	# 1. FIRST remove default samplers to prevent conflicts
	cat("  Removing default samplers...\n")
	mcmcConf$removeSamplers(c("precision", "rho", "beta", "site_effect_sd", "site_effect", "intercept"))

	if (use_legacy_covariate) {
		mcmcConf$removeSamplers("legacy_effect")
	}

	# 2. THEN add specialized samplers - IMPROVED SAMPLER STRATEGY
	cat("  Adding improved samplers for key parameters...\n")

	# Use regular slice samplers for single parameters, AF_slice for blocks
	mcmcConf$addSampler(target = "precision", type = "slice")        # Regular slice sampler for single parameter
	mcmcConf$addSampler(target = "rho", type = "slice")              # Regular slice sampler for single parameter
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")   # Regular slice sampler for single parameter
	mcmcConf$addSampler(target = "intercept", type = "slice")        # Regular slice sampler for single parameter

	if (use_legacy_covariate) {
		mcmcConf$addSampler(target = "legacy_effect", type = "slice")  # Regular slice sampler for single parameter
	}

	# Add specialized samplers for all beta parameters (seasonal + environmental)
		for (i in 1:constants$N.beta) {
			mcmcConf$addSampler(target = paste0("beta[", i, "]"), type = "slice")
			cat("    Added slice sampler for beta[", i, "]\n")
		}

	# 3. Better site effects sampling strategy (matching working approach)
	if (constants$N.site > 1) {
		cat("  Adding improved block sampler for site effects...\n")
		# Use adaptive block sampler for better mixing of correlated parameters
		mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, "]"), type = "AF_slice")
		cat("  Added adaptive block sampler for site_effect[1:", constants$N.site, "]\n")
	} else {
		# Single site effect - use adaptive slice sampler
		cat("  Adding adaptive slice sampler for single site effect...\n")
		mcmcConf$addSampler(target = "site_effect[1]", type = "AF_slice")
	}

	cat("MCMC configured successfully - ALL advanced features RESTORED!\n")

	# Build and compile MCMC
		cat("Building and compiling MCMC...\n")
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)

	cat("MCMC configured successfully\n")

	# Run MCMC with convergence-based sampling (same as main modeling script)
	cat("Running MCMC with convergence-based sampling for restart...\n")
	cat("Running MCMC with convergence-based sampling\n")
	cat("  Initial iterations:", init_iter, "burnin:", burnin, "\n")
	cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
	cat("  Target ESS per chain:", min_eff_size_perchain, "\n")
	cat("  Minimum total iterations:", min_total_iterations, "\n")
	
	cat("Running MCMC with convergence-based sampling for restart\n")
	cat("  Initial iterations:", init_iter, "burnin:", burnin, "\n")
	cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
	cat("  Target ESS per chain:", min_eff_size_perchain, "\n")
	cat("  Minimum total iterations:", min_total_iterations, "\n")

	if (!is.null(restart_inits)) {
		cat("  🔄 RESTART MODE: Using improved initial values from previous chains\n")
		cat("    - Extreme values detected:", sum(restart_inits$extreme_flags), "\n")
		cat("    - Fallback values used:", sum(restart_inits$fallback_used), "\n")
	} else {
		cat("  🆕 FRESH MODE: Starting with fresh initial values\n")
	}

	# Run initial iterations for restart
	cat("  Running initial iterations (", init_iter, " iterations) for restart...\n")
	compiled$run(niter = init_iter, thin = thin, nburnin = 0)
	cat("  Initial iterations completed\n")

	# Get initial samples and check convergence
	initial_samples <- as.matrix(compiled$mvSamples)
	cat("  Initial samples collected, checking convergence...\n")
	cat("  Initial samples dimensions:", dim(initial_samples), "\n")
	
			# Create output directory for checkpoints (same as main modeling script)
	cat("  Creating output directory for checkpoints...\n")
		model_output_dir <- file.path(here("data", "model_outputs", "logit_beta_fixed_priors", model_name))
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	    # Create model_id for consistent naming using package function
    model_id <- create_model_id(model_name, species, min.date, max.date, use_legacy_covariate, "beta_regression")
	
	# Check if we need to continue sampling for convergence
	continue <- TRUE
	loop_counter <- 0
	total_iterations <- init_iter
	
	# Try to check convergence, with fallback if it fails
	tryCatch({
		continue <- check_continue_with_restart(initial_samples, min_eff_size = min_eff_size_perchain)
		continue <- continue$continue  # Extract the continue flag
	}, error = function(e) {
		cat("  WARNING: Convergence check failed, defaulting to continue sampling\n")
		cat("  Error:", e$message, "\n")
		continue <- TRUE
	})
	
	# Store all samples as we go
	all_samples <- initial_samples
	cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
	
	    # Save initial samples as checkpoint using package function
    save_checkpoint_safe(all_samples, total_iterations, 0, model_output_dir, model_id, chain_no, "initial", "beta_regression")
	
	    # Also save a simple progress file using package function
    progress_file <- create_progress_file(model_output_dir, model_id, chain_no, init_iter, "beta_regression")
	
	while ((continue || total_iterations < min_total_iterations) && loop_counter < max_loops) {
		if (continue) {
			cat("  Effective sample size too low; running for another", iter_per_chunk, "iterations\n")
		} else {
			cat("  Minimum iterations not reached; running for another", iter_per_chunk, "iterations\n")
		}
		cat("  Loop", loop_counter + 1, "of", max_loops, "\n")
		
		# Continue sampling without resetting
		compiled$run(niter = iter_per_chunk, thin = thin, nburnin = 0)
		total_iterations <- total_iterations + iter_per_chunk
		
		# Get updated samples and accumulate them
		# CRITICAL: NIMBLE's mvSamples resets between runs, so current_samples contains ONLY the latest iteration
		current_samples <- as.matrix(compiled$mvSamples)
		current_sample_count <- nrow(current_samples)
		previous_sample_count <- nrow(all_samples)
		
		cat("  Current iteration samples:", current_sample_count, "\n")
		cat("  Previous accumulated samples:", previous_sample_count, "\n")
		
		# CORRECTED LOGIC: Always append all current samples since mvSamples resets
		if (current_sample_count > 0) {
			all_samples <- rbind(all_samples, current_samples)
			cat("  ✓ Added", current_sample_count, "new samples, total accumulated:", nrow(all_samples), "\n")
		} else {
			cat("  WARNING: No samples in current iteration\n")
		}
		
		        # Save checkpoint after each loop using package function
        save_checkpoint_safe(all_samples, total_iterations, loop_counter + 1, model_output_dir, model_id, chain_no, paste0("loop", loop_counter + 1), "beta_regression")
		
		        # Update progress file using package function
        update_progress_file(progress_file, total_iterations, loop_counter + 1)
		
		# Check if we need to continue
		continue <- TRUE
		tryCatch({
			continue_result <- check_continue_with_restart(all_samples, min_eff_size = min_eff_size_perchain)
			continue <- continue_result$continue  # Extract the continue flag
		}, error = function(e) {
			cat("  WARNING: Convergence check failed in loop, defaulting to continue sampling\n")
			cat("  Error:", e$message, "\n")
			continue <- TRUE
		})
		loop_counter <- loop_counter + 1
		
		cat("  Total iterations so far:", total_iterations, "\n")
		cat("  Convergence check result:", ifelse(continue, "CONTINUE", "CONVERGED"), "\n")
		cat("  Current accumulated sample size:", nrow(all_samples), "\n")
		cat("  Progress: ", round(loop_counter/max_loops * 100, 1), "% of max loops completed\n")
	}
	
	if (loop_counter >= max_loops) {
		cat("  WARNING: Exceeded maximum loops (", max_loops, "). Stopping sampling.\n")
	} else if (total_iterations >= min_total_iterations) {
		if (continue) {
			cat("  WARNING: Minimum iterations reached but convergence not achieved\n")
		} else {
			cat("  SUCCESS: Convergence reached after", total_iterations, "total iterations\n")
		}
	} else {
		cat("  WARNING: Stopped before minimum iterations due to max loops\n")
	}
	
	    # Update final progress status using package function
		final_status <- if(loop_counter >= max_loops) "Completed (max loops)" else 
			if(total_iterations >= min_total_iterations && !continue) "Converged" else 
				"Completed (min iterations)"
    update_progress_file(progress_file, total_iterations, loop_counter)
		cat("  ✓ Final progress status updated\n")
	
	# Get final samples (use accumulated samples)
	samples <- all_samples

	# Get final samples (samples already collected above)
	# samples variable is already defined from compiled$mvSamples

	cat("MCMC completed successfully\n")
	cat("Final sample dimensions:", dim(samples), "\n")
	cat("Total iterations run:", total_iterations, "\n")
	cat("Restart iterations:", total_iterations, "\n")

	if (!is.null(restart_inits)) {
		cat("🔄 RESTART MODE RESULTS:\n")
		cat("  - Used restart initial values from previous chains\n")
		cat("  - Extreme values detected:", sum(restart_inits$extreme_flags), "\n")
		cat("  - Fallback values used:", sum(restart_inits$fallback_used), "\n")
	} else {
		cat("🆕 FRESH MODE RESULTS:\n")
		cat("  - Started with fresh initial values\n")
	}

	# Final convergence check
	tryCatch({
		final_ess <- effectiveSize(as.mcmc(samples))
		min_final_ess <- min(final_ess, na.rm = TRUE)
		cat("  Final minimum ESS:", round(min_final_ess, 1), "\n")
		cat("  Convergence achieved:", min_final_ess >= 10, "\n")
	}, error = function(e) {
		cat("  Final ESS check failed:", e$message, "\n")
	})

	cat("=== RESTART SAVING SUMMARY ===\n")
	cat("  Initial iterations:", init_iter, "iterations\n")
	cat("  Additional loops:", loop_counter, "iterations\n")
	cat("  Total accumulated samples:", nrow(samples), "iterations\n")
	cat("  This is a restart run with convergence-based sampling\n")

		# Create output directory for this model (same as main modeling script)
	model_output_dir <- file.path(here("data", "model_outputs", "logit_beta_fixed_priors", model_name))
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Save MCMC samples with absolute path
	samples_file <- file.path(model_output_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))

	# Create the complete chain structure with metadata
	chain_output <- list(
		samples = samples,
		metadata = list(
			rank.name = rank.name,
			species = species,
			model_name = model_name,
			model_id = model_id,
			use_legacy_covariate = use_legacy_covariate,
			scenario = scenario,
			min.date = min.date,
			max.date = max.date,
			niter = total_iterations,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			nimble_code = modelCode,
			model_structure = "stable_beta_regression_with_restart_capability",
			restart_info = if (!is.null(restart_inits)) list(
				restart_used = TRUE,
				extreme_values_detected = sum(restart_inits$extreme_flags),
				fallback_values_used = sum(restart_inits$fallback_used),
				fallback_strategy = RESTART_FALLBACK_STRATEGY,
				restart_iterations = total_iterations
			) else list(restart_used = FALSE, restart_iterations = total_iterations)
		)
	)

		# Save with error handling
		tryCatch({
	saveRDS(chain_output, samples_file)
			cat("✓ SUCCESS: Saved MCMC samples to:", samples_file, "\n")
		}, error = function(e) {
			cat("✗ ERROR: Failed to save samples to", samples_file, "\n")
			cat("  Error:", e$message, "\n")
			# Try to save to current directory as fallback
			fallback_file <- paste0("samples_", model_id, "_chain", chain_no, "_FALLBACK.rds")
			tryCatch({
				saveRDS(chain_output, fallback_file)
				cat("✓ FALLBACK: Saved to current directory:", fallback_file, "\n")
			}, error = function(e2) {
				cat("✗ CRITICAL: Failed to save even to fallback location\n")
				cat("  Fallback error:", e2$message, "\n")
			})
		})

	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== Model fitting completed ===\n")

	if (!is.null(restart_inits)) {
		cat("🔄 RESTART MODE: Successfully used improved initial values\n")
	} else {
		cat("🆕 FRESH MODE: Standard model fitting completed\n")
	}

	cat("  - STABLE: Beta regression with logit bounds and precision parameter\n")
	cat("  - All three model types supported: cycl_only, env_only, env_cov\n")
	cat("  - CONVERGENCE-BASED: Adaptive sampling until reasonable ESS reached\n")
	cat("  - ITERATIVE SAVING: Samples accumulated and saved incrementally\n")
	cat("  - RESTART CAPABLE: Can use initial values from previous chains\n")

	return(list(
		status = "SUCCESS",
		samples = samples,
		file = samples_file,
		model_data = model.dat,
		nimble_code = modelCode,
		restart_info = if (!is.null(restart_inits)) restart_inits else NULL,
		metadata = list(
			rank.name = rank.name,
			species = species,
			model_name = model_name,
			model_id = model_id,
			use_legacy_covariate = use_legacy_covariate,
			has_driver_uncertainty = FALSE,  # Flag to identify regular models (no driver uncertainty)
			scenario = scenario,
			min.date = min.date,
			max.date = max.date,
			niter = total_iterations,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			nimble_code = modelCode,
			model_structure = "standardized_beta_regression_with_consistent_priors"
		)
	))

	}, error = function(e) {
		# Enhanced error handling with restart context
		error_time <- Sys.time()
		error_context <- list(
			timestamp = error_time,
			task_idx = j,
			chain_no = chain_no,
			error_message = if(!is.null(e$message) && e$message != "") e$message else "No error message available",
			error_call = if(!is.null(e$call)) paste(deparse(e$call), collapse=" ") else "No call information",
			error_class = class(e)[1],
			system_info = list(
				r_version = R.version.string,
				working_dir = getwd(),
				available_packages = installed.packages()[,"Package"],
				memory_usage = if(exists("gc")) gc() else "GC not available"
			),
			runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA,
			restart_context = list(
				restart_enabled = RESTART_ENABLED,
				restart_inits_available = !is.null(restart_inits),
				extreme_values_detected = if(!is.null(restart_inits)) sum(restart_inits$extreme_flags) else 0,
				fallback_values_used = if(!is.null(restart_inits)) sum(restart_inits$fallback_used) else 0
			)
		)

		# Create detailed error file with absolute path
		model_name <- if(exists("model_name")) model_name else "unknown"
		error_dir <- file.path(getwd(), "data", "model_outputs", "logit_beta_regression", model_name)
		dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
		error_file <- file.path(error_dir, paste0("chain_", j, "_", chain_no, "_ERROR.txt"))

		# Write comprehensive error report
		error_report <- c(
			paste("ERROR DETAILED REPORT -", format(error_time)),
			paste("Task Index:", error_context$task_idx),
			paste("Model Index:", j),
			paste("Chain Number:", error_context$chain_no),
			paste("Error Message:", error_context$error_message),
			paste("Error Call:", error_context$error_call),
			paste("Error Class:", error_context$error_class),
			paste("Runtime (seconds):", round(error_context$runtime, 2)),
			paste("R Version:", error_context$system_info$r_version),
			paste("Working Directory:", error_context$system_info$working_dir),
			paste("Available Packages:", paste(error_context$system_info$available_packages, collapse=", ")),
			paste("Restart Enabled:", error_context$restart_context$restart_enabled),
			paste("Restart Inits Available:", error_context$restart_context$restart_inits_available),
			paste("Extreme Values Detected:", error_context$restart_context$extreme_values_detected),
			paste("Fallback Values Used:", error_context$restart_context$fallback_values_used),
			"",
			"FULL ERROR OBJECT:",
			capture.output(str(e))
		)

		# Save error report with error handling
		tryCatch({
			writeLines(error_report, error_file)
			cat("✓ ERROR REPORT: Saved detailed error to:", error_file, "\n")
		}, error = function(e) {
			cat("✗ ERROR: Failed to save error report to", error_file, "\n")
			cat("  Error:", e$message, "\n")
			# Try to save to current directory as fallback
			fallback_error_file <- paste0("chain_", j, "_", chain_no, "_ERROR_FALLBACK.txt")
			tryCatch({
				writeLines(error_report, fallback_error_file)
				cat("✓ FALLBACK: Saved error report to current directory:", fallback_error_file, "\n")
			}, error = function(e2) {
				cat("✗ CRITICAL: Failed to save error report even to fallback location\n")
				cat("  Fallback error:", e2$message, "\n")
			})
		})

		# Also log to console with detailed information
		cat("ERROR in Model", j, "Chain", error_context$chain_no, ":\n")
		cat("  Message:", error_context$error_message, "\n")
		cat("  Call:", error_context$error_call, "\n")
		cat("  Class:", error_context$error_class, "\n")
		cat("  Runtime:", round(error_context$runtime, 2), "seconds\n")
		cat("  Restart enabled:", error_context$restart_context$restart_enabled, "\n")
		if (error_context$restart_context$restart_inits_available) {
			cat("  Extreme values detected:", error_context$restart_context$extreme_values_detected, "\n")
			cat("  Fallback values used:", error_context$restart_context$fallback_values_used, "\n")
		}
		cat("  Detailed error saved to:", error_file, "\n")

		# Return detailed error information
		return(list(
			status = "ERROR",
			error = error_context$error_message,
			error_details = error_context,
			error_file = error_file
		))
	})
}

# [The rest of the script remains the same as the original 01_fitModels.R]
# Load parameters, set up parallel execution, etc.

# Test with single model for restart functionality
params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 rep("logical", 2),
										 rep("character", 4)))

rerun_list = readRDS(here("data/summary/unconverged_taxa_list.rds"))
converged_list = readRDS(here("data/summary/converged_taxa_list.rds"))

# RESTART CONFIGURATION: Use full unconverged list for production runs
cat("DEBUG: Starting filtering process...\n")
cat("DEBUG: rerun_list length:", length(rerun_list), "\n")
cat("DEBUG: First few rerun_list items:", head(rerun_list, 3), "\n")
cat("DEBUG: converged_list length:", length(converged_list), "\n")

# Create full model names that match the rerun_list format
params_with_full_names <- params_in %>% ungroup %>% 
	filter(scenario %in% c("Legacy with covariate 2013-2018", "2013-06-01_2018-01-01")) %>%
	mutate(full_model_name = ifelse(
		grepl("_with_legacy_covariate$", model_id),
		model_id,  # Already has the suffix
		paste0(model_id, "_with_legacy_covariate")  # Add the suffix
	))

cat("DEBUG: After scenario filtering, params has", nrow(params_with_full_names), "rows\n")
cat("DEBUG: Sample full_model_names:", head(params_with_full_names$full_model_name, 3), "\n")

# Filter for unconverged models
params_filtered <- params_with_full_names %>% filter(full_model_name %in% rerun_list)
cat("DEBUG: After rerun_list filtering, params has", nrow(params_filtered), "rows\n")

# Remove duplicate model entries
params_filtered <- params_filtered %>% distinct(full_model_name, .keep_all = TRUE)
cat("DEBUG: After deduplication, params has", nrow(params_filtered), "rows\n")

# Filter out already converged models (using full model names)
params_filtered <- params_filtered %>% filter(!full_model_name %in% converged_list)
cat("DEBUG: After converged_list filtering, params has", nrow(params_filtered), "rows\n")

# CRITICAL FIX: Order params to match the shell script order (rerun_list order)
# The shell script uses the order from rerun_list, so we must match that exactly
params_filtered$rerun_order <- match(params_filtered$full_model_name, rerun_list)
params <- params_filtered %>% arrange(rerun_order) %>% select(-rerun_order)

cat("DEBUG: Final params order (matching shell script):\n")
for (i in 1:min(5, nrow(params))) {
  cat("  ", i, ":", params$full_model_name[i], "(", params$species[i], ")\n")
}

# Set valid_models for compatibility with runAndSave_task function
valid_models <- params

cat("RESTART CONFIGURATION: Running", nrow(params), "models with restart capability\n")

cat("RESTART CONFIGURATION: Running", nrow(params), "models with restart capability\n")
cat("Restart functionality:", if(RESTART_ENABLED) "ENABLED" else "DISABLED", "\n")
cat("Model configuration:\n")
if (nrow(params) > 0) {
	print(params[, c("rank.name", "species", "model_name", "model_id")])
} else {
	cat("No models to run\n")
}

# Pre-validate data availability for all models to catch errors early
cat("Pre-validating data availability for all models...\n")

# Check if we have any models to validate
if (nrow(params) == 0) {
	cat("❌ ERROR: No models to validate after filtering!\n")
	cat("This means either:\n")
	cat("  1. No models in rerun_list match the scenario filter\n")
	cat("  2. All models in rerun_list are already in converged_list\n")
	cat("  3. There's a mismatch in model naming between params and rerun_list\n")
	
	# Show what we have in the lists for debugging
	cat("\nDebugging information:\n")
	cat("rerun_list length:", length(rerun_list), "\n")
	if (length(rerun_list) > 0) {
		cat("First few rerun_list items:", paste(head(rerun_list, 5), collapse=", "), "\n")
	}
	cat("converged_list length:", length(converged_list), "\n")
	if (length(converged_list) > 0) {
		cat("First few converged_list items:", paste(head(converged_list, 5), collapse=", "), "\n")
	}
	
	# Show sample from params_in for comparison
	if (nrow(params_in) > 0) {
		cat("Sample model_ids from params_in:", paste(head(params_in$model_id, 5), collapse=", "), "\n")
		cat("Available scenarios:", paste(unique(params_in$scenario), collapse=", "), "\n")
	}
	
	stop("No models available for processing after filtering")
}

validation_errors <- list()

for (i in 1:nrow(params)) {
	rank_name <- params$rank.name[i]
	species_name <- params$species[i]

	# Check if rank exists in data
	if (!(rank_name %in% names(all_ranks))) {
		validation_errors <- c(validation_errors,
			paste("Model", i, ": Rank '", rank_name, "' not found in data. Available ranks:",
				  paste(names(all_ranks), collapse=", ")))
		next
	}

	rank_data <- all_ranks[[rank_name]]

	# Check if species exists in rank data
	if (!(species_name %in% colnames(rank_data))) {
		validation_errors <- c(validation_errors,
			paste("Model", i, ": Species '", species_name, "' not found in rank '", rank_name, "'"))
		next
	}

	# Check if species has valid data
	species_data <- rank_data[[species_name]]
	if (all(is.na(species_data)) || all(species_data == 0) || all(species_data == 1)) {
		validation_errors <- c(validation_errors,
			paste("Model", i, ": Species '", species_name, "' has no valid variation (all NA, 0, or 1)"))
	}
}

# Report validation results
if (length(validation_errors) > 0) {
	cat("⚠️  Validation errors found:\n")
	for (error in validation_errors) {
		cat("  ", error, "\n")
	}
	cat("  These models will be filtered out before processing\n")
} else {
	cat("✓ All models passed validation\n")
}

# Create a function to check if a model has valid data
is_valid_model <- function(rank_name, species_name) {
	# Check if rank exists in data
	if (!(rank_name %in% names(all_ranks))) {
		return(FALSE)
	}

	rank_data <- all_ranks[[rank_name]]

	# Check if species exists in rank data
	if (!(species_name %in% colnames(rank_data))) {
		return(FALSE)
	}

	# Check if species has valid data
	species_data <- rank_data[[species_name]]
	if (all(is.na(species_data)) || all(species_data == 0) || all(species_data == 1)) {
		return(FALSE)
	}

	return(TRUE)
}

# Filter the parameters dataframe
valid_indices <- sapply(1:nrow(params), function(i) {
	is_valid_model(params$rank.name[i], params$species[i])
})

params <- params[valid_indices, ]
filtered_n_models <- nrow(params)

# Report filtering results
if (filtered_n_models < nrow(params_in)) {
	n_filtered <- nrow(params_in) - filtered_n_models
	cat("⚠️  Filtered out", n_filtered, "models with unavailable species/ranks\n")
	cat("  Original models:", nrow(params_in), "\n")
	cat("  Valid models:", filtered_n_models, "\n")
} else {
	cat("✓ All", filtered_n_models, "models have valid data\n")
}

# Additional validation: ensure we have at least one valid model
if (filtered_n_models == 0) {
	cat("❌ ERROR: No valid models remaining after filtering!\n")
	cat("Available ranks in data:", paste(names(all_ranks), collapse=", "), "\n")

	# Show some examples of available species for each rank
	for (rank_name in names(all_ranks)) {
		rank_data <- all_ranks[[rank_name]]
		metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
		available_species <- setdiff(colnames(rank_data), metadata_cols)
		cat("  ", rank_name, "species (first 5):", paste(head(available_species, 5), collapse=", "),
		    if(length(available_species) > 5) paste("... (+", length(available_species) - 5, " more)") else "", "\n")
	}

	stop("No valid models to run - check species names and rank names in parameters")
}

# SEQUENTIAL EXECUTION - Respect command line arguments
if (run_all_models) {
cat("Running restart script sequentially for", filtered_n_models, "model ×", nchains, "chains\n")

    # Create a simple task list for sequential execution
all_tasks <- expand.grid(model_idx = 1:filtered_n_models, chain_no = 1:nchains)
    cat("Total sequential tasks:", nrow(all_tasks), "(", filtered_n_models, "models ×", nchains, "chains)\n")
cat("Task details:\n")
print(all_tasks)
    
    # Run everything sequentially with restart capability
    cat("Starting sequential execution with restart capability at:", format(Sys.time()), "\n")
    cat("RESTART: Starting sequential execution for", filtered_n_models, "models with", nchains, "chains\n")
    cat("Expected runtime: Variable (convergence-based sampling)\n")
    cat("  - Models to restart:", filtered_n_models, "models\n")
    cat("  - Chains per model:", nchains, "(total", filtered_n_models * nchains, "sequential tasks)\n")
    cat("  - Target: ESS >= 10 per parameter\n")
    cat("  - Restart functionality:", if(RESTART_ENABLED) "ENABLED" else "DISABLED", "\n")
    start_time <- Sys.time()
    
    # Run all tasks sequentially
    for (task_idx in 1:nrow(all_tasks)) {
        run_sequential_task(task_idx)
    }
} else {
    # Run specific model and chain
    cat("Running specific model index:", specific_model_idx, "chain:", specific_chain_no, "\n")
    
    # Validate model index
    if (specific_model_idx > filtered_n_models) {
        stop("Model index ", specific_model_idx, " exceeds available models (", filtered_n_models, ")")
    }
    
    # Validate chain number
    if (specific_chain_no < 1 || specific_chain_no > nchains) {
        stop("Chain number ", specific_chain_no, " must be between 1 and ", nchains)
    }
    
    cat("Starting specific execution at:", format(Sys.time()), "\n")
    cat("  - Model index:", specific_model_idx, "\n")
    cat("  - Chain number:", specific_chain_no, "\n")
cat("  - Target: ESS >= 10 per parameter\n")
cat("  - Restart functionality:", if(RESTART_ENABLED) "ENABLED" else "DISABLED", "\n")
start_time <- Sys.time()

    # Run the specific task
    result <- run_scenarios_with_restart(j = specific_model_idx, chain_no = specific_chain_no)
    
    # Handle result
    if (result$status == "SUCCESS") {
        cat("✓ SUCCESS: Model", specific_model_idx, "Chain", specific_chain_no, "completed successfully\n")
    } else {
        cat("❌ FAILED: Model", specific_model_idx, "Chain", specific_chain_no, "failed with status:", result$status, "\n")
    }
}

# Sequential execution function that preserves all functionality from 01_fitModels_betaReg.R
run_sequential_task <- function(task_idx) {
  # Get task details first
    task <- all_tasks[task_idx, ]
    model_idx <- task$model_idx
    chain_no <- task$chain_no
    
  cat("=== Starting Model", model_idx, "Chain", chain_no, "===\n")
  cat("Time:", format(Sys.time()), "\n")
  
  # Initialize error tracking
  error_details <- list()
  start_time <- Sys.time()
  
  tryCatch({
    cat("DEBUG: Inside tryCatch block\n")
    
    cat("Task details - model_idx:", model_idx, "chain_no:", chain_no, "\n")
    
    # Validate model index
    if (model_idx > nrow(params)) {
      stop("Model index ", model_idx, " exceeds available models (", nrow(params), ")")
    }
    
    cat("All checks passed, calling run_scenarios_with_restart...\n")
    
    # Run the model with detailed error context
    result <- run_scenarios_with_restart(j = model_idx, chain_no = chain_no)
    
    # Validate result structure
    if (!is.list(result) || !("status" %in% names(result))) {
      stop("Invalid result structure from run_scenarios_with_restart")
    }
    
    # Save result immediately if successful
    if (result$status == "SUCCESS") {
      # Create output directory early with HPC compatibility
      possible_bases <- c(
        here("data", "model_outputs"),
        file.path(getwd(), "data", "model_outputs"),
        file.path(Sys.getenv("HOME"), "data", "model_outputs"),
        file.path(Sys.getenv("PWD"), "data", "model_outputs"),
        "./data/model_outputs"
      )
      
      model_output_dir <- NULL
      for (base_dir in possible_bases) {
        if (!is.null(base_dir) && base_dir != "" && base_dir != "NULL") {
          test_dir <- file.path(base_dir, "logit_beta_regression", params$model_name[model_idx])
          tryCatch({
            dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
            if (dir.exists(test_dir)) {
              model_output_dir <- test_dir
              cat("  Created output directory:", model_output_dir, "\n")
              break
            }
          }, error = function(e) {
            cat("  Failed to create directory with base:", base_dir, "-", e$message, "\n")
          })
        }
      }
      
      if (is.null(model_output_dir)) {
        # Fallback: create in current directory
        model_output_dir <- file.path("data", "model_outputs", "logit_beta_regression", params$model_name[model_idx])
        dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
        cat("  WARNING: Using fallback output directory:", model_output_dir, "\n")
      }
      
      # Create model_id for consistent naming using package function
      use_legacy_covariate <- grepl("Legacy with covariate", params$scenario[model_idx])
      model_id <- create_model_id(params$model_name[model_idx], params$species[model_idx], 
                       params$min.date[model_idx], params$max.date[model_idx], 
                                 use_legacy_covariate, "beta_regression")
      
      # Save MCMC samples immediately
      samples_file <- file.path(model_output_dir, 
                               paste0("samples_", model_id, "_chain", chain_no, ".rds"))
      
      # Create the complete chain structure with metadata
      # Use the metadata from the result if available, otherwise create it
      if ("metadata" %in% names(result) && !is.null(result$metadata)) {
        # Use the complete metadata from the result
        metadata <- result$metadata
        # Add parallel execution specific fields
        metadata$task_idx <- task_idx
        metadata$completed_at <- Sys.time()
      } else {
        # Fallback: create metadata if not available in result
        metadata <- list(
          rank.name = params$rank.name[model_idx],
          species = params$species[model_idx],
          model_name = params$model_name[model_idx],
          model_id = model_id,
          use_legacy_covariate = grepl("Legacy with covariate", params$scenario[model_idx]),
          has_driver_uncertainty = FALSE,  # Flag to identify regular models (no driver uncertainty)
          scenario = params$scenario[model_idx],
          min.date = params$min.date[model_idx],
          max.date = params$max.date[model_idx],
          niter = nrow(result$samples),
          nburnin = 500,  # Default burnin
          thin = 1,       # Default thin
          task_idx = task_idx,
          completed_at = Sys.time(),
          model_data = result$model_data,  # Include model_data from result
          nimble_code = result$nimble_code,  # Include nimble_code if available
          model_structure = "standardized_beta_regression_with_restart_capability"  # Model structure identifier
        )
      }
      
      chain_output <- list(
        samples = result$samples,
        metadata = metadata
      )
      
      saveRDS(chain_output, samples_file)
      cat("SAVED: Chain", chain_no, "for model", model_idx, "to", samples_file, "\n")
      
      # Also save a simple status file to track progress
      status_file <- paste0("chain_", model_idx, "_", chain_no, "_status.txt")
      writeLines(paste("SUCCESS", Sys.time(), sep = "\t"), status_file)
    }
    
    cat("=== Completed Model", model_idx, "Chain", chain_no, "===\n")
    cat("Time:", format(Sys.time()), "\n")
    cat("Status:", result$status, "\n")
    
    return(result)
    
  }, error = function(e) {
    # Capture comprehensive error information
    error_time <- Sys.time()
    error_details <- list(
      timestamp = error_time,
      task_idx = task_idx,
      model_idx = if(exists("model_idx")) model_idx else NA,
      chain_no = if(exists("chain_no")) chain_no else NA,
      error_message = if(!is.null(e$message) && e$message != "") e$message else "No error message available",
      error_call = if(!is.null(e$call)) paste(deparse(e$call), collapse=" ") else "No call information",
      error_class = class(e)[1],
      system_info = list(
        r_version = R.version.string,
        working_dir = getwd(),
        available_packages = installed.packages()[,"Package"],
        memory_usage = if(exists("gc")) gc() else "GC not available"
      ),
      runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA
    )
    
    # Create detailed error file
    task <- all_tasks[task_idx, ]
    error_file <- paste0("chain_", task$model_idx, "_", task$chain_no, "_ERROR.txt")
    
    # Write comprehensive error report
    error_report <- c(
      paste("ERROR DETAILED REPORT -", format(error_time)),
      paste("Task Index:", error_details$task_idx),
      paste("Model Index:", error_details$model_idx),
      paste("Chain Number:", error_details$chain_no),
      paste("Error Message:", error_details$error_message),
      paste("Error Call:", error_details$error_call),
      paste("Error Class:", error_details$error_class),
      paste("Runtime (seconds):", round(error_details$runtime, 2)),
      paste("R Version:", error_details$system_info$r_version),
      paste("Working Directory:", error_details$system_info$working_dir),
      paste("Available Packages:", paste(error_details$system_info$available_packages, collapse=", ")),
      "",
      "FULL ERROR OBJECT:",
      capture.output(str(e))
    )
    
    writeLines(error_report, error_file)
    
    # Also log to console with detailed information
    cat("ERROR in Model", error_details$model_idx, "Chain", error_details$chain_no, ":\n")
    cat("  Message:", error_details$error_message, "\n")
    cat("  Call:", error_details$error_call, "\n")
    cat("  Class:", error_details$error_class, "\n")
    cat("  Runtime:", round(error_details$runtime, 2), "seconds\n")
    cat("  Detailed error saved to:", error_file, "\n")
    
    # Return detailed error information
    return(list(
        status = "ERROR", 
        error = error_details$error_message,
        error_details = error_details,
        error_file = error_file
    ))
  })
}
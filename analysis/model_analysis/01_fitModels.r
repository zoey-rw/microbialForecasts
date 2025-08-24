#!/usr/bin/env Rscript

# Model fitting script for calibration time periods (2013-2018 only)
# TESTING: Runs 2 chains per model for testing iterative saving functionality
# Ensures at least 1 model runs for focused testing
#
# HPC FIXES APPLIED:
# - Fixed "model_output_dir not found" error by moving directory creation earlier
# - Added robust directory creation with multiple fallback paths for HPC compatibility
# - Added fallback checkpoint saving to current directory if model_output_dir fails
# - Added HPC environment debugging information
# - Added directory creation testing at script start

# Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# Run with at least 2 cores available (one MCMC chain per core for testing)
nchains = 4

# Test directory creation for HPC compatibility
test_hpc_directories <- function() {
  cat("=== Testing HPC Directory Creation ===\n")
  
  test_bases <- c(
    getwd(),
    here("data", "model_outputs"),
    file.path(Sys.getenv("HOME"), "data", "model_outputs"),
    file.path(Sys.getenv("PWD"), "data", "model_outputs"),
    "./data/model_outputs"
  )
  
  cat("Testing base directories:\n")
  for (i in seq_along(test_bases)) {
    base <- test_bases[[i]]
    cat("  ", i, ".", base, ":", ifelse(is.null(base) || base == "" || base == "NULL", "INVALID", "VALID"), "\n")
  }
  
  # Test creating a test directory
  test_dir <- file.path(getwd(), "test_hpc_dir")
  tryCatch({
    dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
    if (dir.exists(test_dir)) {
      cat("  ✓ Successfully created test directory:", test_dir, "\n")
      unlink(test_dir, recursive = TRUE)
      cat("  ✓ Successfully cleaned up test directory\n")
    } else {
      cat("  ✗ Failed to create test directory\n")
    }
  }, error = function(e) {
    cat("  ✗ Error creating test directory:", e$message, "\n")
  })
  
  cat("=== End HPC Directory Test ===\n\n")
}

# Run the test
test_hpc_directories()

#### Run on all groups ----

source("source.R")

# Function to check if MCMC should continue based on effective sample size
check_continue <- function(samples, min_eff_size = 10) {
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
	
	cat("  ESS check - Min ESS:", round(min_ess, 1), "Target:", min_eff_size, "\n")
	cat("  Continue sampling:", continue, "\n")
	
	return(continue)
}

params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 rep("logical", 2),
										 rep("character", 4)))

rerun_list = readRDS(here("data/summary/unconverged_taxa_list.rds"))
converged_list = readRDS(here("data/summary/converged_taxa_list.rds"))

# TEST CONFIGURATION: Focus on cycl_only and genus_bac models for testing
params <- params_in %>% ungroup %>% filter(
	# Test ONLY cycl_only models with genus_bac rank for focused testing
	rank.name == "genus_bac" &
	# Focus ONLY on 2013-2018 period (exclude 2015-2018)
	scenario %in% c("Legacy with covariate 2013-2018", "2013-2018") &
	# Test ONLY cycl_only model type for simplicity
	model_name == "cycl_only" 
)

# Filter out already converged models
params <- params <- params %>% filter(!model_id %in% converged_list)

# Sample 1 model for focused testing
set.seed(123)  # For reproducible sampling
params <- params %>%
	sample_n(size = min(1, n()), replace = FALSE) %>%  # Test with 1 model
	ungroup()

cat("TESTING CONFIGURATION: Running", nrow(params), "cycl_only genus_bac model with", nchains, "chains\n")
cat("Model configuration:\n")
if (nrow(params) > 0) {
	print(params[, c("rank.name", "species", "model_name", "model_id")])
} else {
	cat("No models to run\n")
}

# Create function that uses our working approach for each model
run_scenarios_fixed <- function(j, chain_no) {
	# Initialize error tracking and logging
	start_time <- Sys.time()
	error_context <- list()
	
	tryCatch({
	# Load required libraries in each worker
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)

	cat("=== Starting model fitting ===\n")
	cat("Model index:", j, "Chain:", chain_no, "\n")
	cat("Model parameters:\n")
	print(params[j,])
	cat("=============================\n")
	
	# Debug HPC environment information
	cat("HPC Environment Debug Info:\n")
	cat("  Working directory:", getwd(), "\n")
	cat("  HOME:", Sys.getenv("HOME"), "\n")
	cat("  PWD:", Sys.getenv("PWD"), "\n")
	cat("  Current user:", Sys.getenv("USER"), "\n")
	cat("  R session tempdir:", tempdir(), "\n")
	cat("=============================\n")

		# Validate input parameters
	if (is.null(params) || nrow(params) < j) {
		stop("Params data frame not available or index out of bounds")
	}

		# Extract model parameters with validation
	rank.name <- params$rank.name[[j]]
	species <- params$species[[j]]
	model_id <- params$model_id[[j]]
	model_name <- params$model_name[[j]]
	min.date <- params$min.date[[j]]
	max.date <- params$max.date[[j]]
	scenario <- params$scenario[[j]]
		
		# Validate extracted parameters
		if (is.null(rank.name) || is.na(rank.name) || rank.name == "") {
			stop("Invalid rank.name for model index ", j)
		}
		if (is.null(species) || is.na(species) || species == "") {
			stop("Invalid species for model index ", j)
		}
		if (is.null(model_name) || is.na(model_name) || model_name == "") {
			stop("Invalid model_name for model index ", j)
		}

	# Check if this is a legacy covariate model
	use_legacy_covariate <- grepl("Legacy with covariate", scenario)
	
		# Validate data availability
		if (!exists("all_ranks") || is.null(all_ranks)) {
			stop("Data 'all_ranks' not available in worker environment")
		}
	
		# Get the specific group data with validation
	if (!(rank.name %in% names(all_ranks))) {
			stop("Rank name '", rank.name, "' not found in data. Available ranks: ", 
				 paste(names(all_ranks), collapse=", "))
	}
	rank.df <- all_ranks[[rank.name]]
		
		# Validate rank data structure
		if (!is.data.frame(rank.df) || nrow(rank.df) == 0) {
			stop("Rank data for '", rank.name, "' is empty or not a data frame")
		}
		
		# Check if species column exists
		if (!(species %in% colnames(rank.df))) {
			stop("Species '", species, "' not found in rank '", rank.name, "'. Available columns: ",
				 paste(colnames(rank.df), collapse=", "))
		}
	
	cat("Preparing model data for", rank.name, "\n")
	
	# Extract the specific species and create "other" column BEFORE calling prepBetaRegData
	cat("DEBUG: Extracting species", species, "from", rank.name, "\n")
		
		# Validate required columns exist
		required_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", species)
		missing_cols <- required_cols[!required_cols %in% colnames(rank.df)]
		if (length(missing_cols) > 0) {
			stop("Missing required columns in rank data: ", paste(missing_cols, collapse=", "))
		}
		
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!species) %>%
		mutate(other = 1 - !!sym(species))
	
	cat("DEBUG: rank.df_spec dimensions:", dim(rank.df_spec), "\n")
	cat("DEBUG: rank.df_spec columns:", colnames(rank.df_spec), "\n")
		
		# Validate species data
		species_data <- rank.df_spec[[species]]
		if (all(is.na(species_data)) || all(species_data == 0) || all(species_data == 1)) {
			stop("Species '", species, "' has no valid variation (all NA, 0, or 1)")
		}
	
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
	if (is.matrix(model.dat$y)) {
		cat("  model.dat$y first few rows:\n")
		print(head(model.dat$y, 3))
	}
	cat("  N.core calculated:", nrow(model.dat$y), "\n")
	cat("  Model type:", model_name, "\n")
	
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
		
		# Add environmental predictors with validation
		env_predictors <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
		for (pred in env_predictors) {
			if (pred %in% names(model.dat)) {
				constants[[pred]] <- model.dat[[pred]]
				cat("  Added", pred, "predictor\n")
			}
		}

	# Add legacy covariate if needed
	if (use_legacy_covariate) {
			cat("Adding legacy covariate...\n")
		# Create legacy covariate: 1 for legacy period (2013-2015), 0 for post-2015
		legacy_dates <- model.dat$plot_date >= "2013-06-27" & model.dat$plot_date <= "2015-11-30"
		constants$legacy <- as.numeric(legacy_dates)
		
		# Validate legacy covariate to prevent extreme values
		legacy_sum <- sum(constants$legacy)
		legacy_total <- length(constants$legacy)
		if (legacy_sum == 0 || legacy_sum == legacy_total) {
			cat("WARNING: Legacy covariate is all 0s or all 1s - this may cause numerical issues\n")
		}
		
		cat("Legacy covariate added:", legacy_sum, "legacy observations out of", legacy_total, "\n")
		cat("Legacy proportion:", round(legacy_sum/legacy_total, 3), "\n")
	}
	
	# Model hyperparameters - adjust based on model type
		cat("Setting model hyperparameters...\n")
	if (model_name == "env_cycl") {
		constants$N.beta = 8
		# Use more reasonable prior precision for beta coefficients
		constants$omega <- 0.1 * diag(8)  # Reduced from 0.0001 to 0.1
		constants$zeros <- rep(0, 8)
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
		constants$omega <- 0.1 * diag(6)  # Reduced from 0.0001 to 0.1
		constants$zeros <- rep(0, 6)
	} else {
		constants$N.beta = 2
		constants$omega <- 0.1 * diag(2)  # Reduced from 0.0001 to 0.1
		constants$zeros <- rep(0, 2)
	}
	
	# Scale cyclical predictors
		cat("Scaling cyclical predictors...\n")
	constants$sin_mo = scale(constants$sin_mo, center = F) %>% as.numeric()
	constants$cos_mo = scale(constants$cos_mo, center = F) %>% as.numeric()
	
	cat("Constants prepared successfully\n")
	
	# STANDARDIZED MODEL DEFINITIONS - All models use consistent priors
		cat("Building Nimble model...\n")
	if (model_name == "cycl_only" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}
			for (p in 1:N.plot) {
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
				for (t in plot_index[p]:N.date) {
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
						site_effect[plot_site_num[p]] + 
						legacy_effect * legacy[t] +
						intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}

			# STANDARDIZED PRIORS - Consistent across all models
			site_effect_sd ~ dgamma(2, 6)  # Standardized: moderate prior
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			precision ~ dgamma(2, 0.1)  # Standardized: moderate precision prior
			intercept ~ dnorm(-2, sd = 0.8)  # Standardized: normal prior, moderate variance
			rho ~ dbeta(3, 3)  # Standardized: informative beta prior
			legacy_effect ~ dnorm(0, sd = 0.2)  # Standardized: moderate legacy prior
			
			for (b in 1:2) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Standardized: moderate beta prior
			}
		})
	} else if (model_name == "env_cycl" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations - CONVERTED TO BETA REGRESSION
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}

			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
				
				for (t in plot_index[p]:N.date) {
					# Dynamic linear model with environmental + cyclical predictors
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						site_effect[plot_site_num[p]] +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, plot_start[p]] +
						beta[4] * pC[p, plot_start[p]] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						beta[7] * sin_mo[t] +
						beta[8] * cos_mo[t] +
						legacy_effect * legacy[t] +  # LEGACY COVARIATE
						intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}
			
			# Hierarchical priors for site effects - better mixing with constraints
			site_effect_sd ~ dgamma(2, 8)  # Prior on site effect standard deviation
			for (k in 1:N.site) {
				site_effect[k] ~ T(dnorm(0, sd = site_effect_sd), -5, 5)  # Truncated normal priors
			}
			
			# Parameter-specific priors for better convergence
			# precision: Use appropriate prior for beta regression precision parameter
			precision ~ dgamma(2, 0.1)  # Moderate precision prior
			
			# intercept: Use tighter prior for better convergence
			intercept ~ dt(-2, 0.3, df = 3)  # Centered around expected value, tighter variance

			# rho: Use tighter bounded prior for better convergence
			rho ~ dbeta(3, 3)    # More informative, prevents extreme values
			
			# Legacy effect: Account for research facility bias - TIGHT prior to prevent explosion
			legacy_effect ~ dnorm(0, sd = 0.1)  # Very tight prior centered at zero
			
			# Tight variance for ecological effects to ensure stability
			for (b in 1:8) {
				beta[b] ~ dnorm(0, sd=0.15)
			}
		})
	} else if (model_name == "env_cycl" && !use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations - CONVERTED TO BETA REGRESSION
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}

			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
				
				for (t in plot_index[p]:N.date) {
					# Dynamic linear model with environmental + cyclical predictors
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						site_effect[plot_site_num[p]] +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, plot_start[p]] +
						beta[4] * pC[p, plot_start[p]] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						beta[7] * sin_mo[t] +
						beta[8] * cos_mo[t] +
						intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}
			
			# STANDARDIZED PRIORS - Consistent across all models
			site_effect_sd ~ dgamma(2, 6)  # Standardized: moderate prior
			for (k in 1:N.site) {
				site_effect[k] ~ T(dnorm(0, sd = site_effect_sd), -5, 5)  # Truncated normal priors
			}
			
			precision ~ dgamma(2, 0.1)  # Standardized: moderate precision prior
			intercept ~ dnorm(-2, sd = 0.8)  # Standardized: normal prior, moderate variance
			rho ~ dbeta(3, 3)  # Standardized: informative beta prior
			
			for (b in 1:8) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Standardized: moderate beta prior
			}
		})
	} else if (model_name == "env_cov" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations - CONVERTED TO BETA REGRESSION
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}

			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}

				for (t in plot_index[p]:N.date) {
					# Dynamic linear model with environmental predictors only (no seasonal) + legacy covariate
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						site_effect[plot_site_num[p]] +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, plot_start[p]] +
						beta[4] * pC[p, plot_start[p]] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						legacy_effect * legacy[t] +  # LEGACY COVARIATE
						intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}

			# STANDARDIZED PRIORS - Consistent across all models
			site_effect_sd ~ dgamma(2, 6)  # Standardized: moderate prior
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}

			precision ~ dgamma(2, 0.1)  # Standardized: moderate precision prior
			intercept ~ dnorm(-2, sd = 0.8)  # Standardized: normal prior, moderate variance
			rho ~ dbeta(3, 3)  # Standardized: informative beta prior
			legacy_effect ~ dnorm(0, sd = 0.2)  # Standardized: moderate legacy prior

			for (b in 1:6) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Standardized: moderate beta prior
			}
		})
	} else if (model_name == "env_cov" && !use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations - CONVERTED TO BETA REGRESSION
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}

			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}

				for (t in plot_index[p]:N.date) {
					# Dynamic linear model with environmental predictors only (no seasonal)
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						site_effect[plot_site_num[p]] +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, plot_start[p]] +
						beta[4] * pC[p, plot_start[p]] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}

			# STANDARDIZED PRIORS - Consistent across all models
			site_effect_sd ~ dgamma(2, 6)  # Standardized: moderate prior
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}

			precision ~ dgamma(2, 0.1)  # Standardized: moderate precision prior
			intercept ~ dnorm(-2, sd = 0.8)  # Standardized: normal prior, moderate variance
			rho ~ dbeta(3, 3)  # Standardized: informative beta prior

			for (b in 1:6) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Standardized: moderate beta prior
			}
		})
	} else {
		# Default to cycl_only model without legacy
		modelCode <- nimble::nimbleCode({
			for (i in 1:N.core) {
				y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision, 
								 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
			}
			for (p in 1:N.plot) {
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001, 0.9999)
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
				for (t in plot_index[p]:N.date) {
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
						site_effect[plot_site_num[p]] + intercept
					plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
				}
			}
			
			# STANDARDIZED PRIORS - Consistent across all models
			site_effect_sd ~ dgamma(2, 6)  # Standardized: moderate prior
			for (k in 1:N.site) { 
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			precision ~ dgamma(2, 0.1)  # Standardized: moderate precision prior
			intercept ~ dnorm(-2, sd = 0.8)  # Standardized: normal prior, moderate variance
			rho ~ dbeta(3, 3)  # Standardized: informative beta prior
			
			for (b in 1:2) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Standardized: moderate beta prior
			}
		})
	}
	
	# Create inits
		cat("Creating initial values...\n")
	inits <- createInits(constants)
	
	# Calculate data-informed initial values for better convergence
	y_data <- model.dat$y[, 1]  # Get the species abundance data
	y_mean <- mean(y_data, na.rm = TRUE)
	y_sd <- sd(y_data, na.rm = TRUE)
	
	# Simplified initialization strategy - spread chains out properly
	set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model
	
	# Initialize parameters with proper spacing between chains
	inits$rho <- runif(1, 0.2, 0.8)  # Spread chains out more
	inits$intercept <- rnorm(1, logit(y_mean), 0.5)
	inits$precision <- rgamma(1, 2, 0.1)
	inits$beta <- rnorm(constants$N.beta, 0, 0.1)
	
	# Initialize hierarchical parameters
	inits$site_effect_sd <- max(0.001, min(0.3, y_sd * 0.05))
	
	# Initialize legacy effect if using legacy covariate
	if (use_legacy_covariate) {
		inits$legacy_effect <- rnorm(1, 0, 0.05)  # Start close to zero
	}
	
	cat("Model built successfully\n")
	
	# Build model
		cat("Building Nimble model...\n")
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y=model.dat$y), inits = inits)
	
	# Debug: Check data dimensions
	cat("Data dimensions check:\n")
	cat("  model.dat$y dimensions:", dim(model.dat$y), "\n")
	cat("  model.dat$y class:", class(model.dat$y), "\n")
	cat("  constants$N.core:", constants$N.core, "\n")
	cat("  constants$N.spp:", constants$N.spp, "\n")
	
	# Compile model
		cat("Compiling Nimble model...\n")
	cModel <- compileNimble(Rmodel)
	
	cat("Model compiled successfully\n")
	
	# Configure MCMC with proper sampler management - all models now use precision parameter
		cat("Configuring MCMC...\n")
	monitors <- c("beta","precision","site_effect","site_effect_sd","intercept","rho")
	
	if (use_legacy_covariate) {
		monitors <- c(monitors, "legacy_effect")
	}
	mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)
	
	# Remove default samplers before adding specialized ones
	mcmcConf$removeSamplers(c("precision", "rho", "site_effect_sd", "intercept"))
	
	# Remove legacy_effect if using legacy covariate
	if (use_legacy_covariate) {
		mcmcConf$removeSamplers("legacy_effect")
	}
	
	# Remove beta samplers if they exist
	if (constants$N.beta > 1) {
		mcmcConf$removeSamplers(paste0("beta[1:", constants$N.beta, "]"))
	} else {
		mcmcConf$removeSamplers("beta[1]")
	}
	
	# Remove site effect samplers
	if (constants$N.site > 1) {
		mcmcConf$removeSamplers(paste0("site_effect[1:", constants$N.site, "]"))
	} else {
		mcmcConf$removeSamplers("site_effect[1]")
	}
	
	# All models use precision parameter
	mcmcConf$addSampler(target = "precision", type = "slice")
	cat("  Added slice sampler for precision\n")
	
	# site_effect_sd: Site effect scale (all models)
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
	cat("  Added slice sampler for site_effect_sd\n")
	
	# intercept: Use slice sampler for better location parameter exploration (all models)
	mcmcConf$addSampler(target = "intercept", type = "slice")
	cat("  Added slice sampler for intercept\n")
	
	# rho: Use specialized bounded sampler for temporal persistence (all models)
	mcmcConf$addSampler(target = "rho", type = "slice")
	cat("  Added slice sampler for rho\n")
	
	# legacy_effect: Use slice sampler for better mixing of legacy parameter
	if (use_legacy_covariate) {
		mcmcConf$addSampler(target = "legacy_effect", type = "slice")
		cat("  Added slice sampler for legacy_effect\n")
	}
	
	# Use EITHER block OR individual samplers for site effects, NOT BOTH
	if (constants$N.site > 1) {
		# Use block sampling for correlated site effects
		mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, "]"), type = "AF_slice")
		cat("  Added block sampler for site_effect[1:", constants$N.site, "]\n")
	} else {
		# Individual sampler for single site
		mcmcConf$addSampler(target = "site_effect[1]", type = "slice")
		cat("  Added slice sampler for site_effect[1]\n")
	}
	
	# Use block sampler for beta parameters when multiple exist, individual for single
	if (constants$N.beta > 1) {
		# Use block sampler for correlated beta parameters (more efficient)
		mcmcConf$addSampler(target = paste0("beta[1:", constants$N.beta, "]"), type = "AF_slice")
		cat("  Added block AF_slice sampler for beta[1:", constants$N.beta, "]\n")
	} else {
		# Individual sampler for single beta
		mcmcConf$addSampler(target = "beta[1]", type = "slice")
		cat("  Added slice sampler for beta[1]\n")
	}
	
	# Build and compile MCMC
		cat("Building and compiling MCMC...\n")
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# TESTING: Run MCMC with convergence-based sampling for testing
	cat("Running MCMC with convergence-based sampling...\n")
	burnin <- 100    # Reduced burnin for testing
	thin <- 1        # No thinning to preserve all samples
	iter_per_chunk <- 200   # Reduced iterations per convergence check for testing
	init_iter <- 100  # Reduced initial iterations for testing
	min_eff_size_perchain <- 5  # Reduced ESS target for testing
	max_loops <- 10  # Reduced maximum loops for testing
	max_save_size <- 60000  # Maximum samples to keep in memory
	min_total_iterations <- 500  # Reduced minimum iterations for testing
	
	cat("TESTING: Running MCMC with convergence-based sampling (reduced parameters)\n")
	cat("  Initial iterations:", init_iter, "burnin:", burnin, "\n")
	cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
	cat("  Target ESS per chain:", min_eff_size_perchain, "\n")
	cat("  Minimum total iterations:", min_total_iterations, "\n")
	
	# Run initial iterations with progress reporting and adaptation
	cat("  Running initial iterations (", init_iter, " iterations) for adaptation...\n")
	compiled$run(niter = init_iter, thin = thin, nburnin = 0)
	cat("  Initial iterations completed\n")
	
	# Get initial samples and check convergence
	initial_samples <- as.matrix(compiled$mvSamples)
	cat("  Initial samples collected, checking convergence...\n")
	cat("  Initial samples dimensions:", dim(initial_samples), "\n")
	
	# Create output directories early for checkpoint saving
	# Try multiple possible base directories for HPC compatibility
	possible_bases <- c(
		getwd(),
		here("data", "model_outputs"),
		file.path(Sys.getenv("HOME"), "data", "model_outputs"),
		file.path(Sys.getenv("PWD"), "data", "model_outputs"),
		"./data/model_outputs"
	)
	
	model_output_dir <- NULL
	for (base_dir in possible_bases) {
		if (!is.null(base_dir) && base_dir != "" && base_dir != "NULL") {
			test_dir <- file.path(base_dir, "logit_beta_regression", model_name)
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
		model_output_dir <- file.path("data", "model_outputs", "logit_beta_regression", model_name)
		dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
		cat("  WARNING: Using fallback output directory:", model_output_dir, "\n")
	}
	
	# Create model_id for consistent naming with legacy covariate indicator
	legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
	model_id <- paste(model_name, species, min.date, max.date, legacy_indicator, sep = "_")
	
	# Check if we need to continue sampling for convergence
	# Try to check convergence, with fallback if it fails
	tryCatch({
		continue <- check_continue(initial_samples, min_eff_size = min_eff_size_perchain)
	}, error = function(e) {
		cat("  WARNING: Convergence check failed, defaulting to continue sampling\n")
		cat("  Error:", e$message, "\n")
		continue <- TRUE  # Default to continue if check fails
	})
	loop_counter <- 0
	total_iterations <- init_iter
	
	# Store all samples as we go - FIXED: Use initial samples as starting point
	all_samples <- initial_samples
	cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
	
	# Save initial samples as checkpoint
	checkpoint_file <- file.path(model_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_initial.rds"))
	tryCatch({
		saveRDS(list(samples = all_samples, iterations = total_iterations, loop = 0), checkpoint_file)
		cat("  ✓ Checkpoint saved: Initial samples (", nrow(all_samples), " iterations)\n")
	}, error = function(e) {
		cat("  ✗ Failed to save initial checkpoint:", e$message, "\n")
		cat("  Attempting to save to current directory as fallback...\n")
		# Fallback: save to current directory
		fallback_checkpoint <- paste0("checkpoint_", model_id, "_chain", chain_no, "_initial_FALLBACK.rds")
		tryCatch({
			saveRDS(list(samples = all_samples, iterations = total_iterations, loop = 0), fallback_checkpoint)
			cat("  ✓ Fallback checkpoint saved:", fallback_checkpoint, "\n")
		}, error = function(e2) {
			cat("  ✗ CRITICAL: Failed to save even fallback checkpoint:", e2$message, "\n")
		})
	})
	
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
		
		# Get updated samples and accumulate them - FIXED: Get only new samples
		current_samples <- as.matrix(compiled$mvSamples)
		cat("  Current total samples in compiled object:", nrow(current_samples), "\n")
		cat("  Previous accumulated samples:", nrow(all_samples), "\n")
		
		# Only take the new samples (skip the initial ones we already have)
		if (nrow(current_samples) > nrow(initial_samples)) {
			new_samples <- current_samples[(nrow(initial_samples) + 1):nrow(current_samples), , drop = FALSE]
			all_samples <- rbind(all_samples, new_samples)
			cat("  Updated samples collected:", nrow(new_samples), "new samples,", nrow(all_samples), "total accumulated\n")
		} else {
			cat("  WARNING: No new samples detected, using current samples\n")
			all_samples <- current_samples
		}
		
		# Save checkpoint after each loop
		checkpoint_file <- file.path(model_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_loop", loop_counter + 1, ".rds"))
		tryCatch({
			saveRDS(list(samples = all_samples, iterations = total_iterations, loop = loop_counter + 1), checkpoint_file)
			cat("  ✓ Checkpoint saved: Loop", loop_counter + 1, "(", nrow(all_samples), " iterations)\n")
		}, error = function(e) {
			cat("  ✗ Failed to save checkpoint for loop", loop_counter + 1, ":", e$message, "\n")
			cat("  Attempting to save to current directory as fallback...\n")
			# Fallback: save to current directory
			fallback_checkpoint <- paste0("checkpoint_", model_id, "_chain", chain_no, "_loop", loop_counter + 1, "_FALLBACK.rds")
			tryCatch({
				saveRDS(list(samples = all_samples, iterations = total_iterations, loop = loop_counter + 1), fallback_checkpoint)
				cat("  ✓ Fallback checkpoint saved:", fallback_checkpoint, "\n")
			}, error = function(e2) {
				cat("  ✗ CRITICAL: Failed to save even fallback checkpoint:", e2$message, "\n")
			})
		})
		
		# Check if we need to continue
		tryCatch({
			continue <- check_continue(all_samples, min_eff_size = min_eff_size_perchain)
		}, error = function(e) {
			cat("  WARNING: Convergence check failed in loop, defaulting to continue sampling\n")
			cat("  Error:", e$message, "\n")
			continue <- TRUE  # Default to continue if check fails
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
	
	# Get final samples (use accumulated samples)
	samples <- all_samples
	
	cat("MCMC completed successfully\n")
	cat("Final sample dimensions:", dim(samples), "\n")
	cat("Total iterations run:", total_iterations, "\n")
	cat("Convergence loops:", loop_counter, "\n")
	cat("Final ESS check:\n")
	
	# Final convergence check
	tryCatch({
		final_ess <- effectiveSize(as.mcmc(samples))
		min_final_ess <- min(final_ess, na.rm = TRUE)
		cat("  Final minimum ESS:", round(min_final_ess, 1), "\n")
		cat("  Convergence achieved:", min_final_ess >= min_eff_size_perchain, "\n")
	}, error = function(e) {
		cat("  Final ESS check failed:", e$message, "\n")
	})
	
	cat("=== ITERATIVE SAVING SUMMARY ===\n")
	cat("  Initial samples:", nrow(initial_samples), "iterations\n")
	cat("  Additional loops:", loop_counter, "\n")
	cat("  Total accumulated samples:", nrow(all_samples), "iterations\n")
	cat("  Checkpoints saved:", loop_counter + 1, "files\n")
	cat("  Final sample size:", nrow(samples), "iterations\n")
	
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
			nimble_code = modelCode,  # Save the actual Nimble code used
			model_structure = "standardized_beta_regression_with_consistent_priors"  # Model structure identifier
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
	cat("  - TESTING: cycl_only genus_bac model with iterative saving\n")
	cat("  - Beta regression with precision parameter\n")
	cat("  - Block samplers for efficient parameter sampling\n")
	cat("  - CONVERGENCE-BASED: Adaptive sampling until reasonable ESS reached\n")
	cat("  - ITERATIVE SAVING: Samples accumulated and saved incrementally\n")
	
	return(list(
		status = "SUCCESS", 
		samples = samples, 
		file = samples_file,
		model_data = model.dat,  # Include model_data for parallel execution
		nimble_code = modelCode,  # Include nimble code for parallel execution
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
			model_structure = "standardized_beta_regression_with_consistent_priors"
		)
	))
		
	}, error = function(e) {
		# Capture comprehensive error information
		error_time <- Sys.time()
		error_context <- list(
			timestamp = error_time,
			task_idx = j, # Use j for model index
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
			runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA
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

# Run multiple models for testing
cat("Testing with", nrow(params), "models\n")
cat("Models to test:\n")
print(params[, c("rank.name", "species", "model_name", "model_id")])

# Test with all available models across all ranks
test_models <- nrow(params)  # Test all available models across all ranks
cat("\nTesting", test_models, "models across all ranks to verify convergence fixes\n")

# Set up parallel cluster for Nimble
library(parallel)
library(doParallel)

# Create cluster with exactly the number of cores needed for testing
# For testing: 2 chains need exactly 2 cores
ncores <- test_models * nchains  # Exactly 1 model × 2 chains = 2 cores
cat("Creating cluster with exactly", ncores, "cores for", test_models, "model ×", nchains, "chains\n")
cl <- makeCluster(ncores, type = "PSOCK")

# Register the cluster
registerDoParallel(cl)

# Verify cluster size
actual_workers <- length(cl)
cat("Cluster created with", actual_workers, "workers\n")
if (actual_workers != ncores) {
  cat("WARNING: Expected", ncores, "workers but got", actual_workers, "\n")
}

# Export necessary objects to workers
clusterExport(cl, c("params", "k", "nchains"))

# Load data once before parallel execution
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Export the function and data to workers
clusterExport(cl, c("run_scenarios_fixed", "params", "all_ranks", "check_continue"))

# Run in parallel for faster execution
	cat("TESTING: Starting parallel execution for", test_models, "cycl_only genus_bac model with", nchains, "chains using", ncores, "cores\n")
	cat("Expected runtime: Variable (convergence-based sampling)\n")
	cat("  - Models to test:", test_models, "cycl_only genus_bac model\n")
	cat("  - Chains per model:", nchains, "(total", test_models * nchains, "parallel tasks)\n")
	cat("  - Initial iterations: ~", round(test_models * 4 / 6, 1), "minutes\n")
	cat("  - Additional iterations: Variable based on convergence\n")
	cat("  - Target: ESS >= 10 per parameter\n")
start_time <- Sys.time()

# Run TEST model with 2 chains in parallel
cat("TESTING: Running", test_models, "cycl_only genus_bac model with", nchains, "chains in parallel\n")
cat("This tests the iterative saving and convergence-based sampling functionality\n")

# Create a combined task list: (model_idx, chain_no) pairs
all_tasks <- expand.grid(model_idx = 1:test_models, chain_no = 1:nchains)
cat("Total parallel tasks:", nrow(all_tasks), "(", test_models, "models ×", nchains, "chains)\n")
cat("Task details:\n")
print(all_tasks)
cat("Cluster size:", ncores, "workers\n")
cat("Starting parallel execution with foreach...\n")

# Function to monitor progress in real-time
monitor_progress <- function() {
  cat("\n=== REAL-TIME PROGRESS MONITORING ===\n")
  cat("Press Ctrl+C to stop monitoring\n")
  
  while(TRUE) {
    Sys.sleep(30)  # Check every 30 seconds
    
    completed <- 0
    errors <- 0
    for (model_idx in 1:test_models) {
      for (chain_no in 1:nchains) {
        status_file <- paste0("chain_", model_idx, "_", chain_no, "_status.txt")
        error_file <- paste0("chain_", model_idx, "_", chain_no, "_ERROR.txt")
        
        if (file.exists(status_file)) completed <- completed + 1
        if (file.exists(error_file)) errors <- errors + 1
      }
    }
    
    total <- test_models * nchains
    cat(format(Sys.time()), "- Progress:", completed, "/", total, "completed,", 
        errors, "failed,", total - completed - errors, "running\n")
  }
}

# Start progress monitoring in background (optional)
cat("To monitor progress in real-time, run: monitor_progress()\n")
cat("Or check status files manually:\n")
cat("  - chain_[model]_[chain]_status.txt for completed chains\n")
cat("  - chain_[model]_[chain]_ERROR.txt for failed chains\n")

# Run everything in parallel with incremental saving
cat("Starting parallel execution with incremental saving at:", format(Sys.time()), "\n")

# Create a function that saves results as they complete
runAndSave_task <- function(task_idx) {
  # Initialize error tracking
  error_details <- list()
  start_time <- Sys.time()
  
  tryCatch({
    # Get task details
    task <- all_tasks[task_idx, ]
    model_idx <- task$model_idx
    chain_no <- task$chain_no
    
    cat("Worker: Model", model_idx, "Chain", chain_no, "starting at", format(start_time), "\n")
    
    # Log system information for debugging
    cat("Worker: System info - R version:", R.version.string, "\n")
    cat("Worker: Working directory:", getwd(), "\n")
    cat("Worker: Available packages:", paste(installed.packages()[,"Package"], collapse=", "), "\n")
    
    # Check if required packages are loaded
    required_packages <- c("nimble", "microbialForecast", "here", "tidyverse", "coda")
    missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
    if (length(missing_packages) > 0) {
      stop("Missing required packages: ", paste(missing_packages, collapse=", "))
    }
    
    # Check if data is available
    if (!exists("all_ranks") || is.null(all_ranks)) {
      stop("Data 'all_ranks' not available in worker environment")
    }
    
    # Check if params is available
    if (!exists("params") || is.null(params) || nrow(params) == 0) {
      stop("Parameters 'params' not available or empty in worker environment")
    }
    
    # Validate model index
    if (model_idx > nrow(params)) {
      stop("Model index ", model_idx, " exceeds available models (", nrow(params), ")")
    }
    
    cat("Worker: All checks passed, calling run_scenarios_fixed...\n")
    
    # Run the model with detailed error context
    result <- run_scenarios_fixed(j = model_idx, chain_no = chain_no)
    
    # Validate result structure
    if (!is.list(result) || !("status" %in% names(result))) {
      stop("Invalid result structure from run_scenarios_fixed")
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
      
      # Create model_id for consistent naming
      legacy_indicator <- ifelse(grepl("Legacy with covariate", params$scenario[model_idx]), 
                                "with_legacy_covariate", "without_legacy_covariate")
      model_id <- paste(params$model_name[model_idx], params$species[model_idx], 
                       params$min.date[model_idx], params$max.date[model_idx], 
                       legacy_indicator, sep = "_")
      
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
          model_structure = "standardized_beta_regression_with_consistent_priors"  # Model structure identifier
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
    
    cat("Worker: Model", model_idx, "Chain", chain_no, "completed at", format(Sys.time()), "\n")
    return(list(model_idx = model_idx, chain_no = chain_no, result = result))
    
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
      model_idx = error_details$model_idx, 
      chain_no = error_details$chain_no, 
      result = list(
        status = "ERROR", 
        error = error_details$error_message,
        error_details = error_details,
        error_file = error_file
      )
    ))
  })
}

# Export the function to workers
clusterExport(cl, c("runAndSave_task", "params", "all_tasks"))

# Run everything in parallel with incremental saving
all_results_parallel = foreach(task_idx = 1:nrow(all_tasks), 
                             .packages = c("nimble", "microbialForecast", "here", "tidyverse", "coda"),
                             .export = c("runAndSave_task", "params", "all_tasks")) %dopar% {
  runAndSave_task(task_idx)
}
cat("Parallel execution completed at:", format(Sys.time()), "\n")

# Show progress summary
cat("\n=== PROGRESS SUMMARY ===\n")
cat("Checking which chains have been completed...\n")

# Count completed chains
completed_chains <- 0
error_chains <- 0
for (model_idx in 1:test_models) {
  for (chain_no in 1:nchains) {
    status_file <- paste0("chain_", model_idx, "_", chain_no, "_status.txt")
    error_file <- paste0("chain_", model_idx, "_", chain_no, "_ERROR.txt")
    
    if (file.exists(status_file)) {
      completed_chains <- completed_chains + 1
      cat("✓ Model", model_idx, "Chain", chain_no, "completed\n")
    } else if (file.exists(error_file)) {
      error_chains <- error_chains + 1
      cat("✗ Model", model_idx, "Chain", chain_no, "failed\n")
    } else {
      cat("? Model", model_idx, "Chain", chain_no, "status unknown\n")
    }
  }
}

cat("\nProgress Summary:\n")
cat("  Completed chains:", completed_chains, "/", test_models * nchains, "\n")
cat("  Failed chains:", error_chains, "/", test_models * nchains, "\n")
cat("  Success rate:", round(completed_chains / (test_models * nchains) * 100, 1), "%\n")

# Reorganize results by model
all_results <- list()
for (model_idx in 1:test_models) {
  all_results[[model_idx]] <- list()
  for (chain_no in 1:nchains) {
    # Find the result for this model/chain combination
    task_result <- all_results_parallel[[which(all_tasks$model_idx == model_idx & all_tasks$chain_no == chain_no)]]
    all_results[[model_idx]][[chain_no]] <- task_result$result
  }
}

# Stop the cluster
stopCluster(cl)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("ALL MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Summary of all models
cat("\nSummary of All Models:\n")
for (model_idx in 1:test_models) {
  cat("\nModel", model_idx, ":", params$species[model_idx], "(", params$model_name[model_idx], ")\n")
  
  output.list <- all_results[[model_idx]]
  
  # Status summary for this model
  status_summary <- sapply(output.list, function(x) {
    if (is.list(x) && "status" %in% names(x)) {
      x$status
    } else {
      "ERROR"
    }
  })
  
  cat("  Results:", paste(status_summary, collapse = ", "), "\n")
  
  # Detailed status for this model
  for (i in 1:length(output.list)) {
    if (is.list(output.list[[i]]) && "status" %in% names(output.list[[i]])) {
      if (output.list[[i]]$status == "SUCCESS") {
        cat("    Chain", i, ": SUCCESS - Samples:", dim(output.list[[i]]$samples)[1], "iterations\n")
      } else {
        cat("    Chain", i, ": ERROR -", output.list[[i]]$error, "\n")
      }
    } else {
      cat("    Chain", i, ": ERROR - Unexpected output format\n")
    }
  }
}

# Clean up status files
for (chain_no in 1:nchains) {
  status_file <- paste0("chain_", chain_no, "_status.txt")
  if (file.exists(status_file)) {
    unlink(status_file)
  }
}

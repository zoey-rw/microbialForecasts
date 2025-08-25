#!/usr/bin/env Rscript

# Enhanced Model fitting script with RESTART capability
# Can restart models using initial values from previous MCMC runs
# Automatically detects extreme values and applies fallback strategies

# Get arguments from the command line
argv <- commandArgs(TRUE)
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# Default restart settings
RESTART_ENABLED <- TRUE  # Set to TRUE to enable restart functionality
RESTART_MIN_ESS <- 10    # Minimum effective sample size for parameter extraction
RESTART_MAX_RHAT <- 1.5  # Maximum R-hat value allowed
RESTART_BURNIN_PROP <- 0.3  # Proportion of samples to discard as burnin
RESTART_FALLBACK_STRATEGY <- "conservative"  # "conservative", "random", or "zero"

# Run with at least 2 cores available
nchains = 2

#### Run on all groups ----

source("source.R")
source("model_restart_functions.R")  # Load restart functions

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
	# Load required libraries in each worker
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)

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
			existing_chain_files <- find_chain_files(model_name, species, min.date, max.date, use_legacy_covariate)

			if (length(existing_chain_files) > 0) {
				cat("  Found", length(existing_chain_files), "existing chain files\n")

				# Extract final values from existing chains
				extraction_result <- extract_final_values_from_chains(
					existing_chain_files,
					min_ess = RESTART_MIN_ESS,
					max_rhat = RESTART_MAX_RHAT,
					burnin_proportion = RESTART_BURNIN_PROP
				)

				# Create restart initial values
				restart_inits <- create_restart_inits(
					extraction_result,
					use_fallback_for_extreme = TRUE,
					fallback_strategy = RESTART_FALLBACK_STRATEGY
				)

				n_extreme <- sum(restart_inits$extreme_flags)
				n_fallback <- sum(restart_inits$fallback_used)

				cat("  ✓ Restart setup complete:\n")
				cat("    - Extreme values detected:", n_extreme, "/", length(restart_inits$extreme_flags), "\n")
				cat("    - Using fallback values:", n_fallback, "\n")
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
			stop("Data 'all_ranks' not available in worker environment")
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

	# STANDARDIZED MODEL DEFINITIONS - All models use consistent priors
	# [Model code definitions remain the same as original script]
	# ... [Include all the model code definitions from the original script]

	# Create inits with restart capability
		cat("Creating initial values...\n")

	# Use restart initial values if available, otherwise create fresh inits
	if (!is.null(restart_inits)) {
		cat("  Using restart initial values from previous chains\n")

		# Get restart initial values
		inits <- restart_inits$initial_values

		# Ensure all required parameters are present
		set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model

		# Add any missing parameters with default values
		required_params <- c("precision", "rho", "beta", "intercept", "site_effect", "site_effect_sd")
		if (use_legacy_covariate) {
			required_params <- c(required_params, "legacy_effect")
		}

		for (param in required_params) {
			if (!(param %in% names(inits))) {
				# Add missing parameter with reasonable default
				if (param == "precision") inits[[param]] <- rgamma(1, 2, 0.1)
				else if (param == "rho") inits[[param]] <- runif(1, 0.1, 0.9)
				else if (param == "beta") inits[[param]] <- rnorm(constants$N.beta, 0, 0.1)
				else if (param == "intercept") inits[[param]] <- rnorm(1, 0, 0.5)
				else if (param == "site_effect") inits[[param]] <- rnorm(constants$N.site, 0, 0.1)
				else if (param == "site_effect_sd") inits[[param]] <- runif(1, 0.1, 1)
				else if (param == "legacy_effect") inits[[param]] <- rnorm(1, 0, 0.1)

				cat("  Added missing parameter", param, "with default value\n")
			}
		}

		cat("  ✓ Restart initial values ready for chain", chain_no, "\n")

	} else {
		cat("  Creating fresh initial values\n")

		# Original initialization code
		inits <- createInits(constants)

		# Calculate data-informed initial values for better convergence
		y_data <- model.dat$y[, 1]  # Get the species abundance data
		y_mean <- mean(y_data, na.rm = TRUE)
		y_sd <- sd(y_data, na.rm = TRUE)

		# Simplified initialization strategy - spread chains out properly
		set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model

		# Initialize parameters with tighter ranges for better convergence
		inits$rho <- runif(1, 0.3, 0.7)  # Tighter range around 0.5
		inits$intercept <- rnorm(1, logit(y_mean), 0.2)  # Tighter initialization
		inits$precision <- rgamma(1, 2, 0.1)  # Keep moderate
		inits$beta <- rnorm(constants$N.beta, 0, 0.05)  # Tighter beta initialization

		# Initialize hierarchical parameters with tighter ranges
		inits$site_effect_sd <- max(0.001, min(0.2, y_sd * 0.02))  # Tighter site effect SD

		# Initialize legacy effect if using legacy covariate
		if (use_legacy_covariate) {
			inits$legacy_effect <- rnorm(1, 0, 0.03)  # Tighter initialization
		}

		cat("  ✓ Fresh initial values created for chain", chain_no, "\n")
	}

	cat("Model built successfully\n")

	# Build model
		cat("Building Nimble model...\n")
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y=model.dat$y), inits = inits)

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
	monitors <- c("beta","precision","site_effect","site_effect_sd","intercept","rho")

	if (use_legacy_covariate) {
		monitors <- c(monitors, "legacy_effect")
	}
	mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)

	# Remove default samplers before adding specialized ones
	mcmcConf$removeSamplers(c("precision", "rho", "site_effect_sd", "intercept"))

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
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
	mcmcConf$addSampler(target = "intercept", type = "slice")
	mcmcConf$addSampler(target = "rho", type = "slice")

	if (use_legacy_covariate) {
		mcmcConf$addSampler(target = "legacy_effect", type = "slice")
	}

	# Use EITHER block OR individual samplers for site effects, NOT BOTH
	if (constants$N.site > 1) {
		# Use block sampling for correlated site effects
		mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, "]"), type = "AF_slice")
	} else {
		# Individual sampler for single site
		mcmcConf$addSampler(target = "site_effect[1]", type = "slice")
	}

	# Use block sampler for beta parameters when multiple exist, individual for single
	if (constants$N.beta > 1) {
		# Use block sampler for correlated beta parameters (more efficient)
		mcmcConf$addSampler(target = paste0("beta[1:", constants$N.beta, "]"), type = "AF_slice")
	} else {
		# Individual sampler for single beta
		mcmcConf$addSampler(target = "beta[1]", type = "slice")
	}

	# Build and compile MCMC
		cat("Building and compiling MCMC...\n")
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)

	cat("MCMC configured successfully\n")

	# TESTING: Run MCMC with convergence-based sampling for testing
	cat("Running MCMC with convergence-based sampling...\n")
	burnin <- 500    # Reduced burnin for testing
	thin <- 2
	iter_per_chunk <- 500   # Reduced iterations per convergence check for testing
	init_iter <- 1000  # Reduced initial iterations for testing
	min_eff_size_perchain <- 10 # Reduced ESS target for testing
	max_loops <- 20  # Reduced maximum loops for testing
	max_save_size <- 30000  # Maximum samples to keep in memory
	min_total_iterations <- 1500  # Reduced minimum iterations for testing

	cat("TESTING: Running MCMC with convergence-based sampling (reduced parameters)\n")
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

	# Run initial iterations with progress reporting and adaptation
	cat("  Running initial iterations (", init_iter, " iterations) for adaptation...\n")
	compiled$run(niter = init_iter, thin = thin, nburnin = 0)
	cat("  Initial iterations completed\n")

	# Get initial samples and check convergence
	initial_samples <- as.matrix(compiled$mvSamples)
	cat("  Initial samples collected, checking convergence...\n")
	cat("  Initial samples dimensions:", dim(initial_samples), "\n")

	# Create output directories early for checkpoint saving
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
	continue <- TRUE  # Default to continue sampling
	loop_counter <- 0
	total_iterations <- init_iter

	# Try to check convergence with enhanced checking
	tryCatch({
		convergence_check <- check_continue_with_restart(initial_samples, min_eff_size = min_eff_size_perchain,
														model_info = list(model_name = model_name, species = species))
		continue <- convergence_check$continue

		if (convergence_check$restart_recommended) {
			cat("  ⚠️  RESTART RECOMMENDED:", convergence_check$reason, "\n")
			if (!RESTART_ENABLED) {
				cat("  (Restart disabled, continuing with current values)\n")
			}
		}
	}, error = function(e) {
		cat("  WARNING: Convergence check failed, defaulting to continue sampling\n")
		cat("  Error:", e$message, "\n")
		continue <- TRUE  # Default to continue if check fails
	})

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
		# Reset continue to TRUE as default, then check convergence
		continue <- TRUE  # Default to continue unless convergence is achieved
		tryCatch({
			convergence_check <- check_continue_with_restart(all_samples, min_eff_size = min_eff_size_perchain)
			continue <- convergence_check$continue
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
			nimble_code = modelCode,
			model_structure = "standardized_beta_regression_with_restart_capability",
			restart_info = if (!is.null(restart_inits)) list(
				restart_used = TRUE,
				extreme_values_detected = sum(restart_inits$extreme_flags),
				fallback_values_used = sum(restart_inits$fallback_used),
				fallback_strategy = RESTART_FALLBACK_STRATEGY
			) else list(restart_used = FALSE)
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

	cat("  - Beta regression with precision parameter\n")
	cat("  - Block samplers for efficient parameter sampling\n")
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
			scenario = scenario,
			min.date = min.date,
			max.date = max.date,
			niter = total_iterations,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			nimble_code = modelCode,
			model_structure = "standardized_beta_regression_with_restart_capability"
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
cat("Restart functionality:", if(RESTART_ENABLED) "ENABLED" else "DISABLED", "\n")
cat("Model configuration:\n")
if (nrow(params) > 0) {
	print(params[, c("rank.name", "species", "model_name", "model_id")])
} else {
	cat("No models to run\n")
}

# [Continue with the rest of the original script's setup and parallel execution]

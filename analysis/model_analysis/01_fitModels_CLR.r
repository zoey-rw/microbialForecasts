#!/usr/bin/env Rscript

# FIXED CLR model fitting script that addresses critical MCMC issues
# - Removes duplicate model code blocks
# - Fixes sampler conflicts by removing defaults before adding specialized samplers
# - Simplifies initialization strategy
# - Standardizes priors across model types
# - Includes proper metadata saving for all 6 model types
# - Uses CLR transformation instead of beta regression

# Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# LOCAL TESTING: Run with 2 chains on 2 cores for local testing
nchains = 2

# Set number of cores for parallel execution (LOCAL TESTING: 2 cores)
n_cores = 2

#### Run on all groups ----

# Load required packages first
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Source the main script with correct path
source(here("source.R"))

# Function to check if MCMC should continue based on effective sample size
check_continue <- function(samples, min_eff_size = 10) {
	cat("    DEBUG: check_continue called with", nrow(samples), "samples\n")
	
	# Convert to mcmc object for effectiveSize calculation
	if (!inherits(samples, "mcmc")) {
		cat("    DEBUG: Converting samples to mcmc object\n")
		samples <- as.mcmc(samples)
	}
	
	# Calculate effective sample sizes for all parameters
	cat("    DEBUG: Calculating effective sample sizes...\n")
	eff_sizes <- effectiveSize(samples)
	cat("    DEBUG: ESS calculated for", length(eff_sizes), "parameters\n")
	
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

	# Focus on legacy covariate models for 2013-2018 period
	params <- params_in %>% ungroup %>% filter(
		# Test multiple ranks for generalization testing
		rank.name %in% c("phylum_fun", "class_fun", "order_fun") &
		# Focus on 2013-2018 period (both legacy and non-legacy scenarios)
		scenario %in% c("Legacy with covariate 2013-2018", "2013-2018") &
		# All three model types
		model_name %in% c("cycl_only", "env_cov", "env_cycl")
	)

# Filter out already converged models
params <- params %>% filter(!model_id %in% converged_list)

cat("Running", nrow(params), "models for index", k, "\n")
cat("Model configuration:\n")
if (nrow(params) > 0) {
	print(params[k, c("rank.name", "species", "model_name", "model_id")])
} else {
	cat("No models to run\n")
}

# Create function that uses our working approach for each model
run_scenarios_fixed <- function(j, chain_no) {
	# Load required libraries in each worker
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)

	cat("=== Starting CLR model fitting ===\n")
	cat("Model index:", j, "Chain:", chain_no, "\n")
	cat("Model parameters:\n")
	print(params[j,])
	cat("=============================\n")

	# Get the group data
	if (is.null(params) || nrow(params) < j) {
		stop("Params data frame not available or index out of bounds")
	}

	rank.name <- params$rank.name[[j]]
	species <- params$species[[j]]
	model_id <- params$model_id[[j]]
	model_name <- params$model_name[[j]]
	min.date <- params$min.date[[j]]
	max.date <- params$max.date[[j]]
	scenario <- params$scenario[[j]]

	# Check if this is a legacy covariate model
	use_legacy_covariate <- grepl("Legacy with covariate", scenario)
	
	# Use pre-loaded data (no need to load again in each worker)
	
	# Get the specific group data
	if (!(rank.name %in% names(all_ranks))) {
		stop("Rank name not found in data")
	}
	rank.df <- all_ranks[[rank.name]]
	
	cat("Preparing CLR model data for", rank.name, "\n")
	
	# Use prepCLRData for CLR transformation
	cat("DEBUG: Preparing CLR data for species", species, "from", rank.name, "\n")
	
	tryCatch({
		# Use prepCLRData with the same parameters as working script
		model.dat <- prepCLRData(rank.df = rank.df,
														 min.prev = 3,
														 min.date = min.date,
														 max.date = max.date,
														 s = species)
		cat("✅ prepCLRData successful\n")
	}, error = function(e) {
		cat("ERROR in prepCLRData:", e$message, "\n")
		cat("Error traceback:\n")
		print(e)
		stop("CLR data preparation failed")
	})

	cat("CLR data prepared successfully\n")
	
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
	
	# Prepare constants - CLR uses different structure than beta regression
	constants <- list()
	
	# Add basic dimensions
	constants$N.core <- model.dat$N.core
	constants$N.plot <- model.dat$N.plot
	constants$N.site <- model.dat$N.site
	constants$N.date <- model.dat$N.date
	
	# Add response variable - CLR uses vector, not matrix
	constants$y <- as.vector(model.dat$y)
	
	# Add indexing vectors
	constants$timepoint <- model.dat$timepoint
	constants$plot_num <- model.dat$plot_num
	constants$plot_site_num <- model.dat$plot_site_num
	
	# Add seasonal predictors with optimal scaling for CLR models
	constants$sin_mo <- model.dat$sin_mo
	constants$cos_mo <- model.dat$cos_mo
	
	# Apply optimal scaling for CLR models
	sin_mo_scaled <- scale(constants$sin_mo, center = FALSE, scale = TRUE)
	cos_mo_scaled <- scale(constants$cos_mo, center = FALSE, scale = TRUE)
	
	constants$sin_mo <- as.numeric(sin_mo_scaled)
	constants$cos_mo <- as.numeric(cos_mo_scaled)
	
	# Add environmental predictors
	if ("temp" %in% names(model.dat)) constants$temp <- model.dat$temp
	if ("mois" %in% names(model.dat)) constants$mois <- model.dat$mois
	if ("pH" %in% names(model.dat)) constants$pH <- model.dat$pH
	if ("pC" %in% names(model.dat)) constants$pC <- model.dat$pC
	if ("relEM" %in% names(model.dat)) constants$relEM <- model.dat$relEM
	if ("LAI" %in% names(model.dat)) constants$LAI <- model.dat$LAI

	# Add legacy covariate if needed
	if (use_legacy_covariate) {
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
	if (model_name == "env_cycl") {
		constants$N.beta = 8
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
	} else {
		constants$N.beta = 2
	}
	
	# Add plot estimate matrices for monitoring (CLR models need these for downstream analysis)
	constants$plot_estimates <- matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
	constants$plot_predictions <- matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
	
	cat("Constants prepared successfully\n")
	
	# STANDARDIZED MODEL DEFINITIONS - All models use consistent priors and CLR structure
	if (model_name == "cycl_only" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# OBSERVATION MODEL - CLR transformation
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[i], sd = 1/sqrt(precision))
				mu[i] <- intercept +
					beta[1] * sin_mo[timepoint[i]] +
					beta[2] * cos_mo[timepoint[i]] +
					legacy_effect * legacy[timepoint[i]] +
					site_effect[plot_site_num[plot_num[i]]]
			}
			
			# DETERMINISTIC NODES - Plot-level estimates for monitoring
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_estimates[p, t] <- intercept +
						beta[1] * sin_mo[t] +
						beta[2] * cos_mo[t] +
						legacy_effect * legacy[t] +
						site_effect[plot_site_num[p]]
					
					plot_predictions[p, t] <- plot_estimates[p, t]  # For CLR, predictions = estimates
				}
			}

			# 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects (matching beta regression)
			site_effect_sd ~ dgamma(2, 20)  # More informative: dgamma(2, 20) - stable, proven
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			precision ~ dgamma(0.001, 0.001)  # Jeffreys prior - very weak, proven
			intercept ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			legacy_effect ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			
			for (b in 1:2) {
				beta[b] ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			}
		})
	} else if (model_name == "env_cycl" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# OBSERVATION MODEL - CLR transformation
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[i], sd = 1/sqrt(precision))
				mu[i] <- intercept +
					beta[1] * temp[plot_site_num[plot_num[i]], timepoint[i]] +
					beta[2] * mois[plot_site_num[plot_num[i]], timepoint[i]] +
					beta[3] * pH[plot_num[i], timepoint[i]] +
					beta[4] * pC[plot_num[i], timepoint[i]] +
					beta[5] * relEM[plot_num[i], timepoint[i]] +
					beta[6] * LAI[plot_site_num[plot_num[i]], timepoint[i]] +
					beta[7] * sin_mo[timepoint[i]] +
					beta[8] * cos_mo[timepoint[i]] +
					legacy_effect * legacy[timepoint[i]] +
					site_effect[plot_site_num[plot_num[i]]]
			}
			
			# DETERMINISTIC NODES - Plot-level estimates for monitoring
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_estimates[p, t] <- intercept +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, t] +
						beta[4] * pC[p, t] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						beta[7] * sin_mo[t] +
						beta[8] * cos_mo[t] +
						legacy_effect * legacy[t] +
						site_effect[plot_site_num[p]]
					
					plot_predictions[p, t] <- plot_estimates[p, t]  # For CLR, predictions = estimates
				}
			}
			
			# 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects (matching beta regression)
			site_effect_sd ~ dgamma(2, 20)  # More informative: dgamma(2, 20) - stable, proven
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			precision ~ dgamma(0.001, 0.001)  # Jeffreys prior - very weak, proven
			intercept ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			legacy_effect ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			
			for (b in 1:8) {
				beta[b] ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			}
		})
	} else if (model_name == "env_cov" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# OBSERVATION MODEL - CLR transformation
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[i], sd = 1/sqrt(precision))
				mu[i] <- intercept +
					beta[1] * temp[plot_site_num[plot_num[i]], timepoint[i]] +
					beta[2] * mois[plot_site_num[plot_num[i]], timepoint[i]] +
					beta[3] * pH[plot_num[i], timepoint[i]] +
					beta[4] * pC[plot_num[i], timepoint[i]] +
					beta[5] * relEM[plot_num[i], timepoint[i]] +
					beta[6] * LAI[plot_site_num[plot_num[i]], timepoint[i]] +
					legacy_effect * legacy[timepoint[i]] +
					site_effect[plot_site_num[plot_num[i]]]
			}
			
			# DETERMINISTIC NODES - Plot-level estimates for monitoring
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_estimates[p, t] <- intercept +
						beta[1] * temp[plot_site_num[p], t] +
						beta[2] * mois[plot_site_num[p], t] +
						beta[3] * pH[p, t] +
						beta[4] * pC[p, t] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						legacy_effect * legacy[t] +
						site_effect[plot_site_num[p]]
					
					plot_predictions[p, t] <- plot_estimates[p, t]  # For CLR, predictions = estimates
				}
			}

			# 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects (matching beta regression)
			site_effect_sd ~ dgamma(2, 20)  # More informative: dgamma(2, 20) - stable, proven
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}

			precision ~ dgamma(0.001, 0.001)  # Jeffreys prior - very weak, proven
			intercept ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			legacy_effect ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven

			for (b in 1:6) {
				beta[b] ~ dnorm(0, sd = 10)  # Very wide normal - very weak, proven
			}
		})
	} else {
		stop("Unsupported model combination: ", model_name, " with use_legacy_covariate=", use_legacy_covariate, 
		     ". Only models with legacy covariates are supported.")
	}
	
	# Create inits - CLR specific initialization
	cat("Creating CLR-specific initial values...\n")
	y_data <- constants$y  # Get the CLR-transformed data
	y_mean <- mean(y_data, na.rm = TRUE)
	y_sd <- sd(y_data, na.rm = TRUE)
	
	# FIXED: Simplified initialization strategy - spread chains out properly
	set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model
	
	# Initialize parameters with proper spacing between chains (matching beta regression approach)
	inits <- list(
		intercept = rnorm(1, y_mean, 0.5),
		precision = 50,  # Start with moderate precision (matching beta regression)
		beta = rnorm(constants$N.beta, 0, 0.1)  # Start beta parameters close to zero
	)
	
	# Initialize hierarchical parameters (matching beta regression approach)
	inits$site_effect_sd <- 0.5  # Start with moderate site effect SD (matching beta regression)
	
	# Initialize legacy effect if using legacy covariate
	if (use_legacy_covariate) {
		inits$legacy_effect <- 0  # Start legacy effect at 0 (matching beta regression)
	}
	
	# Add site effects starting near zero (matching beta regression approach)
	inits$site_effect <- rnorm(constants$N.site, 0, 0.1)
	
	# Initialize plot estimate matrices for monitoring
	inits$plot_estimates <- matrix(rnorm(constants$N.plot * constants$N.date, y_mean, y_sd * 0.1), 
	                               nrow = constants$N.plot, ncol = constants$N.date)
	inits$plot_predictions <- inits$plot_estimates  # For CLR, predictions = estimates
	
	cat("Model built successfully\n")
	
	# Build model
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y=constants$y), inits = inits)
	
	# Debug: Check data dimensions
	cat("Data dimensions check:\n")
	cat("  constants$y dimensions:", length(constants$y), "\n")
	cat("  constants$y class:", class(constants$y), "\n")
	cat("  constants$N.core:", constants$N.core, "\n")
	cat("  constants$N.beta:", constants$N.beta, "\n")
	
	# Compile model
	cModel <- compileNimble(Rmodel)
	
	cat("Model compiled successfully\n")
	
	# FIXED: Configure MCMC with proper sampler management - all models now use precision parameter
# All models now use CLR regression with precision parameter
monitors <- c("beta","precision","site_effect","site_effect_sd","intercept")

if (use_legacy_covariate) {
    monitors <- c(monitors, "legacy_effect")
}

# Use monitors2 for latent variables at different interval for efficiency (matching beta regression approach)
monitored_latent_params <- c(
    # Monitor latent process variables at different interval for efficiency
    "plot_estimates", "plot_predictions"
)

mcmcConf <- configureMCMC(
    model = cModel, 
    monitors = monitors, 
    monitors2 = monitored_latent_params,  # Use monitors2 for latent variables
    thin = 1,
    thin2 = 20,  # Sample latent variables every 20th iteration for efficiency
    useConjugacy = FALSE
)
	
	# FIXED: Remove default samplers FIRST before adding specialized ones
	cat("  Removing default samplers...\n")
	
	# Remove ALL default samplers to prevent conflicts - all models use precision now
	mcmcConf$removeSamplers(c("precision", "site_effect_sd", "intercept"))
	
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
	
	# FIXED: Add specialized samplers AFTER removing defaults (no conflicts)
	cat("  Adding specialized samplers...\n")
	
	# All models now use precision parameter
	mcmcConf$addSampler(target = "precision", type = "slice")
	cat("  Added slice sampler for precision\n")
	
	# site_effect_sd: Site effect scale (all models)
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
	cat("  Added slice sampler for site_effect_sd\n")
	
	# intercept: Use slice sampler for better location parameter exploration (all models)
	mcmcConf$addSampler(target = "intercept", type = "slice")
	cat("  Added slice sampler for intercept\n")
	
	# legacy_effect: Use slice sampler for better mixing of legacy parameter
	if (use_legacy_covariate) {
		mcmcConf$addSampler(target = "legacy_effect", type = "slice")
		cat("  Added slice sampler for legacy_effect\n")
	}
	
	# FIXED: Use EITHER block OR individual samplers for site effects, NOT BOTH
	if (constants$N.site > 1) {
		# Use block sampling for correlated site effects
		mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, "]"), type = "AF_slice")
		cat("  Added block sampler for site_effect[1:", constants$N.site, "]\n")
	} else {
		# Individual sampler for single site
		mcmcConf$addSampler(target = "site_effect[1]", type = "slice")
		cat("  Added slice sampler for site_effect[1]\n")
	}
	
	# FIXED: Use block sampler for beta parameters when multiple exist, individual for single
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
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# PRODUCTION: Run MCMC with convergence-based sampling until production ESS is reached
	burnin <- 2000   # Adequate burnin for complex hierarchical models
	thin <- 1        # No thinning to preserve all samples
	iter_per_chunk <- 2000   # Larger chunks for efficiency
	init_iter <- 2000 # Initial iterations for adaptation and burnin
	min_eff_size_perchain <- 100  # Production standard ESS per chain
	max_loops <- 50  # Allow adequate sampling for convergence
	max_save_size <- 100000  # Handle larger sample sizes
	min_total_iterations <- 10000  # Minimum iterations needed for complex models
	
	cat("Running PRODUCTION MCMC with convergence-based sampling\n")
	cat("  Burn-in iterations:", burnin, "\n")
	cat("  Initial iterations:", init_iter, "\n")
	cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
	cat("  Target ESS per chain:", min_eff_size_perchain, " (PRODUCTION STANDARD)\n")
	cat("  Minimum total iterations:", min_total_iterations, "\n")
	
	# Run burn-in iterations for proper mixing
	cat("  Running burn-in iterations (", burnin, " iterations) for proper mixing...\n")
	compiled$run(niter = burnin, thin = thin, nburnin = 0)
	cat("  Burn-in completed\n")
	
	# Run initial iterations with progress reporting and adaptation
	cat("  Running initial iterations (", init_iter, " iterations) for adaptation...\n")
	compiled$run(niter = init_iter, thin = thin, nburnin = 0)
	cat("  Initial iterations completed\n")
	
	# Get initial samples and check convergence
	initial_samples <- as.matrix(compiled$mvSamples)
	cat("  Initial samples collected, checking convergence...\n")
	
	# Check if we need to continue sampling for convergence
	cat("  Checking initial convergence with", nrow(initial_samples), "samples...\n")
	
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
	
	cat("  Initial convergence check result:", ifelse(continue, "CONTINUE", "CONVERGED"), "\n")
	
	# Create output directories for checkpoints
	model_output_dir <- file.path(getwd(), "data", "model_outputs", "CLR_regression", model_name)
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create model_id for consistent naming with legacy covariate indicator
	legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
	model_id <- paste(model_name, species, min.date, max.date, legacy_indicator, sep = "_")
	
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
	
	# Create the complete chain structure with enhanced metadata
	chain_output <- list(
		samples = samples,
		samples2 = samples,  # CLR models use same samples for both parameter and plot predictions
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
			model_structure = "standardized_CLR_regression_with_consistent_priors"  # Model structure identifier
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
	
	cat("Saved MCMC samples to:", samples_file, "\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== CLR model fitting completed successfully ===\n")
	cat("  - PRODUCTION STANDARD: ESS ≥ 100 for all parameters\n")
	cat("  - CLR regression with unified precision parameter\n")
	cat("  - Proper legacy covariate handling (present when needed, absent when not)\n")
	cat("  - Block samplers for efficient parameter sampling\n")
	cat("  - PRODUCTION CONVERGENCE: Adaptive sampling until ESS ≥ 100 reached\n")
	cat("  - Production-ready: Adequate burn-in, proper mixing, reliable results\n")
	
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
			model_structure = "standardized_CLR_regression_with_consistent_priors"
		)
	))
}

# Run multiple models for testing
cat("Testing with", nrow(params), "models\n")
cat("Models to test:\n")
print(params[, c("rank.name", "species", "model_name", "model_id")])

# Test with more models to verify generalization across ranks and model types
test_models <- 6  # Test 6 models for comprehensive verification
cat("\nTesting", test_models, "models to verify generalization across ranks and model types\n")

# Set up parallel cluster for LOCAL TESTING
library(parallel)
library(doParallel)

# Create cluster with appropriate number of cores for LOCAL TESTING
ncores <- min(nchains, n_cores)  # Use the explicitly set n_cores for local testing
cl <- makeCluster(ncores, type = "PSOCK")

# Register the cluster
registerDoParallel(cl)

# Export necessary objects to workers
clusterExport(cl, c("params", "k", "nchains"))

# Load data once before parallel execution
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Define valid_models for testing (subset of params for local testing)
test_models <- 6  # Test with 6 models for local verification
valid_models <- params[1:test_models, ]
cat("Testing", test_models, "models to verify generalization across ranks and model types\n")

# Export the function and data to workers
clusterExport(cl, c("run_scenarios_fixed", "params", "all_ranks", "check_continue"))

# Create a combined task list: (model_idx, chain_no) pairs
all_tasks <- expand.grid(model_idx = 1:test_models, chain_no = 1:nchains)
cat("Total parallel tasks:", nrow(all_tasks), "(", test_models, "models ×", nchains, "chains)\n")
cat("Task details:\n")
print(all_tasks)
cat("Cluster size:", n_cores, "workers\n")

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

# Create a function that saves results as they complete
runAndSave_task <- function(task_idx) {
    cat("DEBUG: runAndSave_task called with task_idx =", task_idx, "at", Sys.time(), "\n")
    
    # Test if we can access the function
    if (!exists("run_scenarios_fixed")) {
        cat("ERROR: run_scenarios_fixed function not found in worker\n")
        return(list(error = "run_scenarios_fixed function not found"))
    }
    
    cat("DEBUG: run_scenarios_fixed function exists in worker\n")
    
    # Initialize error tracking
    error_details <- list()
    start_time <- Sys.time()
    
    tryCatch({
        # Get task details
        task <- all_tasks[task_idx, ]
        model_idx <- task$model_idx
        chain_no <- task$chain_no
        
        cat("DEBUG: Worker processing Model", model_idx, "Chain", chain_no, "starting at", format(start_time), "\n")
        
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
        
        # Check if valid_models is available
        if (!exists("valid_models") || is.null(valid_models) || nrow(valid_models) == 0) {
            stop("Parameters 'valid_models' not available or empty in worker environment")
        }
        
        # Validate model index
        if (model_idx > nrow(valid_models)) {
            stop("Model index ", model_idx, " exceeds available models (", nrow(valid_models), ")")
        }
        
        cat("Worker: All checks passed, calling run_scenarios_fixed...\n")
        
        # Run the model with detailed error context
        cat("Worker: About to call run_scenarios_fixed with j =", model_idx, "chain_no =", chain_no, "\n")
        result <- tryCatch({
            run_scenarios_fixed(j = model_idx, chain_no = chain_no)
        }, error = function(e) {
            cat("Worker: ERROR in run_scenarios_fixed:", e$message, "\n")
            cat("Worker: Error call:", if(!is.null(e$call)) paste(deparse(e$call), collapse=" ") else "No call info", "\n")
            stop(e)
        })
        
        cat("Worker: run_scenarios_fixed completed successfully\n")
        cat("Worker: Result class:", class(result), "\n")
        cat("Worker: Result names:", paste(names(result), collapse=", "), "\n")
        if ("status" %in% names(result)) {
            cat("Worker: Result status:", result$status, "\n")
        }
        
        # Validate result structure
        if (!is.list(result) || !("status" %in% names(result))) {
            stop("Invalid result structure from run_scenarios_fixed")
        }
        
        # Save result immediately if successful
        if (result$status == "SUCCESS") {
            # Create output directory early with HPC compatibility
            possible_bases <- c(
                here("data", "model_outputs"),
                file.path(here(), "data", "model_outputs"),
                file.path(Sys.getenv("HOME"), "data", "model_outputs"),
                file.path(here(), "data", "model_outputs")  # Duplicate for consistency
            )
            
            model_output_dir <- NULL
            for (base_dir in possible_bases) {
                if (!is.null(base_dir) && base_dir != "" && base_dir != "NULL") {
                    test_dir <- file.path(base_dir, "CLR_regression", valid_models$model_name[model_idx])
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
                # Fallback: create using here() for consistency
                fallback_output_dir <- here("data", "model_outputs", "CLR_regression", valid_models$model_name[model_idx])
                dir.create(fallback_output_dir, showWarnings = FALSE, recursive = TRUE)
                cat("  WARNING: Using fallback output directory:", fallback_output_dir, "\n")
                model_output_dir <- fallback_output_dir
            }
            
            # Create model_id for consistent naming using helper function
            model_id <- create_model_id(
                valid_models$model_name[model_idx], 
                valid_models$species[model_idx],
                valid_models$min.date[model_idx], 
                valid_models$max.date[model_idx],
                grepl("Legacy with covariate", valid_models$scenario[model_idx]),
                "clr"
            )
            
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
                    rank.name = valid_models$rank.name[model_idx],
                    species = valid_models$species[model_idx],
                    model_name = valid_models$model_name[model_idx],
                    model_id = model_id,
                    use_legacy_covariate = grepl("Legacy with covariate", valid_models$scenario[model_idx]),
                    scenario = valid_models$scenario[model_idx],
                    min.date = valid_models$min.date[model_idx],
                    max.date = valid_models$max.date[model_idx],
                    niter = nrow(result$samples),
                    nburnin = 500,  # Default burnin
                    thin = 1,       # Default thin
                    task_idx = task_idx,
                    completed_at = Sys.time(),
                    model_data = result$model_data,  # Include model_data from result
                    nimble_code = result$nimble_code,  # Include nimble_code if available
                    model_structure = "standardized_CLR_regression_with_consistent_priors"  # Model structure identifier
                )
            }
            
            chain_output <- list(
                samples = result$samples,
                samples2 = result$samples,  # CLR models use same samples for both parameter and plot predictions
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
        error_file <- paste0("chain_", task$model_idx, "_task_", task_idx, "_ERROR.txt")
        
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

# Set start time for runtime calculation
start_time <- Sys.time()

# Define the global model_output_dir variable
model_output_dir <- here("data", "model_outputs", "CLR_regression")

# Create cluster for parallel execution - LOCAL TESTING with 2 cores
cat("Creating cluster with", n_cores, "cores for", nrow(valid_models), "models ×", nchains, "chains\n")
cat("LOCAL TESTING: 2 cores allocated for local testing\n")

# Create cluster with explicit error handling
tryCatch({
    cl <- makeCluster(n_cores, type = "PSOCK")
    cat("✓ Cluster created successfully with", length(cl), "workers\n")
    
    # Test cluster functionality
    cat("Testing cluster functionality...\n")
    test_result <- tryCatch({
        clusterEvalQ(cl, Sys.getpid())
    }, error = function(e) {
        cat("✗ Cluster test failed:", e$message, "\n")
        return(NULL)
    })
    
    if (!is.null(test_result)) {
        cat("✓ Cluster test successful. Worker PIDs:", unlist(test_result), "\n")
    } else {
        cat("✗ Cluster test failed\n")
        cl <- NULL
    }
    
}, error = function(e) {
    cat("✗ ERROR creating cluster:", e$message, "\n")
    cat("Falling back to sequential execution...\n")
    cl <- NULL
})

# Register parallel backend
cat("Registering parallel backend...\n")
registerDoParallel(cl)
cat("Parallel backend registered. Checking registration...\n")
cat("getDoParWorkers():", getDoParWorkers(), "\n")
cat("getDoParName():", getDoParName(), "\n")

# Export ALL necessary variables and functions to workers
cat("Exporting all necessary variables to workers...\n")
clusterExport(cl, c(
    "runAndSave_task", 
    "run_scenarios_fixed", 
    "valid_models", 
    "all_tasks", 
    "model_output_dir",
    "all_ranks",
    "params",
    "n_cores",
    "nchains",
    # Export helper functions
    "load_required_packages",
    "create_directories_safe",
    "create_model_id",
    "save_checkpoint_safe",
    "create_progress_file",
    "update_progress_file"
))

# Also export the here function explicitly
clusterExport(cl, "here")

# Run everything in parallel with incremental saving
cat("Starting foreach loop with", nrow(all_tasks), "tasks\n")

# Quick validation that workers can access functions and data
test_export <- tryCatch({
    clusterEvalQ(cl, exists("runAndSave_task"))
}, error = function(e) {
    cat("✗ ERROR testing function export:", e$message, "\n")
    return(rep(FALSE, length(cl)))
})

test_data <- tryCatch({
    clusterEvalQ(cl, exists("valid_models"))
}, error = function(e) {
    cat("✗ ERROR testing data export:", e$message, "\n")
    return(rep(FALSE, length(cl)))
})

cat("✓ Function export check:", all(unlist(test_export)), "\n")
cat("✓ Data export check:", all(unlist(test_data)), "\n")

all_results_parallel = foreach(task_idx = 1:nrow(all_tasks), 
                               .packages = c()) %dopar% {  # Packages loaded by load_required_packages() in worker
                                   cat("DEBUG: Worker starting task", task_idx, "at", Sys.time(), "\n")
                                   
                                   # Load required packages in worker
                                   load_required_packages()
                                   
                                   # Test if we can access the function
                                   if (!exists("runAndSave_task")) {
                                       cat("ERROR: runAndSave_task function not found in worker\n")
                                       return(list(error = "Function not found"))
                                   }
                                   
                                   # Test if we can access the data
                                   if (!exists("valid_models")) {
                                       cat("ERROR: valid_models data not found in worker\n")
                                       return(list(error = "Data not found"))
                                   }
                                   
                                   # Execute the function with error handling
                                   result <- tryCatch({
                                       runAndSave_task(task_idx)
                                   }, error = function(e) {
                                       cat("ERROR in runAndSave_task:", e$message, "\n")
                                       return(list(error = e$message))
                                   })
                                   
                                   cat("DEBUG: Worker completed task", task_idx, "at", Sys.time(), "\n")
                                   
                                   # Extract the actual status from the nested result structure
                                   actual_status <- "UNKNOWN"
                                   if ("result" %in% names(result) && "status" %in% names(result$result)) {
                                       actual_status <- result$result$status
                                   } else if ("status" %in% names(result)) {
                                       actual_status <- result$status
                                   }
                                   
                                   cat("DEBUG: Result status:", actual_status, "\n")
                                   cat("DEBUG: Result structure:", paste(names(result), collapse=", "), "\n")
                                   if ("result" %in% names(result)) {
                                       cat("DEBUG: Nested result structure:", paste(names(result$result), collapse=", "), "\n")
                                   }
                                   
                                   return(result)
                               }

cat("Foreach loop completed. Results length:", length(all_results_parallel), "\n")
cat("Parallel execution completed at:", format(Sys.time()), "\n")

# Stop cluster
stopCluster(cl)

# Show progress summary
cat("\n=== PROGRESS SUMMARY ===\n")
cat("Checking which chains have been completed...\n")

# Count completed chains from parallel results
completed_chains <- 0
error_chains <- 0
for (i in 1:length(all_results_parallel)) {
    result <- all_results_parallel[[i]]
    if ("error" %in% names(result)) {
        error_chains <- error_chains + 1
        cat("✗ Task", i, "failed with error:", result$error, "\n")
    } else if ("result" %in% names(result) && "status" %in% names(result$result)) {
        if (result$result$status == "SUCCESS") {
            completed_chains <- completed_chains + 1
            cat("✓ Model", result$model_idx, "Chain", result$chain_no, "completed\n")
        } else {
            error_chains <- error_chains + 1
            cat("✗ Model", result$model_idx, "Chain", result$chain_no, "failed with status:", result$result$status, "\n")
        }
    } else if ("status" %in% names(result)) {
        if (result$status == "SUCCESS") {
            completed_chains <- completed_chains + 1
            cat("✓ Model", result$model_idx, "Chain", result$chain_no, "completed\n")
        } else {
            error_chains <- error_chains + 1
            cat("✗ Model", result$model_idx, "Chain", result$chain_no, "failed with status:", result$status, "\n")
        }
    } else {
        error_chains <- error_chains + 1
        cat("? Task", i, "has unknown result structure\n")
        cat("  Available names:", paste(names(result), collapse=", "), "\n")
        if ("result" %in% names(result)) {
            cat("  Nested result names:", paste(names(result$result), collapse=", "), "\n")
        }
    }
}

cat("\nProgress Summary:\n")
cat("  Completed chains:", completed_chains, "/", nrow(valid_models) * nchains, "\n")
cat("  Failed chains:", error_chains, "/", nrow(valid_models) * nchains, "\n")
cat("  Success rate:", round(completed_chains / (nrow(valid_models) * nchains) * 100, 1), "%\n")

# Reorganize results by model
all_results <- list()
for (model_idx in 1:nrow(valid_models)) {
    all_results[[model_idx]] <- list()
    for (chain_no in 1:nchains) {
        # Find the result for this model/chain combination
        task_idx <- which(all_tasks$model_idx == model_idx & all_tasks$chain_no == chain_no)
        if (length(task_idx) > 0) {
            task_result <- all_results_parallel[[task_idx]]
            if ("result" %in% names(task_result)) {
                all_results[[model_idx]][[chain_no]] <- task_result$result
            } else {
                all_results[[model_idx]][[chain_no]] <- task_result
            }
        } else {
            all_results[[model_idx]][[chain_no]] <- list(status = "NOT_FOUND")
        }
    }
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("ALL CLR MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Summary of all models
cat("\nSummary of All CLR Models:\n")
for (model_idx in 1:nrow(valid_models)) {
    cat("\nModel", model_idx, ":", valid_models$species[model_idx], "(", valid_models$model_name[model_idx], ")\n")
    
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
                if ("samples" %in% names(output.list[[i]]) && !is.null(output.list[[i]]$samples)) {
                    sample_dim <- dim(output.list[[i]]$samples)
                    if (length(sample_dim) >= 1) {
                        cat("    Chain", i, ": SUCCESS - Samples:", sample_dim[1], "iterations\n")
                    } else {
                        cat("    Chain", i, ": SUCCESS - Samples: unknown dimensions\n")
                    }
                } else {
                    cat("    Chain", i, ": SUCCESS - No samples data\n")
                }
            } else {
                error_msg <- if ("error" %in% names(output.list[[i]])) output.list[[i]]$error else "Unknown error"
                cat("    Chain", i, ": ERROR -", error_msg, "\n")
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

cat("✓ Proper CLR transformation for compositional data\n")
cat("✓ Precision parameter instead of separate core_sd/sigma parameters\n")
cat("✓ 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects (matching beta regression)\n")
cat("  - precision ~ dgamma(0.001, 0.001) - Jeffreys prior - very weak, proven\n")
cat("  - intercept ~ dnorm(0, sd = 10) - Very wide normal - very weak, proven\n")
cat("  - legacy_effect ~ dnorm(0, sd = 10) - Very wide normal - very weak, proven\n")
cat("  - site_effect_sd ~ dgamma(2, 20) - More informative: dgamma(2, 20) - stable, proven\n")
cat("  - beta ~ dnorm(0, sd = 10) - Very wide normal - very weak, proven\n")
cat("✓ Plot estimate monitoring (plot_estimates, plot_predictions) for comprehensive analysis\n")
cat("✓ Individual slice samplers for beta parameters\n")
cat("✓ Convergence-based sampling with iterative saving\n")
cat("✓ LOCAL TESTING: 2 cores, 2 chains for faster testing\n")
cat("Check output files in:", here("data", "model_outputs", "CLR_regression"), "/[model_name]/\n")



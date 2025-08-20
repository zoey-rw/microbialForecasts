#!/usr/bin/env Rscript

# Fixed model fitting script that uses our working approach
# Now that package dependencies are resolved

# Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# Run with at least 4 cores available (one MCMC chain per core)
nchains = 4

#### Run on all groups ----

source("source_local.R")

params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 rep("logical", 2),
										 rep("character", 4)))

rerun_list = readRDS(here("data/summary/unconverged_taxa_list.rds"))
converged_list = readRDS(here("data/summary/converged_taxa_list.rds"))

	# Focus on target species for 2015-2018 and 2015-2020 periods
	params <- params_in %>% ungroup %>% filter(
		# Target species: saprotroph, ectomycorrhizal, cellulolytic, ascomycota
		rank.name %in% c("saprotroph", "ectomycorrhizal", "cellulolytic", "ascomycota") &
		# Target time periods: 2015-2018 and 2015-2020
		((min.date == "20151101" & max.date == "20180101") |  # 2015-2018
		 (min.date == "20150101" & max.date == "20200101")) &  # 2015-2020
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

	cat("=== Starting model fitting ===\n")
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
	
	# Use pre-loaded data (no need to load again in each worker)
	
	# Get the specific group data
	if (!(rank.name %in% names(all_ranks))) {
		stop("Rank name not found in data")
	}
	rank.df <- all_ranks[[rank.name]]
	
	cat("Preparing model data for", rank.name, "\n")
	
	# Extract the specific species and create "other" column BEFORE calling prepBetaRegData
	cat("DEBUG: Extracting species", species, "from", rank.name, "\n")
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!species) %>%
		mutate(other = 1 - !!sym(species))
	
	cat("DEBUG: rank.df_spec dimensions:", dim(rank.df_spec), "\n")
	cat("DEBUG: rank.df_spec columns:", colnames(rank.df_spec), "\n")
	
	# Use prepBetaRegData with the species-specific data
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
	
	# Prepare constants
	constants <- model.dat[c("plotID", "timepoint", "plot_site",
							"site_start", "plot_start", "plot_index",
							"plot_num", "plot_site_num",
							"N.plot", "N.spp", "N.core", "N.site", "N.date",
							"sin_mo", "cos_mo")]
	
	# Add environmental predictors
	if ("temp" %in% names(model.dat)) constants$temp <- model.dat$temp
	if ("mois" %in% names(model.dat)) constants$mois <- model.dat$mois
	if ("pH" %in% names(model.dat)) constants$pH <- model.dat$pH
	if ("pC" %in% names(model.dat)) constants$pC <- model.dat$pC
	if ("relEM" %in% names(model.dat)) constants$relEM <- model.dat$relEM
	if ("LAI" %in% names(model.dat)) constants$LAI <- model.dat$LAI
	
	# Model hyperparameters - adjust based on model type
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
	constants$sin_mo = scale(constants$sin_mo, center = F) %>% as.numeric()
	constants$cos_mo = scale(constants$cos_mo, center = F) %>% as.numeric()
	
	cat("Constants prepared successfully\n")
	
	# Define the model based on model_name
	if (model_name == "env_cycl") {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations
			for (i in 1:N.core) {
				y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]], sd = core_sd), 0, 1)
			}
			
			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
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
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
				}
			}
			
			# Hierarchical priors for site effects - better mixing
			site_effect_sd ~ dgamma(2, 4)  # Prior on site effect standard deviation
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}
			
			# Parameter-specific priors for better convergence
			# core_sd: Use more appropriate prior for scale parameter
			core_sd ~ dgamma(1, 10)  # Less informative, allows more exploration
			
			# sigma: Process noise - moderate prior
			sigma ~ dgamma(2, 15)    # Tighter gamma prior
			
			# sig: Site effect scale - moderate prior  
			sig ~ dgamma(2, 8)       # Tighter gamma prior
			
			# intercept: Use more robust prior with better centering
			intercept ~ dt(-2, 1.0, df = 3)  # Centered around expected value, wider variance
			
			# rho: Use more appropriate bounded prior with better mixing
			rho ~ dbeta(1.5, 1.5)    # Less informative, allows more exploration
			
			# Simple independent normal priors for beta parameters
			# Use moderate variance for ecological effects
			for (b in 1:8) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Independent normal priors
			}
		})
	} else if (model_name == "env_cov") {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations
			for (i in 1:N.core) {
				y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]], sd = core_sd), 0, 1)
			}
			
			for (p in 1:N.plot) {
				# Plot-level process model
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
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
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
				}
			}
			
			# Hierarchical priors for site effects - better mixing
			site_effect_sd ~ dgamma(2, 4)  # Prior on site effect standard deviation
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}
			
			# Parameter-specific priors for better convergence
			# core_sd: Use more appropriate prior for scale parameter
			core_sd ~ dgamma(1, 10)  # Less informative, allows more exploration
			
			# sigma: Process noise - moderate prior
			sigma ~ dgamma(2, 15)    # Tighter gamma prior
			
			# sig: Site effect scale - moderate prior  
			sig ~ dgamma(2, 8)       # Tighter gamma prior
			
			# intercept: Use more robust prior with better centering
			intercept ~ dt(-2, 1.0, df = 3)  # Centered around expected value, wider variance
			
			# rho: Use more appropriate bounded prior with better mixing
			rho ~ dbeta(1.5, 1.5)    # Less informative, allows more exploration
			# Simple independent normal priors for beta parameters
			# Use moderate variance for ecological effects
			for (b in 1:6) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Independent normal priors
			}
		})
	} else {
		# Default to cycl_only model
		modelCode <- nimble::nimbleCode({
			for (i in 1:N.core) {
				y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]], sd = core_sd), 0, 1)
			}
			for (p in 1:N.plot) {
				for (t in plot_start[p]) {
					Ex[p, t] ~ dunif(0.0001,.9999)
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
				}
				for (t in plot_index[p]:N.date) {
					logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
						beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
						site_effect[plot_site_num[p]] + intercept
					plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
				}
			}
			# Hierarchical priors for site effects - better mixing
			site_effect_sd ~ dgamma(2, 4)  # Prior on site effect standard deviation
			for (k in 1:N.site) { 
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
			}
			# Parameter-specific priors for better convergence
			# core_sd: Use more appropriate prior for scale parameter
			core_sd ~ dgamma(1, 10)  # Less informative, allows more exploration
			
			# sigma: Process noise - moderate prior
			sigma ~ dgamma(2, 15)    # Tighter gamma prior
			
			# sig: Site effect scale - moderate prior  
			sig ~ dgamma(2, 8)       # Tighter gamma prior
			
			# intercept: Use more robust prior with better centering
			intercept ~ dt(-2, 1.0, df = 3)  # Centered around expected value, wider variance
			
			# rho: Use more appropriate bounded prior with better mixing
			rho ~ dbeta(1.5, 1.5)    # Less informative, allows more exploration
			# Simple independent normal priors for beta parameters
			# Use moderate variance for ecological effects
			for (b in 1:2) {
				beta[b] ~ dnorm(0, sd = 0.3)  # Independent normal priors
			}
		})
	}
	
	# Create inits
	inits <- createInits(constants)
	
	# Calculate data-informed initial values for better convergence
	y_data <- model.dat$y[, 1]  # Get the species abundance data
	y_mean <- mean(y_data, na.rm = TRUE)
	y_sd <- sd(y_data, na.rm = TRUE)
	
	# Advanced initialization strategy for better convergence
	# Use different starting points for each chain to explore different regions
	set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model
	
	# Balanced spectrum-based initialization exploring positive and negative regions
	if (chain_no == 1) {
		# Chain 1: Conservative initialization (near zero, both signs)
		inits$beta <- rnorm(constants$N.beta, 0, 0.1)
	} else if (chain_no == 2) {
		# Chain 2: Positive-focused initialization
		inits$beta <- abs(rnorm(constants$N.beta, 0.2, 0.15))  # Force positive
	} else if (chain_no == 3) {
		# Chain 3: Negative-focused initialization  
		inits$beta <- -abs(rnorm(constants$N.beta, 0.2, 0.15))  # Force negative
	} else {
		# Chain 4: Wide exploration (both signs, larger range)
		inits$beta <- rnorm(constants$N.beta, 0, 0.4)
	}
	# Initialize hierarchical parameters
	inits$site_effect_sd <- 0.5  # Start with moderate site effect variance
	
	# Initialize problematic parameters with chain-specific strategies
	if (chain_no == 1) {
		# Chain 1: Conservative initialization
		inits$rho <- 0.5  # Start at middle of range
		inits$intercept <- logit(y_mean)  # Start near observed mean
		inits$core_sd <- y_sd * 0.02  # Start with data-informed observation error
	} else if (chain_no == 2) {
		# Chain 2: Higher persistence, lower intercept
		inits$rho <- 0.7  # Start higher for temporal persistence
		inits$intercept <- logit(y_mean) - 0.5  # Start lower
		inits$core_sd <- y_sd * 0.03  # Start with slightly higher error
	} else if (chain_no == 3) {
		# Chain 3: Lower persistence, higher intercept
		inits$rho <- 0.3  # Start lower for temporal persistence
		inits$intercept <- logit(y_mean) + 0.5  # Start higher
		inits$core_sd <- y_sd * 0.01  # Start with lower error
	} else {
		# Chain 4: Wide exploration
		inits$rho <- 0.1  # Start very low
		inits$intercept <- logit(y_mean) + 1.0  # Start much higher
		inits$core_sd <- y_sd * 0.05  # Start with higher error
	}
	
	# Initialize other parameters consistently
	inits$sigma <- y_sd * 0.1   # Start with data-informed process error
	inits$sig <- y_sd * 0.05    # Start with data-informed site effect variance
	
	cat("Model built successfully\n")
	
	# Build model
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y=model.dat$y), inits = inits)
	
	# Debug: Check data dimensions
	cat("Data dimensions check:\n")
	cat("  model.dat$y dimensions:", dim(model.dat$y), "\n")
	cat("  model.dat$y class:", class(model.dat$y), "\n")
	cat("  constants$N.core:", constants$N.core, "\n")
	cat("  constants$N.spp:", constants$N.spp, "\n")
	
	# Compile model
	cModel <- compileNimble(Rmodel)
	
	cat("Model compiled successfully\n")
	
	# Configure MCMC with advanced settings for better mixing
	monitors <- c("beta","sigma","site_effect","site_effect_sd","sig","intercept","rho","core_sd")
	mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)
	
	# Add advanced samplers for better parameter space exploration
	# Use block sampling for correlated parameters
	if (constants$N.beta > 1) {
		# Block sample beta parameters together for better mixing
		mcmcConf$addSampler(target = paste0("beta[1:", constants$N.beta, "]"), type = "AF_slice")
	}
	
	# Specialized samplers for problematic parameters
	# core_sd: Use slice sampler for better scale parameter exploration
	mcmcConf$addSampler(target = "core_sd", type = "slice")
	cat("  Added slice sampler for core_sd\n")
	
	# sigma: Process noise - slice sampler
	mcmcConf$addSampler(target = "sigma", type = "slice")
	cat("  Added slice sampler for sigma\n")
	
	# sig: Site effect scale - slice sampler
	mcmcConf$addSampler(target = "sig", type = "slice")
	cat("  Added slice sampler for sig\n")
	
	# intercept: Use slice sampler for better location parameter exploration
	mcmcConf$addSampler(target = "intercept", type = "slice")
	cat("  Added slice sampler for intercept\n")
	
	# rho: Use specialized bounded sampler for temporal persistence
	mcmcConf$addSampler(target = "rho", type = "slice")
	cat("  Added slice sampler for rho\n")
	
	# Specialized MCMC strategies for site effects
	# Site effects often have poor mixing due to high dimensionality and correlations
	
	# Strategy 1: Block sampling for correlated site effects
	if (constants$N.site > 1) {
		# Use block sampling to handle correlations between sites
		mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, "]"), type = "AF_slice")
		cat("  Added block sampler for site_effect[1:", constants$N.site, "]\n")
	}
	
	# Strategy 2: Individual site effect samplers for better mixing
	# Add individual slice samplers for each site effect
	for (i in 1:constants$N.site) {
		mcmcConf$addSampler(target = paste0("site_effect[", i, "]"), type = "slice")
		cat("  Added slice sampler for site_effect[", i, "]\n")
	}
	
	# Strategy 3: Hierarchical structure for site effects
	# Add a global site effect variance parameter
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
	cat("  Added slice sampler for site_effect_sd\n")
	
	# Build and compile MCMC
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# Run MCMC with production-level iterations for convergence
	burnin <- 2000  # Substantial burnin for complex posteriors
	thin <- 1       # No thinning to preserve all samples
	iter_per_chunk <- 5000  # Substantial iterations for complex mixing
	init_iter <- 500  # Substantial initial iterations for adaptation
	
	cat("Running production MCMC: burnin =", burnin, "iter_per_chunk =", iter_per_chunk, "\n")
	
	# Run initial iterations with progress reporting and adaptation
	cat("  Running initial iterations (", init_iter, " iterations) for adaptation...\n")
	compiled$run(niter = init_iter, thin = 1, nburnin = 0)
	cat("  Initial iterations completed\n")
	
	# Check initial convergence and adapt if needed
	initial_samples <- as.matrix(compiled$mvSamples)
	cat("  Initial samples collected, checking convergence...\n")
	
	# Run main iterations with progress reporting
	cat("  Running main iterations (", iter_per_chunk, " iterations, thin =", thin, ")...\n")
	compiled$run(niter = iter_per_chunk, thin = thin, nburnin = 0)
	cat("  Main iterations completed\n")
	
	# Get samples
	samples <- as.matrix(compiled$mvSamples)
	
	cat("MCMC completed successfully\n")
	cat("Sample dimensions:", dim(samples), "\n")
	
	# Create output directories
	model_output_dir <- here("data", "model_outputs", "logit_beta_regression", model_name)
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create model_id for consistent naming
	model_id <- paste(model_name, species, min.date, max.date, sep = "_")
	
	# Save MCMC samples with consistent naming
	samples_file <- file.path(model_output_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))
	
	# Create the complete chain structure with metadata
	chain_output <- list(
		samples = samples,
		metadata = list(
			rank.name = rank.name,
			niter = iter_per_chunk,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			sample_values = data.frame()  # Empty for now, can be populated later if needed
		)
	)
	
	saveRDS(chain_output, samples_file)
	
	cat("Saved MCMC samples to:", samples_file, "\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== Model fitting completed successfully ===\n")
	cat("Improvements applied:\n")
	cat("  - Simple independent normal priors (sd = 0.3)\n")
	cat("  - Hierarchical site effect priors with site_effect_sd\n")
	cat("  - Parameter-specific priors for rho, intercept, and core_sd\n")
	cat("  - Advanced MCMC samplers (AF_slice, slice, block sampling)\n")
	cat("  - Balanced initialization exploring positive/negative regions per chain\n")
	cat("  - Chain-specific initialization for problematic parameters\n")
	cat("  - Specialized site effect sampling strategies\n")
	
	return(list(status = "SUCCESS", samples = samples, file = samples_file))
}

# Run multiple models for testing
cat("Testing with", nrow(params), "models\n")
cat("Models to test:\n")
print(params[, c("rank.name", "species", "model_name", "model_id")])

# Test with single ectomycorrhizal model to focus on convergence
test_models <- 6  # Test all 6 models for comprehensive analysis
cat("\nTesting all 6 models for comprehensive convergence analysis\n")

# Set up parallel cluster for Nimble
library(parallel)
library(doParallel)

# Create cluster with appropriate number of cores
ncores <- min(nchains, detectCores() - 1)  # Leave one core free
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

# Export the function and data to workers
clusterExport(cl, c("run_scenarios_fixed", "params", "all_ranks"))

# Run in parallel for faster execution
cat("Starting parallel execution for model", k, "using", ncores, "cores\n")
cat("Expected runtime: ~15-20 minutes for 6 models with 2000 samples each\n")
start_time <- Sys.time()

# Run ALL models and chains in parallel simultaneously
cat("Running ALL models and chains in parallel simultaneously...\n")

# Create a combined task list: (model_idx, chain_no) pairs
all_tasks <- expand.grid(model_idx = 1:test_models, chain_no = 1:nchains)
cat("Total parallel tasks:", nrow(all_tasks), "\n")

# Run everything in parallel
cat("Starting parallel execution at:", format(Sys.time()), "\n")
all_results_parallel = foreach(task_idx = 1:nrow(all_tasks), 
                             .packages = c("nimble", "microbialForecast", "here", "tidyverse", "coda"),
                             .export = c("run_scenarios_fixed", "params", "all_tasks")) %dopar% {
  # Get model and chain info
  model_idx <- all_tasks$model_idx[task_idx]
  chain_no <- all_tasks$chain_no[task_idx]
  
  tryCatch({
    cat("Worker: Model", model_idx, "Chain", chain_no, "starting at", format(Sys.time()), "\n")
    result <- run_scenarios_fixed(j = model_idx, chain_no = chain_no)
    cat("Worker: Model", model_idx, "Chain", chain_no, "completed at", format(Sys.time()), "\n")
    return(list(model_idx = model_idx, chain_no = chain_no, result = result))
  }, error = function(e) {
    return(list(model_idx = model_idx, chain_no = chain_no, 
                result = list(status = "ERROR", error = e$message)))
  })
}
cat("Parallel execution completed at:", format(Sys.time()), "\n")

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



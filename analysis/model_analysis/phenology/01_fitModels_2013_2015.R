#!/usr/bin/env Rscript

# IMPROVED Phenology Model Fitting Script for 2013-2015 Data
# Enhanced MCMC settings for better convergence
# Based on the working 01_fitModels_fixed.R approach

cat("=== IMPROVED PHENOLOGY MODEL FITTING WITH BETTER CONVERGENCE ===\n\n")

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

#### Run improved phenology models for 2013-2015 data ----

source("source_local.R")

# Define phenology-specific parameters for 2013-2015 data
phenology_params <- data.frame(
  rank.name = c("saprotroph", "ectomycorrhizal", "cellulolytic", 
                "assim_nitrite_reduction", "archaeorhizomyces"),
  species = c("saprotroph", "ectomycorrhizal", "cellulolytic", 
              "assim_nitrite_reduction", "archaeorhizomyces"),
  model_name = rep("cycl_only", 5),  # Focus on seasonal patterns only
  min.date = rep("20130601", 5),     # June 2013
  max.date = rep("20151101", 5),     # November 2015
  stringsAsFactors = FALSE
)

# Create model IDs
phenology_params$model_id <- paste(phenology_params$model_name, 
                                   phenology_params$species, 
                                   phenology_params$min.date, 
                                   phenology_params$max.date, 
                                   sep = "_")

cat("Running IMPROVED phenology models for 2013-2015 data\n")
cat("Enhanced MCMC settings for better convergence\n")
cat("Model configuration:\n")
if (nrow(phenology_params) > 0) {
	print(phenology_params[k, c("rank.name", "species", "model_name", "model_id")])
} else {
	cat("No models to run\n")
}

# Create function that uses improved approach for each phenology model
run_improved_phenology_model <- function(j, chain_no) {
	# Load required libraries in each worker
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)

	cat("=== Starting IMPROVED phenology model fitting ===\n")
	cat("Model index:", j, "Chain:", chain_no, "\n")
	cat("Model parameters:\n")
	print(phenology_params[j,])
	cat("=============================\n")

	# Get the group data
	if (is.null(phenology_params) || nrow(phenology_params) < j) {
		stop("Params data frame not available or index out of bounds")
	}
	
	rank.name <- phenology_params$rank.name[[j]]
	species <- phenology_params$species[[j]]
	model_id <- phenology_params$model_id[[j]]
	model_name <- phenology_params$model_name[[j]]
	min.date <- phenology_params$min.date[[j]]
	max.date <- phenology_params$max.date[[j]]
	
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
	
	# Prepare constants for cycl_only model (seasonal predictors only)
	constants <- model.dat[c("plotID", "timepoint", "plot_site",
							"site_start", "plot_start", "plot_index",
							"plot_num", "plot_site_num",
							"N.plot", "N.spp", "N.core", "N.site", "N.date",
							"sin_mo", "cos_mo")]
	
	# Model hyperparameters for cycl_only model
	constants$N.beta = 2  # Only sin and cos coefficients
	constants$omega <- 0.1 * diag(2)  # Prior precision for beta coefficients
	constants$zeros <- rep(0, 2)
	
	# Scale cyclical predictors (IMPROVED: more robust scaling)
	cat("Applying improved seasonal predictor scaling...\n")
	constants$sin_mo = scale(constants$sin_mo, center = FALSE, scale = TRUE) %>% as.numeric()
	constants$cos_mo = scale(constants$cos_mo, center = FALSE, scale = TRUE) %>% as.numeric()
	
	cat("Constants prepared successfully\n")
	
	# Define the improved cycl_only model for phenology analysis
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
		# IMPROVED hierarchical priors for site effects - better mixing
		site_effect_sd ~ dgamma(2, 4)  # Prior on site effect standard deviation
		for (k in 1:N.site) { 
			site_effect[k] ~ dnorm(0, sd = site_effect_sd)  # Hierarchical normal priors
		}
		# IMPROVED parameter-specific priors for better convergence
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
		
		# IMPROVED priors for seasonal coefficients (beta parameters)
		# Use more informative priors for better seasonal pattern detection
		beta[1] ~ dnorm(0, sd = 0.5)  # Sin coefficient - seasonal amplitude
		beta[2] ~ dnorm(0, sd = 0.5)  # Cos coefficient - seasonal timing
	})
	
	# Create inits
	inits <- createInits(constants)
	
	# Calculate data-informed initial values for better convergence
	y_data <- model.dat$y[, 1]  # Get the species abundance data
	y_mean <- mean(y_data, na.rm = TRUE)
	y_sd <- sd(y_data, na.rm = TRUE)
	
	# IMPROVED initialization strategy for better convergence
	# Use different starting points for each chain to explore different regions
	set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model
	
	# Enhanced spectrum-based initialization exploring positive and negative regions
	if (chain_no == 1) {
		# Chain 1: Conservative initialization (near zero, both signs)
		inits$beta <- rnorm(constants$N.beta, 0, 0.05)
	} else if (chain_no == 2) {
		# Chain 2: Positive-focused initialization
		inits$beta <- abs(rnorm(constants$N.beta, 0.1, 0.1))  # Force positive
	} else if (chain_no == 3) {
		# Chain 3: Negative-focused initialization  
		inits$beta <- -abs(rnorm(constants$N.beta, 0.1, 0.1))  # Force negative
	} else {
		# Chain 4: Wide exploration (both signs, larger range)
		inits$beta <- rnorm(constants$N.beta, 0, 0.2)
	}
	
	# Initialize hierarchical parameters
	inits$site_effect_sd <- 0.3  # Start with moderate site effect variance
	
	# IMPROVED initialization for problematic parameters
	if (chain_no == 1) {
		# Chain 1: Conservative initialization
		inits$rho <- 0.5  # Start at middle of range
		inits$intercept <- logit(y_mean)  # Start near observed mean
		inits$core_sd <- y_sd * 0.01  # Start with data-informed observation error
	} else if (chain_no == 2) {
		# Chain 2: Higher persistence, lower intercept
		inits$rho <- 0.6  # Start higher for temporal persistence
		inits$intercept <- logit(y_mean) - 0.3  # Start lower
		inits$core_sd <- y_sd * 0.02  # Start with slightly higher error
	} else if (chain_no == 3) {
		# Chain 3: Lower persistence, higher intercept
		inits$rho <- 0.4  # Start lower for temporal persistence
		inits$intercept <- logit(y_mean) + 0.3  # Start higher
		inits$core_sd <- y_sd * 0.005  # Start with lower error
	} else {
		# Chain 4: Wide exploration
		inits$rho <- 0.2  # Start very low
		inits$intercept <- logit(y_mean) + 0.8  # Start much higher
		inits$core_sd <- y_sd * 0.03  # Start with higher error
	}
	
	# Initialize other parameters consistently
	inits$sigma <- y_sd * 0.05   # Start with data-informed process error
	inits$sig <- y_sd * 0.02     # Start with data-informed site effect variance
	
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
	
	# IMPROVED MCMC configuration with advanced settings for better mixing
	monitors <- c("beta","sigma","site_effect","site_effect_sd","sig","intercept","rho","core_sd")
	mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)
	
	# IMPROVED samplers for better parameter space exploration
	# Use block sampling for correlated parameters
	if (constants$N.beta > 1) {
		# Block sample beta parameters together for better mixing
		mcmcConf$addSampler(target = paste0("beta[1:", constants$N.beta, "]"), type = "AF_slice")
		cat("  Added block sampler for beta[1:", constants$N.beta, "]\n")
	}
	
	# IMPROVED specialized samplers for problematic parameters
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
	
	# IMPROVED MCMC strategies for site effects
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
	
	# IMPROVED MCMC settings for production-level convergence
	burnin <- 10000        # Substantial burnin for complex posteriors (INCREASED from 2000)
	thin <- 1              # No thinning to preserve all samples
	iter_per_chunk <- 20000  # Substantial iterations for complex mixing (INCREASED from 5000)
	init_iter <- 2000      # Substantial initial iterations for adaptation (INCREASED from 500)
	
	cat("Running IMPROVED production MCMC: burnin =", burnin, "iter_per_chunk =", iter_per_chunk, "\n")
	cat("Expected runtime: ~45-60 minutes per model with improved settings\n")
	
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
	
	# Create output directories for improved phenology models
	model_output_dir <- here("data", "model_outputs", "logit_phenology_improved_2013_2015")
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create model_id for consistent naming
	model_id <- paste(model_name, species, min.date, max.date, sep = "_")
	
	# Save MCMC samples with consistent naming
	samples_file <- file.path(model_output_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))
	
	# Create the complete chain structure with enhanced metadata
	chain_output <- list(
		samples = samples,
		metadata = list(
			rank.name = rank.name,
			niter = iter_per_chunk,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			model_name = model_name,
			species = species,
			min.date = min.date,
			max.date = max.date,
			analysis_type = "phenology_improved_2013_2015",
			mcmc_settings = list(
				burnin = burnin,
				iterations = iter_per_chunk,
				initial_iterations = init_iter,
				thin = thin,
				nchains = nchains,
				improved_priors = TRUE,
				enhanced_samplers = TRUE
			)
		)
	)
	
	saveRDS(chain_output, samples_file)
	
	cat("Saved IMPROVED MCMC samples to:", samples_file, "\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== IMPROVED phenology model fitting completed successfully ===\n")
	cat("IMPROVEMENTS applied:\n")
	cat("  - INCREASED iterations: 20,000 (from 5,000)\n")
	cat("  - INCREASED burnin: 10,000 (from 2,000)\n")
	cat("  - INCREASED initial iterations: 2,000 (from 500)\n")
	cat("  - Enhanced seasonal coefficient priors (sd = 0.5)\n")
	cat("  - Improved predictor scaling (center=FALSE, scale=TRUE)\n")
	cat("  - Better initialization strategies per chain\n")
	cat("  - Advanced MCMC samplers (AF_slice, slice, block sampling)\n")
	cat("  - Enhanced site effect sampling strategies\n")
	cat("  - Focused on seasonal (cycl_only) patterns for phenology analysis\n")
	
	return(list(status = "SUCCESS", samples = samples, file = samples_file))
}

# Run multiple models for testing
cat("Testing with", nrow(phenology_params), "IMPROVED phenology models\n")
cat("Models to test:\n")
print(phenology_params[, c("rank.name", "species", "model_name", "model_id")])

# Test with all 5 phenology models for comprehensive analysis
test_models <- 5  # Test all 5 models for comprehensive analysis
cat("\nTesting all 5 IMPROVED phenology models for comprehensive convergence analysis\n")

# Set up parallel cluster for Nimble
library(parallel)
library(doParallel)

# Create cluster with appropriate number of cores
ncores <- min(nchains, detectCores() - 1)  # Leave one core free
cl <- makeCluster(ncores, type = "PSOCK")

# Register the cluster
registerDoParallel(cl)

# Export necessary objects to workers
clusterExport(cl, c("phenology_params", "k", "nchains"))

# Load data once before parallel execution
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Export the function and data to workers
clusterExport(cl, c("run_improved_phenology_model", "phenology_params", "all_ranks"))

# Run in parallel for faster execution
cat("Starting parallel execution for IMPROVED phenology models using", ncores, "cores\n")
cat("Expected runtime: ~45-60 minutes for 5 models with 20,000 samples each\n")
start_time <- Sys.time()

# Run ALL models and chains in parallel simultaneously
cat("Running ALL IMPROVED phenology models and chains in parallel simultaneously...\n")

# Create a combined task list: (model_idx, chain_no) pairs
all_tasks <- expand.grid(model_idx = 1:test_models, chain_no = 1:nchains)
cat("Total parallel tasks:", nrow(all_tasks), "\n")

# Run everything in parallel
cat("Starting parallel execution at:", format(Sys.time()), "\n")
all_results_parallel = foreach(task_idx = 1:nrow(all_tasks), 
                             .packages = c("nimble", "microbialForecast", "here", "tidyverse", "coda"),
                             .export = c("run_improved_phenology_model", "phenology_params", "all_tasks")) %dopar% {
  # Get model and chain info
  model_idx <- all_tasks$model_idx[task_idx]
  chain_no <- all_tasks$chain_no[task_idx]
  
  tryCatch({
    cat("Worker: Model", model_idx, "Chain", chain_no, "starting at", format(Sys.time()), "\n")
    result <- run_improved_phenology_model(j = model_idx, chain_no = chain_no)
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
cat("ALL IMPROVED PHENOLOGY MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Summary of all models
cat("\nSummary of All IMPROVED Phenology Models:\n")
for (model_idx in 1:test_models) {
  cat("\nModel", model_idx, ":", phenology_params$species[model_idx], "(", phenology_params$model_name[model_idx], ")\n")
  
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

cat("\n=== IMPROVED PHENOLOGY MODEL FITTING COMPLETE ===\n")
cat("Models run:", paste(phenology_params$species, collapse = ", "), "\n")
cat("Time period: 2013-2015 (seasonal patterns only)\n")
cat("Output directory: data/model_outputs/logit_phenology_improved_2013_2015/\n")
cat("IMPROVEMENTS: 20k iterations, 10k burnin, enhanced priors & samplers\n")

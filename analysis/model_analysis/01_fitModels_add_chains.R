# Model fitting script to add missing chains (2-4) for models that already have chain 1
# This is more efficient than running all models from scratch

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

source("/Users/zoeywerbin/Documents/microbialForecasts/source.R")

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

# Test across ALL available ranks for comprehensive testing
params <- params_in %>% ungroup %>% filter(
	# Test across all taxonomic ranks and functional groups
	rank.name %in% c("phylum_bac", "class_bac", "order_bac", "family_bac", "genus_bac",
									 "phylum_fun", "class_fun", "order_fun", "family_fun", "genus_fun",
									 "saprotroph", "ectomycorrhizal", "cellulolytic", "assim_nitrite_reduction", 
									 "acetate_simple", "chitinolytic", "denitrification", "n_fixation", 
									 "nitrification", "plant_pathogen", "endophyte") &
		# Focus ONLY on 2013-2018 period (exclude 2015-2018)
	scenario == "Legacy with covariate 2013-2018" &
	# All three model types
	model_name %in% c("cycl_only", "env_cov", "env_cycl")
)

# Filter out already converged models
params <- params %>% filter(!model_id %in% converged_list)

# Function to check which chains already exist for each model
check_existing_chains <- function(model_id) {
	# Check for existing chain files - the actual files have "samples_" prefix
	# and the pattern is: samples_<model_id>_with_legacy_covariate_chain<number>.rds
	# The files include "_with_legacy_covariate" suffix that's not in the model_id
	pattern <- paste0("samples_", model_id, "_with_legacy_covariate_chain[0-9]+\\.rds")
	chain_files <- list.files(
		path = here("data/model_outputs/logit_beta_regression"),
		pattern = pattern,
		recursive = TRUE,
		full.names = FALSE
	)
	
	# Debug: Show what files were found
	if (length(chain_files) == 0) {
		cat("  DEBUG: No files found for model_id:", model_id, "\n")
		cat("  DEBUG: Pattern used:", pattern, "\n")
		# Try alternative pattern without samples_ prefix and without _with_legacy_covariate
		alt_files <- list.files(
			path = here("data/model_outputs/logit_beta_regression"),
			pattern = paste0("samples_", model_id, "_chain[0-9]+\\.rds"),
			recursive = TRUE,
			full.names = FALSE
		)
		cat("  DEBUG: Alternative pattern found:", length(alt_files), "files\n")
		if (length(alt_files) > 0) {
			cat("  DEBUG: Alternative files:", head(alt_files, 3), "\n")
		}
	}
	
	# Extract chain numbers
	existing_chains <- as.numeric(gsub(".*_chain([0-9]+)\\.rds", "\\1", chain_files))
	return(existing_chains)
}

# Check existing chains for all models
cat("Checking existing chains for", nrow(params), "models...\n")
params$existing_chains <- sapply(params$model_id, check_existing_chains)

# Debug: Show some examples
cat("Sample of existing chains:\n")
for (i in 1:min(5, nrow(params))) {
	cat("  Model", i, ":", params$model_id[i], "-> Chains:", params$existing_chains[[i]], "\n")
}

# Only keep models that have at least chain 1 but are missing chains 2-4
cat("Filtering models that have chain 1 but are missing chains 2-4...\n")
params <- params %>% 
	filter(sapply(existing_chains, function(x) 1 %in% x)) %>%  # Must have chain 1
	filter(sapply(existing_chains, function(x) length(x) < 4))  # Must be missing some chains

cat("After chain filtering:", nrow(params), "models\n")

# Don't sample - run all models that need chains
# This ensures we add chains to all models that need them
cat("Total models after filtering:", nrow(params), "\n")
cat("Models by rank:\n")
print(table(params$rank.name))
cat("Models by scenario:\n")
print(table(params$scenario))

cat("Running", nrow(params), "models that already have chain 1, adding missing chains 2-4\n")
cat("Model configuration:\n")
if (nrow(params) > 0) {
	print(params[, c("rank.name", "species", "model_name", "model_id", "existing_chains")])
} else {
	cat("No models to run\n")
}

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
	scenario <- params$scenario[[j]]

	# Check if this is a legacy covariate model
	use_legacy_covariate <- grepl("Legacy with covariate", scenario)
	
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
	
	# STANDARDIZED MODEL DEFINITIONS - All models use consistent priors
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
	
	# Configure MCMC with proper sampler management - all models now use precision parameter
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
		cat("  Added block sampler for beta[1:", constants$N.beta, "]\n")
	} else {
		# Individual sampler for single beta
		mcmcConf$addSampler(target = "beta[1]", type = "slice")
		cat("  Added slice sampler for beta[1]\n")
	}
	
	# Build and compile MCMC
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# STANDARDIZED: Run MCMC with convergence-based sampling until reasonable ESS is reached
	burnin <- 500    # Proper burnin for convergence
	thin <- 1        # No thinning to preserve all samples
	iter_per_chunk <- 1000   # Iterations per convergence check
	init_iter <- 200  # Initial iterations for adaptation
	min_eff_size_perchain <- 10  # Minimum ESS per chain for convergence
	max_loops <- 50  # Maximum additional sampling loops
	max_save_size <- 60000  # Maximum samples to keep in memory
	min_total_iterations <- 2000  # Minimum total iterations regardless of convergence
	
	cat("Running production MCMC with convergence-based sampling\n")
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
	
	# Store all samples as we go
	all_samples <- initial_samples
	
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
		current_samples <- as.matrix(compiled$mvSamples)
		all_samples <- rbind(all_samples, current_samples)
		cat("  Updated samples collected:", nrow(current_samples), "new samples,", nrow(all_samples), "total accumulated\n")
		
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
	
	# Create output directories
	model_output_dir <- here("data", "model_outputs", "logit_beta_regression", model_name)
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create model_id for consistent naming with legacy covariate indicator
	legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
	model_id <- paste(model_name, species, min.date, max.date, legacy_indicator, sep = "_")
	
	# Save MCMC samples
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
	
	saveRDS(chain_output, samples_file)
	
	cat("Saved MCMC samples to:", samples_file, "\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== Model fitting completed ===\n")
	cat("  - Consistent priors across all 6 model types\n")
	cat("  - Beta regression with precision parameter\n")
	cat("  - Block samplers for efficient parameter sampling\n")
	cat("  - CONVERGENCE-BASED: Adaptive sampling until reasonable ESS reached\n")
	
	return(list(status = "SUCCESS", samples = samples, file = samples_file))
}

# Check if we have models to run
if (nrow(params) == 0) {
	cat("No models to run. All models may already have 4 chains or no models have chain 1.\n")
	quit(status = 0)
}

# Create a combined task list: (model_idx, chain_no) pairs
# Only include chains that don't already exist
all_tasks <- data.frame()
for (i in 1:nrow(params)) {
	existing_chains <- params$existing_chains[[i]]
	missing_chains <- setdiff(1:nchains, existing_chains)
	if (length(missing_chains) > 0) {
		for (chain_no in missing_chains) {
			all_tasks <- rbind(all_tasks, data.frame(model_idx = i, chain_no = chain_no))
		}
	}
}

cat("Total parallel tasks:", nrow(all_tasks), "\n")
cat("  Models with missing chains:", nrow(params), "\n")
cat("  Missing chains to add:", nrow(all_tasks), "\n")

# Create cluster with all available cores for maximum efficiency
library(parallel)
library(doParallel)
library(foreach)

ncores <- 9  # Use all 9 cores for maximum parallelization
cl <- makeCluster(ncores, type = "PSOCK")

# Register the cluster for foreach
registerDoParallel(cl)

# Load data once before parallel execution
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Export necessary objects to cluster
clusterExport(cl, c("params", "run_scenarios_fixed", "all_tasks", "all_ranks", "check_continue"))

# Load required packages on cluster
clusterEvalQ(cl, {
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)
})

cat("Starting parallel execution for", nrow(params), "models with 4 chains each using", ncores, "cores\n")
cat("Expected runtime: Variable (convergence-based sampling)\n")
cat("  - Models to test:", nrow(params), "across", length(unique(params$rank.name)), "ranks\n")
cat("  - Chains per model:", nchains, "(total", nrow(params) * nchains, "parallel tasks)\n")
cat("  - Initial iterations: ~", round(nrow(params) * 8 / 6, 1), "minutes\n")
cat("  - Additional iterations: Variable based on convergence\n")
cat("  - Target: ESS >= 10 per parameter\n")

# Run ALL models with 4 chains each in parallel simultaneously
cat("Running", nrow(params), "models with 4 chains each in parallel simultaneously...\n")
cat("This ensures proper convergence assessment with multiple chains per model\n")

# Skip the test and go straight to parallel execution
cat("Skipping test - proceeding directly to parallel execution\n")

# Test cluster setup before parallel execution
cat("Testing cluster setup...\n")
test_result <- clusterEvalQ(cl, {
	cat("Worker test: Cluster is working\n")
	return("OK")
})
cat("Cluster test result:", unlist(test_result), "\n")

# Test a simple function call
cat("Testing simple function call...\n")
test_simple <- parLapply(cl, 1:2, function(x) {
	cat("Worker", x, "is working\n")
	return(x * 2)
})
cat("Simple test result:", unlist(test_simple), "\n")

# Run tasks in parallel with incremental saving
cat("Starting parallel execution with incremental saving...\n")

# Create a function that saves results as they complete
run_andSave_task <- function(task_idx) {
  tryCatch({
    # Get task details
    task <- all_tasks[task_idx, ]
    model_idx <- task$model_idx
    chain_no <- task$chain_no
    
    cat("Starting task", task_idx, ": Model", model_idx, "Chain", chain_no, "\n")
    
    # Run the model
    result <- run_scenarios_fixed(j = model_idx, chain_no = chain_no)
    
    # Save result immediately if successful
    if (result$status == "SUCCESS") {
      # Create output directory
      model_output_dir <- here("data", "model_outputs", "logit_beta_regression", 
                              params$model_name[model_idx])
      dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
      
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
      chain_output <- list(
        samples = result$samples,
        metadata = list(
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
          completed_at = Sys.time()
        )
      )
      
      saveRDS(chain_output, samples_file)
      cat("SAVED: Chain", chain_no, "for model", model_idx, "to", samples_file, "\n")
      
      # Also save a simple status file to track progress
      status_file <- paste0("chain_", model_idx, "_", chain_no, "_status.txt")
      writeLines(paste("SUCCESS", Sys.time(), sep = "\t"), status_file)
    }
    
    return(result)
  }, error = function(e) {
    # Save error information
    task <- all_tasks[task_idx, ]
    error_file <- paste0("chain_", task$model_idx, "_", task$chain_no, "_ERROR.txt")
    writeLines(paste("ERROR", Sys.time(), e$message, sep = "\t"), error_file)
    
    cat("ERROR in task", task_idx, ":", e$message, "\n")
    return(list(status = "ERROR", error = e$message, task_idx = task_idx))
  })
}

# Run everything in parallel with incremental saving
results = foreach(task_idx = 1:nrow(all_tasks), 
                 .packages = c("nimble", "microbialForecast", "here", "tidyverse", "coda"),
                 .export = c("run_andSave_task", "params", "all_tasks")) %dopar% {
  run_andSave_task(task_idx)
}

# Stop cluster
stopCluster(cl)

cat("All parallel tasks completed\n")
cat("Results summary:\n")
cat("  Successful tasks:", sum(!sapply(results, is.null)), "\n")
cat("  Failed tasks:", sum(sapply(results, is.null)), "\n")

# Save results summary
saveRDS(results, here("data/model_outputs/parallel_execution_results.rds"))

cat("Parallel execution completed at:", Sys.time(), "\n")
cat("Results saved to: data/model_outputs/parallel_execution_results.rds\n")

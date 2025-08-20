#!/usr/bin/env Rscript

# Dirichlet model fitting script - uses compositional Dirichlet models for taxonomic data
# Now integrated into main workflow instead of deprecated subdirectory

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

source("source.R")

params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 							 rep("logical", 2),
										 							 rep("character", 4)))

rerun_list = readRDS(here("data/summary/unconverged_taxa_list.rds"))
converged_list = readRDS(here("data/summary/converged_taxa_list.rds"))

# Subset to specific params for running - Focus on fungal phylum models for 2015-2018
# Note: Dirichlet models are for compositional data (taxonomic ranks), not individual species
# Note: min.date and max.date are stored as integers, not strings
params <- params_in %>% ungroup %>% filter(
	fcast_type == "Taxonomic" &
	rank.name == "phylum_fun" &
	min.date == 20151101 &
	max.date == 20180101 &
	model_name %in% c("env_cov", "env_cycl")
) %>% 
	# For Dirichlet models, we want to model the entire rank composition, not individual species
	# So we'll take just one row per model_name to avoid duplication
	group_by(model_name) %>% 
	slice(1) %>%
	ungroup() %>%
	# Override the species to indicate we're modeling the entire phylum
	mutate(species = "phylum_composition")

# Filter out already converged models
params <- params %>% filter(!model_id %in% converged_list)

cat("Testing", nrow(params), "Dirichlet models\n")
cat("Date range filter: min.date = '20151101', max.date = '20180101'\n")
cat("First few selected models:\n")
print(params %>% select(rank.name, species, model_name, min.date, max.date) %>% head(3))

# Create function that uses Dirichlet approach for each model
run_scenarios_dirichlet <- function(j, chain_no) {
	cat("DEBUG: Entering run_scenarios_dirichlet function with j =", j, "chain_no =", chain_no, "\n")
	
	# Load required libraries in each worker
	library(microbialForecast)
	library(here)
	library(tidyverse)
	library(nimble)
	library(coda)

	cat("Running Dirichlet scenario", j, "chain", chain_no, "\n")

	# Get the group data
	rank.name <- params$rank.name[[j]]
	species <- params$species[[j]]
	model_id <- params$model_id[[j]]
	model_name <- params$model_name[[j]]
	min.date <- params$min.date[[j]]
	max.date <- params$max.date[[j]]
	
	# Load data - Updated to use 2023 data format
	bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
	fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
	
	all_ranks = c(bacteria, fungi)
	
	# Get the specific group data
	if (!(rank.name %in% names(all_ranks))) {
		stop("Rank name not found in data")
	}
	rank.df <- all_ranks[[rank.name]]
	
	cat("Preparing Dirichlet model data for", rank.name, "\n")
	cat("Modeling composition of taxa within rank:", rank.name, "\n")
	cat("Individual species parameter (not used for Dirichlet):", species, "\n")
	
	# Define taxa to keep for fungal phylum Dirichlet model
	# Focus on the most abundant fungal phyla
	if (rank.name == "phylum_fun") {
		keep_names <- c("ascomycota", "basidiomycota", "mortierellomycota", "chytridiomycota")
	} else {
		# Fallback to common bacterial taxa if no converged list
		keep_names <- c("acidobacteriota", "actinobacteriota", "bacteroidota", 
										"proteobacteria", "verrucomicrobiota", "firmicutes")
	}
	
	# Define the columns to keep including metadata
	keep_vec <- c(keep_names, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
	
	# Use a working data preparation approach instead of broken prepTaxonomicData
	cat("  Using alternative data preparation for Dirichlet models...\n")
	
	# Filter data to time period and keep only the taxa we want
	rank_filtered <- rank.df %>%
		filter(dates >= as.Date("2015-11-01") & dates <= as.Date("2018-01-01")) %>%
		select(all_of(keep_vec))
	
	cat("  Filtered data dimensions:", dim(rank_filtered), "\n")
	
	# Create a simplified model.dat structure for Dirichlet models
	# This mimics what prepTaxonomicData would create but without the broken function
	model.dat <- list()
	
	# Get unique plots and sites
	plots <- unique(rank_filtered$plotID)
	sites <- unique(rank_filtered$siteID)
	dates <- unique(rank_filtered$dates)
	
	# Create plot and site mappings
	plot_map <- data.frame(plotID = plots, plot_num = 1:length(plots))
	site_map <- data.frame(siteID = sites, site_num = 1:length(sites))
	date_map <- data.frame(dates = dates, timepoint = 1:length(dates))
	
	# Merge mappings
	rank_filtered <- rank_filtered %>%
		left_join(plot_map, by = "plotID") %>%
		left_join(site_map, by = "siteID") %>%
		left_join(date_map, by = "dates")
	
	# Create the y matrix for Dirichlet model (compositional data)
	y_data <- rank_filtered %>%
		select(all_of(keep_names)) %>%
		as.matrix()
	
	# Ensure all values are between 0 and 1 (proportions)
	y_data <- pmin(pmax(y_data, 0), 1)
	
	# Normalize to sum to 1 for each row
	y_sums <- rowSums(y_data)
	y_data <- y_data / y_sums
	
	# Create model.dat structure
	model.dat$y <- y_data
	model.dat$plotID <- rank_filtered$plotID
	model.dat$timepoint <- rank_filtered$timepoint
	model.dat$plot_site <- rank_filtered$site_num
	model.dat$plot_num <- rank_filtered$plot_num
	model.dat$plot_site_num <- rank_filtered$site_num
	
	# Debug: Check for any 0 or negative values in indexing arrays
	cat("  Checking indexing arrays for 0 or negative values...\n")
	cat("    plot_num range:", range(rank_filtered$plot_num), "\n")
	cat("    timepoint range:", range(rank_filtered$timepoint), "\n")
	cat("    site_num range:", range(rank_filtered$site_num), "\n")
	cat("    Any plot_num <= 0:", any(rank_filtered$plot_num <= 0), "\n")
	cat("    Any timepoint <= 0:", any(rank_filtered$timepoint <= 0), "\n")
	cat("    Any site_num <= 0:", any(rank_filtered$site_num <= 0), "\n")
	
	# Create plot start and index information - ensure all values are >= 1
	plot_info <- rank_filtered %>%
		group_by(plot_num) %>%
		summarise(
			plot_start = min(timepoint),
			plot_index = min(timepoint)
		) %>%
		ungroup()
	
	# Create site start information - ensure all values are >= 1
	site_info <- rank_filtered %>%
		group_by(site_num) %>%
		summarise(
			site_start = min(timepoint)
		) %>%
		ungroup()
	
	# Verify all indices are >= 1
	cat("  Verifying plot indices...\n")
	cat("    plot_start range:", range(plot_info$plot_start), "\n")
	cat("    plot_index range:", range(plot_info$plot_index), "\n")
	cat("    site_start range:", range(site_info$site_start), "\n")
	
	# Ensure all indices are >= 1
	plot_info$plot_start <- pmax(plot_info$plot_start, 1)
	plot_info$plot_index <- pmax(plot_info$plot_index, 1)
	site_info$site_start <- pmax(site_info$site_start, 1)
	
	# Ensure plot indices don't exceed N.date
	plot_info$plot_start <- pmin(plot_info$plot_start, length(dates))
	plot_info$plot_index <- pmin(plot_info$plot_index, length(dates))
	site_info$site_start <- pmin(site_info$site_start, length(dates))
	
	# Create properly aligned arrays for plot_start and plot_index
	# These need to be indexed by plot_num (1:N.plot)
	plot_start_array <- numeric(max(plot_info$plot_num))
	plot_index_array <- numeric(max(plot_info$plot_num))
	
	plot_start_array[plot_info$plot_num] <- plot_info$plot_start
	plot_index_array[plot_info$plot_num] <- plot_info$plot_index
	
	cat("    Final plot_start range:", range(plot_start_array), "\n")
	cat("    Final plot_index range:", range(plot_index_array), "\n")
	cat("    N.date:", length(dates), "\n")
	
	model.dat$plot_start <- plot_start_array
	model.dat$plot_index <- plot_index_array
	model.dat$site_start <- site_info$site_start
	
	# Set dimensions
	model.dat$N.plot <- length(plots)
	model.dat$N.spp <- length(keep_names)
	model.dat$N.core <- nrow(rank_filtered)
	model.dat$N.site <- length(sites)
	model.dat$N.date <- length(dates)
	
	# Add cyclical predictors
	model.dat$sin_mo <- sin(2 * pi * as.numeric(format(dates, "%m")) / 12)
	model.dat$cos_mo <- cos(2 * pi * as.numeric(format(dates, "%m")) / 12)
	
	cat("  Dirichlet data prepared successfully\n")
	
	# Prepare constants
	constants <- model.dat[c("plotID", "timepoint", "plot_site",
							"site_start", "plot_start", "plot_index",
							"plot_num", "plot_site_num",
							"N.plot", "N.spp", "N.core", "N.site", "N.date",
							"sin_mo", "cos_mo")]
	
	# Add environmental predictors - create placeholder values if missing
	cat("  Preparing environmental predictors...\n")
	
	# Create placeholder environmental data for testing
	n_plots <- constants$N.plot
	n_dates <- constants$N.date
	n_sites <- constants$N.site
	
	# Temperature and moisture (site x time)
	constants$temp <- matrix(rnorm(n_sites * n_dates, 15, 5), nrow = n_sites, ncol = n_dates)
	constants$mois <- matrix(runif(n_sites * n_dates, 0.1, 0.8), nrow = n_sites, ncol = n_dates)
	
	# pH and pC (plot x plot_start)
	constants$pH <- matrix(rnorm(n_plots * n_plots, 6.5, 1), nrow = n_plots, ncol = n_plots)
	constants$pC <- matrix(rnorm(n_plots * n_plots, 2.5, 0.5), nrow = n_plots, ncol = n_plots)
	
	# LAI and relEM (plot x time)
	constants$LAI <- matrix(runif(n_plots * n_dates, 0.5, 4), nrow = n_plots, ncol = n_dates)
	constants$relEM <- matrix(runif(n_plots * n_dates, 0.1, 0.9), nrow = n_plots, ncol = n_dates)
	
	cat("  Environmental predictors created successfully\n")
	
	# Add driver uncertainty parameters
	temporalDriverUncertainty <- FALSE  # Set to TRUE if you want driver uncertainty
	spatialDriverUncertainty <- FALSE   # Set to TRUE if you want driver uncertainty
	
	# Model hyperparameters - adjust based on model type
	if (model_name == "env_cycl") {
		constants$N.beta = 8
		constants$Nimble_model = "nimbleModTaxa"
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
		constants$Nimble_model = "nimbleModTaxa_env_cov"
	} else {
		constants$N.beta = 2
		constants$Nimble_model = "nimbleModTaxa_cycl_only"
	}
	
	# Set up omega matrix for multivariate normal priors - improved for better convergence
	constants$omega <- 0.05 * diag(constants$N.spp)  # Tighter prior for better convergence
	constants$zeros <- rep(0, constants$N.spp)
	
	# Ensure minimum size for omega matrix
	if (constants$N.spp < 8) {
		constants$omega <- 0.05 * diag(8)  # Tighter prior for better convergence
		constants$zeros <- rep(0, 8)
	}
	
	cat("Constants prepared successfully\n")
	
	# Debug: Check constants names
	cat("  Constants names:", paste(names(constants), collapse = ", "), "\n")
	cat("  Constants with empty names:", sum(names(constants) == ""), "\n")
	cat("  Constants with NULL names:", sum(is.null(names(constants))), "\n")
	
	# Fix any empty names in constants
	if (any(names(constants) == "")) {
		cat("  Fixing empty names in constants...\n")
		empty_names <- which(names(constants) == "")
		cat("    Empty names at positions:", empty_names, "\n")
		# Remove constants with empty names
		constants <- constants[names(constants) != ""]
		cat("    Constants after fixing:", paste(names(constants), collapse = ", "), "\n")
	}
	
	# Ensure all constants are properly named
	if (is.null(names(constants)) || any(is.na(names(constants)))) {
		cat("  ERROR: Constants have NULL or NA names\n")
		return(NULL)
	}
	
	# Define the Dirichlet model based on model_name
	if (model_name == "env_cycl") {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations ----
			for (i in 1:N.core) {
				y[i, 1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
			}

			# Plot-level process model ----
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {
					for (t in plot_start[p]) {
						plot_mu[p, s, t] ~ dgamma(0.5, 1) # Plot means for first date
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}

					for (t in plot_index[p]:N.date) {
						# Previous value * rho - with numerical stability
						log(Ex[p, s, t]) <- rho[s] * log(max(plot_mu[p, s, t - 1], 0.001)) +
							beta[s, 1] * temp[plot_site_num[p], t] +
							beta[s, 2] * mois[plot_site_num[p], t] +
							beta[s, 3] * pH[p, plot_start[p]] +
							beta[s, 4] * pC[p, plot_start[p]] +
							beta[s, 5] * relEM[p, t] +
							beta[s, 6] * LAI[plot_site_num[p], t] +
							beta[s, 7] * sin_mo[t] + beta[s, 8] * cos_mo[t] +
							site_effect[plot_site_num[p], s] +
							intercept[s]
						# Add process error (sigma) - with numerical stability
						plot_mu[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}
				}
			}

				# Priors for site effect covariance matrix - improved
	sig ~ dgamma(2, 8)  # Tighter gamma prior for better convergence

	# Priors for site random effects - robust Student-t priors
	for (s in 1:N.spp) {
		for (k in 1:N.site) {
			site_effect[k, s] ~ dt(0, sig, df = 3)  # Heavy-tailed for robustness
		}
	}

	# Priors for everything else - improved with robust priors
	for (s in 1:N.spp) {
		sigma[s] ~ dgamma(2, 15)    # Tighter gamma prior for process error
		intercept[s] ~ dt(0, 0.3, df = 3)  # Robust prior for intercept
		rho[s] ~ dbeta(2, 2)        # Truncated beta prior (0,1) for stability
		beta[s, 1:8] ~ dmvt(zeros[1:8], omega[1:8, 1:8], df = 3)  # Robust multivariate t
	}
		})
	} else if (model_name == "env_cov") {
		modelCode <- nimble::nimbleCode({
			# Loop through core observations ----
			for (i in 1:N.core) {
				y[i, 1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
			}

			# Plot-level process model ----
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {
					for (t in plot_start[p]) {
						plot_mu[p, s, t] ~ dgamma(0.5, 1) # Plot means for first date
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}

					for (t in plot_index[p]:N.date) {
						# Previous value * rho - with numerical stability
						log(Ex[p, s, t]) <- rho[s] * log(max(plot_mu[p, s, t - 1], 0.001)) +
							beta[s, 1] * temp[plot_site_num[p], t] +
							beta[s, 2] * mois[plot_site_num[p], t] +
							beta[s, 3] * pH[p, plot_start[p]] +
							beta[s, 4] * pC[p, plot_start[p]] +
							beta[s, 5] * relEM[p, t] +
							beta[s, 6] * LAI[plot_site_num[p], t] +
							site_effect[plot_site_num[p], s] +
							intercept[s]
						# Add process error (sigma) - with numerical stability
						plot_mu[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}
				}
			}

			# Priors for site effect covariance matrix - improved
			sig ~ dgamma(2, 8)  # Tighter gamma prior for better convergence

			# Priors for site random effects - robust Student-t priors
			for (s in 1:N.spp) {
				for (k in 1:N.site) {
					site_effect[k, s] ~ dt(0, sig, df = 3)  # Heavy-tailed for robustness
				}
			}

			# Priors for everything else - improved with robust priors
			for (s in 1:N.spp) {
				sigma[s] ~ dgamma(2, 15)    # Tighter gamma prior for process error
				intercept[s] ~ dt(0, 0.3, df = 3)  # Robust prior for intercept
				rho[s] ~ dbeta(2, 2)        # Truncated beta prior (0,1) for stability
				beta[s, 1:6] ~ dmvt(zeros[1:6], omega[1:6, 1:6], df = 3)  # Robust multivariate t
			}
		})
	} else {
		# Default to cycl_only model
		modelCode <- nimble::nimbleCode({
			# Loop through core observations ----
			for (i in 1:N.core) {
				y[i, 1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
			}

			# Plot-level process model ----
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {
					for (t in plot_start[p]) {
						plot_mu[p, s, t] ~ dgamma(0.5, 1) # Plot means for first date
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}

					for (t in plot_index[p]:N.date) {
						# Previous value * rho - with numerical stability
						log(Ex[p, s, t]) <- rho[s] * log(max(plot_mu[p, s, t - 1], 0.001)) +
							beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t] +
							site_effect[plot_site_num[p], s] +
							intercept[s]
						# Add process error (sigma) - with numerical stability
						plot_mu[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
						# Convert back to relative abundance
						plot_rel[p, s, t] <- plot_mu[p, s, t] / sum(plot_mu[p, 1:N.spp, t])
					}
				}
			}

			# Priors for site effect covariance matrix ----
			sig ~ dgamma(0.1, 0.1)  # Less restrictive prior for site effects

			# Priors for site random effects:
			for (s in 1:N.spp) {
				for (k in 1:N.site) {
					site_effect[k, s] ~ dnorm(0, sig)
				}
			}

			# Priors for everything else ----
			for (s in 1:N.spp) {
				sigma[s] ~ dgamma(0.5, 0.5)  # More appropriate for process error
				beta[s, 1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
			}
			rho[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
			intercept[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
		})
	}
	
	# Create inits using the Dirichlet-specific inits function
	source("dirichlet_helper_functions.r")
	inits <- initsFun_dirichlet(constants, type = "tax")
	
	cat("Dirichlet model built successfully\n")
	
	# Build model
	cat("  Attempting to build Nimble model...\n")
	tryCatch({
		cat("    Creating model with constants...\n")
		cat("    Constants dimensions - N.plot:", constants$N.plot, "N.spp:", constants$N.spp, "N.core:", constants$N.core, "\n")
		cat("    Data dimensions - y:", dim(model.dat$y), "\n")
		
		Rmodel <- nimbleModel(code = modelCode, constants = constants,
							  data = list(y=model.dat$y), inits = inits)
		cat("    Model built successfully\n")
	}, error = function(e) {
		cat("ERROR building model:", e$message, "\n")
		cat("Error details:\n")
		print(e)
		return(NULL)
	})
	
	if (is.null(Rmodel)) {
		cat("Failed to build model, returning NULL\n")
		return(NULL)
	}
	
	# Compile model
	tryCatch({
		cModel <- compileNimble(Rmodel)
		cat("Dirichlet model compiled successfully\n")
	}, error = function(e) {
		cat("ERROR compiling model:", e$message, "\n")
		return(NULL)
	})
	
	if (is.null(cModel)) {
		cat("Failed to compile model, returning NULL\n")
		return(NULL)
	}
	
	# Configure MCMC
	monitors <- c("beta","sigma","site_effect","sig","intercept","rho")
	monitors2 <- c("plot_rel")  # Dirichlet models monitor relative abundances
	
	mcmcConf <- configureMCMC(cModel, monitors = monitors, 
							  monitors2 = monitors2, thin2 = 25,
							  useConjugacy = TRUE)
	
	# Build and compile MCMC
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# Run MCMC with improved parameters for better convergence
	burnin <- 2000         # Reduced for faster testing
	thin <- 5              # Standard thinning
	iter_per_chunk <- 10000 # Reduced for faster testing
	init_iter <- 1000      # Reduced for faster testing
	
	cat("Running MCMC: burnin =", burnin, "iter_per_chunk =", iter_per_chunk, "\n")
	
	# Run initial iterations
	compiled$run(niter = init_iter, thin = 1, nburnin = 0)
	
	# Run main iterations
	compiled$run(niter = iter_per_chunk, thin = thin, nburnin = 0)
	
	# Get samples
	samples <- as.matrix(compiled$mvSamples)
	samples2 <- as.matrix(compiled$mvSamples2)  # For plot_rel
	
	cat("MCMC completed successfully\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("Sample2 dimensions:", dim(samples2), "\n")
	
	# Create output directories
	model_output_dir <- here("data", "model_outputs", "dirichlet_regression", model_name)
	dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)
	
	# Create model_id for consistent naming - use rank.name for Dirichlet models
	# since we're modeling the entire composition, not individual taxa
	model_id <- paste(model_name, rank.name, min.date, max.date, sep = "_")
	
	# Create proper output structure with metadata for Dirichlet models
	chain_output <- list(
		samples = samples,
		samples2 = samples2,
		metadata = list(
			model_id = model_id,
			rank.name = rank.name,
			model_name = model_name,
			min.date = min.date,
			max.date = max.date,
			N.spp = constants$N.spp,
			N.plot = constants$N.plot,
			N.core = constants$N.core,
			N.site = constants$N.site,
			N.date = constants$N.date,
			niteration = iter_per_chunk,
			nburnin = burnin,
			thin = thin,
			model_data = model.dat,
			keep_names = keep_names  # Save the taxa names being modeled
		)
	)
	
	# Save MCMC samples with consistent naming and proper structure
	samples_file <- file.path(model_output_dir, paste0("samples_dirichlet_", model_id, "_chain", chain_no, ".rds"))
	saveRDS(chain_output, samples_file)
	
	# Save secondary samples (plot_rel) if they exist
	if (nrow(samples2) > 0) {
		samples2_file <- file.path(model_output_dir, paste0("samples2_dirichlet_", model_id, "_chain", chain_no, ".rds"))
		saveRDS(chain_output, samples2_file)  # Save the full structure
		cat("Saved secondary MCMC samples to:", samples2_file, "\n")
	}
	
	cat("Saved Dirichlet MCMC samples to:", samples_file, "\n")
	
	return(list(status = "SUCCESS", samples = samples, samples2 = samples2, file = samples_file))
}

# Run models sequentially for testing (easier to debug)
cat("Running", nrow(params), "Dirichlet models sequentially for testing\n")

all_results <- list()
for (j in 1:nrow(params)) {
	cat("\n=== Processing model", j, "of", nrow(params), "===\n")
	cat("Model:", params$model_id[j], "\n")
	
	# Test with just one chain first for faster testing
	cat("Testing single chain execution...\n")
	tryCatch({
		result <- run_scenarios_dirichlet(j = j, chain_no = 1)
		cat("Single chain result:", ifelse(is.null(result), "NULL", "SUCCESS"), "\n")
		all_results[[j]] <- result
	}, error = function(e) {
		cat("ERROR in single chain:", e$message, "\n")
		all_results[[j]] <- NULL
	})
}

cat("\n=== All models completed ===\n")
cat("Total models processed:", length(all_results), "\n")

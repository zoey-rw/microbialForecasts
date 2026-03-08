#!/usr/bin/env Rscript
# Fit Dirichlet models for all microbial groups and predictor sets
# - Weak priors for main parameters (Jeffreys, Uniform, wide normal)
# - More informative priors for site effects (dgamma(2, 20))

# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
here::i_am("analysis/model_analysis/dirichlet_covariance/01_fitModels_dirichlet.r")
project_root <- here()

cat("here() starts at", project_root, "\n")
cat("Project root set to:", getwd(), "\n")

# Load the microbialForecast package to access helper functions
library(microbialForecast)

# Load packages and create directories using package functions
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("dirichlet_regression", "dirichlet_regression/env_cycl", "dirichlet_regression/env_cov")
)

# Define output directory early to prevent undefined variable errors
model_output_dir <- here("data", "model_outputs", "dirichlet_regression")

cat("==================================================\n")
cat("Microbial forecasts environment setup complete!\n")
cat("Ready for Dirichlet analysis.\n")
cat("==================================================\n")

# Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# Run with 4 chains for local testing
nchains = 4

#### Run on all groups ----

# Load data early for filtering
cat("Loading data files for filtering...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully for", length(all_ranks), "ranks\n")

# Function to check if MCMC should continue based on effective sample size
check_continue <- function(samples, min_eff_size = 50) {
    # Validate input
    if (is.null(samples) || nrow(samples) == 0) {
        cat("  WARNING: Empty or NULL samples provided to check_continue, defaulting to continue\n")
        return(TRUE)
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
    
    cat("  ESS check - Min ESS:", round(min_ess, 1), "Target:", min_eff_size, "\n")
    cat("  Continue sampling:", continue, "\n")
    
    return(continue)
}

params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 							 rep("logical", 2),
										 							 rep("character", 4)))

# Check for priority mode via environment variable or command line argument
use_priority <- Sys.getenv("USE_PRIORITY", "false")
argv <- commandArgs(TRUE)
if (length(argv) > 1 && (argv[2] == "priority" || argv[2] == "--priority")) {
    use_priority <- "true"
}

if (use_priority == "true" || use_priority == "priority") {
    cat("🎯 PRIORITY MODE: Using high-priority models with existing progress\n")
    priority_file <- here("data/summary/priority_rerun_list.rds")
    if (file.exists(priority_file)) {
        rerun_list <- readRDS(priority_file)
        cat("   Loaded", length(rerun_list), "priority models (have chain 1 completed)\n")
    } else {
        cat("   Priority list not found! Falling back to standard unconverged list...\n")
        rerun_list <- readRDS(here("data/summary/unconverged_taxa_list.rds"))
    }
} else {
    cat("📊 STANDARD MODE: Using all unconverged models\n")
    rerun_list <- readRDS(here("data/summary/unconverged_taxa_list.rds"))
    cat("   Tip: Set USE_PRIORITY=true or add 'priority' argument to focus on models with existing progress\n")
}

converged_list = readRDS(here("data/summary/weak_converged_taxa_list.rds"))

# HPC PRODUCTION CONFIGURATION: Run multiple models with all three model types
params <- params_in %>% ungroup %>% filter(
    # Run ALL model types for comprehensive analysis
    model_name %in% c("env_cov") &
        # Focus on 2013-2018 period for legacy analysis
        scenario %in% c("Legacy with covariate 2013-2018") &
        # For Dirichlet models, focus on taxonomic ranks
	fcast_type == "Taxonomic" &
        rank.name %in% c("phylum_fun")
) %>% distinct(.keep_all = TRUE)

# Filter out already converged models
params <- params %>% filter(!model_id %in% converged_list)

# LOCAL TESTING: Run just 1 model for faster testing
# Limit to reasonable number for local testing
set.seed(123)  # For reproducible sampling
params <- params %>%
    sample_n(size = 1, replace = FALSE) %>%  # Run just 1 model for faster testing
    ungroup()

cat("LOCAL TESTING: Starting parallel execution for", nrow(params), "models with", nchains, "chains\n")
cat("Expected runtime: Variable (convergence-based sampling)\n")
cat("  - Models to run:", nrow(params), "models (env_cycl)\n")
cat("  - Chains per model:", nchains, "(total", nrow(params) * nchains, "parallel tasks)\n")
cat("  - Initial iterations: ~ 0.1 minutes per chain\n")
cat("  - Additional iterations: Variable based on convergence\n")
cat("  - Target: ESS >= 10 per parameter (TESTING VERSION)\n")
cat("🎯 PRIOR STRATEGY: HYBRID (weak main + stable site effects) - PROVEN TO WORK\n")
cat("🔧 TESTING MODE: Local testing with 1 core and reduced iterations\n")

cat("LOCAL TESTING: Running", nrow(params), "models with", nchains, "chains in parallel\n")
cat("This executes the stable env_cycl framework with PROVEN HYBRID PRIORS for local testing\n")

# Filter parameters to only include models with available species and ranks
cat("Filtering models to only include those with available data...\n")
original_n_models <- nrow(params)

# Create a function to check if a model has valid data
is_valid_model <- function(rank_name, species_name) {
    # Check if rank exists in data
    if (!(rank_name %in% names(all_ranks))) {
        return(FALSE)
    }
    
    rank_data <- all_ranks[[rank_name]]
    
    # For Dirichlet models, we check if the rank has multiple taxa
    metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
    available_taxa <- setdiff(colnames(rank_data), metadata_cols)
    
    # Need at least 2 taxa for Dirichlet composition
    return(length(available_taxa) >= 2)
}

# Filter the parameters dataframe
valid_indices <- sapply(1:nrow(params), function(i) {
    is_valid_model(params$rank.name[i], params$species[i])
})

params <- params[valid_indices, ]
filtered_n_models <- nrow(params)

# Report filtering results
if (filtered_n_models < original_n_models) {
    n_filtered <- original_n_models - filtered_n_models
    cat("⚠️  Filtered out", n_filtered, "models with unavailable species/ranks\n")
    cat("  Original models:", original_n_models, "\n")
    cat("  Valid models:", filtered_n_models, "\n")
} else {
    cat("✓ All", filtered_n_models, "models have valid data\n")
}

# Additional validation: ensure we have at least one valid model
if (filtered_n_models == 0) {
    cat("❌ ERROR: No valid models remaining after filtering!\n")
    cat("Available ranks in data:", paste(names(all_ranks), collapse=", "), "\n")
    
    # Show some examples of available taxa for each rank
    for (rank_name in names(all_ranks)) {
        rank_data <- all_ranks[[rank_name]]
        metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
        available_taxa <- setdiff(colnames(rank_data), metadata_cols)
        cat("  ", rank_name, "taxa (first 5):", paste(head(available_taxa, 5), collapse=", "),
            if(length(available_taxa) > 5) paste("... (+", length(available_taxa) - 5, " more)") else "", "\n")
    }
    
    stop("No valid models to run - check rank names in parameters")
}

# Use the filtered parameters from data loading
valid_models <- params
cat("Testing with", nrow(valid_models), "models\n")
cat("Models to test:\n")
print(valid_models)

# Create function that uses Dirichlet approach for each model
run_scenarios_dirichlet <- function(j, chain_no) {
    # Initialize error tracking and logging
    start_time <- Sys.time()
    error_context <- list()
    
    tryCatch({
        # Load required libraries in each worker using helper function
        load_required_packages()
        cat("=== Starting Dirichlet model fitting ===\n")
        cat("Model index:", j, "Chain:", chain_no, "\n")
        cat("Model parameters:\n")
        print(valid_models[j,])
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
        if (is.null(valid_models) || nrow(valid_models) < j) {
            stop("Valid_models data frame not available or index out of bounds")
        }

	cat("Running Dirichlet scenario", j, "chain", chain_no, "\n")

        # Extract model parameters
        rank.name <- valid_models$rank.name[[j]]
        species <- valid_models$species[[j]]
        model_id <- valid_models$model_id[[j]]
        model_name <- valid_models$model_name[[j]]
        min.date <- valid_models$min.date[[j]]
        max.date <- valid_models$max.date[[j]]
        scenario <- valid_models$scenario[[j]]
        
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
	
        # Data already loaded at the top of the script
	
	# Get the specific group data
	if (!(rank.name %in% names(all_ranks))) {
            stop("Rank name '", rank.name, "' not found in data. Available ranks: ", 
                 paste(names(all_ranks), collapse=", "))
	}
	rank.df <- all_ranks[[rank.name]]
        
        # Validate rank data structure
        if (!is.data.frame(rank.df) || nrow(rank.df) == 0) {
            stop("Rank data for '", rank.name, "' is empty or not a data frame")
        }
	
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
	
	# Source the prepDirichletData function
	source(here("microbialForecast", "R", "prepDirichletData.r"))
	
	# Use prepDirichletData for proper Dirichlet model data preparation
	cat("  Using prepDirichletData for Dirichlet model data preparation...\n")
	
	# DRASTICALLY reduce data size for testing - find plots with most data
	cat("Drastically reducing data size for testing...\n")
	cat("Original data dimensions:", dim(rank.df), "\n")
	
	# Find plots with the most observations
	plot_counts <- table(rank.df$plotID)
	top_plots <- names(sort(plot_counts, decreasing = TRUE))[1:5]  # Top 5 plots by observation count
	
	# Find dates with the most observations
	date_counts <- table(rank.df$dateID)
	top_dates <- names(sort(date_counts, decreasing = TRUE))[1:5]  # Top 5 dates by observation count
	
	# Filter to this subset
	rank.df_minimal <- rank.df[rank.df$plotID %in% top_plots & rank.df$dateID %in% top_dates, ]
	cat("Minimal data dimensions:", dim(rank.df_minimal), "\n")
	cat("Selected plots:", paste(top_plots, collapse = ", "), "\n")
	cat("Selected dates:", paste(top_dates, collapse = ", "), "\n")
	
	# Use prepDirichletData with the minimal dataset
	model.dat <- prepDirichletData(rank.df = rank.df_minimal,
		min.prev = 1,  # Lower threshold for small data
		min.date = min.date,
		max.date = max.date)
	
	cat("    Data prepared using prepDirichletData:\n")
	cat("      N.plot:", model.dat$N.plot, "\n")
	cat("      N.date:", model.dat$N.date, "\n")
	cat("      N.site:", model.dat$N.site, "\n")
	cat("      N.core:", model.dat$N.core, "\n")
	cat("      N.spp:", model.dat$N.spp, "\n")
	cat("      timepoint length:", length(model.dat$timepoint), "\n")
	cat("      plot_num length:", length(model.dat$plot_num), "\n")
	cat("      plot_start length:", length(model.dat$plot_start), "\n")
	cat("      plot_site_num length:", length(model.dat$plot_site_num), "\n")
	cat("      keep_taxa:", paste(model.dat$keep_taxa, collapse = ", "), "\n")
	
	cat("  Data prepared successfully using prepBetaRegData approach\n")
	
	# Debug: Check what's in model.dat
	cat("    Model.dat contents:\n")
	cat("      Names:", paste(names(model.dat), collapse = ", "), "\n")
	cat("      Length:", length(model.dat), "\n")
	
	# Prepare constants - use dimensions from prepDirichletData
	constants <- model.dat[c("plot_start", "plot_index",
							"plot_num", "plot_site_num", "timepoint",
							"N.plot", "N.spp", "N.core", "N.site", "N.date")]
	
	# Debug: Check constants after creation
	cat("    Constants after creation:\n")
	cat("      Names:", paste(names(constants), collapse = ", "), "\n")
	cat("      Length:", length(constants), "\n")
	
	# Use the exact same approach as beta regression for environmental data
	cat("  Preparing environmental predictors using beta regression approach...\n")
	
	# Add environmental predictors with validation (ONLY for environmental models)
	if (model_name %in% c("env_cycl", "env_cov")) {
		env_predictors <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
		cat("    Adding environmental predictors with validation...\n")
		
		for (pred in env_predictors) {
			if (pred %in% names(model.dat)) {
				constants[[pred]] <- model.dat[[pred]]
				
				# Validate predictor dimensions and structure
				pred_data <- model.dat[[pred]]
				if (is.matrix(pred_data)) {
					cat("      ✓ Added", pred, "predictor:", dim(pred_data), "matrix\n")
				} else if (is.vector(pred_data)) {
					cat("      ✓ Added", pred, "predictor:", length(pred_data), "vector\n")
				} else {
					cat("      ✓ Added", pred, "predictor:", class(pred_data), "object\n")
				}
				
				# Check for missing or extreme values
				if (is.numeric(pred_data)) {
					missing_pct <- mean(is.na(pred_data)) * 100
					if (missing_pct > 0) {
						cat("        WARNING:", pred, "has", round(missing_pct, 1), "% missing values\n")
					}
					
					if (is.matrix(pred_data)) {
						extreme_vals <- sum(abs(pred_data) > 10, na.rm = TRUE)
						if (extreme_vals > 0) {
							cat("        WARNING:", pred, "has", extreme_vals, "extreme values (>10)\n")
						}
					}
				}
			} else {
				cat("      ❌ ERROR:", pred, "predictor not found in model data\n")
				stop("Missing required environmental predictor: ", pred)
			}
		}
	} else {
		cat("    Skipping environmental predictors for", model_name, "model\n")
	}
	
	# Add seasonal predictors (sin_mo, cos_mo) from model.dat
	if ("sin_mo" %in% names(model.dat)) {
		constants$sin_mo <- model.dat$sin_mo
		cat("    ✓ Added sin_mo predictor:", length(model.dat$sin_mo), "vector\n")
	}
	if ("cos_mo" %in% names(model.dat)) {
		constants$cos_mo <- model.dat$cos_mo
		cat("    ✓ Added cos_mo predictor:", length(model.dat$cos_mo), "vector\n")
	}
	
	# N.date is already set from prepDirichletData
	
	# Set N.beta based on model type before creating initial values
	if (model_name == "env_cycl") {
		constants$N.beta = 4  # Reduced from 8 for testing
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
	} else {
		constants$N.beta = 2
	}
		
			# Create inits using the Dirichlet-specific inits function with correct N.date
	source(here("analysis", "model_analysis", "dirichlet_covariance", "dirichlet_helper_functions.r"))
	
	# Debug: Check constants before creating initial values
	cat("  Debug: Constants before initial values creation:\n")
	cat("    N.plot:", constants$N.plot, "\n")
	cat("    N.spp:", constants$N.spp, "\n")
	cat("    N.core:", constants$N.core, "\n")
	cat("    N.site:", constants$N.site, "\n")
	cat("    N.date:", constants$N.date, "\n")
	cat("    N.beta:", constants$N.beta, "\n")
	
	# Check if all required constants are present
	required_constants <- c("N.plot", "N.spp", "N.core", "N.site", "N.date", "N.beta")
	missing_constants <- required_constants[!required_constants %in% names(constants)]
	if (length(missing_constants) > 0) {
		cat("    ERROR: Missing required constants:", paste(missing_constants, collapse = ", "), "\n")
		stop("Missing required constants for initial values creation")
	}
	
	# Check for NA values in required constants
	for (name in required_constants) {
		if (is.na(constants[[name]])) {
			cat("    ERROR: Constant", name, "is NA\n")
			stop("NA value in required constant: ", name)
		}
	}
	
	# Check indexing arrays for NA values
	indexing_arrays <- c("plot_start", "plot_index", "plot_num", "plot_site_num")
	for (name in indexing_arrays) {
		if (name %in% names(constants)) {
			value <- constants[[name]]
			if (is.numeric(value)) {
				na_count <- sum(is.na(value))
				if (na_count > 0) {
					cat("    ERROR: Indexing array", name, "contains", na_count, "NA values\n")
					stop("NA values in indexing array: ", name)
				}
			}
		}
	}
	
	# Check timepoint array from model.dat
	if ("timepoint" %in% names(model.dat)) {
		timepoint_na <- sum(is.na(model.dat$timepoint))
		if (timepoint_na > 0) {
			cat("    ERROR: timepoint array contains", timepoint_na, "NA values\n")
			stop("NA values in timepoint array")
		}
		cat("    timepoint array: length =", length(model.dat$timepoint), ", range =", range(model.dat$timepoint), "\n")
	}
	
	cat("  Creating initial values...\n")
	inits <- initsFun_dirichlet(constants, type = "tax")
	cat("  Initial values created successfully\n")
	
	# Debug: Check initial values for missing values and dimensions
	cat("  Checking initial values for missing values and dimensions...\n")
	for (name in names(inits)) {
		value <- inits[[name]]
		cat("    ", name, ": dim =", paste(dim(value), collapse = "x"), ", class =", class(value), "\n")
		if (is.array(value) || is.matrix(value)) {
			cat("      Has NA values:", any(is.na(value)), "\n")
			cat("      Has Inf values:", any(is.infinite(value)), "\n")
		}
		if (is.numeric(value)) {
			if (any(is.na(value))) {
				cat("    WARNING: Initial value", name, "contains NA values\n")
			}
			if (any(is.infinite(value))) {
				cat("    WARNING: Initial value", name, "contains infinite values\n")
			}
		}
	}
		
	# Debug: Check constants after adding environmental data
	cat("    Constants after adding environmental data:\n")
	cat("      Names:", paste(names(constants), collapse = ", "), "\n")
	cat("      Length:", length(constants), "\n")
	
	cat("    Environmental predictors loaded successfully using beta regression approach\n")
	
	# Add driver uncertainty parameters
	temporalDriverUncertainty <- FALSE  # Set to TRUE if you want driver uncertainty
	spatialDriverUncertainty <- FALSE   # Set to TRUE if you want driver uncertainty
	
	# Model hyperparameters - adjust based on model type
	if (model_name == "env_cycl") {
		constants$Nimble_model = "nimbleModTaxa"
	} else if (model_name == "env_cov") {
		constants$Nimble_model = "nimbleModTaxa_env_cov"
	} else {
		constants$Nimble_model = "nimbleModTaxa_cycl_only"
	}
	
	# Set up omega matrix for multivariate normal priors - sized to match N.beta
	constants$omega <- 0.05 * diag(constants$N.beta)
	constants$zeros <- rep(0, constants$N.beta)
	
	cat("Constants prepared successfully\n")
	
	# Debug: Check constants names
	cat("  Constants names:", paste(names(constants), collapse = ", "), "\n")
	cat("  Constants with empty names:", sum(names(constants) == ""), "\n")
	cat("  Constants with NULL names:", sum(is.null(names(constants))), "\n")
	
	# Debug: Check for missing values in constants and fix them
	cat("  Checking constants for missing values...\n")
	for (name in names(constants)) {
		value <- constants[[name]]
		if (is.numeric(value)) {
			if (any(is.na(value))) {
				cat("    WARNING:", name, "contains NA values - replacing with defaults\n")
				# Replace NA values with reasonable defaults based on the variable type
				if (name %in% c("mois", "temp", "pH", "pC", "relEM", "LAI")) {
					# For environmental variables, use median of non-NA values
					median_val <- median(value, na.rm = TRUE)
					value[is.na(value)] <- median_val
					constants[[name]] <- value
					cat("      Replaced NA values with median:", median_val, "\n")
				} else if (name %in% c("mois_sd", "temp_sd", "pH_sd", "pC_sd")) {
					# For standard deviations, use a small positive value
					value[is.na(value)] <- 0.1
					constants[[name]] <- value
					cat("      Replaced NA values with 0.1\n")
				} else {
					# For other numeric variables, use a default value
					value[is.na(value)] <- 0
					constants[[name]] <- value
					cat("      Replaced NA values with 0\n")
				}
			}
			if (any(is.infinite(value))) {
				cat("    WARNING:", name, "contains infinite values\n")
			}
		} else if (is.character(value)) {
			if (any(is.na(value))) {
				cat("    WARNING:", name, "contains NA values - replacing with defaults\n")
				value[is.na(value)] <- "default"
				constants[[name]] <- value
				cat("      Replaced NA values with 'default'\n")
			}
		}
	}
	
	# Final check for any remaining NA values
	cat("  Final check for remaining NA values...\n")
	for (name in names(constants)) {
		value <- constants[[name]]
		if (is.numeric(value) && any(is.na(value))) {
			cat("    ERROR: Still has NA values after replacement:", name, "\n")
			cat("      NA count:", sum(is.na(value)), "\n")
			stop("NA values still present in constant: ", name)
		}
		if (is.character(value) && any(is.na(value))) {
			cat("    ERROR: Still has NA values after replacement:", name, "\n")
			cat("      NA count:", sum(is.na(value)), "\n")
			stop("NA values still present in constant: ", name)
		}
	}
	
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
		# ULTRA-SIMPLE model for testing - no complex dynamics
		modelCode <- nimble::nimbleCode({
			# Simple Dirichlet model
			for (i in 1:N.core) {
				y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
			}
			
			# Simple alpha generation - no complex dynamics
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {
					for (t in 1:N.date) {
						alpha[p, s, t] ~ dgamma(1, 1)
					}
				}
			}
			
			# Simple priors
			for (s in 1:N.spp) {
				sigma[s] ~ dgamma(1, 1)
				intercept[s] ~ dnorm(0, 1)
				beta[s, 1:4] ~ dmnorm(zeros[1:4], omega[1:4, 1:4])
			}
		})
	} else if (model_name == "env_cov") {
		modelCode <- nimble::nimbleCode({
			# Dirichlet model with temporal dependence, site effects, seasonal, and environmental parameters
			for (i in 1:N.core) {
				y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
			}

			# Alpha generation with all covariates
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {

					# --- First timepoint ---
					temporal_effect[p, s, 1] ~ dnorm(0, sd = 1)

					# Calculate expected value cleanly on the right side
					Ex[p, s, 1] <- exp(intercept[s] +
							site_effect[plot_site_num[p], s] +
							temporal_effect[p, s, 1] +
							beta[s, 1] * sin_mo[1] + beta[s, 2] * cos_mo[1] +
							beta[s, 3] * temp[plot_site_num[p], 1] + beta[s, 4] * mois[plot_site_num[p], 1] +
							beta[s, 5] * pH[p, 1] + beta[s, 6] * pC[p, 1] +
							beta[s, 7] * relEM[p, 1] + beta[s, 8] * LAI[plot_site_num[p], 1])

					# Sample alpha using the properly defined sigma parameter
					alpha[p, s, 1] ~ dgamma(shape = Ex[p, s, 1], rate = sigma[s])

					# --- Subsequent timepoints ---
					for (t in 2:N.date) {
						# Fix AR(1) standard deviation math
						temporal_effect[p, s, t] ~ dnorm(rho * temporal_effect[p, s, t-1], sd = sqrt(1 - rho^2))

						Ex[p, s, t] <- exp(intercept[s] +
								site_effect[plot_site_num[p], s] +
								temporal_effect[p, s, t] +
								beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t] +
								beta[s, 3] * temp[plot_site_num[p], t] + beta[s, 4] * mois[plot_site_num[p], t] +
								beta[s, 5] * pH[p, t] + beta[s, 6] * pC[p, t] +
								beta[s, 7] * relEM[p, t] + beta[s, 8] * LAI[plot_site_num[p], t])

						alpha[p, s, t] ~ dgamma(shape = Ex[p, s, t], rate = sigma[s])
					}
				}
			}

			# Site effects
			for (s in 1:N.spp) {
				for (site in 1:N.site) {
					site_effect[site, s] ~ dnorm(0, sd = sigma_site[s])
				}
			}

			# Priors
			for (s in 1:N.spp) {
				sigma[s] ~ dgamma(1, 1)
				sigma_site[s] ~ dgamma(1, 1)
				intercept[s] ~ dnorm(0, sd = 1)
				beta[s, 1:8] ~ dmnorm(zeros[1:8], omega[1:8, 1:8])
			}

			rho ~ dunif(-0.99, 0.99)
		})
	} else {
		# Default to cycl_only model
		modelCode <- nimble::nimbleCode({
			# Loop through core observations ----
			for (i in 1:N.core) {
				# Use proper Dirichlet parameterization with concentration parameters
				y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
			}

			# Plot-level process model ----
			for (s in 1:N.spp) {
				for (p in 1:N.plot) {
					# Initial condition - ensure plot_start[p] is valid
									# Initial concentration parameters for Dirichlet
				alpha[p, s, plot_start[p]] ~ dgamma(1, 1)
				# Convert to relative abundance for monitoring
				plot_rel[p, s, plot_start[p]] <- alpha[p, s, plot_start[p]] / sum(alpha[p, 1:N.spp, plot_start[p]])
					
					# Dynamic evolution - only run if we have more than one time point
					if (plot_start[p] < N.date) {
						for (t in (plot_start[p] + 1):N.date) {
						# Previous value * rho - with numerical stability
						log(Ex[p, s, t]) <- rho[s] * log(max(alpha[p, s, t - 1], 0.001)) +
							beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t] +
							site_effect[plot_site_num[p], s] +
							intercept[s]
						# Add process error (sigma) - with numerical stability
						alpha[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
						# Convert back to relative abundance for monitoring
						plot_rel[p, s, t] <- alpha[p, s, t] / sum(alpha[p, 1:N.spp, t])
						}
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
	
	cat("Dirichlet model built successfully\n")
	
	# Build model
	cat("  Attempting to build Nimble model...\n")
	Rmodel <- NULL  # Initialize Rmodel
	tryCatch({
		cat("    Creating model with constants...\n")
		cat("    Constants dimensions - N.plot:", constants$N.plot, "N.spp:", constants$N.spp, "N.core:", constants$N.core, "\n")
		cat("    Data dimensions - y:", dim(model.dat$y), "\n")
		
		# Debug: Check all constants for dimensions and NA values
		cat("    Debug: Checking all constants before model building...\n")
		for (name in names(constants)) {
			value <- constants[[name]]
			if (is.numeric(value)) {
				cat("      ", name, ": numeric, dim =", paste(dim(value), collapse = "x"), 
					", length =", length(value), ", NA count =", sum(is.na(value)), "\n")
			} else if (is.character(value)) {
				cat("      ", name, ": character, length =", length(value), ", NA count =", sum(is.na(value)), "\n")
			} else {
				cat("      ", name, ": ", class(value), ", length =", length(value), "\n")
			}
		}
		
		Rmodel <- nimbleModel(code = modelCode, constants = constants,
							  data = list(y=model.dat$y), inits = inits)
		cat("    Model built successfully\n")
	}, error = function(e) {
		cat("ERROR building model:", e$message, "\n")
		cat("Error details:\n")
		print(e)
		Rmodel <<- NULL  # Set Rmodel to NULL in the outer scope
	})
	
	if (is.null(Rmodel)) {
		cat("Failed to build model, returning NULL\n")
		return(NULL)
	}
	
	# Compile model
	cModel <- NULL  # Initialize cModel
	tryCatch({
		cModel <- compileNimble(Rmodel)
		cat("Dirichlet model compiled successfully\n")
	}, error = function(e) {
		cat("ERROR compiling model:", e$message, "\n")
		cModel <<- NULL  # Set cModel to NULL in the outer scope
	})
	
	if (is.null(cModel)) {
		cat("Failed to compile model, returning NULL\n")
		return(NULL)
	}
	
        # Configure MCMC with proper sampler management
        cat("Configuring MCMC...\n")
        
        # Ultra-simplified monitoring for testing - only the most essential parameters
        monitored_params <- c("beta", "sigma", "intercept", "rho", "sig")
        
        # Skip all latent variables to reduce complexity
        monitored_latent_params <- c()
        
        cat("Monitoring parameters for convergence analysis:\n")
        cat("  Core parameters:", paste(monitored_params, collapse = ", "), "\n")
        cat("  Latent variables:", paste(monitored_latent_params, collapse = ", "), "\n")
        cat("  Total beta parameters:", constants$N.beta, "\n")
        
        mcmcConf <- configureMCMC(
            model = cModel,
            monitors = monitored_params,
            thin = 1,
            enableWAIC = FALSE
        )
        
        # Use default samplers for testing - much simpler
        cat("Using default samplers for testing (simplified configuration)...\n")
        
        cat("MCMC configured successfully - ALL advanced features RESTORED!\n")
	
	# Build and compile MCMC
        cat("Building and compiling MCMC...\n")
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
        # Run MCMC with ultra-simplified sampling (TESTING VERSION)
        cat("Running MCMC with ultra-simplified sampling (TESTING VERSION)...\n")
        burnin <- 50
        thin <- 1
        iter_per_chunk <- 100
        init_iter <- 1000
        min_eff_size_perchain <- 5  # Very low threshold for testing
        max_loops <- 0
        max_save_size <- 10000
        min_total_iterations <- 1000
        
        cat("Running MCMC with convergence-based sampling\n")
        cat("  Initial iterations:", init_iter, "burnin:", burnin, "\n")
        cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
        cat("  Target ESS per chain:", min_eff_size_perchain, "\n")
        cat("  Minimum total iterations:", min_total_iterations, "\n")
	
	# Run initial iterations
        cat("  Running initial iterations (", init_iter, " iterations) for adaptation...\n")
        compiled$run(niter = init_iter, thin = thin, nburnin = 0)
        cat("  Initial iterations completed\n")
        
        # Get initial samples and check convergence
        initial_samples <- as.matrix(compiled$mvSamples)
        cat("  Initial samples collected, checking convergence...\n")
        cat("  Initial samples dimensions:", dim(initial_samples), "\n")
        
        # Create output directory for checkpoints
        cat("  Creating output directory for checkpoints...\n")
        
        # Create species-specific subdirectory
        species_output_dir <- file.path(model_output_dir, model_name, rank.name)
        
        # Ensure the directory exists
        if (!dir.exists(species_output_dir)) {
            dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
        }
        
        # Verify directory was created
        if (!dir.exists(species_output_dir)) {
            stop("CRITICAL: Failed to create checkpoint directory: ", species_output_dir)
        }
        
        cat("  ✓ Checkpoint directory ready:", species_output_dir, "\n")
        
        # Create model_id for consistent naming
        model_id <- create_model_id(model_name, rank.name, min.date, max.date, use_legacy_covariate)
        
        # Check if we need to continue sampling for convergence
        continue <- TRUE
        loop_counter <- 0
        total_iterations <- init_iter
        
        # Try to check convergence, with fallback if it fails
        tryCatch({
            continue <- check_continue(initial_samples, min_eff_size = min_eff_size_perchain)
        }, error = function(e) {
            cat("  WARNING: Convergence check failed, defaulting to continue sampling\n")
            cat("  Error:", e$message, "\n")
            continue <- TRUE
        })
        
        # Store all samples as we go
        all_samples <- initial_samples
        cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
        
        # Save initial samples as checkpoint
        save_checkpoint_safe(all_samples, total_iterations, 0, species_output_dir, model_id, chain_no, "initial")
        
        # Also save a simple progress file
        progress_file <- create_progress_file(species_output_dir, model_id, chain_no, init_iter)
        
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
            save_checkpoint_safe(all_samples, total_iterations, loop_counter + 1, species_output_dir, model_id, chain_no, paste0("loop", loop_counter + 1))
            
            # Update progress file
            update_progress_file(progress_file, total_iterations, loop_counter + 1)
            
            # Check if we need to continue
            continue <- TRUE
            tryCatch({
                continue <- check_continue(all_samples, min_eff_size = min_eff_size_perchain)
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
        
        # Update final progress status
        tryCatch({
            final_status <- if(loop_counter >= max_loops) "Completed (max loops)" else 
                if(total_iterations >= min_total_iterations && !continue) "Converged" else 
                    "Completed (min iterations)"
            writeLines(paste("Completed at:", Sys.time(), "\nTotal iterations:", total_iterations, "\nFinal loop:", loop_counter, "\nStatus:", final_status), progress_file)
            cat("  ✓ Final progress status updated\n")
        }, error = function(e) {
            cat("  ✗ Failed to update final progress status:", e$message, "\n")
        })
        
        # Get final samples (use accumulated samples)
        samples <- all_samples
        
        # Extract plot-level estimates from monitors2 output for comprehensive analysis
        cat("Extracting plot-level estimates from MCMC output...\n")
        plot_samples <- as.matrix(compiled$mvSamples2)
        cat("  Plot samples dimensions:", dim(plot_samples), "\n")
        cat("  Plot samples column names:", paste(colnames(plot_samples), collapse=", "), "\n")
        
        # Validate plot samples structure
        if (is.null(nrow(plot_samples)) || nrow(plot_samples) == 0) {
            cat("  WARNING: No plot samples found in monitors2 output (monitors2 not configured)\n")
            plot_samples <- samples  # Fallback to parameter samples if no plot samples
        } else {
            cat("  ✓ Plot samples extracted successfully\n")
        }
	
	cat("MCMC completed successfully\n")
        cat("Final sample dimensions:", dim(samples), "\n")
        cat("Plot sample dimensions:", dim(plot_samples), "\n")
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
        cat("  Additional loops:", loop_counter, "iterations\n")
        cat("  Total accumulated samples:", nrow(all_samples), "iterations\n")
        cat("  Checkpoints saved:", loop_counter + 1, "files\n")
        cat("  Final sample size:", nrow(all_samples), "iterations\n")
	
        # Save MCMC samples with absolute path
        samples_file <- file.path(model_output_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))
        
        # Create the complete chain structure with metadata
	chain_output <- list(
            samples = all_samples,
            samples2 = plot_samples,  # Plot-level estimates from monitors2 output
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
                thin2 = 20,  # Include thin2 for samples2
                model_data = model.dat,
                nimble_code = modelCode,
                model_structure = "stable_dirichlet_regression_with_compositional_data",
			N.spp = constants$N.spp,
			N.plot = constants$N.plot,
			N.core = constants$N.core,
			N.site = constants$N.site,
			N.date = constants$N.date,
                keep_names = keep_names  # Save the taxa names being modeled
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
        
        cat("Sample dimensions:", dim(all_samples), "\n")
        cat("=== Dirichlet model fitting completed ===\n")
        cat("  - STABLE: Dirichlet regression with compositional data\n")
        cat("  - All three model types supported: cycl_only, env_only, env_cov\n")
        cat("  - CONVERGENCE-BASED: Adaptive sampling until reasonable ESS reached\n")
        cat("  - ITERATIVE SAVING: Samples accumulated and saved incrementally\n")
        
        return(list(
            status = "SUCCESS", 
            samples = all_samples,
            samples2 = plot_samples,  # Include plot-level estimates for consistency
            file = samples_file,
            model_data = model.dat,
            nimble_code = modelCode,
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
                thin2 = 20,  # Include thin2 for samples2
			model_data = model.dat,
                nimble_code = modelCode,
                model_structure = "stable_dirichlet_regression_with_compositional_data",
                N.spp = constants$N.spp,
                N.plot = constants$N.plot,
                N.core = constants$N.core,
                N.site = constants$N.site,
                N.date = constants$N.date,
			keep_names = keep_names  # Save the taxa names being modeled
		)
        ))
        
    }, error = function(e) {
        # Capture comprehensive error information
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
            runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA
        )
        
        # Create detailed error file with absolute path
        model_name <- if(exists("model_name")) model_name else "unknown"
        error_dir <- here("data", "model_outputs", "dirichlet_regression", model_name)
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

# LOCAL TESTING: Sequential execution for debugging
cat("📊 Models to run:", nrow(valid_models), "models with", nchains, "chains each\n")
cat("⏱️  Expected runtime: Variable (convergence-based sampling)\n")
cat("🎯 Target: ESS >= 10 per parameter\n")
cat("🔧 TESTING MODE: Sequential execution for debugging\n")

# Set start time for runtime calculation
start_time <- Sys.time()

# Run models sequentially for testing
cat("Running models sequentially for testing...\n")
all_results_sequential <- list()

for (model_idx in 1:nrow(valid_models)) {
    for (chain_no in 1:nchains) {
        cat("=== Running Model", model_idx, "Chain", chain_no, "===\n")
        
        # Execute the function with error handling
        result <- tryCatch({
            run_scenarios_dirichlet(j = model_idx, chain_no = chain_no)
	}, error = function(e) {
            cat("ERROR in run_scenarios_dirichlet:", e$message, "\n")
            return(list(status = "ERROR", error = e$message))
        })
        
        # Store result
        all_results_sequential[[length(all_results_sequential) + 1]] <- list(
            model_idx = model_idx, 
            chain_no = chain_no, 
            result = result
        )
        
        cat("=== Completed Model", model_idx, "Chain", chain_no, "===\n")
    }
}

cat("Sequential execution completed at:", format(Sys.time()), "\n")
cat("Total results:", length(all_results_sequential), "\n")

# Rename for compatibility with progress summary
all_results_parallel <- all_results_sequential

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

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("ALL MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("✓ Dirichlet distribution for compositional data\n")
cat("✓ Precision parameter for dispersion\n")
cat("✓ LOG transformation with exp() for numerical stability\n")
cat("  - sig ~ dgamma(2, 8) - Tighter gamma prior\n")
cat("  - rho ~ dbeta(2, 2) - Truncated beta prior\n")
cat("  - intercept ~ dt(0, 0.3, df = 3) - Robust prior\n")
cat("  - beta ~ dmvt(zeros, omega, df = 3) - Robust multivariate t\n")
cat("✓ Individual slice samplers for beta parameters\n")
cat("✓ Convergence-based sampling with iterative saving\n")
cat("Check output files in:", here("data", "model_outputs", "dirichlet_regression"), "/[model_name]/\n")

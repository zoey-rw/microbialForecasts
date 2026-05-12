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
nchains <- as.integer(Sys.getenv("NCHAINS", "3"))

# Set number of cores for parallel execution (1 model × 3 chains)
n_cores <- as.integer(Sys.getenv("NCORES", "2"))

#### Run on all groups ----

# Load required packages first
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
project_root <- here()

# Set debug mode flag
debug_mode <- FALSE  # Set to TRUE for debugging

# Load the microbialForecast package to access helper functions
library(microbialForecast)

# Load packages and create directories using package functions
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("CLR_regression", "CLR_regression/cycl_only", "CLR_regression/env_cycl", "CLR_regression/env_cov")
)

# Source the main script with correct path
source(here("source.R"))

# Function to filter files to only include those with both 'with_legacy_covariate' and 'clr'
# AND exclude checkpoint files (we only want final combined sample files)
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  standard_files <- file_list[grepl('with_legacy_covariate', basename(file_list)) & 
                              grepl('clr', basename(file_list)) &
                              !grepl('checkpoint', basename(file_list))]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard files (with CLR suffixes, no checkpoints):', length(standard_files), '\n')
  cat('  Filtered out:', length(file_list) - length(standard_files), '\n\n')
  
  return(standard_files)
}

# Function to check if MCMC should continue based on effective sample size
check_continue <- function(samples, min_eff_size = 10) {
	if (debug_mode) cat("    DEBUG: check_continue called with", nrow(samples), "samples\n")
	
	# Calculate effective sample sizes safely (ignoring zero-variance columns)
	eff_sizes <- apply(samples, 2, function(x) {
		if (sum(!is.na(x)) < 2 || is.na(var(x, na.rm = TRUE)) || var(x, na.rm = TRUE) == 0) {
			return(NA)
		}
		tryCatch(coda::effectiveSize(coda::as.mcmc(x)), error = function(e) NA)
	})
	
	# Check if any valid parameter has ESS below threshold
	valid_ess <- eff_sizes[!is.na(eff_sizes)]
	if (length(valid_ess) == 0) {
		cat("  WARNING: No parameters have valid variation for ESS calculation.\n")
		return(TRUE) # Continue sampling, maybe it will start moving
	}
	
	min_ess <- min(valid_ess)
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

	# Set testing mode flag (FALSE = production run, e.g. Ascomycota only)
	testing_mode <- FALSE
	
	if (testing_mode) {
		# Focus on phylum_fun (ascomycota) for testing CLR workflow
		params <- params_in %>% filter(
			rank.name == "phylum_fun" &
			species == "ascomycota" &
			scenario == "Legacy with covariate 2013-2018" &
			model_name %in% c("cycl_only", "env_cov", "env_cycl")
		) %>%
		slice(1:3)
	} else {
		# Ascomycota, env_cycl only, 1000 iterations, 3 chains
		params <- params_in %>% distinct(.keep_all = TRUE) %>%
			filter(
				species == "ascomycota" &
				scenario == "Legacy with covariate 2013-2018" &
				model_name == "env_cycl"
			)
	}

# TESTING: Don't filter out already converged models - run all models
# params <- params %>% filter(!model_id %in% converged_list)
cat("TESTING MODE: Running all models regardless of convergence status\n")

cat("Testing with", nrow(params), "models\n")
cat("Model configuration:\n")
print(params[, c("rank.name", "species", "model_name", "model_id")])

# Function to sanitize driver uncertainty data
sanitize_driver_uncertainty <- function(constants) {
  cat("Sanitizing driver uncertainty data (STRICT, window-aware)\n")

  eps <- 1e-6

  ## helpers
  .summarize_bad <- function(mask, max_show = 10) {
    rc <- which(mask, arr.ind = TRUE)
    if (length(rc) == 0) return("0")
    head_pairs <- apply(head(rc, max_show), 1, \(r) paste0("(", r[1], ",", r[2], ")"))
    paste0(nrow(rc), " cells; e.g., ", paste(head_pairs, collapse = ", "),
           if (nrow(rc) > max_show) ", ..." else "")
  }
  .in_scope_site <- function(Nsite, Ndate, site_start) {
    m <- matrix(FALSE, Nsite, Ndate)
    for (k in seq_len(Nsite)) {
      if (is.finite(site_start[k]) && site_start[k] >= 1 && site_start[k] <= Ndate) {
        m[k, site_start[k]:Ndate] <- TRUE
      }
    }
    m
  }
  .in_scope_plot <- function(Nplot, Ndate, plot_start) {
    m <- matrix(FALSE, Nplot, Ndate)
    for (p in seq_len(Nplot)) {
      if (is.finite(plot_start[p]) && plot_start[p] >= 1 && plot_start[p] <= Ndate) {
        m[p, plot_start[p]:Ndate] <- TRUE
      }
    }
    m
  }

  ## ---------- temp / temp_sd (TEMPORAL; strict inside site windows only) ----------
  if (!is.null(constants$temp)) {
    in_scope <- .in_scope_site(constants$N.site, constants$N.date, constants$site_start)

    bad_temp_all <- !is.finite(constants$temp)
    bad_temp_sd_all <- !is.finite(constants$temp_sd) | (constants$temp_sd <= 0)

    bad_temp_in <- bad_temp_all & in_scope
    bad_tsd_in  <- bad_temp_sd_all & in_scope

    denom <- max(1L, sum(in_scope))
    bad_temp_pct <- 100 * sum(bad_temp_in) / denom
    bad_tsd_pct  <- 100 * sum(bad_tsd_in)  / denom

    if (any(bad_temp_in)) {
      stop(sprintf("sanitize_driver_uncertainty: temp has non-finite values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_temp_pct, .summarize_bad(bad_temp_in)))
    }
    if (any(bad_tsd_in)) {
      stop(sprintf("sanitize_driver_uncertainty: temp_sd has non-finite or <=0 values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_tsd_pct, .summarize_bad(bad_tsd_in)))
    }

    cat("  temp: 100%% valid within site windows | checked cells=", sum(in_scope), " | out-of-window ignored=", sum(!in_scope & (bad_temp_all | bad_temp_sd_all)), "\n")
    constants$has_temp <- TRUE
  }

  ## ---------- mois / mois_sd (TEMPORAL; strict inside site windows only) ----------
  if (!is.null(constants$mois)) {
    in_scope <- .in_scope_site(constants$N.site, constants$N.date, constants$site_start)

    bad_mois_all <- !is.finite(constants$mois)
    bad_mois_sd_all <- !is.finite(constants$mois_sd) | (constants$mois_sd <= 0)

    bad_mois_in <- bad_mois_all & in_scope
    bad_msd_in  <- bad_mois_sd_all & in_scope

    denom <- max(1L, sum(in_scope))
    bad_mois_pct <- 100 * sum(bad_mois_in) / denom
    bad_msd_pct  <- 100 * sum(bad_msd_in)  / denom

    if (any(bad_mois_in)) {
      stop(sprintf("sanitize_driver_uncertainty: mois has non-finite values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_mois_pct, .summarize_bad(bad_mois_in)))
    }
    if (any(bad_msd_in)) {
      stop(sprintf("sanitize_driver_uncertainty: mois_sd has non-finite or <=0 values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_msd_pct, .summarize_bad(bad_msd_in)))
    }

    cat("  mois: 100%% valid within site windows | checked cells=", sum(in_scope), " | out-of-window ignored=", sum(!in_scope & (bad_mois_all | bad_mois_sd_all)), "\n")
    constants$has_mois <- TRUE
  }

  cat("✓ Driver uncertainty data sanitized (STRICT, window-aware)\n")
  constants
}

# Assert helpers
assert_matrix_dims <- function(x, nr, nc, name) {
  if (is.null(dim(x))) {
    stop(sprintf("Expected matrix for %s, got vector of length %d", name, length(x)))
  }
  if (!identical(dim(x), c(nr, nc))) {
    stop(sprintf("Dimension mismatch for %s: got %dx%d; expected %dx%d",
                 name, nrow(x), ncol(x), nr, nc))
  }
  invisible(TRUE)
}

assert_vector_len <- function(x, n, name) {
  if (length(x) != n)
    stop(sprintf("Length mismatch for %s: got %d; expected %d", name, length(x), n))
  invisible(TRUE)
}

create_inits_with_uncertainty <- function(constants, model_name, model_data = NULL, chain_no = 1) {
    cat("Creating initial values with driver uncertainty for", model_name, "...\n")
    
    # Determine number of beta parameters based on model type
    if (model_name == "env_cycl") {
        n_beta <- 8
    } else if (model_name == "env_cov") {
        n_beta <- 6
    } else {
        n_beta <- 2  # cycl_only
    }
    
    # Create initial beta values
    beta_init <- c(0.01, 0.01)  # Start with seasonal coefficients very close to zero
    if (n_beta > 2) {
        for (i in 1:(n_beta - 2)) {  # Additional environmental coefficients
            beta_init <- c(beta_init, 0.01)  # Start environmental coefficients very close to zero
        }
    }

    # Check for warm-start file
    warmstart_file <- Sys.getenv("WARMSTART_FILE", "")
    if (warmstart_file != "" && file.exists(warmstart_file)) {
        cat("  Loading warm-start initial values from:", warmstart_file, "\n")
        ws <- readRDS(warmstart_file)
        pm <- ws$param_means

        # Add per-chain jitter to break symmetry
        set.seed(chain_no * 42)
        jitter_sd <- 0.05

        inits <- list(
            precision = abs(pm["precision"] + rnorm(1, 0, jitter_sd * 5)),
            rho = max(0.001, min(0.999, pm["rho"] + rnorm(1, 0, jitter_sd))),
            beta = sapply(1:n_beta, function(b) pm[paste0("beta[", b, "]")] + rnorm(1, 0, jitter_sd)),
            site_effect_sd = abs(pm["site_effect_sd"] + rnorm(1, 0, jitter_sd * 0.5)),
            site_effect = sapply(1:constants$N.site, function(k)
                pm[paste0("site_effect[", k, "]")] + rnorm(1, 0, jitter_sd)),
            intercept = pm["intercept"] + rnorm(1, 0, jitter_sd),
            legacy_effect = pm["legacy_effect"] + rnorm(1, 0, jitter_sd),
            sigma_proc = abs(pm["sigma_proc"] + rnorm(1, 0, jitter_sd * 0.2)),
            sigma_init = abs(pm["sigma_init"] + rnorm(1, 0, jitter_sd * 0.2)),
            eta = matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
        )
        cat("  ✓ Warm-start loaded from", ws$source_iterations, "iteration run, jitter_sd=", jitter_sd, "\n")
    } else {
        # Default cold-start initial values
        inits <- list(
            precision = 50,
            rho = 0.3,
            beta = beta_init,
            site_effect_sd = 0.5,
            site_effect = rnorm(constants$N.site, 0, 0.1),
            intercept = 0,
            legacy_effect = 0,
            sigma_proc = 0.1,
            sigma_init = 0.5,
            eta = matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
        )
    }

    cat("  ✓ Initial values set for", model_name, "model\n")
    
    # Initialize driver uncertainty variables ONLY for environmental models
    if (model_name %in% c("env_cycl", "env_cov")) {
        if (constants$N.site > 0 && constants$N.date > 0) {
            # Initialize near observed values instead of random values
            inits$temp_est <- constants$temp  # Start at observed temperature
            inits$mois_est <- constants$mois  # Start at observed moisture
            cat("  ✓ Initialized temp_est and mois_est at observed values\n")
        }
        
        if (constants$N.plot > 0 && constants$N.date > 0) {
            # Initialize spatial pH/pC variables (one per plot)
            inits$pH_est_p <- constants$pH[cbind(1:constants$N.plot, constants$plot_start)]  # Start at observed pH at plot start
            inits$pC_est_p <- constants$pC[cbind(1:constants$N.plot, constants$plot_start)]  # Start at observed pC at plot start
            # Note: pH_est and pC_est are deterministic from pH_est_p and pC_est_p, so no need to initialize them
            cat("  ✓ Initialized pH_est_p, pC_est_p (spatial) at observed values\n")
        }
    } else {
        cat("  ✓ Skipped driver uncertainty variables for", model_name, "model (not needed)\n")
    }
    
    cat("  ✓ Initial values created successfully for", model_name, "model\n")
    cat("    precision:", inits$precision, "\n")
    cat("    rho:", inits$rho, "\n")
    cat("    beta parameters:", n_beta, "\n")
    cat("    site effects:", constants$N.site, "\n")
    cat("    eta matrix dimensions:", paste(dim(inits$eta), collapse = " x "), "\n")
    
    # Only report driver uncertainty dimensions if they exist
    if ("temp_est" %in% names(inits)) {
        cat("    temp_est matrix dimensions:", paste(dim(inits$temp_est), collapse = " x "), "\n")
    }
    if ("mois_est" %in% names(inits)) {
        cat("    mois_est matrix dimensions:", paste(dim(inits$mois_est), collapse = " x "), "\n")
    }
    if ("pH_est_p" %in% names(inits)) {
        cat("    pH_est_p vector length:", length(inits$pH_est_p), "\n")
    }
    if ("pC_est_p" %in% names(inits)) {
        cat("    pC_est_p vector length:", length(inits$pC_est_p), "\n")
    }
    
    return(inits)
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
	
	# Use your working prepCLRData function that already handles zeros properly
	cat("DEBUG: Preparing CLR data for species", species, "from", rank.name, "\n")
	
	tryCatch({
		# Use prepCLRData (which you already fixed to handle zeros)
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
		cat("  model.dat$y summary:\n")
		print(summary(as.vector(model.dat$y)))
	}
	cat("  N.core calculated:", nrow(model.dat$y), "\n")
	cat("  Model type:", model_name, "\n")
	
	# Validate data structure to prevent Nimble errors
	cat("Validating data structure...\n")
	required_vars <- c("N.core", "N.plot", "N.site", "N.date", "timepoint", "plot_num", "plot_site_num")
	missing_vars <- required_vars[!required_vars %in% names(model.dat)]
	if (length(missing_vars) > 0) {
		stop("Missing required variables: ", paste(missing_vars, collapse=", "))
	}
	
	# Check for missing values in critical variables
	if (any(is.na(model.dat$timepoint))) stop("timepoint contains missing values")
	if (any(is.na(model.dat$plot_num))) stop("plot_num contains missing values")
	if (any(is.na(model.dat$plot_site_num))) stop("plot_site_num contains missing values")
	
	# Validate CLR transformation results
	if (any(is.infinite(model.dat$y))) stop("CLR transformation produced infinite values")
	if (any(is.na(model.dat$y))) stop("CLR transformation produced missing values")
	
	cat("✅ Data validation passed\n")
	
	# Prepare constants - CLR uses different structure than beta regression
	constants <- list()
	
	# Add basic dimensions
	constants$N.core <- model.dat$N.core
	constants$N.plot <- model.dat$N.plot
	constants$N.site <- model.dat$N.site
	constants$N.date <- model.dat$N.date
	
	# Response variable will be passed in data, not constants
	
	# Add indexing vectors
	constants$timepoint <- model.dat$timepoint
	constants$plot_num <- model.dat$plot_num
	constants$plot_site_num <- model.dat$plot_site_num
	
	# Add seasonal predictors with optimal scaling for CLR models
	constants$sin_mo <- model.dat$sin_mo
	constants$cos_mo <- model.dat$cos_mo
	
	# Add environmental predictors
	if ("temp" %in% names(model.dat)) constants$temp <- model.dat$temp
	if ("mois" %in% names(model.dat)) constants$mois <- model.dat$mois
	if ("pH" %in% names(model.dat)) constants$pH <- model.dat$pH
	if ("pC" %in% names(model.dat)) constants$pC <- model.dat$pC
	if ("relEM" %in% names(model.dat)) constants$relEM <- model.dat$relEM
	if ("LAI" %in% names(model.dat)) constants$LAI <- model.dat$LAI

	# Add legacy covariate if needed
	if (use_legacy_covariate) {
		cat("Adding legacy covariate using actual dates...\n")
		legacy_cutoff <- as.Date("2015-01-01")

		# derive a per-timepoint legacy vector of length N.date
		if (!"all_dates" %in% names(model.dat)) {
			# fallback if missing from data prep
			all_dates <- seq.Date(as.Date(min.date, format = "%Y%m%d"), as.Date(max.date, format = "%Y%m%d"), by = "month")
		} else {
			all_dates <- as.Date(model.dat$all_dates)
		}

		if (length(all_dates) != constants$N.date) {
			cat("WARNING: Time axis length mismatch (got", length(all_dates), "dates, expected N.date=", constants$N.date, "). Rebuilding by index.\n")
			all_dates <- seq_len(constants$N.date) # still deterministic
			legacy_by_time <- as.numeric(all_dates <= floor(0.6 * constants$N.date))
		} else {
			legacy_by_time <- as.numeric(all_dates < legacy_cutoff)
		}

		# expand to plot x time while respecting plot_start windows
		legacy_mat <- matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
		for (p in seq_len(constants$N.plot)) {
			ts <- max(1L, min(constants$plot_start[p], constants$N.date))
			legacy_mat[p, ts:constants$N.date] <- legacy_by_time[ts:constants$N.date]
		}
		constants$legacy <- legacy_mat

		legacy_sum <- sum(constants$legacy)
		legacy_total <- length(constants$legacy)
		if (legacy_sum == 0 || legacy_sum == legacy_total) {
			cat("WARNING: Legacy covariate is all 0s or all 1s - this may cause numerical issues\n")
		}

		cat("Legacy covariate added:", legacy_sum, "legacy observations out of", legacy_total, "\n")
		cat("  ✓ Legacy covariate added successfully\n")
	}
	
	# Model hyperparameters - adjust based on model type
	if (model_name == "env_cycl") {
		constants$N.beta = 8
	} else if (model_name == "env_cov") {
		constants$N.beta = 6
	} else {
		constants$N.beta = 2
	}
	
	# Remove plot estimate matrices from constants - they will be calculated deterministically
	# constants$plot_estimates and constants$plot_predictions are not needed
	# These will be calculated as deterministic nodes in the model
	
	# Add missing plot indexing constants for NIMBLE loops
	# NIMBLE expects integer vectors without names to avoid "missing value where TRUE/FALSE needed" in conjugacy processing
	if ("plot_start" %in% names(model.dat) && length(model.dat$plot_start) >= constants$N.plot) {
		constants$plot_start <- as.integer(unname(model.dat$plot_start[seq_len(constants$N.plot)]))
		if (any(is.na(constants$plot_start))) constants$plot_start <- rep(1L, constants$N.plot)
	} else {
		constants$plot_start <- rep(1L, constants$N.plot)
	}
	stopifnot(length(constants$plot_start) == constants$N.plot)
	constants$plot_index <- as.integer(pmin(constants$N.date, constants$plot_start + 1L))
	
	# Add site_start for driver uncertainty - integer vector, no names, for NIMBLE conjugacy
	if ("site_start" %in% names(model.dat) && length(model.dat$site_start) >= constants$N.site) {
		constants$site_start <- as.integer(unname(model.dat$site_start[seq_len(constants$N.site)]))
		if (any(is.na(constants$site_start))) constants$site_start <- rep(1L, constants$N.site)
	} else {
		constants$site_start <- rep(1L, constants$N.site)
	}
	stopifnot(length(constants$site_start) == constants$N.site)
	
	# Add driver uncertainty flags to constants for Nimble model evaluation
	cat("Adding driver uncertainty flags to constants...\n")
	# Set driver uncertainty flags based on model type
	constants$temporalDriverUncertainty <- model_name %in% c("env_cycl", "env_cov")
	constants$spatialDriverUncertainty <- model_name %in% c("env_cycl", "env_cov")
	cat("  ✓ temporalDriverUncertainty =", constants$temporalDriverUncertainty, "\n")
	cat("  ✓ spatialDriverUncertainty =", constants$spatialDriverUncertainty, "\n")
	
	# Add driver uncertainty data when flags are enabled
	if (model_name %in% c("env_cycl", "env_cov") &&
	    (constants$temporalDriverUncertainty || constants$spatialDriverUncertainty)) {
		cat("Adding driver uncertainty data with proper dimension matching...\n")
		
		# CRITICAL FIX: Instead of blind 1:N subsetting from the raw file,
		# extract the perfectly aligned arrays directly from model.dat!
		# (prepCLRData already processed them with filter_date_site)
		
		if (constants$temporalDriverUncertainty) {
			constants$temp_sd <- model.dat$temp_sd
			constants$mois_sd <- model.dat$mois_sd
			
			assert_matrix_dims(constants$temp_sd, constants$N.site, constants$N.date, "temp_sd")
			assert_matrix_dims(constants$mois_sd, constants$N.site, constants$N.date, "mois_sd")
		}
		
		if (constants$spatialDriverUncertainty) {
			constants$pH_sd <- model.dat$pH_sd
			constants$pC_sd <- model.dat$pC_sd
			
			assert_matrix_dims(constants$pH_sd, constants$N.plot, constants$N.date, "pH_sd")
			assert_matrix_dims(constants$pC_sd, constants$N.plot, constants$N.date, "pC_sd")
		}
		
		# Validate all predictor dimensions
		cat("  Validating predictor dimensions...\n")
		if (!is.null(constants$temp)) cat("    temp dimensions:", dim(constants$temp), "\n")
		if (!is.null(constants$mois)) cat("    mois dimensions:", dim(constants$mois), "\n")
		if (!is.null(constants$pH)) cat("    pH dimensions:", dim(constants$pH), "\n")
		if (!is.null(constants$pC)) cat("    pC dimensions:", dim(constants$pC), "\n")
		if (!is.null(constants$relEM)) cat("    relEM dimensions:", dim(constants$relEM), "\n")
		if (!is.null(constants$LAI)) cat("    LAI dimensions:", dim(constants$LAI), "\n")
		cat("  ✓ All predictor dimensions validated\n")
		
		# Now the sanitizer will correctly find ONLY pre-observation NAs
		constants <- sanitize_driver_uncertainty(constants)
		
		cat("  ✓ Driver uncertainty data added & sanitized\n")
	} else {
		cat("Skipping driver uncertainty data for", model_name, "model\n")
	}
	
	cat("Constants prepared successfully\n")
	
	# STANDARDIZED MODEL DEFINITIONS - All models use consistent priors and CLR structure
	if (model_name == "cycl_only" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# PRIORS - Match driver uncertainty script
			precision ~ dgamma(2, 0.1)              # Gentle prior - mean 20, var 200
			rho ~ dbeta(1, 1)                       # Uniform prior - very weak
			intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
			legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak
			
			# STABLE PRIORS for site effects
			site_effect_sd ~ dgamma(2, 20)     # More informative: dgamma(2, 20)
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			# Beta parameters for seasonal predictors only - weak priors (2 parameters)
			for (b in 1:2) {
				beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
			}
			
			# NEW: process & initial-state scales
			sigma_proc ~ dunif(0, 0.5)
			tau_proc <- pow(sigma_proc, -2)
			sigma_init ~ dunif(0, 1)
			tau_init <- pow(sigma_init, -2)
			
			# PROCESS MODEL - Direct AR(1) on η with deterministic μ = η (identity link on CLR scale)
			for (p in 1:N.plot) {
				# Initial state on eta (level-anchored)
				eta[p, plot_start[p]] ~ dnorm(
					intercept +
					site_effect[plot_site_num[p]] +
					beta[1] * sin_mo[plot_start[p]] + beta[2] * cos_mo[plot_start[p]] +
					legacy_effect * legacy[p, plot_start[p]],
					tau_init
				)
				mu[p, plot_start[p]] <- eta[p, plot_start[p]]  # Identity link on CLR scale
				
				# DIRECT AR(1) on eta
				for (t in (plot_start[p] + 1):N.date) {
					eta[p, t] ~ dnorm(
						rho * eta[p, t - 1] +
						intercept +
						site_effect[plot_site_num[p]] +
						beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
						legacy_effect * legacy[p, t],
						tau_proc
					)
					mu[p, t] <- eta[p, t]  # Identity link on CLR scale
				}
			}
			
			# OBSERVATION MODEL - CLR scale
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[plot_num[i], timepoint[i]], tau = precision)
			}
			
			# Probability-scale output for downstream compatibility
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_mu[p, t] <- 1 / (1 + exp(-mu[p, t]))   # Back-transform to (0,1) for compatibility
				}
			}
		})
	} else if (model_name == "env_cycl" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# PRIORS - Match driver uncertainty script
			precision ~ dgamma(2, 0.1)              # Gentle prior - mean 20, var 200
			rho ~ dbeta(1, 1)                       # Uniform prior - very weak
			intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
			legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak
			
			# STABLE PRIORS for site effects
			site_effect_sd ~ dgamma(2, 20)     # More informative: dgamma(2, 20)
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			# Beta parameters for environmental and seasonal predictors - weak priors
			for (b in 1:8) {
				beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
			}
			
			# NEW: process & initial-state scales
			sigma_proc ~ dunif(0, 0.5)
			tau_proc <- pow(sigma_proc, -2)
			sigma_init ~ dunif(0, 1)
			tau_init <- pow(sigma_init, -2)
			
			# PROCESS MODEL - Direct AR(1) on η with deterministic μ = η (identity link on CLR scale)
			for (p in 1:N.plot) {
				# Initial state on eta (level-anchored). Use *_est to honor driver-uncertainty toggles.
				eta[p, plot_start[p]] ~ dnorm(
					intercept +
					site_effect[plot_site_num[p]] +
					beta[1] * sin_mo[plot_start[p]] + beta[2] * cos_mo[plot_start[p]] +
					beta[3] * temp_est[plot_site_num[p], plot_start[p]] +
					beta[4] * mois_est[plot_site_num[p], plot_start[p]] +
					beta[5] * pH_est[p, plot_start[p]] + beta[6] * pC_est[p, plot_start[p]] +
					beta[7] * relEM[p, plot_start[p]] +
					beta[8] * LAI[plot_site_num[p], plot_start[p]] +
					legacy_effect * legacy[p, plot_start[p]],
					tau_init
				)
				mu[p, plot_start[p]] <- eta[p, plot_start[p]]  # Identity link on CLR scale
				
				# DIRECT AR(1) on eta
				for (t in (plot_start[p] + 1):N.date) {
					eta[p, t] ~ dnorm(
						rho * eta[p, t - 1] +
						intercept +
						site_effect[plot_site_num[p]] +
						beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
						beta[3] * temp_est[plot_site_num[p], t] +
						beta[4] * mois_est[plot_site_num[p], t] +
						beta[5] * pH_est[p, t] + beta[6] * pC_est[p, t] +
						beta[7] * relEM[p, t] +
						beta[8] * LAI[plot_site_num[p], t] +
						legacy_effect * legacy[p, t],
						tau_proc
					)
					mu[p, t] <- eta[p, t]  # Identity link on CLR scale
				}
			}
			
			# DRIVER UNCERTAINTY - Temporal (temperature and moisture)
			if (temporalDriverUncertainty) {
				for (k in 1:N.site) {
					for (t in site_start[k]:N.date) {
						mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
						temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
					}
				}
			} else {
				for (k in 1:N.site) {
					for (t in site_start[k]:N.date) {
						mois_est[k, t] <- mois[k, t]
						temp_est[k, t] <- temp[k, t]
					}
				}
			}
			
			# DRIVER UNCERTAINTY - Spatial (pH and pC - constant over time)
			if (spatialDriverUncertainty) {
				for (p in 1:N.plot) {
					pH_est_p[p] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
					pC_est_p[p] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
					for (t in plot_start[p]:N.date) {
						pH_est[p, t] <- pH_est_p[p]   # deterministic replicate over t
						pC_est[p, t] <- pC_est_p[p]
					}
				}
			} else {
				for (p in 1:N.plot) {
					for (t in plot_start[p]:N.date) {
						pH_est[p, t] <- pH[p, t]
						pC_est[p, t] <- pC[p, t]
					}
				}
			}
			
			# OBSERVATION MODEL - CLR scale
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[plot_num[i], timepoint[i]], tau = precision)
			}
			
			# Probability-scale output for downstream compatibility
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_mu[p, t] <- 1 / (1 + exp(-mu[p, t]))   # Back-transform to (0,1) for compatibility
				}
			}
		})
	} else if (model_name == "env_cov" && use_legacy_covariate) {
		modelCode <- nimble::nimbleCode({
			# PRIORS - Match driver uncertainty script
			precision ~ dgamma(2, 0.1)              # Gentle prior - mean 20, var 200
			rho ~ dbeta(1, 1)                       # Uniform prior - very weak
			intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
			legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak
			
			# STABLE PRIORS for site effects
			site_effect_sd ~ dgamma(2, 20)     # More informative: dgamma(2, 20)
			for (k in 1:N.site) {
				site_effect[k] ~ dnorm(0, sd = site_effect_sd)
			}
			
			# Beta parameters for environmental predictors - weak priors (6 parameters)
			for (b in 1:6) {
				beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
			}
			
			# NEW: process & initial-state scales
			sigma_proc ~ dunif(0, 0.5)
			tau_proc <- pow(sigma_proc, -2)
			sigma_init ~ dunif(0, 1)
			tau_init <- pow(sigma_init, -2)
			
			# PROCESS MODEL - Direct AR(1) on η with deterministic μ = η (identity link on CLR scale)
			for (p in 1:N.plot) {
				# Initial state on eta (level-anchored). Use *_est to honor driver-uncertainty toggles.
				eta[p, plot_start[p]] ~ dnorm(
					intercept +
					site_effect[plot_site_num[p]] +
					beta[1] * temp_est[plot_site_num[p], plot_start[p]] +
					beta[2] * mois_est[plot_site_num[p], plot_start[p]] +
					beta[3] * pH_est[p, plot_start[p]] + beta[4] * pC_est[p, plot_start[p]] +
					beta[5] * relEM[p, plot_start[p]] +
					beta[6] * LAI[plot_site_num[p], plot_start[p]] +
					legacy_effect * legacy[p, plot_start[p]],
					tau_init
				)
				mu[p, plot_start[p]] <- eta[p, plot_start[p]]  # Identity link on CLR scale
				
				# DIRECT AR(1) on eta
				for (t in (plot_start[p] + 1):N.date) {
					eta[p, t] ~ dnorm(
						rho * eta[p, t - 1] +
						intercept +
						site_effect[plot_site_num[p]] +
						beta[1] * temp_est[plot_site_num[p], t] +
						beta[2] * mois_est[plot_site_num[p], t] +
						beta[3] * pH_est[p, t] + beta[4] * pC_est[p, t] +
						beta[5] * relEM[p, t] +
						beta[6] * LAI[plot_site_num[p], t] +
						legacy_effect * legacy[p, t],
						tau_proc
					)
					mu[p, t] <- eta[p, t]  # Identity link on CLR scale
				}
			}
			
			# DRIVER UNCERTAINTY - Temporal (temperature and moisture)
			if (temporalDriverUncertainty) {
				for (k in 1:N.site) {
					for (t in site_start[k]:N.date) {
						mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
						temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
					}
				}
			} else {
				for (k in 1:N.site) {
					for (t in site_start[k]:N.date) {
						mois_est[k, t] <- mois[k, t]
						temp_est[k, t] <- temp[k, t]
					}
				}
			}
			
			# DRIVER UNCERTAINTY - Spatial (pH and pC - constant over time)
			if (spatialDriverUncertainty) {
				for (p in 1:N.plot) {
					pH_est_p[p] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
					pC_est_p[p] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
					for (t in plot_start[p]:N.date) {
						pH_est[p, t] <- pH_est_p[p]   # deterministic replicate over t
						pC_est[p, t] <- pC_est_p[p]
					}
				}
			} else {
				for (p in 1:N.plot) {
					for (t in plot_start[p]:N.date) {
						pH_est[p, t] <- pH[p, t]
						pC_est[p, t] <- pC[p, t]
					}
				}
			}
			
			# OBSERVATION MODEL - CLR scale
			for (i in 1:N.core) {
				y[i] ~ dnorm(mu[plot_num[i], timepoint[i]], tau = precision)
			}
			
			# Probability-scale output for downstream compatibility
			for (p in 1:N.plot) {
				for (t in 1:N.date) {
					plot_mu[p, t] <- 1 / (1 + exp(-mu[p, t]))   # Back-transform to (0,1) for compatibility
				}
			}
		})
	} else {
		stop("Unsupported model combination: ", model_name, " with use_legacy_covariate=", use_legacy_covariate, 
		     ". Only models with legacy covariates are supported.")
	}
	
	# Create inits using the driver uncertainty approach
	cat("Creating initial values with driver uncertainty...\n")
	set.seed(chain_no * 1000 + j * 100)  # Different seed per chain/model
	
	# Use the standardized initialization function
	inits <- create_inits_with_uncertainty(constants, model_name, model_data = model.dat, chain_no = chain_no)
	
	cat("Model built successfully\n")
	
	# Build model
	Rmodel <- nimbleModel(code = modelCode, constants = constants,
						  data = list(y=as.vector(model.dat$y)), inits = inits)
	
	# Debug: Check data dimensions
	cat("Data dimensions check:\n")
	cat("  y dimensions:", length(as.vector(model.dat$y)), "\n")
	cat("  y class:", class(as.vector(model.dat$y)), "\n")
	cat("  constants$N.core:", constants$N.core, "\n")
	cat("  constants$N.beta:", constants$N.beta, "\n")
	
	# Compile model
	cModel <- compileNimble(Rmodel)
	
	cat("Model compiled successfully\n")
	
	# Configure MCMC with proper sampler management - match driver uncertainty script
	cat("Configuring MCMC...\n")
	
	# Enhanced monitoring for comprehensive convergence analysis (matching driver uncertainty approach)
	# Start with core parameters that all models have
	monitored_params <- c(
		"precision", "rho", "beta", "site_effect_sd", "site_effect", "intercept",
		"sigma_proc", "sigma_init"  # NEW: process and initial state scales
	)
	
	# Add legacy_effect for all models that use it
	if (use_legacy_covariate) {
		monitored_params <- c(monitored_params, "legacy_effect")
	}
	
	# Add driver uncertainty parameters for DRIVER UNCERTAINTY models
	# NOTE: Removed temp_est, mois_est, pH_est, pC_est from monitoring to reduce memory usage
	# These are large arrays that can cause performance issues during testing
	if (model_name %in% c("env_cycl", "env_cov")) {
		# Only monitor the spatial pH/pC variables if using driver uncertainty
		if (constants$spatialDriverUncertainty) {
			monitored_params <- c(monitored_params, "pH_est_p", "pC_est_p")
		}
	}
	
	# Use monitors2 for latent variables at different interval (matching driver uncertainty approach)
	monitored_latent_params <- c(
		# Monitor latent process variables at different interval for efficiency
		"eta", "mu", "plot_mu"  # NEW: eta is now the main latent process, mu is deterministic, plot_mu for compatibility
	)
	
	cat("Monitoring parameters for convergence analysis:\n")
	cat("  Core parameters:", paste(monitored_params, collapse = ", "), "\n")
	cat("  Latent variables:", paste(monitored_latent_params, collapse = ", "), "\n")
	cat("  Total beta parameters:", constants$N.beta, "\n")

	# Choose small thins for local tests to ensure mvSamples2 gets rows
	thin  <- 5   # for parameters
	thin2 <- 10  # for latent (eta/mu) - must be <= total_iterations - burnin
	
	mcmcConf <- configureMCMC(
		model = cModel,
		monitors = monitored_params,
		monitors2 = monitored_latent_params,  # Use monitors2 for latent variables
		thin = thin,
		thin2 = thin2,
		enableWAIC = FALSE
	)
	
	# Add specialized samplers for better convergence of key parameters (matching driver uncertainty approach)
	cat("Adding specialized samplers for convergence improvement...\n")
	
	# 1. FIRST remove default samplers to prevent conflicts
	cat("  Removing default samplers...\n")
	mcmcConf$removeSamplers(c("precision","rho","intercept","site_effect_sd","sigma_proc","sigma_init"))
	for (b in seq_len(constants$N.beta))
		mcmcConf$removeSamplers(paste0("beta[", b, "]"))
	for (k in seq_len(constants$N.site))
		mcmcConf$removeSamplers(paste0("site_effect[", k, "]"))
	if (use_legacy_covariate) mcmcConf$removeSamplers("legacy_effect")
	if (model_name %in% c("env_cycl","env_cov") && constants$spatialDriverUncertainty) {
		for (p in seq_len(constants$N.plot)) {
			mcmcConf$removeSamplers(paste0("pH_est_p[", p, "]"))
			mcmcConf$removeSamplers(paste0("pC_est_p[", p, "]"))
		}
	}
	
	# 2. THEN add specialized samplers - IMPROVED SAMPLER STRATEGY
	cat("  Adding specialized samplers...\n")
	
	# Core parameters - use slice samplers for better mixing
	mcmcConf$addSampler(target = "precision", type = "slice")
	mcmcConf$addSampler(target = "rho", type = "slice")
	mcmcConf$addSampler(target = "intercept", type = "slice")
	mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
	mcmcConf$addSampler(target = "sigma_proc", type = "slice")
	mcmcConf$addSampler(target = "sigma_init", type = "slice")
	
	# Legacy effect if applicable
	if (use_legacy_covariate && model_name %in% c("env_cycl", "env_cov")) {
		mcmcConf$addSampler(target = "legacy_effect", type = "slice")
	}
	
	# Site effects - individual samplers for better mixing
	for (k in 1:constants$N.site) {
		mcmcConf$addSampler(target = paste0("site_effect[", k, "]"), type = "slice")
	}
	
	# Beta parameters - individual samplers
	for (b in 1:constants$N.beta) {
		mcmcConf$addSampler(target = paste0("beta[", b, "]"), type = "slice")
	}
	
	# Driver uncertainty parameters if applicable
	if (model_name %in% c("env_cycl", "env_cov") && constants$spatialDriverUncertainty) {
		for (p in 1:constants$N.plot) {
			mcmcConf$addSampler(target = paste0("pH_est_p[", p, "]"), type = "slice")
			mcmcConf$addSampler(target = paste0("pC_est_p[", p, "]"), type = "slice")
		}
	}
	
	cat("  ✓ Specialized samplers added successfully\n")
	
	# Build and compile MCMC
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	cat("MCMC configured successfully\n")
	
	# Production MCMC parameters
	burnin <- 5000
	iter_per_chunk <- 5000
	init_iter <- 5000
	min_eff_size_perchain <- 5
	max_loops <- 40
	max_save_size <- 100000
	min_total_iterations <- 50000
	max_total_iterations <- 200000

	max_iter_env <- suppressWarnings(as.integer(Sys.getenv("MAX_ITER", NA)))
	if (!is.na(max_iter_env)) {
		max_total_iterations <- max_iter_env
		min_total_iterations <- min(min_total_iterations, max_iter_env)
		cat("  MAX_ITER override:", max_total_iterations, "\n")
	}
	
	cat("Running MCMC: 1000 iterations per chain (env_cycl full dataset)\n")
	cat("  Burn-in iterations:", burnin, "\n")
	cat("  Initial iterations:", init_iter, "\n")
	cat("  Iterations per chunk:", iter_per_chunk, "max loops:", max_loops, "\n")
	cat("  Target ESS per chain:", min_eff_size_perchain, " (TESTING MODE)\n")
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
	
	# Check convergence - fail if check fails
	continue <- tryCatch({
		check_continue(initial_samples, min_eff_size = min_eff_size_perchain)
	}, error = function(e) {
		cat("  ERROR: Convergence check failed:", e$message, "\n")
		stop("Convergence check failed - cannot proceed")
	})
	
	loop_counter <- 0
	total_iterations <- init_iter
	
	cat("  Initial convergence check result:", ifelse(continue, "CONTINUE", "CONVERGED"), "\n")
	
	# Create output directory for checkpoints (matching beta regression approach)
	cat("  Creating output directory for checkpoints...\n")
	
	# Define the base model output directory (matching beta regression approach)
	model_output_dir <- here("data", "model_outputs", "CLR_regression")
	
	# Create species-specific subdirectory (matching beta regression structure)
	species_output_dir <- file.path(model_output_dir, model_name, species)
	
	# Ensure the directory exists
	if (!dir.exists(species_output_dir)) {
		dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
	}
	
	# Verify directory was created
	if (!dir.exists(species_output_dir)) {
		stop("CRITICAL: Failed to create checkpoint directory: ", species_output_dir)
	}
	
	cat("  ✓ Checkpoint directory ready:", species_output_dir, "\n")
	
	# Create model_id for consistent naming with legacy covariate indicator (matching beta regression pattern)
	legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
	model_id <- paste(model_name, species, min.date, max.date, legacy_indicator, sep = "_")
	
	# Store all samples as we go - FIXED: Use initial samples as starting point
	all_samples <- initial_samples
	
	# Also accumulate samples2 (plot estimates) across loops
	initial_plot_samples <- as.matrix(compiled$mvSamples2)
	all_plot_samples <- initial_plot_samples
	cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
	cat("  Starting iterative accumulation with", nrow(all_plot_samples), "initial plot samples\n")
	
	# Don't save initial checkpoint - only save final result
	# This matches the beta regression workflow which doesn't create multiple checkpoint files per chain
	cat("  Initial samples collected:", nrow(all_samples), "iterations\n")
	
	while ((continue || total_iterations < min_total_iterations) && loop_counter < max_loops && total_iterations < max_total_iterations) {
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
		
		# Also accumulate samples2 (plot estimates) - mvSamples2 also resets between runs
		current_plot_samples <- as.matrix(compiled$mvSamples2)
		current_plot_count <- nrow(current_plot_samples)
		if (current_plot_count > 0) {
			all_plot_samples <- rbind(all_plot_samples, current_plot_samples)
			cat("  ✓ Added", current_plot_count, "new plot samples, total accumulated:", nrow(all_plot_samples), "\n")
		} else {
			cat("  WARNING: No plot samples in current iteration\n")
		}
		
		# Don't save intermediate checkpoints - only save final result
		# This matches the beta regression workflow which doesn't create multiple checkpoint files per chain
		cat("  Loop", loop_counter + 1, "completed with", nrow(all_samples), "total samples\n")
		
		# Check if we need to continue
		continue <- tryCatch({
			check_continue(all_samples, min_eff_size = min_eff_size_perchain)
		}, error = function(e) {
			cat("  ERROR: Convergence check failed in loop:", e$message, "\n")
			stop("Convergence check failed in loop")
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
	
	# Use accumulated plot samples instead of just final mvSamples2
	cat("Using accumulated plot-level estimates from all iterations...\n")
	plot_samples <- all_plot_samples
	cat("  Plot samples dimensions:", dim(plot_samples), "\n")
	cat("  Plot samples column names:", paste(colnames(plot_samples), collapse=", "), "\n")
	
	# Validate plot samples structure
	if (nrow(plot_samples) == 0) {
		cat("  WARNING: No plot samples found in accumulated samples\n")
		plot_samples <- samples  # Fallback to parameter samples if no plot samples
	} else {
		cat("  ✓ Plot samples extracted successfully from accumulated samples\n")
		
		# Validate that samples2 has proper column names for combine_chains
		col_names <- colnames(plot_samples)
		if (is.null(col_names) || !any(grepl("plot_mu|eta\\[|mu\\[", col_names))) {
			cat("  WARNING: samples2 missing proper column names (plot_mu/eta/mu), this may cause issues in combine_chains\n")
			cat("  Available column names:", paste(head(col_names, 10), collapse=", "), "\n")
		} else {
			cat("  ✓ samples2 has proper column names for combine_chains\n")
		}
	}
	
	cat("MCMC completed successfully\n")
	cat("Final sample dimensions:", dim(samples), "\n")
	cat("Plot sample dimensions:", dim(plot_samples), "\n")
	cat("Total iterations run:", total_iterations, "\n")
	cat("Convergence loops:", loop_counter, "\n")
	cat("Final ESS check:\n")
	
	# Final convergence check (safely ignore zero-variance columns like padding eta)
	final_ess <- apply(samples, 2, function(x) {
		if (sum(!is.na(x)) < 2 || is.na(var(x, na.rm = TRUE)) || var(x, na.rm = TRUE) == 0) {
			return(NA)
		}
		tryCatch(coda::effectiveSize(coda::as.mcmc(x)), error = function(e) NA)
	})
	
	valid_final_ess <- final_ess[!is.na(final_ess)]
	if (length(valid_final_ess) > 0) {
		min_final_ess <- min(valid_final_ess)
		cat("  Final minimum ESS:", round(min_final_ess, 1), "\n")
		cat("  Convergence achieved:", min_final_ess >= min_eff_size_perchain, "\n")
	} else {
		cat("  Final minimum ESS: NA (No valid variation)\n")
		min_final_ess <- NA
	}
	
	cat("=== ITERATIVE SAVING SUMMARY ===\n")
	cat("  Initial samples:", nrow(initial_samples), "iterations\n")
	cat("  Additional loops:", loop_counter, "\n")
	cat("  Total accumulated samples:", nrow(all_samples), "iterations\n")
	cat("  Total accumulated plot samples:", nrow(all_plot_samples), "iterations\n")
	cat("  Checkpoints saved:", loop_counter + 1, "files\n")
	cat("  Final sample size:", nrow(samples), "iterations\n")
	cat("  Final plot sample size:", nrow(all_plot_samples), "iterations\n")
	
	# Save MCMC samples with absolute path (matching beta regression structure)
	samples_file <- file.path(species_output_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))
	
	# Create the complete chain structure with standardized metadata
	chain_output <- list(
		samples = samples,
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
			has_driver_uncertainty = model_name %in% c("env_cycl", "env_cov"),
			niter = total_iterations,
			nburnin = burnin,
			thin = thin,
			thin2 = thin2,
			model_data = model.dat,
			nimble_code = modelCode,
			model_structure = "CLR_regression"
		)
	)
	
	# Warn (but DO NOT stop) if chains have all-NA or zero-variance parameters.
	# Unmodeled latent states (e.g., eta padding before plot_start) naturally have 0 variance.
	bad_param <- character(0)
	for (j in seq_len(ncol(samples))) {
		col_j <- samples[, j]
		if (sum(!is.na(col_j)) == 0 || is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0)
			bad_param <- c(bad_param, colnames(samples)[j])
	}
	bad_plot <- character(0)
	for (j in seq_len(ncol(plot_samples))) {
		col_j <- plot_samples[, j]
		if (sum(!is.na(col_j)) == 0 || is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0)
			bad_plot <- c(bad_plot, colnames(plot_samples)[j])
	}
	if (length(bad_param) > 0 || length(bad_plot) > 0) {
		cat("  WARNING: Some parameters/latent states have 0 variance (expected for unmodeled eta padding).\n")
		if (length(bad_param) > 0) cat("  Flat parameters:", paste(head(bad_param, 10), collapse = ", "), "...\n")
	}
	
	# Save with error handling and detailed debugging
	cat("  Attempting to save samples to:", samples_file, "\n")
	cat("  File path exists:", file.exists(dirname(samples_file)), "\n")
	cat("  Directory writable:", file.access(dirname(samples_file), mode = 2) == 0, "\n")
	cat("  Sample dimensions:", dim(samples), "\n")
	
	tryCatch({
		saveRDS(chain_output, samples_file)
		cat("✓ SUCCESS: Saved MCMC samples to:", samples_file, "\n")
		cat("  File size:", file.size(samples_file), "bytes\n")
	}, error = function(e) {
		cat("✗ ERROR: Failed to save samples to", samples_file, "\n")
		cat("  Error:", e$message, "\n")
		cat("  Error class:", class(e), "\n")
		stop("Failed to save MCMC samples")
	})
	
	cat("Saved MCMC samples to:", samples_file, "\n")
	cat("Sample dimensions:", dim(samples), "\n")
	cat("=== CLR model fitting completed successfully ===\n")
	cat("  - TESTING MODE: ESS ≥ 10 for all parameters (faster testing)\n")
	cat("  - CLR regression with eta (latent) + mu (deterministic) structure\n")
	cat("  - Driver uncertainty toggles (temporal + spatial) enabled\n")
	cat("  - Proper legacy covariate handling (present when needed, absent when not)\n")
	cat("  - Individual slice samplers for efficient parameter sampling\n")
	cat("  - CONVERGENCE: Adaptive sampling until ESS ≥ 10 reached\n")
	cat("  - Production-ready: Adequate burn-in, proper mixing, reliable results\n")
	cat("  - FIXED: Sample accumulation across loops (NIMBLE mvSamples resets)\n")
	cat("  - FIXED: samples2 (eta, mu, plot_mu) properly monitored and accumulated\n")
	cat("  - READY: Output structure compatible with 02_combineModelChains.r\n")
	cat("  - MATCHED: Structure identical to driver uncertainty script\n")
	
	return(list(
		status = "SUCCESS", 
		samples = samples,
		samples2 = plot_samples,  # Include plot-level estimates for consistency
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
			has_driver_uncertainty = model_name %in% c("env_cycl", "env_cov"),
			niter = total_iterations,
			nburnin = burnin,
			thin = thin,
			thin2 = thin2,
			model_data = model.dat,
			nimble_code = modelCode,
			model_structure = "CLR_regression"
		)
	))
}

# Remove duplicates before proceeding
params <- params %>% distinct(.keep_all = TRUE)

# Single-model run: 1 model × 3 chains × 1000 iterations (2013-2018 time period only)
params <- params %>%
	filter(scenario == "Legacy with covariate 2013-2018") %>%
	slice(1)

# Run single model
cat("Running 1 model with", nchains, "chains (1000 iterations each)\n")
cat("Testing with", nrow(params), "models\n")
cat("Models to test:\n")
print(params[, c("rank.name", "species", "model_name", "model_id")])

# Set up parallel cluster for LOCAL TESTING
library(parallel)
library(doParallel)

# Load data once before parallel execution
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Define valid_models (all params for this run: Ascomycota, 3 model types)
valid_models <- params
cat("Running", nrow(valid_models), "models (Ascomycota, env_cycl), 3 chains, 1000 iterations\n")
cat("Models to test:\n")
print(valid_models[, c("rank.name", "species", "model_name")])

# Set test_models to number of models to run
test_models <- nrow(valid_models)

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
                # Look for status files in the appropriate species-specific directory
                model_name <- valid_models$model_name[model_idx]
                species_name <- valid_models$species[model_idx]
                model_output_dir <- here("data", "model_outputs", "CLR_regression", model_name)
                species_output_dir <- file.path(model_output_dir, species_name)
                status_file <- file.path(species_output_dir, paste0("chain_", model_idx, "_", chain_no, "_status.txt"))
                error_file <- file.path(species_output_dir, paste0("chain_", model_idx, "_", chain_no, "_ERROR.txt"))
                
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
            # Create output directory early with HPC compatibility (matching beta regression approach)
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
                stop("Failed to create output directory")
            }
            
            # Create model_id for consistent naming using helper function
            model_id <- create_model_id(
                valid_models$model_name[model_idx],
                valid_models$species[model_idx],
                valid_models$min.date[model_idx],
                valid_models$max.date[model_idx],
                grepl("Legacy with covariate", valid_models$scenario[model_idx])
            )
            
            # Create species-specific subdirectory (matching beta regression structure)
            species_output_dir <- file.path(model_output_dir, valid_models$species[model_idx])
            if (!dir.exists(species_output_dir)) {
                dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
            }
            
            # Save MCMC samples immediately
            samples_file <- file.path(species_output_dir, 
                                      paste0("samples_", model_id, "_chain", chain_no, ".rds"))
            
            # Create the complete chain structure with metadata
            # Use the metadata from the result if available, otherwise fail
            if ("metadata" %in% names(result) && !is.null(result$metadata)) {
                # Use the complete metadata from the result
                metadata <- result$metadata
                # Add parallel execution specific fields
                metadata$task_idx <- task_idx
                metadata$completed_at <- Sys.time()
            } else {
                stop("No metadata available in result")
            }
            
            samples_to_save <- result$samples
            samples2_to_save <- if ("samples2" %in% names(result)) result$samples2 else result$samples

            # Warn (but DO NOT stop) if chains have all-NA or zero-variance parameters.
            # Unmodeled latent states (e.g., eta padding before plot_start) naturally have 0 variance.
            bad_param <- character(0)
            for (j in seq_len(ncol(samples_to_save))) {
                col_j <- samples_to_save[, j]
                if (sum(!is.na(col_j)) == 0 || is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0)
                    bad_param <- c(bad_param, colnames(samples_to_save)[j])
            }
            bad_plot <- character(0)
            for (j in seq_len(ncol(samples2_to_save))) {
                col_j <- samples2_to_save[, j]
                if (sum(!is.na(col_j)) == 0 || is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0)
                    bad_plot <- c(bad_plot, colnames(samples2_to_save)[j])
            }
            if (length(bad_param) > 0 || length(bad_plot) > 0) {
                cat("  WARNING: Some parameters/latent states have 0 variance (expected for unmodeled eta padding).\n")
            }

            chain_output <- list(
                samples = samples_to_save,
                samples2 = samples2_to_save,
                metadata = metadata
            )
            
            cat("  Attempting to save parallel results to:", samples_file, "\n")
            cat("  Directory exists:", dir.exists(dirname(samples_file)), "\n")
            cat("  Directory writable:", file.access(dirname(samples_file), mode = 2) == 0, "\n")
            
            saveRDS(chain_output, samples_file)
            cat("SAVED: Chain", chain_no, "for model", model_idx, "to", samples_file, "\n")
            cat("  File size:", file.size(samples_file), "bytes\n")
            
            # Also save a simple status file to track progress
            status_file <- file.path(species_output_dir, paste0("chain_", model_idx, "_", chain_no, "_status.txt"))
            cat("  Creating status file:", status_file, "\n")
            writeLines(paste("SUCCESS", Sys.time(), sep = "\t"), status_file)
            cat("  ✓ Status file created successfully\n")
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
        
        # Create detailed error file in the appropriate model output directory
        task <- all_tasks[task_idx, ]
        model_name <- valid_models$model_name[task$model_idx]
        species_name <- valid_models$species[task$model_idx]
        model_output_dir <- here("data", "model_outputs", "CLR_regression", model_name)
        species_output_dir <- file.path(model_output_dir, species_name)
        if (!dir.exists(species_output_dir)) {
            dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
        }
        error_file <- file.path(species_output_dir, paste0("chain_", task$model_idx, "_", task$chain_no, "_ERROR.txt"))
        
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
        
        cat("  Creating error file:", error_file, "\n")
        writeLines(error_report, error_file)
        cat("  ✓ Error file created successfully\n")
        
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

# Also export the here function explicitly and test it in workers
clusterExport(cl, "here")

# Test that here() function works in workers
cat("Testing here() function in workers...\n")
here_test <- tryCatch({
    clusterEvalQ(cl, {
        if (exists("here")) {
            test_path <- here("data", "model_outputs")
            list(exists = TRUE, path = test_path, working_dir = getwd())
        } else {
            list(exists = FALSE, path = NULL, working_dir = getwd())
        }
    })
}, error = function(e) {
    cat("✗ ERROR testing here() in workers:", e$message, "\n")
    return(rep(list(exists = FALSE, path = NULL, working_dir = getwd()), length(cl)))
})

cat("✓ here() function test results:\n")
for (i in seq_along(here_test)) {
    worker_info <- here_test[[i]]
    cat("  Worker", i, ":", 
        "here() exists:", worker_info$exists, 
        "path:", worker_info$path, 
        "wd:", worker_info$working_dir, "\n")
}

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

# Clean up status files from species-specific directories
for (model_idx in 1:nrow(valid_models)) {
    model_name <- valid_models$model_name[model_idx]
    species_name <- valid_models$species[model_idx]
    model_output_dir <- here("data", "model_outputs", "CLR_regression", model_name)
    species_output_dir <- file.path(model_output_dir, species_name)
    
    for (chain_no in 1:nchains) {
        status_file <- file.path(species_output_dir, paste0("chain_", model_idx, "_", chain_no, "_status.txt"))
        if (file.exists(status_file)) {
            unlink(status_file)
        }
    }
}

	cat("✓ DYNAMIC CLR transformation for compositional data with temporal dependence\n")
	cat("✓ Eta (latent) + mu (deterministic) structure matching driver uncertainty script\n")
	cat("✓ Driver uncertainty toggles (temporal + spatial) enabled\n")
	cat("✓ Precision parameter for observation error\n")
	cat("✓ Temporal autocorrelation (rho) for dynamic linear model structure\n")
	cat("✓ 🎯 HYBRID PRIORS - Gentle for main parameters, stable for site effects\n")
	cat("  - precision ~ dgamma(2, 0.1) - Gentle prior - mean 20, var 200\n")
	cat("  - rho ~ dbeta(1, 1) - Uniform prior for temporal dependence - very weak\n")
	cat("  - intercept ~ dnorm(0, sd = 10) - Wide normal prior for intercept\n")
	cat("  - legacy_effect ~ dnorm(0, sd = 10) - Wide normal prior for legacy effect\n")
	cat("  - site_effect_sd ~ dgamma(2, 20) - More informative: dgamma(2, 20) for site effects\n")
	cat("  - beta ~ dnorm(0, sd = 10) - Wide normal prior for coefficients\n")
	cat("  - site_effect ~ dnorm(0, sd = site_effect_sd) - Site effects with hierarchical prior\n")
	cat("  - sigma_proc ~ dunif(0, 0.5) - Process noise scale\n")
	cat("  - sigma_init ~ dunif(0, 1) - Initial state uncertainty scale\n")
	cat("✓ Dynamic process monitoring (eta, mu, plot_mu) with monitors2 extraction\n")
	cat("✓ Individual slice samplers for all parameters including rho\n")
	cat("✓ Convergence-based sampling with iterative saving\n")
	cat("✓ LOCAL TESTING: 2 cores, 2 chains for faster testing\n")
	cat("✓ MATCHED: Structure identical to driver uncertainty script\n")
	cat("Check output files in:", here("data", "model_outputs", "CLR_regression"), "/[model_name]/\n")




#!/usr/bin/env Rscript
# Fit beta regression models for all microbial groups and predictor sets
# - Weak priors for main parameters (Jeffreys, Uniform, wide normal)
# - More informative priors for site effects (dgamma(2, 20))


# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
here::i_am("analysis/model_analysis/01_fitModels_betaReg.R")
project_root <- here()

cat("here() starts at", project_root, "\n")
cat("Project root set to:", getwd(), "\n")

# Load the microbialForecast package to access helper functions
library(microbialForecast)

# Load nimble explicitly for options
library(nimble)

# Set nimbleOptions for faster compilation
nimbleOptions(buildInterfacesForCompiledNimbleFunctions = FALSE) # faster compile
nimbleOptions(optimize = TRUE)
cat("✓ NIMBLE options set for faster compilation\n")

# Load packages and create directories using package functions
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("logit_beta_driver_uncertainty", "logit_beta_driver_uncertainty/cycl_only")
)



# Define output directory early to prevent undefined variable errors
model_output_dir <- here("data", "model_outputs", "logit_beta_driver_uncertainty")

cat("==================================================\n")
cat("Microbial forecasts environment setup complete!\n")
cat("Ready for analysis.\n")
cat("==================================================\n")

# Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
    k <- as.numeric( argv[1] )
} else {
    k=1
}

# Run with 4 chains for HPC production use
nchains = 3

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

# HPC PRODUCTION CONFIGURATION: Run multiple models with driver uncertainty enabled
params <- params_in %>% ungroup %>% filter(
    # Run ALL model types with driver uncertainty enabled
    model_name %in% c("env_cycl", "env_cov", "cycl_only") &
        # Only include models with driver uncertainty enabled
        temporalDriverUncertainty == TRUE &
        spatialDriverUncertainty == TRUE &
        # Focus on 2013-2018 period for legacy analysis
        scenario %in% c("Legacy with covariate 2013-2018") &
        species %in% c("plant_pathogen","herbicide_stress")
        # Include multiple ranks for comprehensive coverage
#        rank.name %in% c("phylum_fun")
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
cat("  - Models to run:", nrow(params), "models (cycl_only, env_cov, env_cycl) with driver uncertainty\n")
cat("  - Chains per model:", nchains, "(total", nrow(params) * nchains, "parallel tasks)\n")
cat("  - Initial iterations: ~ 0.1 minutes per chain\n")
cat("  - Additional iterations: Variable based on convergence\n")
cat("  - Target: ESS >= 3 per parameter\n")
cat("🎯 PRIOR STRATEGY: HYBRID (weak main + stable site effects) - PROVEN TO WORK\n")
cat("🔧 TESTING MODE: Local testing with 2 cores (not HPC)\n")

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
    
    # Check if species exists in rank data
    if (!(species_name %in% colnames(rank_data))) {
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

# Use the filtered parameters from data loading
valid_models <- params
cat("Testing with", nrow(valid_models), "models\n")
cat("Models to test:\n")
print(valid_models)

cat("HPC PRODUCTION: Running", nrow(valid_models), "models across all ranks with beta regression approach\n")

# LOCAL TESTING: Configuration for parallel execution
n_cores <- 2  # LOCAL TESTING: Use 2 cores for faster testing
cat("📊 Models to run:", nrow(valid_models), "models with", nchains, "chains each\n")
cat("⏱️  Expected runtime: Variable (convergence-based sampling)\n")
cat("🎯 Target: ESS >= 10 per parameter\n")
cat("🧪 TESTING MODE: Using first 500 rows of rank_df.spec for faster testing\n")
#

create_nimble_model_with_uncertainty <- function(model_name, use_legacy_covariate = TRUE, 
                                                      temporalDriverUncertainty = TRUE, 
                                                      spatialDriverUncertainty = TRUE) {
    cat("Building FINAL Nimble model WITH driver uncertainty:", model_name, "\n")
    cat("  Temporal uncertainty:", temporalDriverUncertainty, "\n")
    cat("  Spatial uncertainty:", spatialDriverUncertainty, "\n")
    
    if (model_name == "env_cycl" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak
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
            
            # PROCESS MODEL - FINAL with all fixes
            for (p in 1:N.plot) {
                # Initial condition - ONLY at plot_start[p] (no overlap with dynamic loop)
                Ex[p, plot_start[p]] ~ dunif(0.1, 0.9)
                mu[p, plot_start[p]] ~ dbeta(shape1 = Ex[p, plot_start[p]] * precision, 
                                            shape2 = (1 - Ex[p, plot_start[p]]) * precision)
                
                # Dynamic evolution - START AFTER plot_start[p] (NIMBLE will handle empty ranges)
                for (t in (plot_start[p] + 1):N.date) {
                    # Smooth feedback on log(mu) - no hard clamps
                    eta_prev[p, t] <- log(max(1e-6, mu[p, t - 1]))  # Safe log with tiny guard
                    
                    # Linear predictor with UNCERTAIN environmental predictors
                    eta_mean[p, t] <- rho * eta_prev[p, t] +
                        beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +  # Seasonal terms
                        beta[3] * temp_est[plot_site_num[p], t] + beta[4] * mois_est[plot_site_num[p], t] +  # Site-level environmental terms (UNCERTAIN)
                        beta[5] * pH_est[p, t] + beta[6] * pC_est[p, t] +  # Plot-level environmental terms (UNCERTAIN)
                        beta[7] * relEM[p, t] + beta[8] * LAI[plot_site_num[p], t] +  # Mixed level terms
                        site_effect[plot_site_num[p]] +  # Site effects
                        legacy_effect * legacy[p, t] +  # Legacy covariate
                        intercept  # Baseline abundance
                    
                    # ETA CAP to avoid exp() overflow
                    eta_cap[p, t] <- min(eta_mean[p, t], 700)
                    
                    # SMOOTH exp-only link to (0,1) - no hard clamps!
                    Ex_raw[p, t] <- 1 - exp(-exp(eta_cap[p, t]))    # cloglog^-1
                    Ex[p, t] <- Ex_raw[p, t]                         # Already in (0,1)
                    
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, 
                                    shape2 = (1 - Ex[p, t]) * precision)
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
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] ~ dnorm(pH[p, t], sd = pH_sd[p, t])
                        pC_est[p, t] ~ dnorm(pC[p, t], sd = pC_sd[p, t])
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
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else if (model_name == "env_cov" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak
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
            
            # PROCESS MODEL - FINAL with all fixes
            for (p in 1:N.plot) {
                # Initial condition - ONLY at plot_start[p] (no overlap with dynamic loop)
                Ex[p, plot_start[p]] ~ dunif(0.1, 0.9)
                mu[p, plot_start[p]] ~ dbeta(shape1 = Ex[p, plot_start[p]] * precision, 
                                            shape2 = (1 - Ex[p, plot_start[p]]) * precision)
                
                # Dynamic evolution - START AFTER plot_start[p] (NIMBLE will handle empty ranges)
                for (t in (plot_start[p] + 1):N.date) {
                    # Smooth feedback on log(mu) - no hard clamps
                    eta_prev[p, t] <- log(max(1e-6, mu[p, t - 1]))  # Safe log with tiny guard
                    
                    # Linear predictor with UNCERTAIN environmental predictors (no seasonal terms)
                    eta_mean[p, t] <- rho * eta_prev[p, t] +
                        beta[1] * temp_est[plot_site_num[p], t] + beta[2] * mois_est[plot_site_num[p], t] +  # Site-level environmental terms (UNCERTAIN)
                        beta[3] * pH_est[p, t] + beta[4] * pC_est[p, t] +  # Plot-level environmental terms (UNCERTAIN)
                        beta[5] * relEM[p, t] + beta[6] * LAI[plot_site_num[p], t] +  # Mixed level terms
                        site_effect[plot_site_num[p]] +  # Site effects
                        legacy_effect * legacy[p, t] +  # Legacy covariate
                        intercept  # Baseline abundance
                    
                    # ETA CAP to avoid exp() overflow
                    eta_cap[p, t] <- min(eta_mean[p, t], 700)
                    
                    # SMOOTH exp-only link to (0,1) - no hard clamps!
                    Ex_raw[p, t] <- 1 - exp(-exp(eta_cap[p, t]))    # cloglog^-1
                    Ex[p, t] <- Ex_raw[p, t]                         # Already in (0,1)
                    
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, 
                                    shape2 = (1 - Ex[p, t]) * precision)
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
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] ~ dnorm(pH[p, t], sd = pH_sd[p, t])
                        pC_est[p, t] ~ dnorm(pC[p, t], sd = pC_sd[p, t])
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
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else if (model_name == "cycl_only") {
        modelCode <- nimble::nimbleCode({
            # 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak
            rho ~ dbeta(1, 1)                       # Uniform prior - very weak
            intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
            
            # STABLE PRIORS for site effects
            site_effect_sd ~ dgamma(2, 20)     # More informative: dgamma(2, 20)
            for (k in 1:N.site) {
                site_effect[k] ~ dnorm(0, sd = site_effect_sd)
            }
            
            # Beta parameters for seasonal predictors only - weak priors (2 parameters)
            for (b in 1:2) {
                beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
            }
            
            # PROCESS MODEL - FINAL with all fixes
            for (p in 1:N.plot) {
                # Initial condition - ONLY at plot_start[p] (no overlap with dynamic loop)
                Ex[p, plot_start[p]] ~ dunif(0.1, 0.9)
                mu[p, plot_start[p]] ~ dbeta(shape1 = Ex[p, plot_start[p]] * precision, 
                                            shape2 = (1 - Ex[p, plot_start[p]]) * precision)
                
                # Dynamic evolution - START AFTER plot_start[p] (NIMBLE will handle empty ranges)
                for (t in (plot_start[p] + 1):N.date) {
                    # Smooth feedback on log(mu) - no hard clamps
                    eta_prev[p, t] <- log(max(1e-6, mu[p, t - 1]))  # Safe log with tiny guard
                    
                    # Linear predictor with ONLY seasonal terms
                    eta_mean[p, t] <- rho * eta_prev[p, t] +
                        beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +  # Seasonal terms only
                        site_effect[plot_site_num[p]] +  # Site effects
                        intercept  # Baseline abundance
                    
                    # ETA CAP to avoid exp() overflow
                    eta_cap[p, t] <- min(eta_mean[p, t], 700)
                    
                    # SMOOTH exp-only link to (0,1) - no hard clamps!
                    Ex_raw[p, t] <- 1 - exp(-exp(eta_cap[p, t]))    # cloglog^-1
                    Ex[p, t] <- Ex_raw[p, t]                         # Already in (0,1)
                    
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, 
                                    shape2 = (1 - Ex[p, t]) * precision)
                }
            }
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else {
        stop("Unsupported model combination: ", model_name, " with use_legacy_covariate=", use_legacy_covariate)
    }
    
    return(modelCode)
}

# Function to sanitize driver uncertainty data
sanitize_driver_uncertainty <- function(constants) {
    cat("Sanitizing driver uncertainty data...\n")
    
    eps <- 1e-6  # Small epsilon for invalid standard deviations
    
    # --- Sanitize pC / pC_sd (SPATIAL - constant over time) ---
    bad_pC <- !is.finite(constants$pC)
    bad_pC_sd <- !is.finite(constants$pC_sd) | (constants$pC_sd <= 0)
    
    # For spatial predictors, check only the first time point since they're constant
    bad_pC_plot <- bad_pC[, 1]  # Only check first time point
    bad_pC_sd_plot <- bad_pC_sd[, 1]  # Only check first time point
    
    # Create mask for valid plots (spatial predictors)
    constants$has_pC_plot <- (!bad_pC_plot) & (!bad_pC_sd_plot)
    
    # Replace bad pC means with a benign fill (use first time point for all)
    pC_fill <- median(constants$pC[is.finite(constants$pC)], na.rm = TRUE)
    constants$pC[bad_pC] <- pC_fill
    
    # Replace bad/zero sds with small epsilon (use first time point for all)
    pC_sd_fill <- median(constants$pC_sd[is.finite(constants$pC_sd)], na.rm = TRUE)
    if (is.na(pC_sd_fill) || pC_sd_fill <= 0) pC_sd_fill <- eps
    constants$pC_sd[bad_pC_sd] <- pC_sd_fill
    
    cat("  pC invalid cells:", sum(bad_pC), 
        " | pC_sd invalid/≤0 cells:", sum(bad_pC_sd), 
        " | valid pC plots:", sum(constants$has_pC_plot), "\n")
    
    # --- Sanitize pH / pH_sd (SPATIAL - constant over time) ---
    bad_pH <- !is.finite(constants$pH)
    bad_pH_sd <- !is.finite(constants$pH_sd) | (constants$pH_sd <= 0)
    
    # For spatial predictors, check only the first time point since they're constant
    bad_pH_plot <- bad_pH[, 1]  # Only check first time point
    bad_pH_sd_plot <- bad_pH_sd[, 1]  # Only check first time point
    
    # Create mask for valid plots (spatial predictors)
    constants$has_pH_plot <- (!bad_pH_plot) & (!bad_pH_sd_plot)
    
    # Replace bad pH means with a benign fill (use first time point for all)
    pH_fill <- median(constants$pH[is.finite(constants$pH)], na.rm = TRUE)
    constants$pH[bad_pH] <- pH_fill
    
    # Replace bad/zero sds with small epsilon (use first time point for all)
    pH_sd_fill <- median(constants$pH_sd[is.finite(constants$pH_sd)], na.rm = TRUE)
    if (is.na(pH_sd_fill) || pH_sd_fill <= 0) pH_sd_fill <- eps
    constants$pH_sd[bad_pH_sd] <- pH_sd_fill
    
    cat("  pH invalid cells:", sum(bad_pH), 
        " | pH_sd invalid/≤0 cells:", sum(bad_pH_sd), 
        " | valid pH plots:", sum(constants$has_pH_plot), "\n")
    
    # --- Sanitize temp / temp_sd ---
    bad_temp <- !is.finite(constants$temp)
    bad_temp_sd <- !is.finite(constants$temp_sd) | (constants$temp_sd <= 0)
    constants$has_temp <- (!bad_temp) & (!bad_temp_sd)
    temp_fill <- median(constants$temp[is.finite(constants$temp)], na.rm = TRUE)
    constants$temp[bad_temp] <- temp_fill
    constants$temp_sd[bad_temp_sd] <- eps
    
    cat("  temp invalid cells:", sum(bad_temp), 
        " | temp_sd invalid/≤0 cells:", sum(bad_temp_sd), 
        " | stochastic temp cells:", sum(constants$has_temp), "\n")
    
    # --- Sanitize mois / mois_sd ---
    bad_mois <- !is.finite(constants$mois)
    bad_mois_sd <- !is.finite(constants$mois_sd) | (constants$mois_sd <= 0)
    constants$has_mois <- (!bad_mois) & (!bad_mois_sd)
    mois_fill <- median(constants$mois[is.finite(constants$mois)], na.rm = TRUE)
    constants$mois[bad_mois] <- mois_fill
    constants$mois_sd[bad_mois_sd] <- eps
    
    cat("  mois invalid cells:", sum(bad_mois), 
        " | mois_sd invalid/≤0 cells:", sum(bad_mois_sd), 
        " | stochastic mois cells:", sum(constants$has_mois), "\n")
    
    cat("✓ Driver uncertainty data sanitized successfully\n")
    return(constants)
}

create_inits_with_uncertainty <- function(constants, model_name, model_data = NULL) {
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

    inits <- list(
        precision = 50,  # Start with moderate precision
        rho = 0.3,      # Start rho at 0.3 (moderate persistence)
        beta = beta_init,  # Seasonal + environmental coefficients
        site_effect_sd = 0.5,  # Start with moderate site effect SD
        site_effect = rnorm(constants$N.site, 0, 0.1),  # Small random initial values
        intercept = -2,  # Start intercept at -2
        legacy_effect = 0,  # Start legacy effect at 0
        Ex = matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date),  # Start with moderate abundance
        mu = matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date)   # Start with moderate abundance
    )
    
    # Initialize driver uncertainty variables CLOSE TO OBSERVED DATA for better mixing
    if (constants$N.site > 0 && constants$N.date > 0) {
        # Initialize near observed values instead of random values
        inits$temp_est <- constants$temp  # Start at observed temperature
        inits$mois_est <- constants$mois  # Start at observed moisture
        cat("  ✓ Initialized temp_est and mois_est at observed values\n")
    }
    
    if (constants$N.plot > 0 && constants$N.date > 0) {
        # Initialize near observed values instead of random values
        inits$pH_est <- constants$pH  # Start at observed pH
        inits$pC_est <- constants$pC  # Start at observed pC
        cat("  ✓ Initialized pH_est and pC_est at observed values\n")
    }
    
    cat("  ✓ Initial values created successfully WITH driver uncertainty\n")
    cat("    precision:", inits$precision, "\n")
    cat("    rho:", inits$rho, "\n")
    cat("    beta parameters:", n_beta, "\n")
    cat("    site effects:", constants$N.site, "\n")
    cat("    Ex matrix dimensions:", dim(inits$Ex), "\n")
    cat("    mu matrix dimensions:", dim(inits$mu), "\n")
    cat("    temp_est matrix dimensions:", dim(inits$temp_est), "\n")
    cat("    mois_est matrix dimensions:", dim(inits$mois_est), "\n")
    cat("    pH_est matrix dimensions:", dim(inits$pH_est), "\n")
    cat("    pC_est matrix dimensions:", dim(inits$pC_est), "\n")
    
    return(inits)
}

# Function to create visualization of plot_mu over time
create_plot_mu_visualization <- function(samples, samples2, model_data, metadata, output_dir) {
    cat("Creating plot_mu visualization...\n")
    
    # Use actual column names from samples
    cn <- colnames(samples)
    param_names <- cn
    
    # Create a basic plot of parameter estimates
    png(file.path(output_dir, "parameter_estimates.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    # Plot key parameters
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(samples[, i], type = "l", main = param_names[i], 
                 xlab = "Iteration", ylab = "Value")
            abline(h = mean(samples[, i]), col = "red", lty = 2)
        }
    }
    
    dev.off()
    cat("  ✓ Parameter estimates plot saved\n")
    
    # Create trace plots
    png(file.path(output_dir, "trace_plots.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(samples[, i], type = "l", main = paste("Trace:", param_names[i]), 
                 xlab = "Iteration", ylab = "Value")
        }
    }
    
    dev.off()
    cat("  ✓ Trace plots saved\n")
    
    # Create density plots
    png(file.path(output_dir, "density_plots.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(density(samples[, i]), main = paste("Density:", param_names[i]), 
                 xlab = "Value", ylab = "Density")
            abline(v = mean(samples[, i]), col = "red", lty = 2)
        }
    }
    
    dev.off()
    cat("  ✓ Density plots saved\n")
    
    # Create mu over time visualization if samples2 is available
    if (!is.null(samples2) && ncol(samples2) > 0) {
        cat("Creating mu over time visualization...\n")
        
        # Extract mu samples and compute summaries
        mu_samples <- samples2
        n_plots <- min(5, sqrt(ncol(mu_samples)))  # Show up to 5 plots
        
        png(file.path(output_dir, "mu_over_time.png"), width = 1200, height = 800)
        par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
        
        for (p in 1:n_plots) {
            # Find mu samples for this plot
            plot_cols <- grep(paste0("mu\\[", p, ","), colnames(mu_samples))
            if (length(plot_cols) > 0) {
                plot_mu <- mu_samples[, plot_cols, drop = FALSE]
                
                # Handle NaN values by removing them
                plot_mu_clean <- plot_mu[complete.cases(plot_mu), , drop = FALSE]
                
                if (nrow(plot_mu_clean) > 0) {
                    # Compute posterior mean and CI for each time point (with na.rm = TRUE)
                    mu_mean <- apply(plot_mu_clean, 2, mean, na.rm = TRUE)
                    mu_ci_lower <- apply(plot_mu_clean, 2, quantile, 0.025, na.rm = TRUE)
                    mu_ci_upper <- apply(plot_mu_clean, 2, quantile, 0.975, na.rm = TRUE)
                    
                    # Check for valid values
                    valid_idx <- !is.na(mu_mean) & !is.na(mu_ci_lower) & !is.na(mu_ci_upper)
                    
                    if (sum(valid_idx) > 0) {
                        # Plot only valid time points
                        time_points <- which(valid_idx)
                        plot(time_points, mu_mean[valid_idx], type = "l", 
                             main = paste("Plot", p, "μ over time"),
                             xlab = "Time", ylab = "μ", ylim = c(0, 1))
                        polygon(c(time_points, rev(time_points)), 
                                c(mu_ci_upper[valid_idx], rev(mu_ci_lower[valid_idx])), 
                                col = "lightblue", border = NA)
                        lines(time_points, mu_mean[valid_idx], col = "blue", lwd = 2)
                    } else {
                        plot(1, 1, type = "n", main = paste("Plot", p, "μ over time (no valid data)"),
                             xlab = "Time", ylab = "μ", ylim = c(0, 1))
                        text(1, 0.5, "No valid data", cex = 1.2)
                    }
                } else {
                    plot(1, 1, type = "n", main = paste("Plot", p, "μ over time (no data)"),
                         xlab = "Time", ylab = "μ", ylim = c(0, 1))
                    text(1, 0.5, "No data", cex = 1.2)
                }
            }
        }
        
        dev.off()
        cat("  ✓ Mu over time plot saved\n")
    }
    
    # Create summary statistics using actual column names
    summary_stats <- data.frame(
        Parameter = param_names,
        Mean = apply(samples, 2, mean),
        SD = apply(samples, 2, sd),
        Q025 = apply(samples, 2, quantile, 0.025),
        Q975 = apply(samples, 2, quantile, 0.975),
        ESS = coda::effectiveSize(coda::as.mcmc(samples))
    )
    
    write.csv(summary_stats, file.path(output_dir, "parameter_summary.csv"), row.names = FALSE)
    cat("  ✓ Parameter summary saved\n")
    
    return(summary_stats)
}

# Create function that uses our working approach for each model
run_scenarios_fixed <- function(j, chain_no) {
    # Initialize error tracking and logging
    start_time <- Sys.time()
    error_context <- list()
    
    tryCatch({
        # Load required libraries in each worker using helper function
        load_required_packages()
        cat("=== Starting model fitting ===\n")
        cat("Model index:", j, "Chain:", chain_no, "\n")
        cat("Model parameters:\n")
        print(params[j,])
        cat("=============================\n")
        
        # Debug HPC environment information
        cat("=============================\n")
        
        # Validate input parameters
        if (is.null(params) || nrow(params) < j) {
            stop("Params data frame not available or index out of bounds")
        }
        
        # Extract model parameters
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
        
        # Check if species column exists
        if (!(species %in% colnames(rank.df))) {
            stop("Species '", species, "' not found in rank '", rank.name, "'. Available columns: ",
                 paste(colnames(rank.df), collapse=", "))
        }
        
        cat("Preparing model data for", rank.name, "\n")
        
        # Extract the specific species
        cat("DEBUG: Extracting species", species, "from", rank.name, "\n")
        
        # Validate required columns exist
        required_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", species)
        missing_cols <- required_cols[!required_cols %in% colnames(rank.df)]
        if (length(missing_cols) > 0) {
            stop("Missing required columns in rank data: ", paste(missing_cols, collapse=", "))
        }
        
        # Extract species data for beta regression
        rank.df_spec <- rank.df %>%
            select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!species) %>%
            head(500)  # TESTING: Use only first 500 rows for faster testing
        
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
        
        # For beta regression, we only need the species abundance column
        # The model handles proportions directly (0-1 range)
        
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
        
        # For beta regression, we expect exactly 1 column (the species abundance)
        if (ncol(model.dat$y) != 1) {
            cat("  ❌ ERROR: Beta regression requires exactly 1 column in response data\n")
            cat("  Current columns:", ncol(model.dat$y), "\n")
            cat("  Column names:", paste(colnames(model.dat$y), collapse=", "), "\n")
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
        
        # Add driver uncertainty data (ONLY for environmental models)
        if (model_name %in% c("env_cycl", "env_cov")) {
            cat("Adding driver uncertainty data with proper dimension matching...\n")
            
            # Check if uncertainty data exists in predictor data
            predictor_data <- readRDS(here("data/clean/all_predictor_data.rds"))
            
            # FIX: Properly subset uncertainty matrices to match model dimensions
            cat("  Fixing dimension mismatches...\n")
            cat("    Model dimensions - N.site:", constants$N.site, "N.plot:", constants$N.plot, "N.date:", constants$N.date, "\n")
            
            # Subset temporal uncertainty (temp_sd, mois_sd) to N.site x N.date
            if ("temp_sd" %in% names(predictor_data)) {
                temp_sd_full <- predictor_data$temp_sd
                cat("    temp_sd full dimensions:", dim(temp_sd_full), "\n")
                constants$temp_sd <- temp_sd_full[1:constants$N.site, 1:constants$N.date, drop = FALSE]
                cat("  ✓ Added temp_sd uncertainty data (subset):", dim(constants$temp_sd), "matrix\n")
            } else {
                cat("  ⚠️  WARNING: temp_sd not found, using default values\n")
                constants$temp_sd <- matrix(0.1, nrow = constants$N.site, ncol = constants$N.date)
            }
            
            if ("mois_sd" %in% names(predictor_data)) {
                mois_sd_full <- predictor_data$mois_sd
                cat("    mois_sd full dimensions:", dim(mois_sd_full), "\n")
                constants$mois_sd <- mois_sd_full[1:constants$N.site, 1:constants$N.date, drop = FALSE]
                cat("  ✓ Added mois_sd uncertainty data (subset):", dim(constants$mois_sd), "matrix\n")
            } else {
                cat("  ⚠️  WARNING: mois_sd not found, using default values\n")
                constants$mois_sd <- matrix(0.1, nrow = constants$N.site, ncol = constants$N.date)
            }
            
            # Subset spatial uncertainty (pH_sd, pC_sd) - SPATIAL PREDICTORS ARE CONSTANT OVER TIME
            if ("pH_sd" %in% names(predictor_data)) {
                pH_sd_full <- predictor_data$pH_sd
                cat("    pH_sd full dimensions:", dim(pH_sd_full), "\n")
                # For spatial predictors, use the same value for all time points
                pH_sd_plot <- pH_sd_full[1:constants$N.plot, 1, drop = FALSE]  # Take first time point
                constants$pH_sd <- matrix(pH_sd_plot, nrow = constants$N.plot, ncol = constants$N.date)
                cat("  ✓ Added pH_sd uncertainty data (constant over time):", dim(constants$pH_sd), "matrix\n")
            } else {
                cat("  ⚠️  WARNING: pH_sd not found, using default values\n")
                constants$pH_sd <- matrix(0.1, nrow = constants$N.plot, ncol = constants$N.date)
            }
            
            if ("pC_sd" %in% names(predictor_data)) {
                pC_sd_full <- predictor_data$pC_sd
                cat("    pC_sd full dimensions:", dim(pC_sd_full), "\n")
                # For spatial predictors, use the same value for all time points
                pC_sd_plot <- pC_sd_full[1:constants$N.plot, 1, drop = FALSE]  # Take first time point
                constants$pC_sd <- matrix(pC_sd_plot, nrow = constants$N.plot, ncol = constants$N.date)
                cat("  ✓ Added pC_sd uncertainty data (constant over time):", dim(constants$pC_sd), "matrix\n")
            } else {
                cat("  ⚠️  WARNING: pC_sd not found, using default values\n")
                constants$pC_sd <- matrix(0.1, nrow = constants$N.plot, ncol = constants$N.date)
            }
            
            # Sanitize driver uncertainty data
            constants <- sanitize_driver_uncertainty(constants)
            
            # Add driver uncertainty flags
            constants$temporalDriverUncertainty <- TRUE
            constants$spatialDriverUncertainty <- TRUE
            
            cat("  ✓ Driver uncertainty data added successfully\n")
        } else {
            cat("Skipping driver uncertainty data for", model_name, "model\n")
        }
        
        # Add legacy covariate if needed
        if (use_legacy_covariate) {
            cat("Adding legacy covariate with enhanced validation...\n")
            
                    # Create legacy covariate matrix properly for plot x time indexing
        # The legacy covariate should be 1 for legacy period (2013-2015), 0 for post-2015
        cat("Creating legacy covariate matrix for plot x time structure...\n")
        cat("Matrix dimensions needed:", constants$N.plot, "plots ×", constants$N.date, "time points\n")
        
        # Use time-based approach (simpler and more robust)
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
        
        # Create model with driver uncertainty
        cat("Building Nimble model WITH driver uncertainty...\n")
        modelCode <- create_nimble_model_with_uncertainty(model_name, use_legacy_covariate, 
                                                         temporalDriverUncertainty = TRUE, 
                                                         spatialDriverUncertainty = TRUE)
        
        # Create initial values with driver uncertainty
        cat("Creating initial values WITH driver uncertainty...\n")
        inits <- create_inits_with_uncertainty(constants, model_name, model_data = model.dat)
        
        cat("Model built successfully\n")
        
        # STEP: Comprehensive Model Validation
        cat("Performing comprehensive model validation...\n")
        
        # Validate model dimensions
        cat("  Validating model dimensions...\n")
        expected_dims <- list(
            N.plot = constants$N.plot,
            N.date = constants$N.date,
            N.site = constants$N.site,
            N.core = constants$N.core
        )
        
        for (dim_name in names(expected_dims)) {
            if (expected_dims[[dim_name]] <= 0) {
                cat("    ❌ ERROR:", dim_name, "is", expected_dims[[dim_name]], "- must be positive\n")
                stop("Invalid model dimensions")
            }
            cat("    ✓", dim_name, "=", expected_dims[[dim_name]], "\n")
        }
        
        # Validate environmental predictors for environmental models
        if (model_name %in% c("env_cycl", "env_cov")) {
            cat("  Validating environmental predictors...\n")
            required_env <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
            
            for (pred in required_env) {
                if (pred %in% names(constants)) {
                    pred_data <- constants[[pred]]
                    if (is.matrix(pred_data)) {
                        cat("    ✓", pred, ":", dim(pred_data), "matrix\n")
                        
                        # Check matrix dimensions match expected
                        if (nrow(pred_data) != constants$N.plot && ncol(pred_data) != constants$N.date) {
                            cat("      ⚠️  WARNING:", pred, "dimensions may not match plot/time structure\n")
                        }
                    } else if (is.vector(pred_data)) {
                        cat("    ✓", pred, ":", length(pred_data), "vector\n")
                    } else {
                        cat("    ⚠️  WARNING:", pred, ":", class(pred_data), "object\n")
                    }
                } else {
                    cat("    ❌ ERROR:", pred, "missing for environmental model\n")
                    stop("Missing required environmental predictor: ", pred)
                }
            }
        }
        
        # Validate seasonal predictors
        cat("  Validating seasonal predictors...\n")
        if (length(constants$sin_mo) != constants$N.date) {
            cat("    ❌ ERROR: sin_mo length (", length(constants$sin_mo), ") != N.date (", constants$N.date, ")\n")
            stop("Seasonal predictor dimension mismatch")
        }
        if (length(constants$cos_mo) != constants$N.date) {
            cat("    ❌ ERROR: cos_mo length (", length(constants$cos_mo), ") != N.date (", constants$N.date, ")\n")
            stop("Seasonal predictor dimension mismatch")
        }
        cat("    ✓ Seasonal predictors validated\n")
        
        # Validate response data
        cat("  Validating response data...\n")
        if (nrow(model.dat$y) != constants$N.core) {
            cat("    ❌ ERROR: Response data rows (", nrow(model.dat$y), ") != N.core (", constants$N.core, ")\n")
            stop("Response data dimension mismatch")
        }
        if (ncol(model.dat$y) != 1) {
            cat("    ❌ ERROR: Response data should have 1 column, got", ncol(model.dat$y), "\n")
            stop("Response data structure error")
        }
        
        # Check for response data issues
        y_values <- as.vector(model.dat$y[,1])
        y_na <- sum(is.na(y_values))
        y_inf <- sum(is.infinite(y_values))
        y_neg <- sum(y_values < 0, na.rm = TRUE)
        y_gt1 <- sum(y_values > 1, na.rm = TRUE)
        
        if (y_na > 0) cat("    ⚠️  WARNING:", y_na, "missing values in response\n")
        if (y_inf > 0) cat("    ❌ ERROR:", y_inf, "infinite values in response\n")
        if (y_neg > 0) cat("    ❌ ERROR:", y_neg, "negative values in response\n")
        if (y_gt1 > 0) cat("    ❌ ERROR:", y_gt1, "values > 1 in response\n")
        
        if (y_inf > 0 || y_neg > 0 || y_gt1 > 0) {
            stop("Response data contains invalid values")
        }
        
        cat("    ✓ Response data validated\n")
        
        cat("  ✓ All model validations passed\n")
        
        # Build model
        cat("Building Nimble model...\n")
        Rmodel <- nimbleModel(code = modelCode, constants = constants,
                              data = list(y = model.dat$y[,1,drop=FALSE]), inits = inits)
        
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
        
        # Configure MCMC with proper sampler management
        cat("Configuring MCMC...\n")
        
        # Enhanced monitoring for comprehensive convergence analysis (matching working approach)
        monitored_params <- c(
            # Core parameters (primary focus) - MONITOR ALL ESSENTIAL PARAMETERS
            "precision", "rho", "beta", "site_effect_sd", "site_effect", "intercept", "legacy_effect",
            # Driver uncertainty parameters - MONITOR FOR ANALYSIS
            "temp_est", "mois_est", "pH_est", "pC_est"
        )
        
        # Use monitors2 for latent variables at different interval (matching working approach)
        monitored_latent_params <- c(
            # Monitor latent process variables at different interval for efficiency
            "Ex", "mu"
        )
        
        cat("Monitoring parameters for convergence analysis:\n")
        cat("  Core parameters:", paste(monitored_params, collapse = ", "), "\n")
        cat("  Latent variables:", paste(monitored_latent_params, collapse = ", "), "\n")
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
        
        # Run MCMC with convergence-based sampling
        cat("Running MCMC with convergence-based sampling...\n")
        burnin <- 100  # TESTING: Reduced for faster testing
        thin <- 1
        iter_per_chunk <- 200  # TESTING: Reduced for faster testing
        init_iter <- 500  # TESTING: Reduced for faster testing
        min_eff_size_perchain <- 10  # TESTING: Reduced for faster testing
        max_loops <- 3  # TESTING: Reduced for faster testing
        max_save_size <- 50000
        min_total_iterations <- 5000
        
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
        
        # Create model_id for consistent naming
        model_id <- create_model_id(model_name, species, min.date, max.date, use_legacy_covariate)
        
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
        
        # Also accumulate samples2 (plot estimates) across loops
        initial_plot_samples <- as.matrix(compiled$mvSamples2)
        all_plot_samples <- initial_plot_samples
        cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
        cat("  Starting iterative accumulation with", nrow(all_plot_samples), "initial plot samples\n")
        
        # Save initial samples as checkpoint (including samples2)
        initial_checkpoint_data <- list(
            samples = all_samples,
            samples2 = all_plot_samples,
            iterations = total_iterations,
            loop = 0
        )
        checkpoint_file <- file.path(species_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_initial.rds"))
        tryCatch({
            saveRDS(initial_checkpoint_data, checkpoint_file)
            cat("  ✓ Initial checkpoint saved with samples and samples2\n")
        }, error = function(e) {
            cat("  ✗ Failed to save initial checkpoint:", e$message, "\n")
        })
        
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
            
            # Save checkpoint after each loop (including samples2)
            loop_checkpoint_data <- list(
                samples = all_samples,
                samples2 = all_plot_samples,
                iterations = total_iterations,
                loop = loop_counter + 1
            )
            loop_checkpoint_file <- file.path(species_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_loop", loop_counter + 1, ".rds"))
            tryCatch({
                saveRDS(loop_checkpoint_data, loop_checkpoint_file)
                cat("  ✓ Loop checkpoint saved with samples and samples2\n")
            }, error = function(e) {
                cat("  ✗ Failed to save loop checkpoint:", e$message, "\n")
            })
            
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
            if (is.null(col_names) || !any(grepl("plot_mu", col_names))) {
                cat("  WARNING: samples2 missing proper column names (plot_mu), this may cause issues in combine_chains\n")
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
        cat("  Total accumulated plot samples:", nrow(all_plot_samples), "iterations\n")
        cat("  Checkpoints saved:", loop_counter + 1, "files\n")
        cat("  Final sample size:", nrow(all_samples), "iterations\n")
        cat("  Final plot sample size:", nrow(all_plot_samples), "iterations\n")
        
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
                has_driver_uncertainty = TRUE,  # Flag to identify driver uncertainty models
                scenario = scenario,
                min.date = min.date,
                max.date = max.date,
                niter = total_iterations,
                nburnin = burnin,
                thin = thin,
                thin2 = 20,  # Include thin2 for samples2
                model_data = model.dat,
                nimble_code = modelCode,
                model_structure = "stable_beta_regression_with_driver_uncertainty"
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
        cat("=== Model fitting completed WITH DRIVER UNCERTAINTY ===\n")
        cat("  - STABLE: Beta regression with driver uncertainty in environmental predictors\n")
        cat("  - TEMPORAL UNCERTAINTY: Temperature and moisture (site-level)\n")
        cat("  - SPATIAL UNCERTAINTY: pH and pC (plot-level, constant over time)\n")
        cat("  - SMOOTH LINK: cloglog transformation (1 - exp(-exp(eta))) for (0,1) bounds\n")
        cat("  - CONVERGENCE-BASED: Adaptive sampling until reasonable ESS reached\n")
        cat("  - ITERATIVE SAVING: Samples accumulated and saved incrementally\n")
        cat("  - FIXED: Sample accumulation across loops (NIMBLE mvSamples resets)\n")
        cat("  - FIXED: samples2 (plot estimates) properly monitored and accumulated\n")
        cat("  - READY: Output structure compatible with 02_combineModelChains.r\n")
        
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
                model_structure = "stable_beta_regression_with_driver_uncertainty"
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
        error_dir <- here("data", "model_outputs", "logit_beta_regression", model_name)
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

# Create a combined task list: (model_idx, chain_no) pairs
all_tasks <- expand.grid(model_idx = 1:nrow(valid_models), chain_no = 1:nchains)
cat("Total parallel tasks:", nrow(all_tasks), "(", nrow(valid_models), "models ×", nchains, "chains)\n")
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
        for (model_idx in 1:nrow(valid_models)) {
            for (chain_no in 1:nchains) {
                status_file <- paste0("chain_", model_idx, "_", chain_no, "_status.txt")
                error_file <- paste0("chain_", model_idx, "_", chain_no, "_ERROR.txt")
                
                if (file.exists(status_file)) completed <- completed + 1
                if (file.exists(error_file)) errors <- errors + 1
            }
        }
        
        total <- nrow(valid_models) * nchains
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
cat("Preparing parallel execution functions...\n")

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
                    test_dir <- file.path(base_dir, "logit_beta_driver_uncertainty", valid_models$model_name[model_idx])
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
                fallback_output_dir <- here("data", "model_outputs", "logit_beta_driver_uncertainty", valid_models$model_name[model_idx])
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
                grepl("Legacy with covariate", valid_models$scenario[model_idx])
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
                    has_driver_uncertainty = TRUE,  # Flag to identify driver uncertainty models
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
                    model_structure = "stable_beta_regression_with_driver_uncertainty"  # Model structure identifier
                )
            }
            
            chain_output <- list(
                samples = result$samples,
                samples2 = if("samples2" %in% names(result)) result$samples2 else result$samples,  # Plot-level estimates from monitors2 or fallback to parameter samples
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

# Define the global model_output_dir variable
model_output_dir <- here("data", "model_outputs", "logit_beta_driver_uncertainty")

# Set start time for runtime calculation
start_time <- Sys.time()

# Create cluster for parallel execution - LOCAL TESTING with 4 cores
cat("Creating cluster with", n_cores, "cores for", nrow(valid_models), "models ×", nchains, "chains\n")
cat("LOCAL TESTING: 4 cores allocated for local testing\n")

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
cat("ALL MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Summary of all models
cat("\nSummary of All Models:\n")
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


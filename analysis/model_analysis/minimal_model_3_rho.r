#!/usr/bin/env Rscript

# Add Rho Parameter - Phase 3 Minimal Model
# Focus: Add temporal persistence (rho parameter) to working minimal model
# Based on: Successful rho implementation from working models
# Test: <500 iterations to verify rho parameter convergence

library(nimble)
library(tidyverse)
library(here)
library(coda)
library(microbialForecast)

cat("=== Add Rho Parameter - Phase 3 Minimal Model ===\n")

# Load data directly like the working model
cat("Loading data files...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully\n")

# Get the same data structure as the working model
params_in = read.csv(here("data/clean/model_input_df.csv"),
                     colClasses = c(rep("character", 4),
                     rep("logical", 2),
                     rep("character", 4)))

# Focus on site effects period - use 2013-2018 range for site indexing
params <- params_in %>% ungroup %>% filter(
  rank.name %in% c("cellulolytic") &          
  scenario %in% c("Legacy with covariate 2013-2018") & 
  model_name %in% c("cycl_only") &            
  min.date == "20130601" &                    
  max.date == "20180101"                      
) %>% slice(1:1)  

cat("Running rho parameter test for", nrow(params), "models\n")

# Get the group data
rank.name <- params$rank.name[[1]]
species <- params$species[[1]]
min.date <- params$min.date[[1]]
max.date <- params$max.date[[1]]

cat("Testing:", species, "from", min.date, "to", max.date, "\n")

# Get the specific group data
rank.df <- all_ranks[[rank.name]]
rank.df_spec <- rank.df %>%
  select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!species) %>%
  mutate(other = 1 - !!sym(species))

# Prepare model data
cat("Preparing model data...\n")
model.dat <- prepBetaRegData(rank.df = rank.df_spec,
                            min.prev = 3,
                            min.date = min.date,
                            max.date = max.date)

cat("Data prepared successfully\n")
cat("  N.core:", nrow(model.dat$y), "\n")
cat("  N.plot:", model.dat$N.plot, "\n")
cat("  N.date:", model.dat$N.date, "\n")
cat("  N.site:", model.dat$N.site, "\n")
cat("  Species mean abundance:", round(mean(model.dat$y[,1], na.rm = TRUE), 3), "\n")

# Verify site indexing
cat("Site indexing verification:\n")
cat("  plot_site_num range:", range(model.dat$plot_site_num), "\n")
cat("  plot_start range:", range(model.dat$plot_start), "\n")
cat("  Any NA in plot_site_num:", any(is.na(model.dat$plot_site_num)), "\n")

# PHASE 3 MINIMAL MODEL: Add rho parameter for temporal persistence
cat("Defining Phase 3 minimal model with rho parameter...\n")
phase3_rho_model <- nimbleCode({
  # PRIORS - Minimal set with rho parameter (removed redundant sigma)
  precision ~ dgamma(2, 0.1)  # Observation precision
  
  # NEW: Rho parameter for temporal persistence (based on working models)
  rho ~ dbeta(5, 5)  # Tighter beta prior centered at 0.5, prevents extreme values
  
  # Site effects - Hierarchical structure (VERIFIED WORKING)
  site_effect_sd ~ dgamma(2, 0.1)
  for(s in 1:N.site) {
    site_effect[s] ~ dnorm(0, sd = site_effect_sd)
  }
  
  # PROCESS MODEL - Site effects + rho temporal persistence (FIXED: Use beta stochasticity)
  for(p in 1:N.plot) {
    # Initial condition - single time point per plot
    for(t in plot_start[p]) {
      Ex[p, t] ~ dunif(0.1, 0.9)
      mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
    }
    
    # Dynamic evolution - site effects + rho temporal persistence (FIXED: Remove normal noise)
    for(t in (plot_start[p] + 1):N.date) {
      # ENHANCED: Add rho parameter for temporal persistence (like working models)
      logit_Ex_prev[p, t] <- logit(Ex[p, t-1])
      logit_Ex_mean[p, t] <- rho * logit_Ex_prev[p, t] + site_effect[plot_site_num[p]]
      Ex[p, t] <- ilogit(logit_Ex_mean[p, t])  # FIXED: Remove normal noise, use deterministic evolution
      mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
    }
  }
  
  # OBSERVATION MODEL - Direct indexing 
  for(i in 1:N.core) {
    y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                     shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
  }
})

# Prepare data for NIMBLE
cat("Preparing NIMBLE data...\n")
nimble_data <- list(
  y = model.dat$y[,1,drop=FALSE]  # Single column for simplicity
)

nimble_constants <- list(
  N.site = model.dat$N.site,
  N.plot = model.dat$N.plot,
  N.date = model.dat$N.date,
  N.core = nrow(model.dat$y),
  plot_site_num = model.dat$plot_site_num,
  plot_start = model.dat$plot_start,
  plot_num = model.dat$plot_num,
  timepoint = model.dat$timepoint
)

# Initial values - Simple and stable with rho parameter
cat("Setting initial values with rho parameter...\n")
initial_values <- list(
  precision = 50,  # IMPROVED: Start with moderate precision (not too high)
  rho = 0.3,      # IMPROVED: Start rho at 0.3 (moderate persistence, not 0.5)
  site_effect_sd = 0.5,  # IMPROVED: Start with moderate site effect SD
  site_effect = rnorm(model.dat$N.site, 0, 0.1),  # IMPROVED: Small random initial values
  Ex = matrix(0.3, nrow = model.dat$N.plot, ncol = model.dat$N.date),  # IMPROVED: Start at 0.3 not 0.5
  mu = matrix(0.3, nrow = model.dat$N.plot, ncol = model.dat$N.date)   # IMPROVED: Start at 0.3 not 0.5
)

# Build and compile model
cat("Building NIMBLE model with rho parameter...\n")
model <- nimbleModel(
  code = phase3_rho_model,
  constants = nimble_constants,
  data = nimble_data,
  inits = initial_values
)

cat("Compiling model...\n")
compiled_model <- compileNimble(model)

# IMPROVED MCMC configuration for rho parameter convergence
cat("Configuring improved MCMC for rho parameter convergence...\n")

# Enhanced monitoring for comprehensive convergence analysis
monitored_params <- c(
  # Core parameters (primary focus) - MONITOR ONLY ESSENTIAL PARAMETERS
  "precision", "rho", "site_effect_sd", "site_effect"
)

# IMPROVED: Use monitors2 for latent variables at different interval
monitored_latent_params <- c(
  # Monitor latent process variables at different interval for efficiency
  "Ex", "mu"
)

cat("Monitoring parameters for convergence analysis:\n")
cat("  Core parameters:", paste(c("precision", "rho", "site_effect_sd"), collapse = ", "), "\n")
cat("  Site effects:", "site_effect[1:", model.dat$N.site, "]\n")
cat("  IMPROVED: Latent variables (Ex, mu) monitored via monitors2 at different interval\n")
cat("  Total parameters monitored:", length(monitored_params), "\n")
cat("  Latent variables monitored via monitors2:", length(monitored_latent_params), "\n")

mcmc_config <- configureMCMC(
  model = compiled_model,
  monitors = monitored_params,
  monitors2 = monitored_latent_params,  # IMPROVED: Use monitors2 for latent variables
  thin = 1,
  thin2 = 10,  # IMPROVED: Sample latent variables every 10th iteration for efficiency
  enableWAIC = FALSE
)

# Add specialized samplers for better convergence of key parameters
cat("Adding specialized samplers for convergence improvement...\n")

# 1. FIRST remove default samplers to prevent conflicts
cat("  Removing default samplers...\n")
mcmc_config$removeSamplers(c("precision", "rho", "site_effect_sd"))

# 2. THEN add specialized samplers - IMPROVED SAMPLER STRATEGY
cat("  Adding improved samplers for key parameters...\n")

# IMPROVED: Use regular slice samplers for single parameters, AF_slice for blocks
mcmc_config$addSampler(target = "precision", type = "slice")        # Regular slice sampler for single parameter
mcmc_config$addSampler(target = "rho", type = "slice")              # Regular slice sampler for single parameter
mcmc_config$addSampler(target = "site_effect_sd", type = "slice")   # Regular slice sampler for single parameter

# 3. IMPROVED: Better site effects sampling strategy
if (model.dat$N.site > 1) {
  cat("  Adding improved block sampler for site effects...\n")
  # Use adaptive block sampler for better mixing of correlated parameters
  mcmc_config$addSampler(target = paste0("site_effect[1:", model.dat$N.site, "]"), type = "AF_slice")
  cat("  Added adaptive block sampler for site_effect[1:", model.dat$N.site, "]\n")
} else {
  # Single site effect - use adaptive slice sampler
  cat("  Adding adaptive slice sampler for single site effect...\n")
  mcmc_config$addSampler(target = "site_effect[1]", type = "AF_slice")
}

# 4. Verify sampler configuration
cat("  Final sampler configuration:\n")
cat("    precision sampler:", mcmc_config$getSamplers()$precision$type, "\n")
cat("    rho sampler:", mcmc_config$getSamplers()$rho$type, "\n")  # NEW: Verify rho sampler
cat("    site_effect_sd sampler:", mcmc_config$getSamplers()$site_effect_sd$type, "\n")
cat("    IMPROVED: Single parameters use slice samplers, site effects use adaptive block sampler (AF_slice)\n")

# Build MCMC
cat("Building MCMC...\n")
mcmc_built <- buildMCMC(mcmc_config)
compiled_mcmc <- compileNimble(mcmc_built, project = compiled_model)

cat("Running MCMC test with rho parameter...\n")
start_time <- Sys.time()

mcmc_samples <- runMCMC(
  compiled_mcmc,
  niter = 400, 
  nburnin = 100,  
  nchains = 3, 
  samplesAsCodaMCMC = TRUE
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("MCMC completed!\n")
cat("Runtime:", round(runtime, 2), "minutes\n")
cat("Number of chains:", length(mcmc_samples), "\n")
cat("Samples per chain:", nrow(mcmc_samples[[1]]), "\n")
cat("Total samples:", length(mcmc_samples) * nrow(mcmc_samples[[1]]), "\n")

# Convergence diagnostics for rho parameter
cat("\n Convergence diagnostics for rho parameter:\n")

# Check MCMC output structure first
cat("Checking MCMC output structure...\n")
cat("  MCMC output type:", class(mcmc_samples), "\n")
cat("  Number of chains:", length(mcmc_samples), "\n")
cat("  MCMC output structure:", str(mcmc_samples), "\n")

# Check if we have valid output
if(length(mcmc_samples) == 0) {
  cat("ERROR: No MCMC output generated!\n")
  stop("MCMC failed to produce output")
}

# Handle different output structures
if(inherits(mcmc_samples, "mcmc.list")) {
  cat("  Output is mcmc.list (correct format for CodaMCMC)\n")
  cat("  Actual number of chains in mcmc.list:", length(mcmc_samples), "\n")
  cat("  First chain type:", class(mcmc_samples[[1]]), "\n")
  cat("  First chain dimensions:", dim(mcmc_samples[[1]]), "\n")
  cat("  First chain columns:", ncol(mcmc_samples[[1]]), "\n")
  
  # Continue with mcmc.list diagnostics...
  
} else if(is.list(mcmc_samples) && "samples" %in% names(mcmc_samples)) {
  cat("  Output is monitors2 structure (samples + samples2)\n")
  cat("  Main samples (monitors):", length(mcmc_samples$samples), "chains\n")
  cat("  Latent samples (monitors2):", length(mcmc_samples$samples2), "chains\n")
  
  # Extract main parameters from samples
  main_samples <- mcmc_samples$samples
  if(inherits(main_samples, "mcmc.list")) {
    cat("  Main samples type: mcmc.list\n")
    cat("  First chain dimensions:", dim(main_samples[[1]]), "\n")
    cat("  First chain columns:", ncol(main_samples[[1]]), "\n")
    
    # Get parameters from main samples
    all_params <- colnames(main_samples[[1]])
    cat("  Main parameters monitored:", length(all_params), "\n")
    if(length(all_params) > 0) {
      cat("  Parameters:", paste(head(all_params, 10), collapse = ", "), 
          ifelse(length(all_params) > 10, "...", ""), "\n")
    }
    
    # Check for missing values in main parameters
    cat("\nChecking main parameters for missing values...\n")
    missing_check <- sapply(all_params, function(param) {
      missing_count <- 0
      for(chain in 1:length(main_samples)) {
        param_data <- main_samples[[chain]][, param]
        missing_count <- missing_count + sum(is.na(param_data) | is.infinite(param_data))
      }
      return(missing_count)
    })
    
    problem_params <- names(missing_check[missing_check > 0])
    if(length(problem_params) > 0) {
      for(param in head(problem_params, 5)) {
        cat("  ", param, ":", missing_check[param], "missing/infinite values\n")
      }
      if(length(problem_params) > 5) {
        cat("  ... and", length(problem_params) - 5, "more parameters with missing values\n")
      }
    } else {
      cat("  No missing/infinite values in main parameters - excellent!\n")
    }
    
    # Parameter summary for main parameters
    valid_params <- names(missing_check[missing_check == 0])
    if(length(valid_params) > 0) {
      cat("\nMain parameter summary (valid parameters only):\n")
      cat("  Valid parameters:", length(valid_params), "out of", length(all_params), "\n")
      
      # Focus on core parameters for summary
      core_params <- intersect(valid_params, c("precision", "rho", "site_effect_sd", 
                                               paste0("site_effect[", 1:min(5, model.dat$N.site), "]")))
      if(length(core_params) > 0) {
        cat("  Core parameter summary:\n")
        # Safely combine samples from all chains
        combined_samples <- tryCatch({
          do.call(rbind, lapply(main_samples, function(x) x[, core_params, drop = FALSE]))
        }, error = function(e) {
          cat("  Error combining samples:", e$message, "\n")
          return(NULL)
        })
        
        if(!is.null(combined_samples)) {
          print(summary(combined_samples))
        } else {
          cat("  Could not generate parameter summary\n")
        }
      }
    } else {
      cat("\nWARNING: No main parameters have valid values for summary!\n")
    }
    
    # Check for convergence issues in main parameters
    cat("\nChecking main parameters for convergence issues...\n")
    first_chain <- main_samples[[1]]
    if(any(is.na(first_chain))) {
      cat("WARNING: NA values detected in main parameters!\n")
    } else {
      cat("No NA values in main parameters - good!\n")
    }
    
    if(any(is.infinite(first_chain))) {
      cat("WARNING: Infinite values detected in main parameters!\n")
    } else {
      cat("No infinite values in main parameters - good!\n")
    }
    
  } else {
    cat("  ERROR: Main samples is not mcmc.list!\n")
    stop("Invalid main samples structure")
  }
  
  # Check latent variables from samples2
  latent_samples <- mcmc_samples$samples2
  if(inherits(latent_samples, "mcmc.list")) {
    cat("\nLatent variables (monitors2) analysis:\n")
    cat("  Latent samples type: mcmc.list\n")
    cat("  First chain dimensions:", dim(latent_samples[[1]]), "\n")
    cat("  First chain columns:", ncol(latent_samples[[1]]), "\n")
    cat("  Note: Latent variables sampled every 10th iteration (thin2=10)\n")
    
    # Check for missing values in latent variables (expected due to seasonal sampling)
    cat("  Expected: Some missing values in latent variables due to seasonal sampling\n")
  } else {
    cat("  ERROR: Latent samples is not mcmc.list!\n")
    stop("Invalid latent samples structure")
  }
  
} else if(is.list(mcmc_samples)) {
  cat("  Output is a list (not mcmc.list)\n")
  cat("  List length:", length(mcmc_samples), "\n")
  
  # Check what's in the list
  for(i in 1:min(3, length(mcmc_samples))) {
    cat("  List element", i, "type:", class(mcmc_samples[[i]]), "\n")
    if(is.matrix(mcmc_samples[[i]]) || is.data.frame(mcmc_samples[[i]])) {
      cat("    Dimensions:", dim(mcmc_samples[[i]]), "\n")
      cat("    Columns:", ncol(mcmc_samples[[i]]), "\n")
    }
  }
  
  # Try to extract parameters from first list element
  if(length(mcmc_samples) > 0) {
    first_element <- mcmc_samples[[1]]
    if(is.matrix(first_element) || is.data.frame(first_element)) {
      all_params <- colnames(first_element)
      cat("  Parameters from first element:", length(all_params), "\n")
      if(length(all_params) > 0) {
        cat("  First few parameters:", paste(head(all_params, 10), collapse = ", "), "\n")
      }
    }
  }
  
} else {
  cat("  ERROR: Unexpected MCMC output structure!\n")
  cat("  Class:", class(mcmc_samples), "\n")
  stop("Invalid MCMC output structure")
}

# Site effects convergence analysis
cat("\nSite effects convergence analysis:\n")
site_params <- grep("site_effect", colnames(main_samples[[1]]), value = TRUE)
if(length(site_params) > 0) {
  cat("Site effects parameters:", length(site_params), "\n")
  
  # Check effective sample sizes for site effects
  cat("\nEffective sample sizes for site effects:\n")
  for(param in head(site_params, 10)) {  # Limit to first 10 for readability
    ess <- effectiveSize(main_samples[[1]][, param])
    cat("  ", param, ":", round(ess, 1), "\n")
  }
  if(length(site_params) > 10) {
    cat("  ... and", length(site_params) - 10, "more site effects\n")
  }
} else {
  cat("  No site effects parameters found\n")
}

# MULTI-CHAIN CONVERGENCE ANALYSIS
cat("\nMulti-chain convergence analysis:\n")

# Gelman-Rubin diagnostics for key parameters
key_params <- intersect(valid_params, c("precision", "rho", "site_effect_sd"))
if(length(key_params) > 0 && length(main_samples) > 1) {
  cat("Gelman-Rubin diagnostics (values < 1.1 indicate convergence):\n")
  
  for(param in key_params) {
    if(param %in% valid_params) {
      # Extract parameter from all chains
      param_chains <- lapply(main_samples, function(x) x[, param])
      param_mcmc <- mcmc.list(lapply(param_chains, mcmc))
      
      # Calculate Gelman-Rubin diagnostic
      tryCatch({
        gelman_result <- gelman.diag(param_mcmc)
        cat("  ", param, ": R-hat =", round(gelman_result$psrf[1], 3), "\n")
      }, error = function(e) {
        cat("  ", param, ": R-hat calculation failed -", e$message, "\n")
      })
    }
  }
} else {
  cat("  Note: Multi-chain diagnostics require multiple chains and valid parameters\n")
}

# Effective sample size analysis
cat("\nEffective sample size analysis:\n")
for(param in key_params) {
  if(param %in% valid_params) {
    # Calculate ESS for each chain
    ess_per_chain <- sapply(main_samples, function(x) effectiveSize(x[, param]))
    total_ess <- sum(ess_per_chain)
    cat("  ", param, ":\n")
    cat("    Chain ESS:", paste(round(ess_per_chain, 1), collapse = ", "), "\n")
    cat("    Total ESS:", round(total_ess, 1), "\n")
  }
}

# Rho parameter specific analysis (NEW)
if("rho" %in% colnames(mcmc_samples[[1]])) {
  cat("\nRho parameter analysis (temporal persistence):\n")
  rho_samples <- mcmc_samples[[1]][, "rho"]
  cat("  Mean:", round(mean(rho_samples), 4), "\n")
  cat("  SD:", round(sd(rho_samples), 4), "\n")
  cat("  Range:", round(range(rho_samples), 4), "\n")
  cat("  Effective sample size:", round(effectiveSize(rho_samples), 1), "\n")
  
  # Interpret rho values
  rho_mean <- mean(rho_samples)
  if(rho_mean > 0.7) {
    cat("  Interpretation: Strong temporal persistence (high autocorrelation)\n")
  } else if(rho_mean > 0.3) {
    cat("  Interpretation: Moderate temporal persistence\n")
  } else if(rho_mean > 0) {
    cat("  Interpretation: Weak temporal persistence\n")
  } else {
    cat("  Interpretation: Negative temporal persistence (unusual)\n")
  }
}

# Precision parameter specific analysis
if("precision" %in% colnames(mcmc_samples[[1]])) {
  cat("\nPrecision parameter analysis:\n")
  precision_samples <- mcmc_samples[[1]][, "precision"]
  cat("  Mean:", round(mean(precision_samples), 4), "\n")
  cat("  SD:", round(sd(precision_samples), 4), "\n")
  cat("  Range:", round(range(precision_samples), 4), "\n")
  cat("  Effective sample size:", round(effectiveSize(precision_samples), 1), "\n")
}

# Site effect SD parameter specific analysis
if("site_effect_sd" %in% colnames(mcmc_samples[[1]])) {
  cat("\nSite effect SD parameter analysis:\n")
  site_sd_samples <- mcmc_samples[[1]][, "site_effect_sd"]
  cat("  Mean:", round(mean(site_sd_samples), 4), "\n")
  cat("  SD:", round(sd(site_sd_samples), 4), "\n")
  cat("  Range:", round(range(site_sd_samples), 4), "\n")
  cat("  Effective sample size:", round(effectiveSize(site_sd_samples), 1), "\n")
}

# Latent variable convergence analysis
cat("\nLatent variable convergence analysis:\n")
latent_params <- c("Ex", "mu")
for(param in latent_params) {
  if(param %in% colnames(mcmc_samples[[1]])) {
    param_samples <- mcmc_samples[[1]][, param]
    cat("  ", param, ":\n")
    cat("    Mean:", round(mean(param_samples), 4), "\n")
    cat("    SD:", round(sd(param_samples), 4), "\n")
    cat("    Range:", round(range(param_samples), 4), "\n")
    cat("    Effective sample size:", round(effectiveSize(param_samples), 1), "\n")
  }
}

# Comprehensive convergence summary
cat("\nConvergence summary:\n")
cat("  Main parameters (monitors):", length(all_params), "\n")
cat("  Core parameters (precision, rho, site_effect_sd):", sum(all_params %in% c("precision", "rho", "site_effect_sd")), "\n")
cat("  Site effects:", sum(grepl("site_effect", all_params)), "\n")
cat("  IMPROVED: Latent variables (Ex, mu)  via monitors2 with thin2=10\n")
cat("  Latent variables in monitors2:", ncol(mcmc_samples$samples2[[1]]), "\n")
cat("  Total MCMC iterations:", nrow(main_samples[[1]]), "\n")
cat("  Burn-in iterations:", 100, "\n")
cat("  Number of chains:", length(main_samples), "\n")

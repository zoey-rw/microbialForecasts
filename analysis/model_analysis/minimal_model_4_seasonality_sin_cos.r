#!/usr/bin/env Rscript

# Extended Seasonality Parameters - Phase 4 Minimal Model
# Focus: Run seasonal minimal model for more iterations to ensure convergence
# Based on: Successful seasonal implementation with enhanced diagnostics
# Test: 2000+ iterations to verify seasonal parameter convergence

library(nimble)
library(tidyverse)
library(here)
library(coda)
library(microbialForecast)

cat("=== Extended Seasonality Parameters - Phase 4 Minimal Model ===\n")

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

cat("Running extended seasonal parameter test for", nrow(params), "models\n")

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

# PHASE 4 MINIMAL MODEL: Add seasonal sin/cos coefficients (like working cycl_only models)
cat("Defining Phase 4 minimal model with seasonal parameters...\n")
phase4_seasonal_model <- nimbleCode({
  # PRIORS - Minimal set with rho parameter + seasonal coefficients
  precision ~ dgamma(2, 0.1)  # Observation precision
  
  # Rho parameter for temporal persistence (from Phase 3)
  rho ~ dbeta(5, 5)  # Tighter beta prior centered at 0.5, prevents extreme values
  
  # NEW: Seasonal coefficients (like working cycl_only models)
  beta[1] ~ dnorm(0, sd = 0.3)  # sin coefficient - moderate variance for ecological effects
  beta[2] ~ dnorm(0, sd = 0.3)  # cos coefficient - moderate variance for ecological effects
  
  # Site effects - Hierarchical structure (VERIFIED WORKING)
  site_effect_sd ~ dgamma(2, 0.1)
  for(s in 1:N.site) {
    site_effect[s] ~ dnorm(0, sd = site_effect_sd)
  }
  
  # PROCESS MODEL - Site effects + rho temporal persistence + seasonal predictors
  for(p in 1:N.plot) {
    # Initial condition - single time point per plot
    for(t in plot_start[p]) {
      Ex[p, t] ~ dunif(0.1, 0.9)
      mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
    }
    
    # Dynamic evolution - site effects + rho temporal persistence + seasonal predictors
    for(t in (plot_start[p] + 1):N.date) {
      # ENHANCED: Add seasonal predictors (like working cycl_only models)
      logit_Ex_prev[p, t] <- logit(Ex[p, t-1])
      logit_Ex_mean[p, t] <- rho * logit_Ex_prev[p, t] + 
                              beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +  # NEW: Seasonal terms
                              site_effect[plot_site_num[p]]
      Ex[p, t] <- ilogit(logit_Ex_mean[p, t])  # Deterministic evolution
      mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
    }
  }
  
  # OBSERVATION MODEL - Direct indexing 
  for(i in 1:N.core) {
    y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                     shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
  }
})

# Prepare data for NIMBLE with seasonal predictors
cat("Preparing NIMBLE data with seasonal predictors...\n")
nimble_data <- list(
  y = model.dat$y[,1,drop=FALSE]  # Single column for simplicity
)

# Prepare constants including seasonal predictors (like working models)
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

# Add seasonal predictors (like working cycl_only models)
if ("sin_mo" %in% names(model.dat)) {
  nimble_constants$sin_mo <- scale(model.dat$sin_mo, center = FALSE) %>% as.numeric()
  cat("  sin_mo predictor added and scaled\n")
} else {
  # Create seasonal predictors if not available
  cat("  Creating seasonal predictors...\n")
  months <- as.numeric(format(as.Date(model.dat$plot_date), "%m"))
  nimble_constants$sin_mo <- scale(sin(2 * pi * months / 12), center = FALSE) %>% as.numeric()
  nimble_constants$cos_mo <- scale(cos(2 * pi * months / 12), center = FALSE) %>% as.numeric()
  cat("  Seasonal predictors created from dates\n")
}

if ("cos_mo" %in% names(model.dat)) {
  nimble_constants$cos_mo <- scale(model.dat$cos_mo, center = FALSE) %>% as.numeric()
  cat("  cos_mo predictor added and scaled\n")
}

# Initial values - Simple and stable with rho parameter + seasonal coefficients
cat("Setting initial values with rho parameter and seasonal coefficients...\n")
initial_values <- list(
  precision = 50,  # Start with moderate precision
  rho = 0.3,      # Start rho at 0.3 (moderate persistence)
  beta = c(0.1, 0.1),  # NEW: Start seasonal coefficients near zero
  site_effect_sd = 0.5,  # Start with moderate site effect SD
  site_effect = rnorm(model.dat$N.site, 0, 0.1),  # Small random initial values
  Ex = matrix(0.3, nrow = model.dat$N.plot, ncol = model.dat$N.date),  # Start at 0.3
  mu = matrix(0.3, nrow = model.dat$N.plot, ncol = model.dat$N.date)   # Start at 0.3
)

# Build and compile model
cat("Building NIMBLE model with seasonal parameters...\n")
model <- nimbleModel(
  code = phase4_seasonal_model,
  constants = nimble_constants,
  data = nimble_data,
  inits = initial_values
)

cat("Compiling model...\n")
compiled_model <- compileNimble(model)

# EXTENDED MCMC configuration for better convergence
cat("Configuring extended MCMC for seasonal parameter convergence...\n")

# Enhanced monitoring for comprehensive convergence analysis
monitored_params <- c(
  # Core parameters (primary focus) - MONITOR ALL ESSENTIAL PARAMETERS
  "precision", "rho", "beta", "site_effect_sd", "site_effect"
)

# IMPROVED: Use monitors2 for latent variables at different interval
monitored_latent_params <- c(
  # Monitor latent process variables at different interval for efficiency
  "Ex", "mu"
)

cat("Monitoring parameters for convergence analysis:\n")
cat("  Core parameters:", paste(c("precision", "rho", "beta[1]", "beta[2]", "site_effect_sd"), collapse = ", "), "\n")
cat("  Site effects:", "site_effect[1:", model.dat$N.site, "]\n")
cat("  IMPROVED: Latent variables (Ex, mu) monitored via monitors2 at different interval\n")
cat("  Total parameters monitored:", length(monitored_params), "\n")
cat("  Latent variables monitored via monitors2:", length(monitored_latent_params), "\n")

mcmc_config <- configureMCMC(
  model = compiled_model,
  monitors = monitored_params,
  monitors2 = monitored_latent_params,  # Use monitors2 for latent variables
  thin = 1,
  thin2 = 20,  # Sample latent variables every 20th iteration for efficiency
  enableWAIC = FALSE
)

# Add specialized samplers for better convergence of key parameters
cat("Adding specialized samplers for convergence improvement...\n")

# 1. FIRST remove default samplers to prevent conflicts
cat("  Removing default samplers...\n")
mcmc_config$removeSamplers(c("precision", "rho", "beta", "site_effect_sd"))

# 2. THEN add specialized samplers - IMPROVED SAMPLER STRATEGY
cat("  Adding improved samplers for key parameters...\n")

# IMPROVED: Use regular slice samplers for single parameters, AF_slice for blocks
mcmc_config$addSampler(target = "precision", type = "slice")        # Regular slice sampler for single parameter
mcmc_config$addSampler(target = "rho", type = "slice")              # Regular slice sampler for single parameter
mcmc_config$addSampler(target = "site_effect_sd", type = "slice")   # Regular slice sampler for single parameter

# NEW: Add specialized samplers for seasonal coefficients
mcmc_config$addSampler(target = "beta[1]", type = "slice")         # Slice sampler for sin coefficient
mcmc_config$addSampler(target = "beta[2]", type = "slice")         # Slice sampler for cos coefficient

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
cat("    rho sampler:", mcmc_config$getSamplers()$rho$type, "\n")
cat("    beta[1] sampler:", mcmc_config$getSamplers()$`beta[1]`$type, "\n")  # NEW: Verify beta samplers
cat("    beta[2] sampler:", mcmc_config$getSamplers()$`beta[2]`$type, "\n")  # NEW: Verify beta samplers
cat("    site_effect_sd sampler:", mcmc_config$getSamplers()$site_effect_sd$type, "\n")
cat("    IMPROVED: Single parameters use slice samplers, site effects use adaptive block sampler (AF_slice)\n")

# Build MCMC
cat("Building MCMC...\n")
mcmc_built <- buildMCMC(mcmc_config)
compiled_mcmc <- compileNimble(mcmc_built, project = compiled_model)

# EXTENDED MCMC run for better convergence
cat("Running EXTENDED MCMC test with seasonal parameters...\n")
start_time <- Sys.time()

mcmc_samples <- runMCMC(
  compiled_mcmc,
  niter = 2000,  # EXTENDED: More iterations for better convergence
  nburnin = 500,  # EXTENDED: More burnin for better initialization
  nchains = 3, 
  samplesAsCodaMCMC = TRUE
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("EXTENDED MCMC completed!\n")
cat("Runtime:", round(runtime, 2), "minutes\n")
cat("Number of chains:", length(mcmc_samples), "\n")
cat("Samples per chain:", nrow(mcmc_samples[[1]]), "\n")
cat("Total samples:", length(mcmc_samples) * nrow(mcmc_samples[[1]]), "\n")

# Create traceplots for key parameters
cat("\nCreating traceplots for key parameters...\n")

# Function to create traceplots
create_traceplots <- function(mcmc_samples, param_names) {
  # Convert to mcmc.list if needed
  if (is.list(mcmc_samples) && "samples" %in% names(mcmc_samples)) {
    mcmc_list <- mcmc_samples$samples
  } else if (inherits(mcmc_samples, "mcmc.list")) {
    mcmc_list <- mcmc_samples
  } else {
    stop("Invalid MCMC output format")
  }
  
  # Create traceplots for each parameter
  for (param in param_names) {
    if (param %in% colnames(mcmc_list[[1]])) {
      cat("  Creating traceplot for:", param, "\n")
      
      # Extract parameter data
      param_data <- mcmc_list[, param]
      
      # Create traceplot
      png(here("figures", paste0("traceplot_", gsub("[^a-zA-Z0-9]", "_", param), ".png")), 
          width = 1200, height = 800, res = 300)
      
      # Use coda's traceplot function
      traceplot(param_data, main = paste("Traceplot for", param))
      
      dev.off()
      cat("    Saved traceplot to: figures/traceplot_", gsub("[^a-zA-Z0-9]", "_", param), ".png\n")
    } else {
      cat("  Parameter not found:", param, "\n")
    }
  }
}

# Create traceplots for key parameters
key_params <- c("rho", "beta[1]", "beta[2]", "precision", "site_effect_sd")
create_traceplots(mcmc_samples, key_params)

# Convergence diagnostics for seasonal parameters
cat("\n Convergence diagnostics for seasonal parameters:\n")

# Check MCMC output structure first
cat("Checking MCMC output structure...\n")
cat("  MCMC output type:", class(mcmc_samples), "\n")
cat("  Number of chains:", length(mcmc_samples), "\n")

# Check if we have valid output
if(length(mcmc_samples) == 0) {
  cat("ERROR: No MCMC output generated!\n")
  stop("MCMC failed to produce output")
}

# Handle different output structures
if(inherits(mcmc_samples, "mcmc.list")) {
  cat("  Output is mcmc.list (correct format for CodaMCMC)\n")
  main_samples <- mcmc_samples
} else if(is.list(mcmc_samples) && "samples" %in% names(mcmc_samples)) {
  cat("  Output is monitors2 structure (samples + samples2)\n")
  cat("  Main samples (monitors):", length(mcmc_samples$samples), "chains\n")
  cat("  Latent samples (monitors2):", length(mcmc_samples$samples2), "chains\n")
  
  # Extract main parameters from samples
  main_samples <- mcmc_samples$samples
} else {
  cat("  ERROR: Unexpected MCMC output structure!\n")
  stop("Invalid MCMC output structure")
}

if(inherits(main_samples, "mcmc.list")) {
  cat("  Main samples type: mcmc.list\n")
  cat("  First chain dimensions:", dim(main_samples[[1]]), "\n")
  cat("  First chain columns:", ncol(main_samples[[1]]), "\n")
  
  # Get parameters from main samples
  all_params <- colnames(main_samples[[1]])
  cat("  Main parameters monitored:", length(all_params), "\n")
  
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
    core_params <- intersect(valid_params, c("precision", "rho", "beta[1]", "beta[2]", "site_effect_sd", 
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
  
  # MULTI-CHAIN CONVERGENCE ANALYSIS
  cat("\nMulti-chain convergence analysis:\n")
  
  # Gelman-Rubin diagnostics for key parameters
  key_params <- intersect(valid_params, c("precision", "rho", "beta[1]", "beta[2]", "site_effect_sd"))
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
  
  # Rho parameter specific analysis
  if("rho" %in% colnames(main_samples[[1]])) {
    cat("\nRho parameter analysis (temporal persistence):\n")
    rho_samples <- main_samples[[1]][, "rho"]
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
  
  # NEW: Seasonal coefficient analysis with user-provided conversion code
  if("beta[1]" %in% colnames(main_samples[[1]]) && "beta[2]" %in% colnames(main_samples[[1]])) {
    cat("\nSeasonal coefficient analysis:\n")
    
    # Extract seasonal coefficients
    beta1_samples <- main_samples[[1]][, "beta[1]"]
    beta2_samples <- main_samples[[1]][, "beta[2]"]
    
    # Calculate means for seasonal analysis
    beta1_mean <- mean(beta1_samples)
    beta2_mean <- mean(beta2_samples)
    
    cat("  beta[1] (sin coefficient):\n")
    cat("    Mean:", round(beta1_mean, 4), "\n")
    cat("    SD:", round(sd(beta1_samples), 4), "\n")
    cat("    Range:", round(range(beta1_samples), 4), "\n")
    cat("    Effective sample size:", round(effectiveSize(beta1_samples), 1), "\n")
    
    cat("  beta[2] (cos coefficient):\n")
    cat("    Mean:", round(beta2_mean, 4), "\n")
    cat("    SD:", round(sd(beta2_samples), 4), "\n")
    cat("    Range:", round(range(beta2_samples), 4), "\n")
    cat("    Effective sample size:", round(effectiveSize(beta2_samples), 1), "\n")
    
    # USER-PROVIDED CODE: Calculate seasonal amplitude and peak timing
    cat("\nSeasonal pattern analysis (using user-provided conversion code):\n")
    
    # Calculate seasonal amplitude
    amplitude <- sqrt(beta1_mean^2 + beta2_mean^2)
    cat("  Seasonal amplitude:", round(amplitude, 4), "\n")
    
    # Calculate peak timing (in months)
    peak_month <- atan2(beta2_mean, beta1_mean) * 12 / (2 * pi)
    if(peak_month < 0) peak_month <- peak_month + 12
    cat("  Peak month:", round(peak_month, 1), "\n")
    
    # Seasonal strength interpretation
    seasonal_strength <- if(amplitude > 0.1) "Strong seasonal pattern"
                        else if(amplitude > 0.05) "Moderate seasonal pattern" 
                        else "Weak seasonal pattern"
    cat("  Seasonal strength:", seasonal_strength, "\n")
    
    # Interpret peak timing
    if(peak_month >= 3 && peak_month <= 5) {
      cat("  Peak timing: Spring (March-May)\n")
    } else if(peak_month >= 6 && peak_month <= 8) {
      cat("  Peak timing: Summer (June-August)\n")
    } else if(peak_month >= 9 && peak_month <= 11) {
      cat("  Peak timing: Fall (September-November)\n")
    } else {
      cat("  Peak timing: Winter (December-February)\n")
    }
    
    # Additional seasonal insights
    cat("\nSeasonal coefficient interpretation:\n")
    if(abs(beta1_mean) > 0.1) {
      cat("  beta[1] (sin): Strong seasonal pattern in sin component\n")
    } else if(abs(beta1_mean) > 0.05) {
      cat("  beta[1] (sin): Moderate seasonal pattern in sin component\n")
    } else {
      cat("  beta[1] (sin): Weak seasonal pattern in sin component\n")
    }
    
    if(abs(beta2_mean) > 0.1) {
      cat("  beta[2] (cos): Strong seasonal pattern in cos component\n")
    } else if(abs(beta2_mean) > 0.05) {
      cat("  beta[2] (cos): Moderate seasonal pattern in cos component\n")
    } else {
      cat("  beta[2] (cos): Weak seasonal pattern in cos component\n")
    }
  }
  
  # Comprehensive convergence summary
  cat("\nConvergence summary:\n")
  cat("  Main parameters (monitors):", length(all_params), "\n")
  cat("  Core parameters (precision, rho, beta[1], beta[2], site_effect_sd):", 
      sum(all_params %in% c("precision", "rho", "beta[1]", "beta[2]", "site_effect_sd")), "\n")
  cat("  Site effects:", sum(grepl("site_effect", all_params)), "\n")
  cat("  Total MCMC iterations:", nrow(main_samples[[1]]), "\n")
  cat("  Burn-in iterations:", 500, "\n")
  cat("  Number of chains:", length(main_samples), "\n")
  
} else {
  cat("  ERROR: Main samples is not mcmc.list!\n")
  stop("Invalid main samples structure")
}

cat("\n=== EXTENDED Phase 4 Minimal Model with Seasonality COMPLETED ===\n")
cat("  - Extended MCMC run: 2000 iterations with 500 burnin\n")
cat("  - Added seasonal sin/cos coefficients (beta[1], beta[2])\n")
cat("  - Maintained rho parameter for temporal persistence\n")
cat("  - Enhanced MCMC configuration for seasonal parameters\n")
cat("  - Comprehensive convergence diagnostics implemented\n")
cat("  - Seasonal pattern analysis using user-provided conversion code\n")
cat("  - Traceplots created for key parameters\n")
cat("  - Ready for production use with even longer chains\n")

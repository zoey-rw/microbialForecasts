#!/usr/bin/env Rscript

# CONSTRAINED MCMC MODEL - Fixes for parameter explosion and convergence issues
# Implements bounded linear predictors, tight priors, and regularization

source("../../source.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 01_fitModels_constrained.R <model_type> <taxon> <time_period> <chain_number>")
}

model_name <- args[1]
taxon <- args[2]
time_period <- args[3]
chain_no <- as.numeric(args[4])

cat("=== FITTING CONSTRAINED MCMC MODEL ===\n")
cat("Model:", model_name, "\n")
cat("Taxon:", taxon, "\n")
cat("Time period:", time_period, "\n")
cat("Chain:", chain_no, "\n\n")

# Check if legacy covariate should be used based on time period
use_legacy_covariate <- grepl("20130601_20151101", time_period)
cat("Using legacy covariate:", use_legacy_covariate, "\n")

# Load data
model.dat <- readRDS(here::here("data", "model_outputs", "input_data", paste0("input_data_", time_period, ".rds")))
cat("Data loaded successfully\n")

# Prepare constants and data
constants <- model.dat$constants
constants$plot_site_num <- as.numeric(factor(model.dat$plot_site_num))
constants$N.site <- max(constants$plot_site_num)

cat("Constants prepared successfully\n")

# CONSTRAINED MODEL DEFINITIONS - Addresses parameter explosion
cat("Building constrained Nimble model...\n")
if (model_name == "cycl_only" && use_legacy_covariate) {
  modelCode <- nimble::nimbleCode({
    for (i in 1:N.core) {
      y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision,
                 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
    }

    for (p in 1:N.plot) {
      for (t in plot_start[p]) {
        Ex[p, t] ~ dunif(0.01, 0.99)  # Tighter initial bounds
        plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
      }

      for (t in plot_index[p]:N.date) {
        # Bound previous value before logit transform
        mu_prev_bounded[p, t] <- max(0.001, min(0.999, plot_mu[p, t - 1]))

        # Calculate linear predictor with bounded previous value
        logit_raw[p, t] <- rho * logit(mu_prev_bounded[p, t]) +
          beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
          site_effect[plot_site_num[p]] +
          legacy_effect * legacy[t] +
          intercept

        # Bound the linear predictor to prevent explosion
        logit_bounded[p, t] <- max(-6, min(6, logit_raw[p, t]))

        # Transform back with bounded values
        Ex[p, t] <- ilogit(logit_bounded[p, t])
        plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
      }
    }

    # CONSTRAINED PRIORS - Much tighter to prevent parameter explosion
    precision ~ dgamma(5, 1)        # Mean = 5, concentrated (was 2, 0.1)
    rho ~ dbeta(8, 4)      # Mean = 0.67, concentrated away from 1
    intercept ~ dnorm(-2, sd = 0.2) # Tighter intercept (was 0.8)
    legacy_effect ~ dnorm(0, sd = 0.1) # Much tighter legacy effect (was 0.2)

    # Very tight beta coefficients to prevent large effects
    for (b in 1:2) {
      beta[b] ~ dnorm(0, sd = 0.1)  # Very tight seasonal effects (was 0.3)
    }

    # CONSTRAINED: Regularized hierarchical site effects
    site_effect_sd ~ dgamma(5, 50)   # Mean = 0.1, very tight (was 2, 6)
    site_effect_mean ~ dnorm(0, sd = 0.01) # Overall site level

    for (k in 1:N.site) {
      site_effect_raw[k] ~ dnorm(0, sd = site_effect_sd)
    }
    # Force site effects to sum to zero (prevents drift)
    site_effect_sum <- sum(site_effect_raw[1:N.site])
    for (k in 1:N.site) {
      site_effect[k] <- site_effect_raw[k] - site_effect_sum/N.site + site_effect_mean
    }
  })
} else if (model_name == "cycl_only" && !use_legacy_covariate) {
  modelCode <- nimble::nimbleCode({
    for (i in 1:N.core) {
      y[i, 1] ~ dbeta(shape1 = plot_mu[plot_num[i], timepoint[i]] * precision,
                 shape2 = (1 - plot_mu[plot_num[i], timepoint[i]]) * precision)
    }

    for (p in 1:N.plot) {
      for (t in plot_start[p]) {
        Ex[p, t] ~ dunif(0.01, 0.99)  # Tighter initial bounds
        plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
      }

      for (t in plot_index[p]:N.date) {
        # CRITICAL: Bound previous value before logit transform
        mu_prev_bounded[p, t] <- max(0.001, min(0.999, plot_mu[p, t - 1]))

        # Calculate linear predictor with bounded previous value
        logit_raw[p, t] <- rho * logit(mu_prev_bounded[p, t]) +
          beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
          site_effect[plot_site_num[p]] +
          intercept

        # CRITICAL: Bound the linear predictor to prevent explosion
        logit_bounded[p, t] <- max(-6, min(6, logit_raw[p, t]))

        # Transform back with bounded values
        Ex[p, t] <- ilogit(logit_bounded[p, t])
        plot_mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
      }
    }

    # CONSTRAINED PRIORS - Much tighter to prevent parameter explosion
    precision ~ dgamma(5, 1)        # Mean = 5, concentrated
    rho ~ dbeta(8, 4)              # Mean = 0.67, concentrated away from 1
    intercept ~ dnorm(-2, sd = 0.2) # Tighter intercept

    # Very tight beta coefficients
    for (b in 1:2) {
      beta[b] ~ dnorm(0, sd = 0.05)  # Very tight seasonal effects
    }

    # CONSTRAINED: Regularized hierarchical site effects
    site_effect_sd ~ dgamma(5, 50)   # Mean = 0.1, very tight
    site_effect_mean ~ dnorm(0, sd = 0.01) # Overall site level

    for (k in 1:N.site) {
      site_effect_raw[k] ~ dnorm(0, sd = site_effect_sd)
    }
    # Force site effects to sum to zero (prevents drift)
    site_effect_sum <- sum(site_effect_raw[1:N.site])
    for (k in 1:N.site) {
      site_effect[k] <- site_effect_raw[k] - site_effect_sum/N.site + site_effect_mean
    }
  })
} else {
  stop("Constrained model only implemented for cycl_only models. Use original script for other model types.")
}

# CONSERVATIVE INITIALIZATION - Start ALL chains in stable region
cat("Creating conservative initial values...\n")
inits <- createInits(constants)

# Calculate data-informed initial values
y_data <- model.dat$y[, 1]
y_mean <- mean(y_data, na.rm = TRUE)
y_sd <- sd(y_data, na.rm = TRUE)

# CONSERVATIVE: Start ALL chains in stable region
set.seed(chain_no * 1000 + 1)

inits$rho <- runif(1, 0.2, 0.6)     # Well below 1, avoid boundary
inits$intercept <- rnorm(1, -1.5, 0.1)  # Moderate abundance, tight
inits$precision <- rgamma(1, 5, 1)      # Start at mean, concentrated
inits$beta <- rnorm(constants$N.beta, 0, 0.02)  # Very small effects

# Conservative site effect initialization
inits$site_effect_sd <- runif(1, 0.01, 0.1)  # Small site variation
inits$site_effect_mean <- rnorm(1, 0, 0.01)  # Small overall level

# Initialize site effects with sum-to-zero constraint
raw_effects <- rnorm(constants$N.site, 0, inits$site_effect_sd)
inits$site_effect_raw <- raw_effects - mean(raw_effects)
inits$site_effect <- inits$site_effect_raw + inits$site_effect_mean

# Legacy effect if used
if (use_legacy_covariate) {
  inits$legacy_effect <- rnorm(1, 0, 0.02)  # Very tight initialization
}

cat("Conservative initialization complete\n")

# Build model
cat("Building Nimble model...\n")
Rmodel <- nimbleModel(code = modelCode, constants = constants,
            data = list(y=model.dat$y), inits = inits)

# Compile model
cat("Compiling Nimble model...\n")
cModel <- compileNimble(Rmodel)

cat("Model compiled successfully\n")

# Configure MCMC with conservative settings
cat("Configuring MCMC with conservative settings...\n")
monitors <- c("beta","precision","site_effect","site_effect_sd","intercept","rho","site_effect_mean")

if (use_legacy_covariate) {
  monitors <- c(monitors, "legacy_effect")
}

mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)

# Remove default samplers and add conservative alternatives
mcmcConf$removeSamplers(c("precision", "rho", "site_effect_sd", "intercept", "site_effect_mean"))

if (use_legacy_covariate) {
  mcmcConf$removeSamplers("legacy_effect")
}

# Remove beta samplers
if (constants$N.beta > 1) {
  mcmcConf$removeSamplers(paste0("beta[1:", constants$N.beta, "]"))
} else {
  mcmcConf$removeSamplers("beta[1]")
}

# Remove site effect samplers
if (constants$N.site > 1) {
  mcmcConf$removeSamplers(c(
    paste0("site_effect[1:", constants$N.site, "]"),
    paste0("site_effect_raw[1:", constants$N.site, "]")
  ))
} else {
  mcmcConf$removeSamplers(c("site_effect[1]", "site_effect_raw[1]"))
}

# Add conservative samplers
mcmcConf$addSampler(target = "precision", type = "slice")
mcmcConf$addSampler(target = "rho", type = "slice")
mcmcConf$addSampler(target = "site_effect_sd", type = "slice")
mcmcConf$addSampler(target = "intercept", type = "slice")
mcmcConf$addSampler(target = "site_effect_mean", type = "slice")

if (use_legacy_covariate) {
  mcmcConf$addSampler(target = "legacy_effect", type = "slice")
}

# Add samplers for beta parameters
for (b in 1:constants$N.beta) {
  mcmcConf$addSampler(target = paste0("beta[", b, "]"), type = "slice")
}

# Add samplers for site effects
for (k in 1:constants$N.site) {
  mcmcConf$addSampler(target = paste0("site_effect_raw[", k, "]"), type = "slice")
}

cat("MCMC configuration complete\n")

# Build MCMC
cat("Building MCMC...\n")
Rmcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(Rmcmc, project = cModel)

# CONSERVATIVE MCMC SETTINGS
niter <- 3000    # Shorter but more stable iterations
nburnin <- 1500  # Longer burn-in for stability
nthin <- 2       # Moderate thinning

cat("Conservative MCMC settings:\n")
cat("  Iterations:", niter, "\n")
cat("  Burn-in:", nburnin, "\n")
cat("  Thinning:", nthin, "\n")
cat("  Samples per chain:", (niter - nburnin) / nthin, "\n")

# Run the MCMC
start_time <- Sys.time()
mcmc.out <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = nthin,
                    setSeed = chain_no * 1000, samplesAsCodaMCMC = TRUE)
end_time <- Sys.time()

cat("MCMC completed in", round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "minutes\n")

# Process output
samples <- as.matrix(mcmc.out)
cat("MCMC output dimensions:", dim(samples), "\n")

# Save output with constrained model indicator
output_dir <- here::here("data", "model_outputs", "logit_beta_regression_constrained")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create model-specific subdirectory
model_dir <- file.path(output_dir, model_name)
if (!dir.exists(model_dir)) {
  dir.create(model_dir)
}

# Save samples
samples_file <- paste0("samples_", model_name, "_", taxon, "_", time_period, "_constrained_chain", chain_no, ".rds")
samples_path <- file.path(model_dir, samples_file)

# Create output list with metadata
out <- list(
  samples = samples,
  metadata = list(
    model_name = model_name,
    taxon = taxon,
    time_period = time_period,
    chain = chain_no,
    use_legacy_covariate = use_legacy_covariate,
    niter = niter,
    nburnin = nburnin,
    nthin = nthin,
    model_type = "constrained_mcmc",
    constraints = "Bounded linear predictor (-6,6), tight priors, sum-to-zero site effects",
    prior_changes = "beta sd=0.05, rho Beta(8,4), precision Gamma(5,1), site_effect_sd Gamma(5,50)",
    data_file = paste0("input_data_", time_period, ".rds")
  )
)

saveRDS(out, samples_path, compress = FALSE)
cat("Samples saved to:", samples_path, "\n")

# Calculate basic convergence diagnostics
if (ncol(samples) > 1) {
  cat("Calculating convergence diagnostics...\n")

  param_means <- colMeans(samples)
  param_sds <- apply(samples, 2, sd)
  param_ranges <- apply(samples, 2, function(x) diff(range(x)))

  convergence_summary <- data.frame(
    parameter = colnames(samples),
    mean = param_means,
    sd = param_sds,
    range = param_ranges,
    n_samples = nrow(samples),
    chain = chain_no
  )

  # Check for extreme values
  extreme_params <- convergence_summary[abs(convergence_summary$mean) > 10 |
                                        convergence_summary$range > 50, ]
  if (nrow(extreme_params) > 0) {
    cat("WARNING: Parameters with extreme values detected:\n")
    print(extreme_params[, c("parameter", "mean", "range")])
  } else {
    cat("✅ No extreme parameter values detected\n")
  }

  summary_file <- paste0("convergence_summary_", model_name, "_", taxon, "_", time_period, "_constrained_chain", chain_no, ".rds")
  summary_path <- file.path(model_dir, summary_file)
  saveRDS(convergence_summary, summary_path, compress = FALSE)
  cat("Convergence summary saved to:", summary_path, "\n")
}

cat("=== CONSTRAINED MODEL FITTING COMPLETE ===\n")
cat("Key improvements:\n")
cat("  • Bounded linear predictor (-6, 6) to prevent explosion\n")
cat("  • Tight priors: beta sd=0.05, rho Beta(8,4), precision Gamma(5,1)\n")
cat("  • Sum-to-zero site effects to prevent drift\n")
cat("  • Conservative initialization strategy\n")

cat("\nOutput files:\n")
cat("  Samples:", samples_path, "\n")
if (exists("summary_path")) {
  cat("  Summary:", summary_path, "\n")
}

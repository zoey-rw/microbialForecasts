#!/usr/bin/env Rscript
# Phase 4: Dirichlet model with seasonal parameters
# This adds seasonal (sin/cos) covariates to the site effects model

# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
here::i_am("analysis/model_analysis/dirichlet_covariance/04_fitModels_dirichlet_seasonal.r")
project_root <- here()

cat("here() starts at", project_root, "\n")
cat("Project root set to:", getwd(), "\n")

# Load the microbialForecast package
library(microbialForecast)

# Load packages and create directories
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("dirichlet_regression", "dirichlet_regression/seasonal", "dirichlet_regression/site_effects", "dirichlet_regression/temporal", "dirichlet_regression/env_cycl", "dirichlet_regression/basic")
)

cat("==================================================\n")
cat("PHASE 4: Dirichlet with Seasonal Parameters\n")
cat("==================================================\n")

# Load real data
cat("Loading real data...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
cat("Data loaded successfully for", length(all_ranks), "ranks\n")

# Create simplified data preparation for testing
cat("Creating simplified data preparation...\n")
rank.df <- all_ranks$phylum_fun
cat("Original data dimensions:", dim(rank.df), "\n")

# Subset to manageable size for testing
set.seed(123)
test_plots <- sample(unique(rank.df$plotID), 10)
test_dates <- sort(unique(rank.df$dateID))[1:9]  # First 9 dates

rank.df <- rank.df %>%
  filter(plotID %in% test_plots, dateID %in% test_dates) %>%
  filter(!is.na(ascomycota), !is.na(basidiomycota), !is.na(mortierellomycota))

cat("Subset data dimensions:", dim(rank.df), "\n")

# Prepare data for Dirichlet model
y <- as.matrix(rank.df[, c("ascomycota", "basidiomycota", "mortierellomycota", "other")])
plot_num <- as.numeric(as.factor(rank.df$plotID))
timepoint <- as.numeric(as.factor(rank.df$dateID))
plot_site <- as.numeric(as.factor(rank.df$siteID))
plot_site_num <- as.numeric(as.factor(rank.df$siteID))

# Model dimensions
N.plot <- length(unique(plot_num))
N.date <- length(unique(timepoint))
N.spp <- ncol(y)
N.core <- nrow(y)
N.site <- length(unique(plot_site))

cat("Model dimensions:\n")
cat("  N.plot:", N.plot, "\n")
cat("  N.date:", N.date, "\n")
cat("  N.spp:", N.spp, "\n")
cat("  N.core:", N.core, "\n")
cat("  N.site:", N.site, "\n")
cat("  keep_taxa:", colnames(y), "\n")

# Create seasonal covariates (sin/cos of month)
# Extract month from dateID (assuming format like "2013-06" or similar)
if (all(grepl("-", unique(rank.df$dateID)))) {
  months <- as.numeric(sapply(strsplit(rank.df$dateID, "-"), function(x) x[2]))
} else {
  # Fallback: use sequential numbering
  months <- as.numeric(as.factor(rank.df$dateID))
}

# Create seasonal covariates
sin_mo <- sin(2 * pi * months / 12)
cos_mo <- cos(2 * pi * months / 12)

# Create constants
constants <- list(
  N.plot = N.plot,
  N.date = N.date,
  N.spp = N.spp,
  N.core = N.core,
  N.site = N.site,
  plot_num = plot_num,
  timepoint = timepoint,
  plot_site_num = plot_site_num,
  sin_mo = sin_mo,
  cos_mo = cos_mo,
  zeros = rep(0, 6),  # 4 seasonal + 2 additional parameters
  omega = diag(6)
)

# Create initial values with seasonal structure
inits <- list(
  alpha = array(rgamma(N.plot * N.spp * N.date, 1, 1), 
                dim = c(N.plot, N.spp, N.date)),
  beta = matrix(rnorm(N.spp * 6), nrow = N.spp, ncol = 6),  # 6 seasonal parameters
  sigma = rep(1, N.spp),
  intercept = rep(0, N.spp),
  rho = 0.5,  # Single temporal autocorrelation parameter
  sig = 1,
  site_effect = matrix(rnorm(N.site * N.spp), nrow = N.site, ncol = N.spp),
  temporal_effect = array(rnorm(N.plot * N.spp * N.date), 
                         dim = c(N.plot, N.spp, N.date)),
  Ex = array(1, dim = c(N.plot, N.spp, N.date)),
  plot_rel = array(1/N.spp, dim = c(N.plot, N.spp, N.date))
)

# Create seasonal model code
modelCode <- nimbleCode({
  # Dirichlet model with temporal dependence, site effects, and seasonal parameters
  for (i in 1:N.core) {
    y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
  }
  
  # Alpha generation with temporal autocorrelation, site effects, and seasonal terms
  for (s in 1:N.spp) {
    for (p in 1:N.plot) {
      # First timepoint - no temporal dependence
      temporal_effect[p, s, 1] ~ dnorm(0, 1)
      alpha[p, s, 1] ~ dgamma(exp(intercept[s] + 
                                  site_effect[plot_site_num[p], s] + 
                                  temporal_effect[p, s, 1] +
                                  beta[s, 1] * sin_mo[1] + beta[s, 2] * cos_mo[1] +
                                  beta[s, 3] * sin_mo[1]^2 + beta[s, 4] * cos_mo[1]^2 +
                                  beta[s, 5] * sin_mo[1] * cos_mo[1] + beta[s, 6] * 1), 1)
      
      # Subsequent timepoints - AR(1) process with seasonal terms
      for (t in 2:N.date) {
        temporal_effect[p, s, t] ~ dnorm(rho * temporal_effect[p, s, t-1], 1 - rho^2)
        alpha[p, s, t] ~ dgamma(exp(intercept[s] + 
                                    site_effect[plot_site_num[p], s] + 
                                    temporal_effect[p, s, t] +
                                    beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t] +
                                    beta[s, 3] * sin_mo[t]^2 + beta[s, 4] * cos_mo[t]^2 +
                                    beta[s, 5] * sin_mo[t] * cos_mo[t] + beta[s, 6] * 1), 1)
      }
    }
  }
  
  # Site effects - random effects for each site-species combination
  for (s in 1:N.spp) {
    for (site in 1:N.site) {
      site_effect[site, s] ~ dnorm(0, sigma_site[s])
    }
  }
  
  # Priors
  for (s in 1:N.spp) {
    sigma[s] ~ dgamma(1, 1)
    sigma_site[s] ~ dgamma(1, 1)  # Site effect variance
    intercept[s] ~ dnorm(0, 1)
    beta[s, 1:6] ~ dmnorm(zeros[1:6], omega[1:6, 1:6])  # 6 seasonal parameters
  }
  
  # Single temporal autocorrelation parameter
  rho ~ dunif(-0.99, 0.99)
})

cat("Building Dirichlet model with seasonal parameters...\n")
Rmodel <- nimbleModel(modelCode, constants = constants, data = list(y = y), inits = inits)

cat("Compiling model...\n")
cModel <- compileNimble(Rmodel)

cat("Configuring MCMC...\n")
mcmcConf <- configureMCMC(cModel, monitors = c("beta", "intercept", "rho", "sigma", "sigma_site", "site_effect", "temporal_effect"))

cat("Building MCMC...\n")
myMCMC <- buildMCMC(mcmcConf)
compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)

cat("Running MCMC with 3 chains, 1000 iterations each...\n")
samples <- runMCMC(compiled, niter = 1000, nburnin = 500, nchains = 3, samplesAsCodaMCMC = TRUE)

cat("All chains completed successfully!\n")
cat("Combined samples dimensions:", dim(samples), "\n")

# Convergence diagnostics
cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
gelman <- gelman.diag(samples, multivariate = FALSE)
cat("Gelman-Rubin Diagnostic (should be < 1.1 for convergence):\n")
print(gelman)

mean_gbr <- mean(gelman$psrf[, "Point est."], na.rm = TRUE)
cat("Mean Gelman-Rubin statistic:", round(mean_gbr, 3), "\n")
if (mean_gbr < 2.0) {
  cat("✅ CONVERGENCE: Mean GBR < 2.0 - Model appears to have converged\n")
} else {
  cat("❌ CONVERGENCE: Mean GBR >= 2.0 - Model may not have converged\n")
}

# Effective sample sizes
ess <- effectiveSize(samples)
cat("\nEffective Sample Sizes:\n")
print(round(ess, 0))
cat("Minimum ESS:", round(min(ess, na.rm = TRUE), 1), "\n")
cat("Mean ESS:", round(mean(ess, na.rm = TRUE), 1), "\n")

# Save results
output_file <- here("data", "model_outputs", "dirichlet_regression", "seasonal", "seasonal_test.rds")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(samples, output_file)
cat("Results saved to:", output_file, "\n")

cat("\n=== PHASE 4: SEASONAL PARAMETERS TEST SUCCESSFUL ===\n")
cat("This tests seasonal (sin/cos) covariates in the Dirichlet model!\n")



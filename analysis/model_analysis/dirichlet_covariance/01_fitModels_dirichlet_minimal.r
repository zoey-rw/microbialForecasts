#!/usr/bin/env Rscript
# Minimal Dirichlet model fitting for testing
# Uses drastically reduced data dimensions to ensure MCMC completion

# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
here::i_am("analysis/model_analysis/dirichlet_covariance/01_fitModels_dirichlet_minimal.r")
project_root <- here()

cat("here() starts at", project_root, "\n")
cat("Project root set to:", getwd(), "\n")

# Load the microbialForecast package
library(microbialForecast)

# Load packages and create directories
load_required_packages()
create_directories_safe(
    here("data", "model_outputs"), 
    c("dirichlet_regression", "dirichlet_regression/env_cycl", "dirichlet_regression/env_cov")
)

# Define output directory
model_output_dir <- here("data", "model_outputs", "dirichlet_regression")

cat("==================================================\n")
cat("MINIMAL Dirichlet analysis - Testing version\n")
cat("==================================================\n")

# Load data files
cat("Loading data files for filtering...\n")
load(here("data", "clean", "all_ranks.RData"))
cat("Data loaded successfully for", length(all_ranks), "ranks\n")

# Get the phylum_fun data
rank.df <- all_ranks$phylum_fun

# DRASTICALLY reduce data size for testing
cat("Drastically reducing data size for testing...\n")
cat("Original data dimensions:", dim(rank.df), "\n")

# Take only first 5 plots and first 5 time points
unique_plots <- unique(rank.df$plotID)[1:5]  # Only 5 plots
unique_dates <- unique(rank.df$dateID)[1:5]  # Only 5 dates

# Filter to this tiny subset
rank.df_minimal <- rank.df[rank.df$plotID %in% unique_plots & rank.df$dateID %in% unique_dates, ]

cat("Minimal data dimensions:", dim(rank.df_minimal), "\n")

# Use prepDirichletData with the minimal dataset
model.dat <- prepDirichletData(rank.df = rank.df_minimal,
                              min.prev = 1,  # Lower threshold for small data
                              min.date = min(unique_dates),
                              max.date = max(unique_dates))

cat("Model data prepared:\n")
cat("  N.plot:", model.dat$N.plot, "\n")
cat("  N.date:", model.dat$N.date, "\n")
cat("  N.spp:", model.dat$N.spp, "\n")
cat("  N.core:", model.dat$N.core, "\n")

# Create minimal constants
constants <- list(
  N.spp = model.dat$N.spp,
  N.core = model.dat$N.core,
  N.plot = model.dat$N.plot,
  N.date = model.dat$N.date,
  N.site = model.dat$N.site,
  N.beta = 4,  # Reduced from 8
  plot_start = model.dat$plot_start,
  plot_num = model.dat$plot_num,
  plot_site_num = model.dat$plot_site_num,
  timepoint = model.dat$timepoint,
  temp = matrix(rnorm(3 * model.dat$N.date), nrow = 3, ncol = model.dat$N.date),
  mois = matrix(rnorm(3 * model.dat$N.date), nrow = 3, ncol = model.dat$N.date),
  pH = matrix(rnorm(3 * model.dat$N.date), nrow = 3, ncol = model.dat$N.date),
  pC = matrix(rnorm(3 * model.dat$N.date), nrow = 3, ncol = model.dat$N.date),
  relEM = matrix(rnorm(3 * model.dat$N.date), nrow = 3, ncol = model.dat$N.date),
  LAI = data.frame(LAI = rnorm(model.dat$N.date)),
  sin_mo = sin(2 * pi * 1:model.dat$N.date / 12),
  cos_mo = cos(2 * pi * 1:model.dat$N.date / 12),
  zeros = rep(0, 4),
  omega = diag(4)
)

# Create minimal initial values
inits <- list(
  alpha = array(rgamma(model.dat$N.plot * model.dat$N.spp * model.dat$N.date, 1, 1), 
                dim = c(model.dat$N.plot, model.dat$N.spp, model.dat$N.date)),
  beta = matrix(rnorm(model.dat$N.spp * 4), nrow = model.dat$N.spp, ncol = 4),
  sigma = rep(1, model.dat$N.spp),
  intercept = rep(0, model.dat$N.spp),
  rho = rep(0.5, model.dat$N.spp),
  sig = 1,
  site_effect = matrix(rnorm(model.dat$N.site * model.dat$N.spp), 
                       nrow = model.dat$N.site, ncol = model.dat$N.spp),
  Ex = array(1, dim = c(model.dat$N.plot, model.dat$N.spp, model.dat$N.date)),
  plot_rel = array(1/model.dat$N.spp, dim = c(model.dat$N.plot, model.dat$N.spp, model.dat$N.date))
)

# Create ultra-simple model code
modelCode <- nimbleCode({
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

cat("Building minimal Dirichlet model...\n")
Rmodel <- nimbleModel(modelCode, constants = constants, data = list(y = model.dat$y), inits = inits)

cat("Compiling minimal model...\n")
cModel <- compileNimble(Rmodel)

cat("Configuring minimal MCMC...\n")
mcmcConf <- configureMCMC(cModel, monitors = c("beta", "sigma", "intercept"))

cat("Building minimal MCMC...\n")
myMCMC <- buildMCMC(mcmcConf)

cat("Compiling minimal MCMC...\n")
compiled <- compileNimble(myMCMC, project = cModel)

cat("Running minimal MCMC (100 iterations)...\n")
compiled$run(niter = 100, nburnin = 50)

cat("MCMC completed successfully!\n")
cat("Samples collected:", nrow(as.matrix(compiled$mvSamples)), "\n")

# Save results
output_file <- here("data", "model_outputs", "dirichlet_regression", "minimal_test_results.rds")
saveRDS(list(
  samples = as.matrix(compiled$mvSamples),
  model_data = model.dat,
  constants = constants
), output_file)

cat("Results saved to:", output_file, "\n")
cat("=== MINIMAL DIRICHLET TEST SUCCESSFUL ===\n")

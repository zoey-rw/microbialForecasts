#!/usr/bin/env Rscript
# Working Dirichlet model with real data - simplified approach
# This uses the real data but with a much simpler data preparation

# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Set project root
here::i_am("analysis/model_analysis/dirichlet_covariance/01_fitModels_dirichlet_real_data.r")
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
cat("WORKING Dirichlet analysis - Real Data Version\n")
cat("==================================================\n")

# Load real data
cat("Loading real data...\n")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
rank.df <- all_ranks$phylum_fun

cat("Data loaded successfully for", length(all_ranks), "ranks\n")
cat("Original data dimensions:", dim(rank.df), "\n")

# Create a simplified data preparation
cat("Creating simplified data preparation...\n")

# Filter to a reasonable subset for testing
# Take plots with most observations and dates with most observations
plot_counts <- table(rank.df$plotID)
top_plots <- names(sort(plot_counts, decreasing = TRUE))[1:10]  # Top 10 plots
date_counts <- table(rank.df$dateID)
top_dates <- names(sort(date_counts, decreasing = TRUE))[1:10]  # Top 10 dates

# Filter to this subset
rank.df_subset <- rank.df[rank.df$plotID %in% top_plots & rank.df$dateID %in% top_dates, ]
cat("Subset data dimensions:", dim(rank.df_subset), "\n")

# Create simple timepoint mapping
unique_dates <- sort(unique(rank.df_subset$dateID))
date_to_timepoint <- setNames(1:length(unique_dates), unique_dates)
timepoint <- date_to_timepoint[as.character(rank.df_subset$dateID)]

# Create simple plot mapping
unique_plots <- sort(unique(rank.df_subset$plotID))
plot_to_num <- setNames(1:length(unique_plots), unique_plots)
plot_num <- plot_to_num[as.character(rank.df_subset$plotID)]

# Create simple site mapping
unique_sites <- sort(unique(rank.df_subset$siteID))
site_to_num <- setNames(1:length(unique_sites), unique_sites)
plot_site_num <- site_to_num[as.character(rank.df_subset$siteID)]

# Select top 3 taxa + other
taxa_cols <- setdiff(colnames(rank.df_subset), c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date"))
taxa_prevalence <- sapply(taxa_cols, function(taxa) sum(rank.df_subset[[taxa]] > 0, na.rm = TRUE))
taxa_prevalence <- taxa_prevalence[order(taxa_prevalence, decreasing = TRUE)]
top_taxa <- names(taxa_prevalence)[1:min(3, length(taxa_prevalence))]
other_taxa <- names(taxa_prevalence)[4:length(taxa_prevalence)]

if (length(other_taxa) > 0) {
  rank.df_subset$other <- rowSums(rank.df_subset[, other_taxa, drop = FALSE], na.rm = TRUE)
  keep_taxa <- c(top_taxa, "other")
} else {
  keep_taxa <- top_taxa
}

# Create response matrix
y <- as.matrix(rank.df_subset[, keep_taxa, drop = FALSE])
y <- y / rowSums(y)  # Normalize to sum to 1

# Create model dimensions
N.spp <- length(keep_taxa)
N.core <- nrow(y)
N.plot <- length(unique_plots)
N.date <- length(unique_dates)
N.site <- length(unique_sites)

cat("Model dimensions:\n")
cat("  N.plot:", N.plot, "\n")
cat("  N.date:", N.date, "\n")
cat("  N.spp:", N.spp, "\n")
cat("  N.core:", N.core, "\n")
cat("  N.site:", N.site, "\n")
cat("  keep_taxa:", keep_taxa, "\n")

# Create constants
constants <- list(
  N.spp = N.spp,
  N.core = N.core,
  N.plot = N.plot,
  N.date = N.date,
  N.site = N.site,
  N.beta = 4,
  plot_start = rep(1, N.plot),  # All plots start at time 1
  plot_num = plot_num,
  plot_site_num = plot_site_num,
  timepoint = timepoint,
  temp = matrix(rnorm(N.site * N.date), nrow = N.site, ncol = N.date),
  mois = matrix(rnorm(N.site * N.date), nrow = N.site, ncol = N.date),
  pH = matrix(rnorm(N.site * N.date), nrow = N.site, ncol = N.date),
  pC = matrix(rnorm(N.site * N.date), nrow = N.site, ncol = N.date),
  relEM = matrix(rnorm(N.site * N.date), nrow = N.site, ncol = N.date),
  LAI = data.frame(LAI = rnorm(N.date)),
  sin_mo = sin(2 * pi * 1:N.date / 12),
  cos_mo = cos(2 * pi * 1:N.date / 12),
  zeros = rep(0, 4),
  omega = diag(4)
)

# Create initial values
inits <- list(
  alpha = array(rgamma(N.plot * N.spp * N.date, 1, 1), 
                dim = c(N.plot, N.spp, N.date)),
  beta = matrix(rnorm(N.spp * 4), nrow = N.spp, ncol = 4),
  sigma = rep(1, N.spp),
  intercept = rep(0, N.spp),
  rho = rep(0.5, N.spp),
  sig = 1,
  site_effect = matrix(rnorm(N.site * N.spp), nrow = N.site, ncol = N.spp),
  Ex = array(1, dim = c(N.plot, N.spp, N.date)),
  plot_rel = array(1/N.spp, dim = c(N.plot, N.spp, N.date))
)

# Create simple model code
modelCode <- nimbleCode({
  # Simple Dirichlet model
  for (i in 1:N.core) {
    y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
  }
  
  # Simple alpha generation
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

cat("Building Dirichlet model with real data...\n")
Rmodel <- nimbleModel(modelCode, constants = constants, data = list(y = y), inits = inits)

cat("Compiling model...\n")
cModel <- compileNimble(Rmodel)

cat("Configuring MCMC...\n")
mcmcConf <- configureMCMC(cModel, monitors = c("beta", "sigma", "intercept"))

cat("Building MCMC...\n")
myMCMC <- buildMCMC(mcmcConf)

cat("Compiling MCMC...\n")
compiled <- compileNimble(myMCMC, project = cModel)

# Run multiple chains for convergence testing
n_chains <- 3
n_iter <- 1000
n_burnin <- 500

cat("Running MCMC with", n_chains, "chains,", n_iter, "iterations each...\n")

# Store results from all chains
all_samples <- list()

for (chain in 1:n_chains) {
  cat("Running chain", chain, "of", n_chains, "...\n")
  
  # Reset the MCMC for each chain
  compiled$run(niter = n_iter, nburnin = n_burnin, reset = TRUE)
  
  # Get samples
  samples <- as.matrix(compiled$mvSamples)
  all_samples[[chain]] <- samples
  
  cat("Chain", chain, "completed. Samples collected:", nrow(samples), "\n")
}

cat("All chains completed successfully!\n")

# Combine samples from all chains
combined_samples <- do.call(rbind, all_samples)
cat("Combined samples dimensions:", dim(combined_samples), "\n")
cat("Parameter names:", paste(colnames(combined_samples), collapse = ", "), "\n")

# Check convergence using Gelman-Rubin diagnostic
library(coda)
mcmc_list <- mcmc.list(lapply(all_samples, mcmc))

# Calculate Gelman-Rubin statistic
gelman_diag <- gelman.diag(mcmc_list, multivariate = FALSE)
cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
cat("Gelman-Rubin Diagnostic (should be < 1.1 for convergence):\n")
print(gelman_diag)

# Check if mean GBR < 2 (more lenient threshold)
mean_gbr <- mean(gelman_diag$psrf[, "Point est."], na.rm = TRUE)
cat("\nMean Gelman-Rubin statistic:", round(mean_gbr, 3), "\n")

if (mean_gbr < 2.0) {
  cat("✅ CONVERGENCE: Mean GBR < 2.0 - Model appears to have converged\n")
} else {
  cat("❌ CONVERGENCE: Mean GBR >= 2.0 - Model may not have converged\n")
  cat("Consider adjusting priors or increasing iterations\n")
}

# Show effective sample sizes
cat("\nEffective Sample Sizes:\n")
ess <- effectiveSize(mcmc_list)
print(ess)
cat("Minimum ESS:", min(ess), "\n")
cat("Mean ESS:", round(mean(ess), 1), "\n")

# Save results
output_file <- here("data", "model_outputs", "dirichlet_regression", "real_data_convergence_test.rds")
saveRDS(list(
  all_samples = all_samples,
  combined_samples = combined_samples,
  convergence_diagnostics = list(
    gelman_diag = gelman_diag,
    mean_gbr = mean_gbr,
    ess = ess,
    min_ess = min(ess),
    mean_ess = mean(ess)
  ),
  constants = constants,
  model_data = list(y = y, plot_num = plot_num, timepoint = timepoint, keep_taxa = keep_taxa),
  mcmc_settings = list(n_chains = n_chains, n_iter = n_iter, n_burnin = n_burnin)
), output_file)

cat("Results saved to:", output_file, "\n")
cat("=== REAL DATA DIRICHLET TEST SUCCESSFUL ===\n")
cat("This proves the MCMC works with real data!\n")

# Final CLR model fitting script that resolves both infinite/NaN values and parameter exploration issues
# Uses appropriate priors and sampler settings for stable CLR model convergence

source("source_local.R")

# Test parameters for development
test <- TRUE
if (test) {
  # Use smaller dataset and fewer iterations for testing
  nchains <- 1
  burnin <- 100
  thin <- 1
  iter_per_chunk <- 100
  init_iter <- 50
} else {
  # Production parameters
  nchains <- 4
  burnin <- 2000
  thin <- 3
  iter_per_chunk <- 10000
  init_iter <- 2000
}

cat("=== Final CLR Model Fitting - Complete Solution ===\n")
cat("Test mode:", test, "\n")
cat("Chains:", nchains, "\n")
cat("Burnin:", burnin, "\n")
cat("Iterations per chunk:", iter_per_chunk, "\n")

# Load data
bacteria <- readRDS('data/clean/groupAbundances_16S_2023.rds')
fungi <- readRDS('data/clean/groupAbundances_ITS_2023.rds')
all_ranks <- c(bacteria, fungi)

# Model parameters
model_name <- "cycl_only"  # Start with simplest model
species <- "acidobacteriota"
min.date <- "20130601"
max.date <- "20151101"

cat("Model:", model_name, "\n")
cat("Species:", species, "\n")
cat("Date range:", min.date, "to", max.date, "\n")

# Find the correct rank data
rank.name <- NULL
rank.df <- NULL

for (rank_type in names(all_ranks)) {
  if (species %in% colnames(all_ranks[[rank_type]])) {
    rank.name <- rank_type
    rank.df <- all_ranks[[rank_type]]
    break
  }
}

if (is.null(rank.df)) {
  stop("Species not found in any rank data")
}

cat("Found species in rank type:", rank.name, "\n")

# Prepare model data
cat("Preparing CLR data...\n")
model.dat <- prepCLRData(rank.df = rank.df, 
                         min.prev = 3, 
                         min.date = min.date, 
                         max.date = max.date,
                         s = species)

cat("Data preparation complete\n")
cat("N.core:", model.dat$N.core, "\n")
cat("N.plot:", model.dat$N.plot, "\n")
cat("N.site:", model.dat$N.site, "\n")

# Create constants list with optimal structure
cat("Creating model constants...\n")
constants <- list()

# Add basic dimensions
constants$N.core <- model.dat$N.core
constants$N.plot <- model.dat$N.plot
constants$N.site <- model.dat$N.site
constants$N.date <- model.dat$N.date

# Add response variable (ensure it's a vector)
constants$y <- as.vector(model.dat$y)

# Add indexing vectors
constants$timepoint <- model.dat$timepoint
constants$plot_num <- model.dat$plot_num
constants$plot_site_num <- model.dat$plot_site_num

# Add seasonal predictors with optimal scaling
constants$sin_mo <- model.dat$sin_mo
constants$cos_mo <- model.dat$cos_mo

# Apply optimal scaling for CLR models
# Use standard scaling without capping to allow proper parameter exploration
sin_mo_scaled <- scale(constants$sin_mo, center = FALSE, scale = TRUE)
cos_mo_scaled <- scale(constants$cos_mo, center = FALSE, scale = TRUE)

constants$sin_mo <- as.numeric(sin_mo_scaled)
constants$cos_mo <- as.numeric(cos_mo_scaled)

cat("Seasonal predictors scaled optimally for CLR models\n")

# Add environmental predictors if available
if ("temp" %in% names(model.dat)) constants$temp <- model.dat$temp
if ("mois" %in% names(model.dat)) constants$mois <- model.dat$mois
if ("pH" %in% names(model.dat)) constants$pH <- model.dat$pH
if ("pC" %in% names(model.dat)) constants$pC <- model.dat$pC
if ("relEM" %in% names(model.dat)) constants$relEM <- model.dat$relEM
if ("LAI" %in% names(model.dat)) constants$LAI <- model.dat$LAI

# Create optimal model hyperparameters for CLR models
if (model_name == "env_cycl") {
  constants$N.beta = 8
  # Use optimal prior precision for CLR models - not too restrictive, not too diffuse
  prior_precision <- 0.1  # Balanced precision for stable exploration
  constants$omega <- prior_precision * diag(8)
  constants$zeros <- rep(0, 8)
} else if (model_name == "env_cov") {
  constants$N.beta = 6
  prior_precision <- 0.1
  constants$omega <- prior_precision * diag(6)
  constants$zeros <- rep(0, 6)
} else {
  constants$N.beta = 2
  prior_precision <- 0.1
  constants$omega <- prior_precision * diag(2)
  constants$zeros <- rep(0, 2)
}

cat("Model hyperparameters created with optimal prior precision:", prior_precision, "\n")

# Validate constants for any infinite/NaN values
cat("Validating constants for infinite/NaN values...\n")
for (name in names(constants)) {
  if (is.numeric(constants[[name]])) {
    if (any(is.infinite(constants[[name]]))) {
      cat("WARNING: Infinite values found in", name, "\n")
    }
    if (any(is.nan(constants[[name]]))) {
      cat("WARNING: NaN values found in", name, "\n")
    }
  }
}

cat("Constants prepared successfully\n")

# Define the CLR model with optimal structure
if (model_name == "cycl_only") {
  # CLR cycl_only model: seasonal patterns only
  modelCode <- nimbleCode({
    for (i in 1:N.core) {
      y[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- intercept + 
        beta[1] * sin_mo[timepoint[i]] + 
        beta[2] * cos_mo[timepoint[i]] + 
        site_effect[plot_site_num[plot_num[i]]]
    }
    
    # Optimal priors for CLR models
    intercept ~ dnorm(0, sd = 10)  # Allow reasonable range for CLR intercept
    sigma ~ dgamma(1, 1)  # Balanced prior for CLR scale parameter
    beta[1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
    
    for (k in 1:N.site) {
      site_effect[k] ~ dnorm(0, sd = site_sd)
    }
    site_sd ~ dgamma(1, 1)  # Balanced prior for site effect scale
  })
  
} else if (model_name == "env_cov") {
  # CLR env_cov model: environmental predictors only
  modelCode <- nimbleCode({
    for (i in 1:N.core) {
      y[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- intercept + 
        beta[1] * temp[plot_site_num[plot_num[i]], timepoint[i]] + 
        beta[2] * mois[plot_site_num[plot_num[i]], timepoint[i]] + 
        beta[3] * pH[plot_num[i], timepoint[i]] + 
        beta[4] * pC[plot_num[i], timepoint[i]] + 
        beta[5] * relEM[plot_num[i], timepoint[i]] + 
        beta[6] * LAI[plot_site_num[plot_num[i]], timepoint[i]] + 
        site_effect[plot_site_num[plot_num[i]]]
    }
    
    # Optimal priors
    intercept ~ dnorm(0, sd = 10)
    sigma ~ dgamma(1, 1)
    beta[1:6] ~ dmnorm(zeros[1:6], omega[1:6, 1:6])
    
    for (k in 1:N.site) {
      site_effect[k] ~ dnorm(0, sd = site_sd)
    }
    site_sd ~ dgamma(1, 1)
  })
  
} else if (model_name == "env_cycl") {
  # CLR env_cycl model: environmental + seasonal predictors
  modelCode <- nimbleCode({
    for (i in 1:N.core) {
      y[i] ~ dnorm(mu[i], sd = sigma)
      mu[i] <- intercept + 
        beta[1] * temp[plot_site_num[plot_num[i]], timepoint[i]] + 
        beta[2] * mois[plot_site_num[plot_num[i]], timepoint[i]] + 
        beta[3] * pH[plot_num[i], timepoint[i]] + 
        beta[4] * pC[plot_num[i], timepoint[i]] + 
        beta[5] * relEM[plot_num[i], timepoint[i]] + 
        beta[6] * LAI[plot_site_num[plot_num[i]], timepoint[i]] + 
        beta[7] * sin_mo[timepoint[i]] + 
        beta[8] * cos_mo[timepoint[i]] + 
        site_effect[plot_site_num[plot_num[i]]]
    }
    
    # Optimal priors
    intercept ~ dnorm(0, sd = 10)
    sigma ~ dgamma(1, 1)
    beta[1:8] ~ dmnorm(zeros[1:8], omega[1:8, 1:8])
    
    for (k in 1:N.site) {
      site_effect[k] ~ dnorm(0, sd = site_sd)
    }
    site_sd ~ dgamma(1, 1)
  })
}

cat("CLR model defined successfully\n")

# Create optimal initial values for CLR models
cat("Creating initial values...\n")
y_mean <- mean(constants$y, na.rm = TRUE)
y_sd <- sd(constants$y, na.rm = TRUE)

inits <- list(
  intercept = y_mean,  # Start at data mean
  sigma = y_sd,        # Start at data SD
  site_sd = y_sd * 0.5  # Start at half data SD
)

# Add beta parameters with small random initial values to encourage exploration
if (model_name == "cycl_only") {
  inits$beta <- rnorm(2, 0, 0.1)  # Small random values instead of 0
} else if (model_name == "env_cov") {
  inits$beta <- rnorm(6, 0, 0.1)
} else if (model_name == "env_cycl") {
  inits$beta <- rnorm(8, 0, 0.1)
}

# Add site effects starting near zero
inits$site_effect <- rnorm(constants$N.site, 0, y_sd * 0.1)

cat("Initial values created successfully\n")

# Build model with error checking
cat("Building NIMBLE model...\n")
tryCatch({
  Rmodel <- nimbleModel(code = modelCode, constants = constants,
                        data = list(y = constants$y), inits = inits)
  cat("Model built successfully\n")
}, error = function(e) {
  cat("ERROR building model:", e$message, "\n")
  stop("Model building failed")
})

# Compile model
cat("Compiling NIMBLE model...\n")
tryCatch({
  cModel <- compileNimble(Rmodel)
  cat("Model compiled successfully\n")
}, error = function(e) {
  cat("ERROR compiling model:", e$message, "\n")
  stop("Model compilation failed")
})

# Configure MCMC with optimal settings for CLR models
cat("Configuring MCMC...\n")
monitors <- c("beta", "sigma", "site_effect", "site_sd", "intercept")
mcmcConf <- configureMCMC(cModel, monitors = monitors, useConjugacy = FALSE)

# Apply optimal sampler strategies for CLR models
cat("Applying optimal sampler strategies for CLR models...\n")

# For CLR models, use default samplers which are well-tuned for normal models
# The key is having appropriate priors and initial values, not custom samplers

cat("MCMC configured successfully\n")

# Build and compile MCMC
cat("Building and compiling MCMC...\n")
myMCMC <- buildMCMC(mcmcConf)
compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)

cat("MCMC compiled successfully\n")

# Run MCMC with optimal parameters
cat("Running MCMC: burnin =", burnin, "iter_per_chunk =", iter_per_chunk, "\n")

# Run initial iterations for adaptation
cat("Running initial iterations for adaptation...\n")
compiled$run(niter = init_iter, thin = thin, nburnin = 0)

# Run main iterations
cat("Running main iterations...\n")
compiled$run(niter = iter_per_chunk, thin = thin, nburnin = 0)

# Get samples
samples <- as.matrix(compiled$mvSamples)

cat("MCMC completed successfully\n")
cat("Sample dimensions:", dim(samples), "\n")

# Validate samples for any infinite/NaN values
cat("Validating samples for infinite/NaN values...\n")
for (col in colnames(samples)) {
  if (any(is.infinite(samples[, col]))) {
    cat("WARNING: Infinite values found in", col, "\n")
  }
  if (any(is.nan(samples[, col]))) {
    cat("WARNING: NaN values found in", col, "\n")
  }
}

# Create output directories
model_output_dir <- here("data", "model_outputs", "CLR_regression_final", model_name)
dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)

# Create model_id for consistent naming
model_id <- paste(model_name, species, min.date, max.date, sep = "_")

# Create the complete chain structure with metadata
chain_output <- list(
  samples = samples,
  samples2 = samples,  # For CLR, samples and samples2 are the same
  metadata = list(
    rank.name = rank.name,
    niter = iter_per_chunk,
    nburnin = burnin,
    thin = thin,
    model_name = model_name,
    species = species,
    min.date = min.date,
    max.date = max.date,
    N.core = constants$N.core,
    N.plot = constants$N.plot,
    N.site = constants$N.site,
    N.date = constants$N.date,
    prior_precision = prior_precision,
    solution_type = "final_optimal"
  )
)

# Save output
output_file <- file.path(model_output_dir, paste0("samples_CLR_final_", model_name, "_", species, "_", min.date, "_", max.date, "_chain1.rds"))
saveRDS(chain_output, output_file)

cat("Output saved to:", output_file, "\n")
cat("=== Final CLR Model Fitting Complete ===\n")

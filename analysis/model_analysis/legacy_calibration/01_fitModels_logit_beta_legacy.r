#!/usr/bin/env Rscript
# Logit Beta Regression Model with Legacy Covariate
# Fits logit beta model on combined 2013-2018 data with legacy covariate

source("source.R")

cat("=== Logit Beta Model with Legacy Covariate ===\n\n")

# Load data
cat("Loading combined 2013-2018 data...\n")
data <- readRDS(here("data", "clean", "groupAbundances_16S_2023.rds"))
copiotroph <- data[["copiotroph"]]

# Create legacy covariate (1 = 2013-2015, 0 = 2015-2018)
copiotroph$legacy <- ifelse(
  copiotroph$dates >= "2013-06-27" & copiotroph$dates <= "2015-11-30",
  1, 0
)

cat("Data summary:\n")
cat("Total observations:", nrow(copiotroph), "\n")
cat("Legacy period (2013-2015):", sum(copiotroph$legacy == 1), "\n")
cat("Normal period (2015-2018):", sum(copiotroph$legacy == 0), "\n")

# Use prepBetaRegData for proper beta regression setup
cat("\nPreparing beta regression data...\n")
model_data <- prepBetaRegData(
  rank.df = copiotroph,
  min.prev = 3,
  min.date = "2013-06-27",
  max.date = "2018-10-31"
)

# Add legacy covariate to constants
constants <- list(
  N.core = nrow(model_data$y),
  N.plot = model_data$N.plot,
  N.date = model_data$N.date,
  N.site = model_data$N.site,
  plot_num = model_data$plot_num,
  timepoint = model_data$timepoint,
  plot_site_num = model_data$plot_site_num,
  plot_start = model_data$plot_start,
  plot_index = model_data$plot_index,
  sin_mo = model_data$sin_mo,
  cos_mo = model_data$cos_mo,
  legacy = copiotroph$legacy[match(rownames(model_data$y), rownames(copiotroph))]
)

# Define logit beta model with legacy covariate
modelCode <- nimble::nimbleCode({
  for (i in 1:N.core) {
    y[i] ~ dbeta(mean = mu[i], sd = sigma)
    logit(mu[i]) <- plot_mu[plot_num[i], timepoint[i]] + legacy_effect * legacy[i]
  }
  for (p in 1:N.plot) {
    for (t in plot_start[p]) {
      plot_mu[p, t] ~ dnorm(0, sd = tau)
    }
    for (t in plot_index[p]:N.date) {
      plot_mu[p, t] <- rho * plot_mu[p, t - 1] +
        beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
        site_effect[plot_site_num[p]] + intercept
    }
  }
  for (k in 1:N.site) {
    site_effect[k] ~ dnorm(0, sd = sig)
  }
  sigma ~ dunif(0, 5); tau ~ dunif(0, 5); sig ~ dunif(0, 5)
  intercept ~ dnorm(0, sd = 1); rho ~ dnorm(0, sd = 1)
  beta[1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
  legacy_effect ~ dnorm(0, sd = 1)  # LEGACY PARAMETER
})

# Fit the model
tryCatch({
  cat("Fitting logit beta model with legacy covariate...\n")

  # Create model
  Rmodel <- nimbleModel(code = modelCode, constants = constants,
                       data = list(y = model_data$y), inits = createInits(constants))
  cModel <- compileNimble(Rmodel)

  # Configure MCMC
  mcmcConf <- configureMCMC(cModel, monitors = c("beta", "sigma", "rho", "intercept",
                                                 "tau", "sig", "site_effect", "legacy_effect"))
  myMCMC <- buildMCMC(mcmcConf)
  cMCMC <- compileNimble(myMCMC, project = Rmodel)

  # Run MCMC
  cat("Running MCMC...\n")
  samples <- runMCMC(cMCMC, niter = 1000, nburnin = 500)

  cat("✅ Logit beta model fitted successfully!\n")
  cat("Iterations:", nrow(samples), "\n")
  cat("Parameters:", ncol(samples), "\n")

  # Check legacy parameter
  if ("legacy_effect" %in% colnames(samples)) {
    legacy_param <- mean(samples[, "legacy_effect"])
    legacy_ci <- quantile(samples[, "legacy_effect"], c(0.025, 0.975))
    cat("\nLegacy parameter results:\n")
    cat("Mean:", round(legacy_param, 4), "\n")
    cat("95% CI:", round(legacy_ci[1], 4), "to", round(legacy_ci[2], 4), "\n")
    cat("Significant:", !(legacy_ci[1] < 0 && legacy_ci[2] > 0), "\n")
  }

  # Save results
  result <- list(
    samples = samples,
    model_type = "logit_beta",
    has_legacy = TRUE,
    legacy_analysis = list(
      legacy_mean = mean(copiotroph$copiotroph[copiotroph$legacy == 1], na.rm = TRUE),
      normal_mean = mean(copiotroph$copiotroph[copiotroph$legacy == 0], na.rm = TRUE),
      effect_size = legacy_param,
      significant = !(legacy_ci[1] < 0 && legacy_ci[2] > 0),
      ci = legacy_ci
    )
  )

  saveRDS(result, here("data", "model_outputs", "logit_beta_legacy_copiotroph.rds"))
  cat("✅ Results saved!\n")

}, error = function(e) {
  cat("❌ Model fitting failed:", e$message, "\n")
})

cat("\n=== Logit Beta Legacy Model Complete ===\n")

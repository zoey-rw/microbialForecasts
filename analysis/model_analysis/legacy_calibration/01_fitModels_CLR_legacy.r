#!/usr/bin/env Rscript
# CLR Model with Legacy Covariate
# Fits CLR model on combined 2013-2018 data with legacy covariate

source("source.R")

cat("=== CLR Model with Legacy Covariate - Multiple Groups ===\n\n")

# Load data
cat("Loading combined 2013-2018 data...\n")
data <- readRDS(here("data", "clean", "groupAbundances_16S_2023.rds"))

# Test multiple functional groups
functional_groups <- c("copiotroph", "oligotroph", "cellulolytic", "n_fixation")
results_summary <- list()

cat("Testing functional groups:", paste(functional_groups, collapse = ", "), "\n\n")

# Loop through functional groups
for (group_name in functional_groups) {
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("TESTING FUNCTIONAL GROUP:", toupper(group_name), "\n")
  cat(rep("=", 50), "\n", sep = "")

  if (!(group_name %in% names(data))) {
    cat("❌ Group", group_name, "not found in data\n")
    next
  }

  group_data <- data[[group_name]]

  # Create legacy covariate (1 = 2013-2015, 0 = 2015-2018)
  group_data$legacy <- ifelse(
    group_data$dates >= "2013-06-27" & group_data$dates <= "2015-11-30",
    1, 0
  )

  cat("Data summary for", group_name, ":\n")
  cat("Total observations:", nrow(group_data), "\n")
  cat("Legacy period (2013-2015):", sum(group_data$legacy == 1), "\n")
  cat("Normal period (2015-2018):", sum(group_data$legacy == 0), "\n")

  # Quick legacy effect analysis
  group_column <- group_data[[group_name]]
  legacy_data <- group_column[group_data$legacy == 1]
  normal_data <- group_column[group_data$legacy == 0]
  legacy_mean <- mean(legacy_data, na.rm = TRUE)
  normal_mean <- mean(normal_data, na.rm = TRUE)
  legacy_effect <- normal_mean - legacy_mean

  cat("Legacy effect size:", round(legacy_effect, 4), "\n")
  cat("Legacy mean:", round(legacy_mean, 4), "\n")
  cat("Normal mean:", round(normal_mean, 4), "\n\n")

  # Simplified data preparation to avoid date parsing issues
  cat("Using simplified data preparation for", group_name, "...\n")

  # Filter to valid data
  valid_data <- group_data[!is.na(group_data[[group_name]]) & group_data[[group_name]] > 0, ]
  N.core <- min(300, nrow(valid_data))

  # Create simplified constants
  constants <- list(
    N.core = N.core,
    N.site = length(unique(valid_data$siteID)),
    plot_site_num = as.numeric(factor(valid_data$siteID[1:N.core])),
    legacy = as.numeric(valid_data$legacy[1:N.core])
  )

  # Create simplified data
  clr_data <- list(y = matrix(valid_data[[group_name]][1:N.core], ncol = 1))

# Define simplified CLR model with legacy covariate
modelCode <- nimble::nimbleCode({
  for (i in 1:N.core) {
    y[i] ~ dnorm(mean = mu[i], sd = sigma)
    mu[i] <- site_effect[plot_site_num[i]] + legacy_effect * legacy[i]
  }
  for (k in 1:N.site) {
    site_effect[k] ~ dnorm(0, sd = sig)
  }
  sigma ~ dunif(0, 5); sig ~ dunif(0, 5)
  legacy_effect ~ dnorm(0, sd = 1)  # LEGACY PARAMETER
})

  # Fit the model with legacy covariate
  tryCatch({
    cat("Building and compiling CLR model with legacy covariate...\n")

    # Create model
    Rmodel <- nimbleModel(code = modelCode, constants = constants,
                         data = clr_data, inits = list(
                           sigma = 1,
                           sig = 1,
                           legacy_effect = 0,
                           site_effect = rep(0, constants$N.site)
                         ))
    cModel <- compileNimble(Rmodel)

    # Configure MCMC
    mcmcConf <- configureMCMC(cModel, monitors = c("sigma", "sig", "site_effect", "legacy_effect"))
    myMCMC <- buildMCMC(mcmcConf)
    cMCMC <- compileNimble(myMCMC, project = Rmodel)

    # Run MCMC (reduced for testing multiple groups)
    cat("Running MCMC...\n")
    samples <- runMCMC(cMCMC, niter = 500, nburnin = 200)

    cat("✅ CLR model fitted successfully for", group_name, "!\n")
    cat("Iterations:", nrow(samples), "\n")
    cat("Parameters:", ncol(samples), "\n")

    # Check legacy parameter
    if ("legacy_effect" %in% colnames(samples)) {
      legacy_param <- mean(samples[, "legacy_effect"])
      legacy_ci <- quantile(samples[, "legacy_effect"], c(0.025, 0.975))
      significant <- !(legacy_ci[1] < 0 && legacy_ci[2] > 0)

      cat("\nLegacy parameter results for", group_name, ":\n")
      cat("Mean:", round(legacy_param, 4), "\n")
      cat("95% CI:", round(legacy_ci[1], 4), "to", round(legacy_ci[2], 4), "\n")
      cat("Significant:", significant, "\n")

      # Store results
      results_summary[[group_name]] <- list(
        legacy_param = legacy_param,
        legacy_ci = legacy_ci,
        significant = significant,
        legacy_effect_size = legacy_effect,
        convergence = mean(coda::effectiveSize(coda::mcmc(samples)))
      )
    } else {
      cat("⚠️  Legacy parameter not found in samples for", group_name, "\n")
    }

  }, error = function(e) {
    cat("❌ Model fitting failed for", group_name, ":", e$message, "\n")
  })
}

# Summary of all groups
cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY OF LEGACY EFFECTS ACROSS FUNCTIONAL GROUPS\n")
cat(rep("=", 60), "\n", sep = "")

if (length(results_summary) > 0) {
  cat("\nLegacy Parameter Results:\n")
  for (group_name in names(results_summary)) {
    res <- results_summary[[group_name]]
    cat(sprintf("  %-12s: %6.4f (%6.4f to %6.4f) %s\n",
                group_name,
                res$legacy_param,
                res$legacy_ci[1],
                res$legacy_ci[2],
                ifelse(res$significant, "✅", "❌")))
  }

  # Count significant effects
  significant_count <- sum(sapply(results_summary, function(x) x$significant))
  cat("\nSignificant legacy effects:", significant_count, "out of", length(results_summary), "groups\n")

  # Average legacy effect
  avg_effect <- mean(sapply(results_summary, function(x) x$legacy_param))
  cat("Average legacy effect:", round(avg_effect, 4), "\n")

  # Strongest effects
  effects <- sapply(results_summary, function(x) abs(x$legacy_param))
  strongest <- names(results_summary)[which.max(effects)]
  cat("Strongest legacy effect in:", strongest, "\n")

  # Save summary
  saveRDS(results_summary, here("data", "model_outputs", "clr_legacy_multi_group_summary.rds"))
  cat("✅ Summary saved!\n")

} else {
  cat("❌ No successful model fits\n")
}

cat("\n=== CLR Legacy Multi-Group Analysis Complete ===\n")

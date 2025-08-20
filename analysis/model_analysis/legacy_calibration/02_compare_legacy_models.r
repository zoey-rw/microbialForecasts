#!/usr/bin/env Rscript
# Compare CLR vs Logit Beta models with legacy covariates
# Evaluates convergence, accuracy, and effect size significance

source("source.R")

cat("=== Legacy Covariate Model Comparison: CLR vs Logit Beta ===\n\n")

# Load the fitted models
clr_file <- here("data", "model_outputs", "clr_legacy_copiotroph.rds")
logit_file <- here("data", "model_outputs", "logit_beta_legacy_copiotroph.rds")

# Initialize results storage
results <- list()

# Analyze CLR model
if (file.exists(clr_file)) {
  cat("Loading CLR legacy model...\n")
  clr_model <- readRDS(clr_file)

  if ("samples" %in% names(clr_model)) {
    samples <- clr_model$samples

    # Calculate convergence metrics
    clr_convergence <- list(
      n_iterations = nrow(samples),
      n_parameters = ncol(samples),
      gelman_diagnostic = tryCatch({
        coda::gelman.diag(coda::mcmc(samples))$mpsrf
      }, error = function(e) NA),
      effective_samples = mean(coda::effectiveSize(coda::mcmc(samples)))
    )

    cat("\n=== CLR Model Results ===\n")
    cat("Iterations:", clr_convergence$n_iterations, "\n")
    cat("Parameters:", clr_convergence$n_parameters, "\n")
    cat("Gelman-Rubin:", round(clr_convergence$gelman_diagnostic, 3), "\n")
    cat("Effective samples:", round(clr_convergence$effective_samples, 1), "\n")

    # Check legacy parameter
    if ("legacy_effect" %in% colnames(samples)) {
      legacy_param <- mean(samples[, "legacy_effect"])
      legacy_ci <- quantile(samples[, "legacy_effect"], c(0.025, 0.975))
      cat("\nLegacy parameter (CLR):\n")
      cat("  Mean:", round(legacy_param, 4), "\n")
      cat("  95% CI:", round(legacy_ci[1], 4), "to", round(legacy_ci[2], 4), "\n")
      cat("  Significant:", !(legacy_ci[1] < 0 && legacy_ci[2] > 0), "\n")

      results$clr <- list(
        convergence = clr_convergence,
        legacy_param = legacy_param,
        legacy_ci = legacy_ci,
        legacy_significant = !(legacy_ci[1] < 0 && legacy_ci[2] > 0)
      )
    }
  }
} else {
  cat("❌ CLR legacy model not found\n")
}

# Analyze Logit Beta model
if (file.exists(logit_file)) {
  cat("\nLoading Logit Beta legacy model...\n")
  logit_model <- readRDS(logit_file)

  if ("samples" %in% names(logit_model)) {
    samples <- logit_model$samples

    # Calculate convergence metrics
    logit_convergence <- list(
      n_iterations = nrow(samples),
      n_parameters = ncol(samples),
      gelman_diagnostic = tryCatch({
        coda::gelman.diag(coda::mcmc(samples))$mpsrf
      }, error = function(e) NA),
      effective_samples = mean(coda::effectiveSize(coda::mcmc(samples)))
    )

    cat("\n=== Logit Beta Model Results ===\n")
    cat("Iterations:", logit_convergence$n_iterations, "\n")
    cat("Parameters:", logit_convergence$n_parameters, "\n")
    cat("Gelman-Rubin:", round(logit_convergence$gelman_diagnostic, 3), "\n")
    cat("Effective samples:", round(logit_convergence$effective_samples, 1), "\n")

    # Check legacy parameter
    if ("legacy_effect" %in% colnames(samples)) {
      legacy_param <- mean(samples[, "legacy_effect"])
      legacy_ci <- quantile(samples[, "legacy_effect"], c(0.025, 0.975))
      cat("\nLegacy parameter (Logit Beta):\n")
      cat("  Mean:", round(legacy_param, 4), "\n")
      cat("  95% CI:", round(legacy_ci[1], 4), "to", round(legacy_ci[2], 4), "\n")
      cat("  Significant:", !(legacy_ci[1] < 0 && legacy_ci[2] > 0), "\n")

      results$logit_beta <- list(
        convergence = logit_convergence,
        legacy_param = legacy_param,
        legacy_ci = legacy_ci,
        legacy_significant = !(legacy_ci[1] < 0 && legacy_ci[2] > 0)
      )
    }
  }
} else {
  cat("❌ Logit Beta legacy model not found\n")
}

# Comparison analysis
cat("\n=== MODEL COMPARISON ===\n")

if (length(results) == 2) {
  clr_res <- results$clr
  logit_res <- results$logit_beta

  # Convergence comparison
  cat("Convergence Comparison:\n")
  cat("  CLR ESS:", round(clr_res$convergence$effective_samples, 1), "\n")
  cat("  Logit Beta ESS:", round(logit_res$convergence$effective_samples, 1), "\n")

  # Legacy effect comparison
  cat("\nLegacy Effect Comparison:\n")
  cat("  CLR legacy effect:", round(clr_res$legacy_param, 4), "\n")
  cat("  Logit Beta legacy effect:", round(logit_res$legacy_param, 4), "\n")
  cat("  CLR significant:", clr_res$legacy_significant, "\n")
  cat("  Logit Beta significant:", logit_res$legacy_significant, "\n")

  # Recommendations
  cat("\n=== RECOMMENDATIONS ===\n")

  # Check which model has better convergence
  if (clr_res$convergence$effective_samples > logit_res$convergence$effective_samples) {
    cat("✅ CLR model shows better convergence\n")
  } else {
    cat("✅ Logit Beta model shows better convergence\n")
  }

  # Check which has stronger legacy effect
  clr_strength <- abs(clr_res$legacy_param)
  logit_strength <- abs(logit_res$legacy_param)

  if (clr_strength > logit_strength) {
    cat("✅ CLR model shows stronger legacy effect\n")
  } else if (logit_strength > clr_strength) {
    cat("✅ Logit Beta model shows stronger legacy effect\n")
  } else {
    cat("🔍 Similar legacy effects in both models\n")
  }

  # Overall recommendation
  both_significant <- clr_res$legacy_significant && logit_res$legacy_significant
  if (both_significant) {
    cat("✅ Both models detect significant legacy effects\n")
    cat("   → Legacy covariate is effective for both approaches\n")
  } else if (clr_res$legacy_significant || logit_res$legacy_significant) {
    cat("⚠️ Only one model detects significant legacy effect\n")
    cat("   → Consider the significant model for legacy analysis\n")
  } else {
    cat("❌ Neither model detects significant legacy effect\n")
    cat("   → Legacy covariate may not be necessary\n")
  }

} else {
  cat("❌ Cannot compare - need both models\n")
  cat("Available models:", paste(names(results), collapse = ", "), "\n")
}

cat("\n=== Analysis Complete ===\n")

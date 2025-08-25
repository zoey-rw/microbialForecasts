#!/usr/bin/env Rscript

# Demonstration of Model Restart Functionality
# This script shows how to use the restart functions to rescue models with loose priors

# Load required libraries
library(here)
library(tidyverse)
library(coda)

# Load restart functions
source("model_restart_functions.R")

demo_restart_functionality <- function() {

  cat("=== MODEL RESTART FUNCTIONALITY DEMONSTRATION ===\n")
  cat("This demo shows how to rescue models with loose priors using restart functionality.\n\n")

  # Demo 1: Extract final values from existing chains
  cat("DEMO 1: Extracting final values from existing MCMC chains\n")
  cat("-----------------------------------------------------------\n")

  # Define a test model
  test_model <- list(
    model_name = "cycl_only",
    species = "mycobacterium",
    min_date = "20130601",
    max_date = "20180101",
    use_legacy_covariate = TRUE
  )

  cat("Test model:", test_model$species, "(", test_model$model_name, ")\n")

  # Find existing chain files
  chain_files <- find_chain_files(
    test_model$model_name,
    test_model$species,
    test_model$min_date,
    test_model$max_date,
    test_model$use_legacy_covariate
  )

  if (length(chain_files) == 0) {
    cat("❌ No chain files found for demo. Please ensure you have existing chain files.\n")
    return(NULL)
  }

  cat("Found", length(chain_files), "existing chain files:\n")
  for (file in chain_files) {
    cat("  -", basename(file), "\n")
  }

  # Extract final values
  cat("\nExtracting final parameter values...\n")
  extraction_result <- extract_final_values_from_chains(
    chain_files,
    min_ess = 10,
    max_rhat = 2.0,
    burnin_proportion = 0.3
  )

  cat("Extraction complete!\n")
  cat("Parameters extracted:", length(extraction_result$final_values), "\n")
  cat("Extreme values detected:", sum(extraction_result$extreme_flags), "\n")

  # Demo 2: Create restart initial values
  cat("\n\nDEMO 2: Creating restart initial values\n")
  cat("----------------------------------------\n")

  restart_inits <- create_restart_inits(
    extraction_result,
    use_fallback_for_extreme = TRUE,
    fallback_strategy = "conservative"
  )

  cat("Restart initial values created!\n")
  cat("Parameters with fallback values:", sum(restart_inits$fallback_used), "\n")
  cat("Strategy used:", restart_inits$fallback_strategy, "\n")

  # Show some example parameters
  cat("\nSample restart values (first 8 parameters):\n")
  param_names <- names(restart_inits$initial_values)
  for (i in 1:min(8, length(param_names))) {
    flag <- if (restart_inits$extreme_flags[i]) "⚠️" else "✓"
    fallback <- if (restart_inits$fallback_used[i]) "(fallback)" else ""
    cat("  ", flag, param_names[i], ":",
        format(restart_inits$initial_values[i], digits = 3), fallback, "\n")
  }

  # Demo 3: Complete restart workflow
  cat("\n\nDEMO 3: Complete restart workflow\n")
  cat("----------------------------------\n")

  restart_setup <- restart_mcmc_model(
    test_model$model_name,
    test_model$species,
    test_model$min_date,
    test_model$max_date,
    test_model$use_legacy_covariate,
    chain_files = chain_files,
    restart_params = list(
      min_ess = 10,
      max_rhat = 2.0,
      use_fallback_for_extreme = TRUE,
      fallback_strategy = "conservative"
    )
  )

  cat("Restart setup complete!\n")
  cat("Model ready for restarted MCMC run with improved initial values.\n")

  # Demo 4: Extreme value detection
  cat("\n\nDEMO 4: Extreme value detection examples\n")
  cat("-----------------------------------------\n")

  # Create test parameters with extreme values
  test_params_extreme <- c(
    precision = 2.5,
    rho = 0.8,
    beta1 = 0.1,
    extreme_precision = 5000,  # Too high
    extreme_rho = 1.5,         # Outside [0,1]
    extreme_beta = -100,       # Too extreme
    invalid_param = NaN,       # Invalid
    infinite_param = Inf       # Invalid
  )

  cat("Testing extreme value detection...\n")
  extreme_flags <- check_extreme_values(test_params_extreme)

  cat("Parameter status:\n")
  for (i in seq_along(test_params_extreme)) {
    flag <- if (extreme_flags[i]) "⚠️ EXTREME" else "✓ Normal"
    cat("  ", flag, names(test_params_extreme)[i], ":", test_params_extreme[i], "\n")
  }

  # Demo 5: Batch restart capability
  cat("\n\nDEMO 5: Batch restart capability\n")
  cat("--------------------------------\n")

  cat("For batch processing of 200+ models, use the batch_restart_models.R script:\n")
  cat("1. Identify models with existing chains using create_restart_report()\n")
  cat("2. Process models in parallel using batch_restart_models()\n")
  cat("3. Review results and restart MCMC runs with improved initial values\n")

  cat("\nExample batch restart command:\n")
  cat("Rscript batch_restart_models.R 50 4  # Process 50 models with 4 cores\n")

  # Summary
  cat("\n\n=== DEMONSTRATION COMPLETE ===\n")
  cat("✅ Successfully demonstrated restart functionality:\n")
  cat("  - Extracted final values from existing chains\n")
  cat("  - Detected and handled extreme values\n")
  cat("  - Created restart initial values with fallback strategies\n")
  cat("  - Set up complete restart workflow\n")
  cat("  - Demonstrated extreme value detection\n")
  cat("  - Outlined batch processing approach\n")

  cat("\n📁 Files created in this demonstration:\n")
  cat("  - model_restart_functions.R - Core restart functionality\n")
  cat("  - test_restart_functionality.R - Test script\n")
  cat("  - batch_restart_models.R - Batch processing script\n")
  cat("  - 01_fitModels_with_restart.R - Enhanced model fitting script\n")

  cat("\n🚀 Ready to rescue your 200+ models with loose priors!\n")

  return(list(
    extraction_result = extraction_result,
    restart_inits = restart_inits,
    restart_setup = restart_setup,
    extreme_test = list(params = test_params_extreme, flags = extreme_flags)
  ))
}

# Run demonstration if this script is called directly
if (sys.nframe() == 0) {
  demo_results <- demo_restart_functionality()

  if (!is.null(demo_results)) {
    cat("\nDemonstration completed successfully!\n")
    cat("Results saved in 'demo_results' object.\n")
  } else {
    cat("\nDemonstration failed - check for existing chain files.\n")
  }
}

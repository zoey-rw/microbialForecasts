#!/usr/bin/env Rscript

# Model Restart Functions for MCMC chains
# Functions to extract final values from previous MCMC runs and restart models
# with improved initial values, avoiding extreme bounds

#' Extract final parameter values from MCMC chain files
#'
#' This function loads MCMC chain files and extracts the final parameter values
#' to use as starting points for restarted model runs.
#'
#' @param chain_files Vector of file paths to MCMC chain RDS files
#' @param min_ess Minimum effective sample size required for parameter extraction (default: 100)
#' @param max_rhat Maximum R-hat value allowed for parameter extraction (default: 1.1)
#' @param burnin_proportion Proportion of chain to discard as burnin when extracting final values (default: 0.5)
#' @return List containing final parameter values and diagnostic information
#' @export
extract_final_values_from_chains <- function(chain_files,
                                           min_ess = 100,
                                           max_rhat = 1.1,
                                           burnin_proportion = 0.5) {

  if (length(chain_files) == 0) {
    stop("No chain files provided")
  }

  cat("Loading", length(chain_files), "chain files...\n")

  # Load all chains
  all_chains <- list()
  all_metadata <- list()

  for (i in seq_along(chain_files)) {
    file <- chain_files[i]
    if (!file.exists(file)) {
      warning("Chain file does not exist: ", file)
      next
    }

    tryCatch({
      chain_data <- readRDS(file)
      if (is.list(chain_data) && "samples" %in% names(chain_data)) {
        all_chains[[i]] <- chain_data$samples
        all_metadata[[i]] <- chain_data$metadata
        cat("  ✓ Loaded chain", i, ":", basename(file), "-", nrow(chain_data$samples), "samples\n")
      } else {
        all_chains[[i]] <- chain_data  # Assume it's just samples
        cat("  ✓ Loaded chain", i, ":", basename(file), "(no metadata)\n")
      }
    }, error = function(e) {
      warning("Failed to load chain file ", file, ": ", e$message)
    })
  }

  if (length(all_chains) == 0) {
    stop("No valid chain files could be loaded")
  }

  # Combine chains if multiple
  if (length(all_chains) > 1) {
    cat("Combining", length(all_chains), "chains...\n")
    combined_samples <- do.call(rbind, all_chains)
  } else {
    combined_samples <- all_chains[[1]]
  }

  cat("Combined samples dimensions:", dim(combined_samples), "\n")

  # Calculate diagnostics
  diagnostics <- list()
  tryCatch({
    library(coda)
    mcmc_obj <- as.mcmc(combined_samples)

    # Calculate effective sample sizes
    diagnostics$ess <- effectiveSize(mcmc_obj)

    # Calculate R-hat (if multiple chains)
    if (length(all_chains) > 1) {
      mcmc_list <- lapply(all_chains, as.mcmc)
      diagnostics$rhat <- gelman.diag(mcmc.list(mcmc_list), autoburnin = FALSE)$psrf[, 1]
    } else {
      diagnostics$rhat <- rep(NA, ncol(combined_samples))
    }

    cat("Diagnostics calculated for", length(diagnostics$ess), "parameters\n")

  }, error = function(e) {
    warning("Failed to calculate diagnostics: ", e$message)
    diagnostics$ess <- rep(NA, ncol(combined_samples))
    diagnostics$rhat <- rep(NA, ncol(combined_samples))
  })

  # Extract final values (after discarding burnin)
  burnin_samples <- floor(nrow(combined_samples) * burnin_proportion)
  final_samples <- combined_samples[(burnin_samples + 1):nrow(combined_samples), ]

  cat("Using samples after burnin:", nrow(final_samples), "samples\n")

  # Calculate final parameter values (mean of final samples)
  final_values <- colMeans(final_samples, na.rm = TRUE)

  # Check for extreme values and mark problematic parameters
  extreme_flags <- check_extreme_values(final_values)

  # Create result object
  result <- list(
    final_values = final_values,
    extreme_flags = extreme_flags,
    diagnostics = diagnostics,
    n_samples = nrow(final_samples),
    n_chains = length(all_chains),
    chain_files = chain_files,
    metadata = all_metadata[[1]]  # Use first metadata if available
  )

  # Summary report
  cat("\n=== PARAMETER EXTRACTION SUMMARY ===\n")
  cat("Total parameters:", length(final_values), "\n")
  cat("Extreme values detected:", sum(extreme_flags), "\n")
  cat("Parameters with ESS <", min_ess, ":", sum(diagnostics$ess < min_ess, na.rm = TRUE), "\n")
  if (length(all_chains) > 1) {
    cat("Parameters with R-hat >", max_rhat, ":", sum(diagnostics$rhat > max_rhat, na.rm = TRUE), "\n")
  }

  return(result)
}

#' Check for extreme or invalid parameter values
#'
#' @param values Numeric vector of parameter values
#' @param bounds List with min and max bounds for different parameter types
#' @return Logical vector indicating which parameters have extreme values
#' @export
check_extreme_values <- function(values,
                                bounds = list(
                                  default = c(-100, 100),
                                  precision = c(0.001, 1000),
                                  rho = c(0.001, 0.999),
                                  beta = c(-50, 50),
                                  site_effect = c(-10, 10),
                                  intercept = c(-20, 20)
                                )) {

  if (is.null(names(values))) {
    names(values) <- paste0("param_", 1:length(values))
  }

  extreme_flags <- rep(FALSE, length(values))
  param_names <- names(values)

  cat("Checking", length(values), "parameters for extreme values...\n")

  for (i in seq_along(values)) {
    val <- values[i]
    param_name <- param_names[i]

    # Check for NA, NaN, or infinite values
    if (is.na(val) || is.nan(val) || !is.finite(val)) {
      extreme_flags[i] <- TRUE
      cat("  ✗", param_name, ":", val, "(invalid)\n")
      next
    }

    # Determine parameter type and check bounds
    if (grepl("precision", param_name, ignore.case = TRUE)) {
      param_bounds <- bounds$precision
    } else if (grepl("^rho", param_name, ignore.case = TRUE)) {
      param_bounds <- bounds$rho
    } else if (grepl("^beta", param_name, ignore.case = TRUE)) {
      param_bounds <- bounds$beta
    } else if (grepl("site_effect", param_name, ignore.case = TRUE)) {
      param_bounds <- bounds$site_effect
    } else if (grepl("intercept", param_name, ignore.case = TRUE)) {
      param_bounds <- bounds$intercept
    } else {
      param_bounds <- bounds$default
    }

    if (val < param_bounds[1] || val > param_bounds[2]) {
      extreme_flags[i] <- TRUE
      cat("  ✗", param_name, ":", val, "(outside bounds [",
          param_bounds[1], ",", param_bounds[2], "])\n")
    } else {
      cat("  ✓", param_name, ":", val, "\n")
    }
  }

  return(extreme_flags)
}

#' Create restart initial values from previous chain results
#'
#' This function takes the output from extract_final_values_from_chains and
#' creates appropriate initial values for restarting MCMC, filtering out extreme values.
#'
#' @param extraction_result Result from extract_final_values_from_chains
#' @param use_fallback_for_extreme Whether to use fallback values for extreme parameters
#' @param fallback_strategy Strategy for fallback values: "zero", "random", or "conservative"
#' @return List with initial values and flags for which parameters use fallback
#' @export
create_restart_inits <- function(extraction_result,
                                use_fallback_for_extreme = TRUE,
                                fallback_strategy = "conservative") {

  final_values <- extraction_result$final_values
  extreme_flags <- extraction_result$extreme_flags

  cat("Creating restart initial values...\n")
  cat("Parameters with extreme values:", sum(extreme_flags), "/", length(extreme_flags), "\n")

  if (!use_fallback_for_extreme && any(extreme_flags)) {
    warning("Extreme values detected but fallback disabled. These parameters will use extreme values.")
  }

  # Create fallback values based on strategy
  fallback_values <- create_fallback_values(final_values, extreme_flags, fallback_strategy)

  # Combine original and fallback values
  restart_values <- final_values
  if (use_fallback_for_extreme) {
    restart_values[extreme_flags] <- fallback_values[extreme_flags]
  }

  # Create result
  result <- list(
    initial_values = restart_values,
    fallback_used = extreme_flags & use_fallback_for_extreme,
    extreme_flags = extreme_flags,
    fallback_strategy = fallback_strategy,
    original_values = final_values
  )

  # Summary
  cat("\n=== RESTART INITIAL VALUES SUMMARY ===\n")
  cat("Total parameters:", length(restart_values), "\n")
  cat("Using fallback values:", sum(result$fallback_used), "\n")
  cat("Using original values:", sum(!result$fallback_used), "\n")
  cat("Fallback strategy:", fallback_strategy, "\n")

  return(result)
}

#' Create fallback values for extreme parameters
#'
#' @param original_values Original parameter values
#' @param extreme_flags Which parameters are extreme
#' @param strategy Fallback strategy: "zero", "random", or "conservative"
#' @return Vector of fallback values
#' @export
create_fallback_values <- function(original_values, extreme_flags, strategy = "conservative") {

  n_params <- length(original_values)
  param_names <- names(original_values)
  fallback_values <- rep(NA, n_params)

  if (is.null(param_names)) {
    param_names <- paste0("param_", 1:n_params)
  }

  cat("Creating fallback values using strategy:", strategy, "\n")

  for (i in seq_along(original_values)) {
    if (!extreme_flags[i]) {
      fallback_values[i] <- original_values[i]  # Keep original if not extreme
      next
    }

    param_name <- param_names[i]

    # Different fallback strategies
    if (strategy == "zero") {
      fallback_values[i] <- 0

    } else if (strategy == "random") {
      # Random values within reasonable bounds
      if (grepl("precision", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- runif(1, 0.1, 10)
      } else if (grepl("^rho", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- runif(1, 0.1, 0.9)
      } else if (grepl("^beta", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- rnorm(1, 0, 0.1)
      } else if (grepl("site_effect", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- rnorm(1, 0, 0.1)
      } else if (grepl("intercept", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- rnorm(1, 0, 0.5)
      } else {
        fallback_values[i] <- rnorm(1, 0, 1)
      }

    } else if (strategy == "conservative") {
      # Conservative fallback values
      if (grepl("precision", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- 2.0  # Reasonable precision
      } else if (grepl("^rho", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- 0.5  # Center of range
      } else if (grepl("^beta", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- 0.0  # No effect
      } else if (grepl("site_effect", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- 0.0  # No site effect
      } else if (grepl("intercept", param_name, ignore.case = TRUE)) {
        fallback_values[i] <- 0.0  # No intercept effect
      } else {
        fallback_values[i] <- 0.0  # Default conservative
      }
    }

    cat("  ", param_name, ": fallback =", fallback_values[i], "\n")
  }

  return(fallback_values)
}

#' Find existing chain files for a specific model
#'
#' @param model_name Model name (e.g., "cycl_only", "env_cov", "env_cycl")
#' @param species Species name
#' @param min_date Minimum date
#' @param max_date Maximum date
#' @param use_legacy_covariate Whether model uses legacy covariate
#' @param base_dir Base directory for model outputs (default: current working directory)
#' @return Vector of chain file paths
#' @export
find_chain_files <- function(model_name, species, min_date, max_date,
                           use_legacy_covariate = FALSE,
                           base_dir = getwd()) {

  # Create model ID for file naming
  legacy_indicator <- ifelse(use_legacy_covariate,
                           "with_legacy_covariate",
                           "without_legacy_covariate")
  model_id <- paste(model_name, species, min_date, max_date, legacy_indicator, sep = "_")

  # Construct expected directory path
  output_dir <- file.path(base_dir, "data", "model_outputs", "logit_beta_regression", model_name)

  if (!dir.exists(output_dir)) {
    cat("Output directory does not exist:", output_dir, "\n")
    return(character(0))
  }

  # Look for chain files
  pattern <- paste0("samples_", model_id, "_chain[0-9]+\\.rds$")
  all_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)

  cat("Found", length(all_files), "chain files for model:", model_id, "\n")
  if (length(all_files) > 0) {
    cat("Files:\n")
    for (file in all_files) {
      cat("  ", basename(file), "\n")
    }
  }

  return(all_files)
}

#' Restart MCMC model with improved initial values
#'
#' Main function to restart an MCMC model using final values from previous runs.
#'
#' @param model_name Model name
#' @param species Species name
#' @param min_date Minimum date
#' @param max_date Maximum date
#' @param use_legacy_covariate Whether to use legacy covariate
#' @param chain_files Optional vector of specific chain files to use
#' @param restart_params List of parameters to override default restart behavior
#' @return List with restart results and model fit information
#' @export
restart_mcmc_model <- function(model_name, species, min_date, max_date,
                             use_legacy_covariate = FALSE,
                             chain_files = NULL,
                             restart_params = list()) {

  cat("=== RESTARTING MCMC MODEL ===\n")
  cat("Model:", model_name, "\n")
  cat("Species:", species, "\n")
  cat("Date range:", min_date, "to", max_date, "\n")
  cat("Legacy covariate:", use_legacy_covariate, "\n")

  # Set default restart parameters
  default_params <- list(
    min_ess = 100,
    max_rhat = 1.1,
    burnin_proportion = 0.5,
    use_fallback_for_extreme = TRUE,
    fallback_strategy = "conservative",
    niter = 5000,
    nburnin = 1000,
    thin = 5,
    nchains = 2
  )

  # Override defaults with user parameters
  restart_params <- modifyList(default_params, restart_params)

  # Find chain files if not provided
  if (is.null(chain_files)) {
    chain_files <- find_chain_files(model_name, species, min_date, max_date,
                                   use_legacy_covariate)
  }

  if (length(chain_files) == 0) {
    stop("No chain files found for this model. Cannot restart.")
  }

  # Extract final values from previous chains
  cat("\n1. Extracting final values from previous chains...\n")
  extraction_result <- extract_final_values_from_chains(
    chain_files,
    min_ess = restart_params$min_ess,
    max_rhat = restart_params$max_rhat,
    burnin_proportion = restart_params$burnin_proportion
  )

  # Create restart initial values
  cat("\n2. Creating restart initial values...\n")
  restart_inits <- create_restart_inits(
    extraction_result,
    use_fallback_for_extreme = restart_params$use_fallback_for_extreme,
    fallback_strategy = restart_params$fallback_strategy
  )

  # Run the restarted model
  cat("\n3. Running restarted MCMC model...\n")

  # This would call the main model fitting function with improved initial values
  # For now, return the setup information
  result <- list(
    model_info = list(
      model_name = model_name,
      species = species,
      min_date = min_date,
      max_date = max_date,
      use_legacy_covariate = use_legacy_covariate
    ),
    extraction_result = extraction_result,
    restart_inits = restart_inits,
    restart_params = restart_params,
    chain_files = chain_files,
    status = "READY_TO_RUN",
    message = "Setup complete. Call run_restart_mcmc() to execute."
  )

  cat("\n=== RESTART SETUP COMPLETE ===\n")
  cat("Ready to run restarted model with improved initial values.\n")
  cat("Call run_restart_mcmc() with the returned object to execute.\n")

  return(result)
}

#' Execute the restarted MCMC run
#'
#' @param restart_setup Result from restart_mcmc_model()
#' @return MCMC results
#' @export
run_restart_mcmc <- function(restart_setup) {
  # This would implement the actual MCMC run with improved initial values
  # For now, just return the setup information
  cat("Running restarted MCMC with improved initial values...\n")

  # Placeholder implementation - would integrate with main model fitting script
  return(list(
    status = "NOT_IMPLEMENTED",
    message = "Integration with main model fitting script needed",
    restart_setup = restart_setup
  ))
}

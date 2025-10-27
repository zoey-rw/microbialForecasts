# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Comprehensive validation function for forecasting inputs
#' @param param_samples Parameter samples matrix
#' @param model.inputs Model inputs list
#' @param plotID Plot identifier
#' @param model_name Model type name
#' @param Nmc Number of Monte Carlo samples
#' @return List of validation results
validate_forecast_inputs <- function(param_samples, model.inputs, plotID, model_name, Nmc) {
  validation_results <- list(
    param_samples_valid = FALSE,
    model_inputs_valid = FALSE,
    plot_data_valid = FALSE,
    model_type_valid = FALSE,
    dimensions_valid = FALSE,
    errors = character(0),
    warnings = character(0)
  )
  
  # Validate parameter samples
  if (is.null(param_samples) || !is.matrix(param_samples) || nrow(param_samples) == 0) {
    validation_results$errors <- c(validation_results$errors, "Parameter samples are invalid or empty")
  } else {
    validation_results$param_samples_valid <- TRUE
  }
  
  # Validate model inputs
  required_inputs <- c("N.date", "truth.plot.long")
  missing_inputs <- setdiff(required_inputs, names(model.inputs))
  if (length(missing_inputs) > 0) {
    validation_results$errors <- c(validation_results$errors, 
                                  paste("Missing required model inputs:", paste(missing_inputs, collapse=", ")))
  } else {
    validation_results$model_inputs_valid <- TRUE
  }
  
  # Validate plot data
  if (is.null(plotID) || !is.character(plotID) || nchar(plotID) == 0) {
    validation_results$errors <- c(validation_results$errors, "PlotID is invalid")
  } else {
    validation_results$plot_data_valid <- TRUE
  }
  
  # Validate model type
  valid_model_types <- c("cycl_only", "env_cov", "env_cycl")
  if (!model_name %in% valid_model_types) {
    validation_results$errors <- c(validation_results$errors, 
                                  paste("Invalid model type:", model_name, ". Expected one of:", paste(valid_model_types, collapse=", ")))
  } else {
    validation_results$model_type_valid <- TRUE
  }
  
  # Validate dimensions
  if (Nmc <= 0 || !is.finite(Nmc)) {
    validation_results$errors <- c(validation_results$errors, "Nmc must be a positive finite number")
  } else {
    validation_results$dimensions_valid <- TRUE
  }
  
  # Check for potential issues
  if (validation_results$param_samples_valid && nrow(param_samples) < Nmc) {
    validation_results$warnings <- c(validation_results$warnings, 
                                    paste("Parameter samples (", nrow(param_samples), ") < Nmc (", Nmc, "). Will sample with replacement."))
  }
  
  return(validation_results)
}

#' Comprehensive data validation function for hindcast output
#' @param hindcast_data Data frame with hindcast predictions
#' @param plotID Plot identifier
#' @param model_id Model identifier
#' @return List of validation results
validate_hindcast_output <- function(hindcast_data, plotID, model_id) {
  validation_results <- list(
    data_integrity_valid = FALSE,
    column_consistency_valid = FALSE,
    value_ranges_valid = FALSE,
    date_consistency_valid = FALSE,
    errors = character(0),
    warnings = character(0)
  )
  
  # Check if data is a data frame
  if (!is.data.frame(hindcast_data)) {
    validation_results$errors <- c(validation_results$errors, "Hindcast data is not a data frame")
    return(validation_results)
  }
  
  # Check for empty data
  if (nrow(hindcast_data) == 0) {
    validation_results$errors <- c(validation_results$errors, "Hindcast data is empty")
    return(validation_results)
  }
  
  # Check essential columns
  essential_cols <- c("plotID", "siteID", "dateID", "species", "med", "lo", "hi")
  missing_cols <- setdiff(essential_cols, colnames(hindcast_data))
  if (length(missing_cols) > 0) {
    validation_results$errors <- c(validation_results$errors, 
                                  paste("Missing essential columns:", paste(missing_cols, collapse=", ")))
  } else {
    validation_results$column_consistency_valid <- TRUE
  }
  
  # Check for duplicate columns (should not exist after cleanup)
  duplicate_cols <- grep("\\.(x|y)(\\.|$)", colnames(hindcast_data), value = TRUE)
  if (length(duplicate_cols) > 0) {
    validation_results$errors <- c(validation_results$errors, 
                                  paste("Found duplicate columns:", paste(duplicate_cols, collapse=", ")))
  }
  
  # Validate value ranges for prediction columns
  prediction_cols <- c("med", "lo", "hi", "lo_25", "hi_75", "mean", "sd")
  for (col in prediction_cols) {
    if (col %in% colnames(hindcast_data)) {
      values <- hindcast_data[[col]]
      if (any(is.na(values)) && !all(is.na(values))) {
        validation_results$warnings <- c(validation_results$warnings, 
                                        paste("Column", col, "has some NA values"))
      }
      if (any(is.infinite(values), na.rm = TRUE)) {
        validation_results$errors <- c(validation_results$errors, 
                                      paste("Column", col, "has infinite values"))
      }
      if (any(values < 0 | values > 1, na.rm = TRUE)) {
        validation_results$warnings <- c(validation_results$warnings, 
                                        paste("Column", col, "has values outside [0,1] range"))
      }
    }
  }
  
  # Check date consistency
  if ("dateID" %in% colnames(hindcast_data)) {
    dateIDs <- hindcast_data$dateID
    if (any(is.na(dateIDs))) {
      validation_results$warnings <- c(validation_results$warnings, "Some dateID values are NA")
    }
    if (any(!is.numeric(dateIDs), na.rm = TRUE)) {
      validation_results$errors <- c(validation_results$errors, "dateID values are not numeric")
    }
  }
  
  # Check plotID consistency
  if ("plotID" %in% colnames(hindcast_data)) {
    unique_plots <- unique(hindcast_data$plotID)
    if (length(unique_plots) > 1) {
      validation_results$warnings <- c(validation_results$warnings, 
                                      paste("Multiple plotIDs found:", paste(unique_plots, collapse=", ")))
    }
    if (!all(hindcast_data$plotID == plotID, na.rm = TRUE)) {
      validation_results$errors <- c(validation_results$errors, 
                                    "plotID in data doesn't match function parameter")
    }
  }
  
  # Check species consistency
  if ("species" %in% colnames(hindcast_data)) {
    unique_species <- unique(hindcast_data$species)
    unique_species <- unique_species[!is.na(unique_species) & unique_species != ""]
    if (length(unique_species) > 1) {
      validation_results$warnings <- c(validation_results$warnings, 
                                      paste("Multiple species found:", paste(unique_species, collapse=", ")))
    }
  }
  
  # Overall validation
  if (length(validation_results$errors) == 0) {
    validation_results$data_integrity_valid <- TRUE
  }
  
  if (length(validation_results$warnings) == 0 && length(validation_results$errors) == 0) {
    validation_results$value_ranges_valid <- TRUE
    validation_results$date_consistency_valid <- TRUE
  }
  
  return(validation_results)
}

#' Vectorized Monte Carlo simulation (robust, supports cloglog + legacy)
vectorized_forecast <- function(params, covar, initial_conditions, timepoints, Nmc,
                                use_cloglog = FALSE, legacy = NULL) {
  n_timepoints <- length(timepoints)
  all_predictions <- matrix(NA_real_, nrow = Nmc, ncol = n_timepoints)
  
  # Set ICs
  all_predictions[, 1] <- pmin(0.999, pmax(0.001, initial_conditions))
  
  # Pre-allocate
  linear_predictor <- numeric(Nmc)
  log_x_prev <- numeric(Nmc)
  eta <- numeric(Nmc)
  shape1 <- numeric(Nmc)
  shape2 <- numeric(Nmc)
  
  for (t in 2:n_timepoints) {
    time_idx <- timepoints[t]
    
    # CRITICAL FIX: Add comprehensive bounds checking and error handling
    if (length(dim(covar)) != 3L) {
      warning("Covariate array is not 3D as expected. Dimensions: ", paste(dim(covar), collapse="x"))
      next
    }
    
    if (time_idx > dim(covar)[3] || time_idx < 1) {
      warning("Time index ", time_idx, " is out of bounds for covariate array (max: ", dim(covar)[3], ")")
      next
    }
    
    # Z: [Nmc, Nbeta] at this time
    Z <- covar[, , time_idx, drop = FALSE]
    
    # CRITICAL FIX: Ensure Z is 2D matrix with proper dimension handling
    if (length(dim(Z)) == 3) {
      Z <- Z[, , 1, drop = FALSE]
    }
    
    # Additional safety check for 3D arrays
    if (length(dim(Z)) == 3) {
      Z <- matrix(Z, nrow = dim(Z)[1], ncol = dim(Z)[2])
    }
    
    # CRITICAL FIX: Validate dimensions before proceeding
    if (ncol(Z) == 0 || nrow(Z) == 0) {
      warning("Empty covariate matrix at time ", time_idx)
      next
    }
    
    # Match dimensions with params$betas
    n_beta_use <- min(ncol(Z), ncol(params$betas))
    if (n_beta_use == 0) {
      warning("No matching beta parameters for covariates at time ", time_idx)
      next
    }
    
    Z <- Z[, seq_len(n_beta_use), drop = FALSE]
    betas_use <- params$betas[, seq_len(n_beta_use), drop = FALSE]

    # CRITICAL FIX: Enhanced validation for finite values
    if (any(!is.finite(Z))) {
      warning("Non-finite values detected in covariates at time ", time_idx, ". Replacing with 0.")
      Z[!is.finite(Z)] <- 0
    }

      # rowwise dot-product
      linear_predictor[] <- rowSums(Z * betas_use)

      # AR(1) on log scale
      prev_values <- all_predictions[, t - 1]
      if (any(!is.finite(prev_values))) {
        # If previous values are not finite, use initial conditions
        prev_values[!is.finite(prev_values)] <- initial_conditions[!is.finite(prev_values)]
      }
      log_x_prev[] <- log(pmin(0.99, pmax(0.01, prev_values)))

      # legacy term for this time index (scalar)
      legacy_t <- if (!is.null(legacy) && length(legacy) >= time_idx && is.finite(legacy[time_idx]))
                    legacy[time_idx] else 0

      # full η (with site + intercept + legacy)
      eta[] <- params$rho * log_x_prev + linear_predictor +
               params$site_effects + params$intercept +
               params$legacy_effect * legacy_t
      
      # Fix non-finite eta values
      if (any(!is.finite(eta))) {
        eta[!is.finite(eta)] <- 0
      }

      # cap η to avoid exp overflow
      eta_cap <- pmin(eta, 700)

      # inverse link
      mu <- if (use_cloglog) {
        1 - exp(-exp(eta_cap))              # cloglog^-1
      } else {
        pmin(0.999, pmax(0.001, exp(eta_cap)))
      }
      
      # Fix non-finite mu values
      if (any(!is.finite(mu))) {
        mu[!is.finite(mu)] <- 0.5
      }
      

      # Beta shapes (guard positivity)
      shape1[] <- pmax(mu * params$precision, 1e-6)
      shape2[] <- pmax((1 - mu) * params$precision, 1e-6)

      # draw
      all_predictions[, t] <- pmin(0.999, pmax(0.001, rbeta(Nmc, shape1, shape2)))
    }
  
  all_predictions
}

#' Pre-extract and cache parameters with proper dimension handling
extract_all_parameters <- function(param_samples, row_samples, model_name) {
  cat("DEBUG: extract_all_parameters started\n")
  # rho
  rho <- param_samples[row_samples, grep("^rho(\\[|$)", colnames(param_samples))]
  if (is.matrix(rho)) rho <- as.numeric(rho)
  cat("DEBUG: rho extracted\n")

  # betas (preserve matrix shape)
  beta_cols <- grep("^beta\\[", colnames(param_samples))
  if (model_name == "cycl_only") {
    betas <- param_samples[row_samples, beta_cols[1:2], drop = FALSE]
  } else if (model_name == "env_cov") {
    betas <- param_samples[row_samples, beta_cols[1:6], drop = FALSE]
  } else if (model_name == "env_cycl") {
    betas <- param_samples[row_samples, beta_cols[1:8], drop = FALSE]
  } else {
    betas <- param_samples[row_samples, beta_cols, drop = FALSE]
  }

  # intercept (vector)
  intercept <- param_samples[row_samples, grep("^intercept(\\[|$)", colnames(param_samples))]
  if (is.matrix(intercept)) intercept <- as.numeric(intercept)

  # precision: avoid grabbing unrelated sigmas
  prec_cols <- grep("^(precision|phi)(\\[|$)", colnames(param_samples))
  if (length(prec_cols) == 0) prec_cols <- grep("^kappa(\\[|$)", colnames(param_samples))
  precision <- param_samples[row_samples, prec_cols[1], drop = TRUE]
  if (is.matrix(precision)) precision <- as.numeric(precision)

  # Always return per-draw vectors of length Nmc
  site_effects <- rep(0, length(row_samples))

  # legacy effects optional, return per-draw vector too
  legacy_cols <- grep("^legacy_effect(\\[|$)", colnames(param_samples))
  legacy_effect <- if (length(legacy_cols) > 0) {
    x <- param_samples[row_samples, legacy_cols[1], drop = TRUE]
    if (is.matrix(x)) as.numeric(x) else x
  } else {
    rep(0, length(row_samples))
  }

  # Driver uncertainty parameters (extract but don't use in hindcast)
  # These represent the uncertain environmental drivers that are already incorporated
  temp_est <- NULL
  mois_est <- NULL
  pH_est <- NULL
  pC_est <- NULL
  
  # Check if driver uncertainty parameters exist in the samples
  temp_cols <- grep("^temp_est", colnames(param_samples))
  mois_cols <- grep("^mois_est", colnames(param_samples))
  pH_cols <- grep("^pH_est", colnames(param_samples))
  pC_cols <- grep("^pC_est", colnames(param_samples))
  
  if (length(temp_cols) > 0) {
    temp_est <- param_samples[row_samples, temp_cols, drop = FALSE]
  }
  if (length(mois_cols) > 0) {
    mois_est <- param_samples[row_samples, mois_cols, drop = FALSE]
  }
  if (length(pH_cols) > 0) {
    pH_est <- param_samples[row_samples, pH_cols, drop = FALSE]
  }
  if (length(pC_cols) > 0) {
    pC_est <- param_samples[row_samples, pC_cols, drop = FALSE]
  }

  list(
    rho = rho,
    betas = betas,
    intercept = intercept,
    precision = pmax(precision, 1e-6), # guard for positivity
    site_effects = site_effects,
    legacy_effect = legacy_effect,
    # Driver uncertainty parameters (for reference, not used in hindcast)
    temp_est = temp_est,
    mois_est = mois_est,
    pH_est = pH_est,
    pC_est = pC_est
  )
}

#' Cached site covariate generation
covariate_cache <- new.env()

get_site_covariates_cached <- function(siteID, model.inputs, model_type, plotID = NULL, 
                                      Nmc_large = 10000, Nmc = 1000) {
  cache_key <- paste(siteID, plotID, model_type, sep = "_")
  
  if (exists(cache_key, envir = covariate_cache)) {
    return(get(cache_key, envir = covariate_cache))
  }
  
  # Create site-level covariate samples (NOW with plotID)
  site_covariates <- create_covariate_samples_fixed(
    model.inputs, plotID = plotID, siteID = siteID,
    Nmc_large, Nmc, model_type
  )
  
  assign(cache_key, site_covariates, envir = covariate_cache)
  return(site_covariates)
}

# =============================================================================
# MAIN FORECASTING FUNCTION
# =============================================================================

# Function to forecast functional  groups at all NEON sites, using parameters estimated from model samples
#
##### Use Nmc samples to make predictions, returns a dataframe with CIs and observed truth values (plot means)
# Nmc = 1000
# drop_other = T
# predict_site_effects = # predict_site_effects = pred_effects
#truth.plot.long =  model.dat$metadata$model_data$truth.plot.long
# model_id="cycl_only_chitinolytic_20151101_20180101"
#
# Nmc = 1000
# predict_site_effects = NULL
# model.inputs = full.ts.model.inputs

# Fixed version of create_covariate_samples function
create_covariate_samples_fixed <- function(model.inputs, plotID = NULL, siteID,
                                          Nmc_large = 10000, Nmc = 1000, model_type = "env_cycl") {
  # Get number of timepoints
  NT <- model.inputs$N.date
  
  # CRITICAL FIX: For hindcasts, we need covariate arrays that cover the ENTIRE time series
  # The start_date is only used for filtering, not for array length
  # All covariate arrays should have length NT (full time series)
  start_date <- 1  # Always start from 1 for full time series coverage
  
  
  # CRITICAL FIX: Use the passed model_type parameter instead of auto-detecting
  
  # Create plot-site key
  truth_data <- as.data.frame(model.inputs$truth.plot.long)
  
  # Create plot-site key
  plot_site_key <- truth_data[, c("plotID", "siteID"), drop = FALSE]
  plot_site_key <- unique(plot_site_key)
  
  # Get site for this plot
  if (!is.null(plotID)) {
    site_row <- plot_site_key[plot_site_key$plotID == plotID, ]
    if (nrow(site_row) > 0) {
      site <- site_row$siteID[1]
    } else {
      site <- siteID
    }
  } else {
    site <- siteID
  }
  
  # Create covariate samples based on model type
  # Access environmental data from site-level data, not plot-level data
  if (site %in% rownames(model.inputs$temp)) {
    # Fix subscript out of bounds error by checking actual dimensions
    temp_cols <- ncol(model.inputs$temp)
    mois_cols <- ncol(model.inputs$mois)
    
    temp_end <- min(NT, temp_cols)
    mois_end <- min(NT, mois_cols)
    
    if (start_date <= temp_end) {
      temp_data <- model.inputs$temp[site, start_date:temp_end]
    } else {
      temp_data <- rep(NA, NT - start_date + 1)
    }
    
    if (start_date <= mois_end) {
      mois_data <- model.inputs$mois[site, start_date:mois_end]
    } else {
      mois_data <- rep(NA, NT - start_date + 1)
    }
  }
  
  # CRITICAL FIX: Handle actual model types (cycl_only, env_cov, env_cycl)
  if (model_type %in% c("env_cov", "env_cycl")) {
    # Full environmental model
    # CRITICAL FIX: Handle cases where environmental data doesn't extend to full time series
    # Get the actual available timepoints for each environmental variable
    temp_available <- min(NT, ncol(model.inputs$temp))
    mois_available <- min(NT, ncol(model.inputs$mois))
    pH_available <- min(NT, ncol(model.inputs$pH))
    pC_available <- min(NT, ncol(model.inputs$pC))
    relEM_available <- min(NT, ncol(model.inputs$relEM))
    sin_mo_available <- min(NT, length(model.inputs$sin_mo))
    cos_mo_available <- min(NT, length(model.inputs$cos_mo))
    
    # Create environmental data arrays, extending with NA for missing timepoints
    # CRITICAL FIX: Always create arrays with full time series length (NT)
    temp_data <- rep(NA, NT)
    mois_data <- rep(NA, NT)
    pH_data <- rep(NA, NT)
    pC_data <- rep(NA, NT)
    relEM_data <- rep(NA, NT)
    LAI_data <- rep(NA, NT)  # CRITICAL FIX: Add LAI array
    sin_mo_data <- rep(NA, NT)
    cos_mo_data <- rep(NA, NT)
    
    
    # Fill available data - use direct indexing since arrays are now full length NT
    if (1 <= temp_available) {
      temp_end <- min(NT, temp_available)
      temp_data[1:temp_end] <- model.inputs$temp[site, 1:temp_end]
    }
    if (1 <= mois_available) {
      mois_end <- min(NT, mois_available)
      mois_data[1:mois_end] <- model.inputs$mois[site, 1:mois_end]
    }
    if (!is.null(model.inputs$pH) && !is.null(plotID) && is.character(plotID) && plotID %in% rownames(model.inputs$pH)) {
      # pH is constant over time - extend to full time series length
      pH_value <- model.inputs$pH[plotID, 1]
      pH_data[] <- pH_value
    }
    if (!is.null(model.inputs$pC) && !is.null(plotID) && is.character(plotID) && plotID %in% rownames(model.inputs$pC)) {
      # pC is constant over time - extend to full time series length
      pC_value <- model.inputs$pC[plotID, 1]
      pC_data[] <- pC_value
    }
    if (1 <= relEM_available && !is.null(plotID) && is.character(plotID) && plotID %in% rownames(model.inputs$relEM)) {
      relEM_end <- min(NT, relEM_available)
      relEM_data[1:relEM_end] <- model.inputs$relEM[plotID, 1:relEM_end]
    }
    if (!is.null(model.inputs$LAI) && !is.null(plotID)) {
      # LAI data has site-level rownames, so extract site ID from plot ID
      siteID_for_LAI <- gsub("_.*", "", plotID)  # Extract site ID from plot ID (e.g., "BART_001" -> "BART")
      if (siteID_for_LAI %in% rownames(model.inputs$LAI)) {
        LAI_available <- min(NT, ncol(model.inputs$LAI))
        if (1 <= LAI_available) {
          LAI_end <- min(NT, LAI_available)
          LAI_data[1:LAI_end] <- model.inputs$LAI[siteID_for_LAI, 1:LAI_end]
        }
      } else {
      }
    }
    if (1 <= sin_mo_available) {
      sin_mo_end <- min(NT, sin_mo_available)
      sin_mo_data[1:sin_mo_end] <- model.inputs$sin_mo[1:sin_mo_end]
    }
    if (1 <= cos_mo_available) {
      cos_mo_end <- min(NT, cos_mo_available)
      cos_mo_data[1:cos_mo_end] <- model.inputs$cos_mo[1:cos_mo_end]
    }
    
    # CRITICAL FIX: Ensure all covariate arrays have the same length (full time series)
    # All arrays should already be length NT, but verify and extend if needed
    expected_length <- NT
    
    # Extend arrays that are shorter than expected
    if (length(temp_data) < expected_length) {
      temp_data <- c(temp_data, rep(NA, expected_length - length(temp_data)))
    }
    if (length(mois_data) < expected_length) {
      mois_data <- c(mois_data, rep(NA, expected_length - length(mois_data)))
    }
    if (length(pH_data) < expected_length) {
      pH_data <- c(pH_data, rep(NA, expected_length - length(pH_data)))
    }
    if (length(pC_data) < expected_length) {
      pC_data <- c(pC_data, rep(NA, expected_length - length(pC_data)))
    }
    if (length(relEM_data) < expected_length) {
      relEM_data <- c(relEM_data, rep(NA, expected_length - length(relEM_data)))
    }
    if (length(LAI_data) < expected_length) {
      LAI_data <- c(LAI_data, rep(NA, expected_length - length(LAI_data)))
    }
    if (length(sin_mo_data) < expected_length) {
      sin_mo_data <- c(sin_mo_data, rep(NA, expected_length - length(sin_mo_data)))
    }
    if (length(cos_mo_data) < expected_length) {
      cos_mo_data <- c(cos_mo_data, rep(NA, expected_length - length(cos_mo_data)))
    }
    
    
    covariate_samples <- list(
      temp = temp_data,
      mois = mois_data,
      pH = pH_data,
      pC = pC_data,
      relEM = relEM_data,
      LAI = LAI_data,  # CRITICAL FIX: Add LAI to covariate samples
      sin_mo = sin_mo_data,
      cos_mo = cos_mo_data,
      temp_sd = if (!is.null(model.inputs$temp_sd)) {
        temp_sd_data <- rep(0.1, NT)
        if (1 <= temp_available) {
          temp_sd_end <- min(NT, temp_available)
          temp_sd_data[1:temp_sd_end] <- model.inputs$temp_sd[site, 1:temp_sd_end]
        }
        temp_sd_data
      } else rep(0.1, NT),
      mois_sd = if (!is.null(model.inputs$mois_sd)) {
        mois_sd_data <- rep(0.1, NT)
        if (1 <= mois_available) {
          mois_sd_end <- min(NT, mois_available)
          mois_sd_data[1:mois_sd_end] <- model.inputs$mois_sd[site, 1:mois_sd_end]
        }
        mois_sd_data
      } else rep(0.1, NT),
      pH_sd = if (!is.null(model.inputs$pH_sd) && !is.null(plotID) && plotID %in% rownames(model.inputs$pH_sd)) {
        # pH_sd is constant over time - extend to full time series length
        pH_sd_value <- model.inputs$pH_sd[plotID, 1]
        rep(pH_sd_value, NT)
      } else rep(0.1, NT),
      pC_sd = if (!is.null(model.inputs$pC_sd) && !is.null(plotID) && plotID %in% rownames(model.inputs$pC_sd)) {
        # pC_sd is constant over time - extend to full time series length
        pC_sd_value <- model.inputs$pC_sd[plotID, 1]
        rep(pC_sd_value, NT)
      } else rep(0.1, NT),
      LAI_sd = if (!is.null(model.inputs$LAI_sd) && !is.null(plotID)) {
        # LAI_sd has site-level rownames, so extract site ID from plot ID
        siteID_for_LAI_sd <- gsub("_.*", "", plotID)
        if (siteID_for_LAI_sd %in% rownames(model.inputs$LAI_sd)) {
          LAI_sd_value <- model.inputs$LAI_sd[siteID_for_LAI_sd, 1]
          rep(LAI_sd_value, NT)
        } else {
          rep(0.1, NT)  # Default value if site not found
        }
      } else rep(0.1, NT)  # CRITICAL FIX: Add LAI_sd
    )
  } else if (model_type == "cycl_only") {
    # Cyclical-only model - only seasonal data
    sin_mo_data <- rep(NA, NT)
    cos_mo_data <- rep(NA, NT)
    
    # Get available seasonal data
    sin_mo_available <- min(NT, length(model.inputs$sin_mo))
    cos_mo_available <- min(NT, length(model.inputs$cos_mo))
    
    if (1 <= sin_mo_available) {
      sin_mo_end <- min(NT, sin_mo_available)
      sin_mo_data[1:sin_mo_end] <- model.inputs$sin_mo[1:sin_mo_end]
    }
    if (1 <= cos_mo_available) {
      cos_mo_end <- min(NT, cos_mo_available)
      cos_mo_data[1:cos_mo_end] <- model.inputs$cos_mo[1:cos_mo_end]
    }
    
    covariate_samples <- list(
      sin_mo = sin_mo_data,
      cos_mo = cos_mo_data
    )
  } else {
    # Unknown model type - error
    stop("Unknown model type: ", model_type, ". Expected: cycl_only, env_cov, or env_cycl")
  }
  
  # CRITICAL FIX: Always create 8 covariates in consistent order regardless of model type
  # This ensures the indexing in the main function works correctly
  # Order: temp, mois, pH, pC, relEM, LAI, sin_mo, cos_mo (indices 1-8)
  N.beta <- 8
  covariate_names <- c("temp", "mois", "pH", "pC", "relEM", "LAI", "sin_mo", "cos_mo")
  
  # CRITICAL FIX: Create 3D array for covariate samples with comprehensive validation
  # Each Monte Carlo sample gets its own covariate realization
  if (length(covariate_samples) > 0) {
    # Create 3D array: [Nmc, N.beta, NT]
    covariate_samples_array <- array(NA, dim = c(Nmc, N.beta, NT))
    
    # CRITICAL FIX: Validate input dimensions before processing
    if (Nmc <= 0 || N.beta <= 0 || NT <= 0) {
      stop("Invalid dimensions: Nmc=", Nmc, ", N.beta=", N.beta, ", NT=", NT)
    }
    
    # Debug: Check what covariates are available
    cat("DEBUG: Available covariates:", names(covariate_samples), "\n")
    cat("DEBUG: Expected covariates:", covariate_names, "\n")
    cat("DEBUG: Array dimensions: [", Nmc, ", ", N.beta, ", ", NT, "]\n")
    
    for (i in 1:Nmc) {
      for (j in 1:N.beta) {
        cov_name <- covariate_names[j]
        if (cov_name %in% names(covariate_samples)) {
          # Get mean and sd for this covariate
          mean_val <- covariate_samples[[cov_name]]
          sd_name <- paste0(cov_name, "_sd")
          sd_val <- if (sd_name %in% names(covariate_samples)) {
            covariate_samples[[sd_name]]
          } else {
            rep(0.1, length(mean_val))  # Default SD if not available
          }
          
          # Ensure we're working with vectors, not lists
          if (is.list(mean_val)) {
            mean_val <- unlist(mean_val)
          }
          if (is.list(sd_val)) {
            sd_val <- unlist(sd_val)
          }
          
          # Ensure sd_val has the same length as mean_val
          if (length(sd_val) != length(mean_val)) {
            sd_val <- rep(sd_val[1], length(mean_val))
          }
          
          # CRITICAL FIX: Enhanced validation and error handling for covariate sampling
          if (any(is.na(mean_val)) || any(is.na(sd_val))) {
            # Use default values for NA entries
            mean_val[is.na(mean_val)] <- 0
            sd_val[is.na(sd_val)] <- 0.1
            warning("NA values detected in ", cov_name, " - replaced with defaults")
          }
          
          # Debug: check for invalid values before rnorm
          if (any(is.infinite(mean_val)) || any(is.infinite(sd_val))) {
            mean_val[is.infinite(mean_val)] <- 0
            sd_val[is.infinite(sd_val)] <- 0.1
            warning("Infinite values detected in ", cov_name, " - replaced with defaults")
          }
          
          if (any(sd_val <= 0)) {
            sd_val[sd_val <= 0] <- 0.1
            warning("Non-positive standard deviation detected in ", cov_name, " - replaced with 0.1")
          }
          
          # CRITICAL FIX: Validate array dimensions before assignment
          if (length(mean_val) != NT || length(sd_val) != NT) {
            warning("Length mismatch for ", cov_name, ": mean_val=", length(mean_val), ", sd_val=", length(sd_val), ", NT=", NT)
            # Pad or truncate to match NT
            if (length(mean_val) < NT) {
              mean_val <- c(mean_val, rep(mean_val[length(mean_val)], NT - length(mean_val)))
            } else {
              mean_val <- mean_val[1:NT]
            }
            if (length(sd_val) < NT) {
              sd_val <- c(sd_val, rep(sd_val[length(sd_val)], NT - length(sd_val)))
            } else {
              sd_val <- sd_val[1:NT]
            }
          }
          
          # Store samples in the 3D array for this specific MC sample
          tryCatch({
            covariate_samples_array[i, j, ] <- rnorm(NT, mean_val, sd_val)
          }, error = function(e) {
            warning("Error sampling ", cov_name, " for MC sample ", i, ": ", e$message)
            covariate_samples_array[i, j, ] <<- rep(0, NT)
          })
        } else {
          # If covariate not available, fill with zeros
          covariate_samples_array[i, j, ] <- rep(0, NT)
        }
      }
    }
    
    # Return the 3D array directly
    return(covariate_samples_array)
  } else {
    # If no covariate samples, return empty array
    return(array(0, dim = c(Nmc, N.beta, NT)))
  }
}

# OPTIMIZATION FUNCTIONS - Add these before the main function


#' @title fg_fcast_beta
#' @description Forecast functional groups at NEON plots and sites, using beta regression output
#' @export
#' @title fcast_logit_beta
#' @description Unified forecasting function for observed and unobserved sites
#' @export
fcast_logit_beta <- function(plotID,
													model.inputs,
													param_samples,
													truth.plot.long,
                             plot_summary = NULL,
													Nmc = 1000,
													predict_site_effects = NULL,
                             rank.name = NULL,
													model_id,
													metadata = NULL,
													...) {

  cat("DEBUG: fcast_logit_beta function started\n")
  
  # Extract samples2 if the combined file object was passed
  if (is.list(param_samples) && "samples2" %in% names(param_samples)) {
    samples2_from_file <- param_samples$samples2
    param_samples_actual <- param_samples$samples
  } else {
    samples2_from_file <- NULL
    param_samples_actual <- param_samples
  }
  
  # Convert param_samples if needed
	if (inherits(param_samples_actual, "mcmc.list")) {
		# Convert each chain to matrix and combine
		param_matrices <- lapply(param_samples_actual, as.matrix)
		param_samples_actual <- do.call(rbind, param_matrices)
	}
	
	# Update param_samples variable to use the actual samples
	param_samples <- param_samples_actual
	
	cat("DEBUG: param_samples converted successfully\n")

	siteID <- substr(plotID, 1, 4)

  # Determine if site is observed or unobserved
  unobserved_sites <- c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", 
                        "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", 
                        "WREF", "YELL")
  is_new_site <- siteID %in% unobserved_sites
  
  # Process truth data
  truth_data <- if (is.data.frame(truth.plot.long)) {
    truth.plot.long
  } else if (is.list(truth.plot.long)) {
    if ("data" %in% names(truth.plot.long)) truth.plot.long$data else truth.plot.long[[1]]
		} else {
    stop("Invalid truth.plot.long format")
		}
	truth_data <- as.data.frame(truth_data)
	
  # Extract taxon name with enhanced validation
  taxon_name <- if (!is.null(metadata) && is.list(metadata) && "species" %in% names(metadata)) {
    metadata$species
  } else if (length(unique(truth_data$species)) > 0) {
    unique(truth_data$species)[1]
	} else {
    "unknown_taxon"
  }
  
  # CRITICAL FIX: Validate taxon name consistency
  # Check if taxon_name matches what's in truth_data to prevent mismatched information
  if (nrow(truth_data) > 0) {
    truth_species <- unique(truth_data$species)
    truth_species <- truth_species[!is.na(truth_species) & truth_species != ""]
    
    if (length(truth_species) > 0 && !taxon_name %in% truth_species) {
      warning("Taxon name mismatch detected: function parameter '", taxon_name, 
              "' does not match truth_data species: ", paste(truth_species, collapse=", "))
      # Use the most common species from truth_data as the authoritative taxon name
      taxon_name <- names(sort(table(truth_data$species), decreasing = TRUE))[1]
      cat("DEBUG: Corrected taxon_name to:", taxon_name, "\n")
    }
  }
  
  # CRITICAL FIX: Use the rank.name parameter directly instead of loading data
  # The rank.name should already be correctly determined by the calling script
  # This eliminates the massive data loading bottleneck that was causing CPU spikes
  
  # Parse model type
  model_name <- if (grepl("env_cycl", model_id)) {
    "env_cycl"
  } else if (grepl("env_cov", model_id)) {
    "env_cov"
  } else if (grepl("cycl_only", model_id)) {
    "cycl_only"
	} else {
    stop("Cannot determine model type from model_id: ", model_id)
  }
  
  # Sample parameters
  cat("DEBUG: About to sample parameters...\n")
  Nmc_large <- min(3000, nrow(param_samples))
  cat("DEBUG: Nmc_large =", Nmc_large, "\n")
  row_samples <- if (Nmc > Nmc_large) {
    sample.int(Nmc_large, Nmc, replace = TRUE)
	} else {
    sample.int(Nmc_large, Nmc, replace = FALSE)
	}
  cat("DEBUG: row_samples created, length =", length(row_samples), "\n")

  # OPTIMIZATION 1: Pre-extract all parameters
  cat("DEBUG: About to call extract_all_parameters...\n")
  params <- extract_all_parameters(param_samples, row_samples, model_name)
  cat("DEBUG: extract_all_parameters completed successfully\n")
  
  # Handle site effects for new vs observed sites
	if (is_new_site) {
    if (!is.null(predict_site_effects) && is.data.frame(predict_site_effects) && siteID %in% predict_site_effects$siteID) {
      # Use predicted site effects
      site_effect_val <- predict_site_effects$fit[predict_site_effects$siteID == siteID]
      params$site_effects[] <- site_effect_val  # Fill the pre-allocated vector
		} else {
			# Sample from site effect variance
			site_effect_sd_cols <- grep("site_effect_sd", colnames(param_samples))
      if (length(site_effect_sd_cols) > 0) {
        site_effect_sd <- mean(param_samples[row_samples, site_effect_sd_cols], na.rm = TRUE)
        params$site_effects[] <- rnorm(Nmc, 0, site_effect_sd)
      }
    }
			} else {
    # For observed sites, extract the specific site effect
		if ("site_num" %in% colnames(truth_data)) {
      site_num <- unique(truth_data[truth_data$siteID == siteID, ]$site_num)[1]
		} else {
      # Fallback: try to infer from site_start order
      site_num <- which(names(model.inputs$site_start) == siteID)
      if (length(site_num) == 0) site_num <- 1
    }
    
    site_param <- paste0("site_effect[", site_num, "]")
		if (site_param %in% colnames(param_samples)) {
      params$site_effects[] <- as.numeric(param_samples[row_samples, site_param])
    }
  }
  
  # OPTIMIZATION 2: Get covariates directly (bypass cache for parallel compatibility)
  cat("DEBUG: About to call create_covariate_samples_fixed...\n")
  covar <- create_covariate_samples_fixed(model.inputs, plotID = plotID, siteID = siteID,
                                         Nmc_large, Nmc, model_name)
  cat("DEBUG: create_covariate_samples_fixed completed successfully\n")
  
  # Select covariates based on model type
	if (model_name == "cycl_only") {
    covar <- covar[, 7:8, ]  # sin_mo, cos_mo
	} else if (model_name == "env_cov") {
    covar <- covar[, 1:6, ]   # temp, mois, pH, pC, relEM, LAI
	} else if (model_name == "env_cycl") {
    covar <- covar[, 1:8, ]   # all 8
  }
  
  # CRITICAL FIX: Enhanced validation with detailed error messages
  if (length(dim(covar)) != 3L) {
    stop("Covariate array must be 3D. Got dimensions: ", paste(dim(covar), collapse="x"))
  }
  if (dim(covar)[1] != Nmc) {
    stop("Covariate array first dimension (", dim(covar)[1], ") must match Nmc (", Nmc, ")")
  }
  if (dim(covar)[3] != model.inputs$N.date) {
    stop("Covariate array third dimension (", dim(covar)[3], ") must match model.inputs$N.date (", model.inputs$N.date, ")")
  }
  
  # Additional validation for parameter dimensions
  if (ncol(params$betas) != dim(covar)[2]) {
    warning("Beta parameter count (", ncol(params$betas), ") doesn't match covariate count (", dim(covar)[2], ")")
  }
  
  cat("DEBUG: Covariate array validation passed. Dimensions: [", dim(covar)[1], ", ", dim(covar)[2], ", ", dim(covar)[3], "]\n")
  
  # Create date key - use trained time map if provided, otherwise fallback
  if (!is.null(metadata) && is.list(metadata) && !is.null(metadata$trained_time_map)) {
    date_key <- metadata$trained_time_map %>%
      dplyr::transmute(
        date_num = trained_date_num,
        dateID = as.numeric(as.character(dateID))  # Ensure numeric type
      )
  } else {
    date_sequence <- seq.Date(from = as.Date("2013-06-01"), by = "month", 
                              length.out = model.inputs$N.date)
    date_key <- data.frame(
      date_num = 1:model.inputs$N.date,
      dateID = as.numeric(format(date_sequence, "%Y%m"))
    )
  }
  
  # Determine plot start date
  if (plotID %in% names(model.inputs$plot_start)) {
    plot_start_date <- model.inputs$plot_start[plotID]
  } else {
    # For plots not in model.inputs$plot_start, infer from truth data
    plot_dates <- truth_data %>% 
      filter(plotID == !!plotID, !is.na(date_num)) %>%
      pull(date_num)
    
    if (length(plot_dates) > 0) {
      plot_start_date <- min(plot_dates, na.rm = TRUE)
    } else {
      # Ultimate fallback - derive Jan 2018 index programmatically
      cal_boundary <- which(date_key$dateID == 201801)
      if (length(cal_boundary) == 0) cal_boundary <- 56
      plot_start_date <- if (is_new_site) 1 else cal_boundary
    }
  }
  
  # Set defaults first
  ic_mean <- 0.5
  ic_sd <- 0.1
  
  if (is_new_site) {
    # Use taxon-specific initial condition
    tryCatch({
      observed_data <- truth_data %>% 
        filter(species == taxon_name, !is.na(truth)) %>% 
        pull(truth)
      if (length(observed_data) > 0) {
        ic_mean <- mean(observed_data)
        ic_sd <- max(min(sd(observed_data), 0.3), 0.05)
      }
    }, error = function(e) {
      cat("Warning: Could not extract IC for unobserved site", plotID, ":", e$message, "\n")
      # Falls back to defaults already set
    })
  } else {
    # For observed sites: EXTRACT FROM PLOT_MU SAMPLES
    tryCatch({
      # First try to get samples2 (raw plot estimates) if available
      samples2_matrix <- NULL
      if (!is.null(samples2_from_file)) {
        # samples2 was extracted from combined file
        samples2_raw <- samples2_from_file
        if (inherits(samples2_raw, "mcmc.list")) {
          samples2_matrix <- as.matrix(do.call(rbind, samples2_raw))
        } else if (is.matrix(samples2_raw)) {
          samples2_matrix <- samples2_raw
        }
      }
      
      # Check if samples2 was passed separately or in plot_summary
      if (is.null(samples2_matrix) && !is.null(plot_summary)) {
        if (is.list(plot_summary) && "plot_mu" %in% names(plot_summary)) {
          samples2_matrix <- plot_summary$plot_mu
        }
      }
      
      if (!is.null(samples2_matrix) && is.matrix(samples2_matrix)) {
        # Extract from samples2 matrix
        
        # Get plot index for this plotID
        plot_indices <- unique(truth.plot.long$plot_num[truth.plot.long$plotID == plotID])
        if (length(plot_indices) > 0) {
          plot_idx <- plot_indices[1]  # Use first plot index if multiple
          
          # Get calibration timepoints
          if (!is.null(metadata) && is.list(metadata) && "model_data" %in% names(metadata)) {
            calibration_timepoints <- unique(metadata$model_data$truth.plot.long$timepoint)
          } else {
            calibration_timepoints <- 1:ncol(samples2_matrix)
          }
          
          # Extract plot_mu values for this plot during calibration period
          valid_timepoints <- calibration_timepoints[calibration_timepoints <= ncol(samples2_matrix)]
          plot_mu_values <- samples2_matrix[plot_idx, valid_timepoints]
          
          # Find last non-NA value
          last_valid_idx <- which(!is.na(plot_mu_values))
          if (length(last_valid_idx) > 0) {
            last_idx <- max(last_valid_idx)
            ic_mean <- plot_mu_values[last_idx]
            
            # Estimate sd from the variability in the last few timepoints
            recent_values <- plot_mu_values[max(1, last_idx-2):last_idx]
            recent_values <- recent_values[!is.na(recent_values)]
            if (length(recent_values) > 1) {
              ic_sd <- max(sd(recent_values), 0.05)
            } else {
              ic_sd <- 0.1
            }
          } else {
            # Fallback if no valid plot_mu values
            ic_mean <- 0.5
            ic_sd <- 0.1
          }
        } else {
          # Fallback if plot not found in truth data
          ic_mean <- 0.5
          ic_sd <- 0.1
        }
      } else if (!is.null(plot_summary)) {
        # FALLBACK: Try old plot_summary format
        if (inherits(plot_summary, "summary.mcmc")) {
          plot_est <- data.frame(
            rowname = rownames(plot_summary$quantiles),
            plot_summary$quantiles,
            stringsAsFactors = FALSE
          )
          plot_est <- plot_est %>%
            mutate(
              plotID = gsub("plot_mu\\[([^,]+),.*", "\\1", rowname),
              timepoint = as.numeric(gsub("plot_mu\\[.*,([^]]+)\\]", "\\1", rowname))
            ) %>%
            filter(grepl("^plot_mu\\[", rowname))
        } else if (is.data.frame(plot_summary)) {
          plot_est <- plot_summary
        } else if (is.list(plot_summary) && length(plot_summary) >= 2) {
          plot_est <- plot_summary[[2]]
        } else {
          plot_est <- NULL
        }
        
        # Extract last calibration timepoint estimate
        if (!is.null(plot_est) && nrow(plot_est) > 0) {
          plot_est_filtered <- plot_est %>% filter(plotID == !!plotID)
          
          # Get calibration timepoints only
          if (!is.null(metadata) && is.list(metadata) && "model_data" %in% names(metadata)) {
            calibration_timepoints <- unique(metadata$model_data$truth.plot.long$timepoint)
            plot_est_calibration <- plot_est_filtered %>% 
              filter(timepoint %in% calibration_timepoints)
          } else {
            plot_est_calibration <- plot_est_filtered
          }
          
          # Find last valid timepoint with quantile data
          possible_median_cols <- c("X50.", "50%", "med", "median")
          valid_col <- NULL
          for (col in possible_median_cols) {
            if (col %in% colnames(plot_est_calibration)) {
              valid_col <- col
              break
            }
          }
          
          if (!is.null(valid_col) && nrow(plot_est_calibration) > 0) {
            valid_timepoints <- plot_est_calibration %>% 
              filter(!is.na(.data[[valid_col]])) %>%
              pull(timepoint)
            
            if (length(valid_timepoints) > 0) {
              last_timepoint <- max(valid_timepoints)
              last_obs <- plot_est_calibration %>% filter(timepoint == last_timepoint)
              
              # Extract median and IQR for IC
              ic_mean_raw <- last_obs[[valid_col]]
              ic_mean <- if (length(ic_mean_raw) > 0) ic_mean_raw[1] else 0.5
              
              # Calculate sd from IQR if available
              if (all(c("25%", "75%") %in% colnames(last_obs))) {
                ic_sd <- (last_obs$`75%`[1] - last_obs$`25%`[1]) / (2 * qnorm(0.75))
              } else {
                ic_sd <- 0.1
              }
              
              if (is.na(ic_sd) || ic_sd <= 0) ic_sd <- 0.1
            } else {
              # Fallback if no valid estimates
              ic_mean <- 0.5
              ic_sd <- 0.1
            }
          } else {
            # No plot summary available, use fallback
            ic_mean <- 0.5
            ic_sd <- 0.1
          }
        } else {
          # No plot summary, use fallback
          ic_mean <- 0.5
          ic_sd <- 0.1
        }
      } else {
        # No plot_summary provided at all
        ic_mean <- 0.5
        ic_sd <- 0.1
      }
    }, error = function(e) {
      cat("Warning: Could not extract IC from plot_summary for", plotID, ":", e$message, "\n")
      # Falls back to defaults already set
    })
  }
  
  ic <- truncnorm::rtruncnorm(Nmc, mean = ic_mean, sd = ic_sd, a = 0, b = 1)
  
  # Define timepoints to forecast
  valid_timepoints <- plot_start_date:model.inputs$N.date
  
  # Determine if this is a driver uncertainty model
  use_cloglog <- if (!is.null(metadata) && is.list(metadata) && "has_driver_uncertainty" %in% names(metadata)) {
    metadata$has_driver_uncertainty
  } else {
    # Fallback: check if model_id contains driver uncertainty indicators
    if (is.character(model_id)) {
      grepl("driver_uncertainty", model_id) || grepl("logit_beta_driver_uncertainty", model_id)
    } else {
      FALSE
    }
  }
  
  
  # Build legacy vector (length NT) from model inputs/metadata
  NT <- model.inputs$N.date
  legacy_vec <- NULL

  if (!is.null(model.inputs$legacy)) {
    # preferred: plot x time matrix, use the row for this plot if available
    if (!is.null(plotID) && plotID %in% rownames(model.inputs$legacy)) {
      legacy_vec <- as.numeric(model.inputs$legacy[plotID, seq_len(NT)])
    } else {
      # fallback: site-level or overall mean across plots
      legacy_vec <- as.numeric(colMeans(model.inputs$legacy[, seq_len(NT), drop = FALSE], na.rm = TRUE))
    }
  } else if (!is.null(metadata) && is.list(metadata) && "legacy_by_time" %in% names(metadata)) {
    legacy_vec <- as.numeric(metadata$legacy_by_time)[seq_len(NT)]
  } else {
    # last-resort heuristic: first 60% legacy = 1, then 0
    legacy_vec <- as.numeric(seq_len(NT) <= floor(0.6 * NT))
  }

  legacy_vec[!is.finite(legacy_vec)] <- 0
  
  # OPTIMIZATION 3: Vectorized Monte Carlo simulation
  all_predictions <- vectorized_forecast(
		params = params,
		covar = covar,
		initial_conditions = ic,
		timepoints = valid_timepoints,
		Nmc = Nmc,
		use_cloglog = use_cloglog,
		legacy = legacy_vec
	)

	# Create output dataframe
  # Check if we have any valid predictions
  if (ncol(all_predictions) == 0 || nrow(all_predictions) == 0) {
    warning("No valid predictions generated for plotID ", plotID, " - returning empty data frame")
    ci <- data.frame(
      lo = numeric(0), lo_25 = numeric(0), med = numeric(0), hi_75 = numeric(0), hi = numeric(0),
      mean = numeric(0), sd = numeric(0), plotID = character(0), siteID = character(0), 
      species = character(0), new_site = logical(0)
    )
  } else {
    ci <- as.data.frame(t(apply(all_predictions, 2, quantile, 
                                c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))) %>%
      mutate(
        mean = apply(all_predictions, 2, mean, na.rm = TRUE),
        sd = apply(all_predictions, 2, sd, na.rm = TRUE),
        plotID = plotID,
        siteID = siteID,
        species = taxon_name,
        new_site = is_new_site
      )
  }
  
  # Only proceed with column assignments if we have data
  if (nrow(ci) > 0) {
    colnames(ci)[1:5] <- c("lo", "lo_25", "med", "hi_75", "hi")
    
    # CRITICAL FIX: Ensure confidence interval columns are numeric, not lists
    # This prevents the list column issue that causes plotting problems
    ci$lo <- as.numeric(ci$lo)
    ci$lo_25 <- as.numeric(ci$lo_25)
    ci$med <- as.numeric(ci$med)
    ci$hi_75 <- as.numeric(ci$hi_75)
    ci$hi <- as.numeric(ci$hi)
    ci$mean <- as.numeric(ci$mean)
    ci$sd <- as.numeric(ci$sd)
  }
  
  # FIX: Replace calibration period with fitted posterior instead of re-simulation
  # Use samples2 matrix if available, otherwise fall back to plot_summary
  plot_mu_samples <- NULL
  if (!is.null(samples2_matrix)) {
    plot_mu_samples <- samples2_matrix
  } else if (!is.null(plot_summary) && is.list(plot_summary) && !is.null(plot_summary$plot_mu)) {
    plot_mu_samples <- plot_summary$plot_mu
  }
  
  if (!is.null(plot_mu_samples) && nrow(ci) > 0) {
    # Determine calibration boundary
    cal_end <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
      as.numeric(metadata$cal_end_dateID)
    } else {
      201801  # legacy fallback
    }
    
    # Find which timepoints are in calibration period
    cal_timepoints <- valid_timepoints[valid_timepoints <= which(date_key$dateID == cal_end)]
    
    if (length(cal_timepoints) > 0 && ncol(plot_mu_samples) >= max(cal_timepoints)) {
      # Extract fitted posterior for calibration period
      cal_mu_samples <- plot_mu_samples[, cal_timepoints, drop = FALSE]
      
      # Calculate quantiles from fitted posterior
      cal_quantiles <- t(apply(cal_mu_samples, 2, quantile, 
                              c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
      cal_mean <- apply(cal_mu_samples, 2, mean, na.rm = TRUE)
      cal_sd <- apply(cal_mu_samples, 2, sd, na.rm = TRUE)
      
      # Replace calibration rows with fitted posterior
      cal_rows <- which(valid_timepoints %in% cal_timepoints)
      if (length(cal_rows) > 0) {
        ci[cal_rows, 1:5] <- cal_quantiles
        ci[cal_rows, "mean"] <- cal_mean
        ci[cal_rows, "sd"] <- cal_sd
      }
    }
  }
  
  # Add date information
  ci$date_num <- valid_timepoints
  
  ci <- left_join(ci, date_key, by = "date_num") %>%
    mutate(
      dates = as.Date(paste0(dateID, "01"), format = "%Y%m%d"),
      dateID = as.numeric(as.character(dateID))  # Ensure consistent numeric type
    )
  
  # Add truth values - ensure dateID format matches
  plot_obs <- truth_data %>% filter(plotID == !!plotID)
  if (nrow(plot_obs) > 0) {
    # Ensure dateID is numeric in both datasets
    truth_subset <- plot_obs %>%
      mutate(dateID = as.numeric(as.character(dateID))) %>%
      group_by(dateID) %>%
      summarise(truth = mean(truth, na.rm = TRUE), .groups = "drop")
    
    # Ensure ci$dateID is also numeric before joining
    ci <- ci %>%
      mutate(dateID = as.numeric(as.character(dateID))) %>%
      left_join(truth_subset, by = "dateID")
  } else {
    ci$truth <- NA
  }
  
  # Note: Raw observation data (sampleID level) is available in model.inputs$sample_values
  # but is not included in the main output to avoid breaking downstream rbind operations
  
  # Add metadata and start date information
  ci <- ci %>% mutate(
    model_id = model_id,
    model_name = model_name,
    fcast_period = if (is_new_site) {
      # For unobserved sites, entire period is hindcast (2013-2020)
      "hindcast"
    } else {
      # For observed sites: Extract calibration end date from metadata
      cal_end <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
        as.numeric(metadata$cal_end_dateID)
      } else {
        # legacy fallback if not provided
        201801
      }
      ifelse(dateID <= cal_end, "calibration", "hindcast")
    },
    plot_start = plot_start_date,  # Already defined earlier in function
    site_start = ifelse(is_new_site, 1, 
                       ifelse(siteID %in% names(model.inputs$site_start), 
                              model.inputs$site_start[siteID], NA))
  )
  
  # CRITICAL FIX: Validate output data before returning to prevent mismatched information
  validation_results <- validate_hindcast_output(ci, plotID, model_id)
  
  if (!validation_results$data_integrity_valid) {
    stop("Data validation failed for plotID ", plotID, ": ", 
         paste(validation_results$errors, collapse="; "))
  }
  
  if (length(validation_results$warnings) > 0) {
    warning("Data validation warnings for plotID ", plotID, ": ", 
            paste(validation_results$warnings, collapse="; "))
  }
  
  cat("DEBUG: Data validation passed for plotID", plotID, "\n")
  
  return(ci)
}

# =============================================================================
# DIAGNOSTIC FUNCTIONS
# =============================================================================

#' Generate diagnostic visualizations for hindcast data
#' @param hindcast_data Data frame with hindcast predictions
#' @param model_id Model identifier
#' @param taxon Taxon name
#' @param out_dir Optional output directory. If not provided, will create default structure
#' @export
generate_hindcast_diagnostics <- function(hindcast_data, model_id, taxon, out_dir = NULL) {
  # Ensure out_dir is provided and valid
  stopifnot(!is.null(out_dir) && nzchar(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

  # Sanity columns
  need <- c("plotID","siteID","dates","med","lo","hi","fcast_period")
  miss <- setdiff(need, names(hindcast_data))
  if (length(miss)) stop("hindcast_data missing: ", paste(miss, collapse=", "))

  # Dates must be Date
  if (!inherits(hindcast_data$dates, "Date")) {
    # handle numeric YYYYMM or POSIX-ish numbers
    if (is.numeric(hindcast_data$dates) && all(hindcast_data$dates > 10000, na.rm = TRUE) &&
        !"dateID" %in% names(hindcast_data)) {
      hindcast_data$dates <- as.Date(hindcast_data$dates, origin = "1970-01-01")
    } else if ("dateID" %in% names(hindcast_data)) {
      hindcast_data$dates <- as.Date(paste0(hindcast_data$dateID, "01"), "%Y%m%d")
    } else {
      hindcast_data$dates <- as.Date(hindcast_data$dates)
    }
  }

  # Pick one plot per site (your original behavior)
  sel_plots <- hindcast_data |>
    dplyr::group_by(siteID) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(plotID)

  # If nothing to plot, leave a breadcrumb and return
  if (!length(sel_plots)) {
    marker <- file.path(out_dir, "_NO_PLOTS_FOUND.txt")
    writeLines("No selected plots in hindcast_data.", marker)
    return(invisible(out_dir))
  }

  # Has truth?
  has_truth <- "truth" %in% names(hindcast_data)

  for (plot_id in sel_plots) {
    site_id <- hindcast_data$siteID[match(plot_id, hindcast_data$plotID)]
    df <- dplyr::filter(hindcast_data, plotID == plot_id)
    
    # Check if we have both site effect types
    has_site_effects <- "predicted_site_effect" %in% names(df)
    if (has_site_effects) {
      modeled_df <- dplyr::filter(df, predicted_site_effect == TRUE)
      random_df <- dplyr::filter(df, predicted_site_effect == FALSE)
    } else {
      modeled_df <- df
      random_df <- NULL
    }

    # Convert to numeric if needed, with better error handling
    convert_numeric <- function(data_df) {
      if (is.list(data_df$med)) {
        data_df$med <- sapply(data_df$med, function(x) {
          if (is.numeric(x) && length(x) > 0 && !is.na(x[1])) x[1] else NA_real_
        })
      } else {
        data_df$med <- as.numeric(data_df$med)
      }
      
      if (is.list(data_df$lo)) {
        data_df$lo <- sapply(data_df$lo, function(x) {
          if (is.numeric(x) && length(x) > 0 && !is.na(x[1])) x[1] else NA_real_
        })
      } else {
        data_df$lo <- as.numeric(data_df$lo)
      }
      
      if (is.list(data_df$hi)) {
        data_df$hi <- sapply(data_df$hi, function(x) {
          if (is.numeric(x) && length(x) > 0 && !is.na(x[1])) x[1] else NA_real_
        })
      } else {
        data_df$hi <- as.numeric(data_df$hi)
      }
      return(data_df)
    }
    
    # Convert both dataframes
    modeled_df <- convert_numeric(modeled_df)
    if (!is.null(random_df)) {
      random_df <- convert_numeric(random_df)
    }
    
    # Filter data for both site effect types
    filter_data <- function(data_df) {
      cal <- dplyr::filter(data_df, fcast_period == "calibration", !is.na(med), !is.na(lo), !is.na(hi), 
                           is.finite(med), is.finite(lo), is.finite(hi))
      hin <- dplyr::filter(data_df, fcast_period == "hindcast",    !is.na(med), !is.na(lo), !is.na(hi),
                           is.finite(med), is.finite(lo), is.finite(hi))
      return(list(cal = cal, hin = hin))
    }
    
    modeled_data <- filter_data(modeled_df)
    random_data <- if (!is.null(random_df)) filter_data(random_df) else list(cal = data.frame(), hin = data.frame())
    
    # Skip this plot if there's insufficient valid data
    total_cal <- nrow(modeled_data$cal) + nrow(random_data$cal)
    total_hin <- nrow(modeled_data$hin) + nrow(random_data$hin)
    if (total_cal == 0 && total_hin == 0) {
      cat("Skipping plot", plot_id, "- no valid data (all values are NA or non-finite)\n")
      next
    }

    # Only create plot if we have at least some valid data
    if (total_cal > 0 || total_hin > 0) {
      p <- ggplot2::ggplot(modeled_df, ggplot2::aes(x = dates))
      
      # Add modeled effects (blue)
      if (nrow(modeled_data$cal) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = modeled_data$cal, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightblue") +
          ggplot2::geom_line(  data = modeled_data$cal, ggplot2::aes(y = med), linewidth = 1, color = "blue")
      }
      if (nrow(modeled_data$hin) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = modeled_data$hin, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightgreen") +
          ggplot2::geom_line(  data = modeled_data$hin, ggplot2::aes(y = med), linewidth = 1, color = "green")
      }
      
      # Add random effects (red/orange) - only if we have both types
      if (!is.null(random_df) && nrow(random_data$cal) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = random_data$cal, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "lightcoral") +
          ggplot2::geom_line(  data = random_data$cal, ggplot2::aes(y = med), linewidth = 1, color = "red", linetype = "dashed")
      }
      if (!is.null(random_df) && nrow(random_data$hin) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = random_data$hin, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "lightyellow") +
          ggplot2::geom_line(  data = random_data$hin, ggplot2::aes(y = med), linewidth = 1, color = "orange", linetype = "dashed")
      }
      
      # Add truth data (same for both site effect types)
      if (has_truth && any(!is.na(modeled_df$truth))) {
        # CRITICAL FIX: Validate and correct truth values to prevent dateID corruption
        truth_values <- modeled_df$truth
        
        # Check if truth values look like dateIDs (large numbers > 10000)
        if (is.numeric(truth_values) && any(truth_values > 10000, na.rm = TRUE)) {
          # This is dateID corruption - set all truth values to NA to prevent plotting
          cat("WARNING: Truth values appear to be corrupted with dateID values (range:", 
              min(truth_values, na.rm = TRUE), "to", max(truth_values, na.rm = TRUE), 
              "). Setting truth values to NA.\n")
          modeled_df$truth <- NA_real_
          has_truth <- FALSE
        } else if (is.numeric(truth_values) && any(truth_values < 0 | truth_values > 1, na.rm = TRUE)) {
          # Truth values should be in [0,1] range
          cat("WARNING: Truth values outside [0,1] range (range:", 
              min(truth_values, na.rm = TRUE), "to", max(truth_values, na.rm = TRUE), 
              "). Setting out-of-range values to NA.\n")
          modeled_df$truth[modeled_df$truth < 0 | modeled_df$truth > 1] <- NA_real_
        }
        
        if (has_truth && any(!is.na(modeled_df$truth))) {
          p <- p +
            ggplot2::geom_point(data = modeled_df, ggplot2::aes(y = truth), color = "red", size = 1.8, alpha = 0.8) +
            ggplot2::geom_line( data = modeled_df, ggplot2::aes(y = truth), color = "red", linewidth = 0.8, alpha = 0.8)
        }
      }

      # boundary line
      if (nrow(modeled_data$cal) && nrow(modeled_data$hin)) {
        boundary <- max(modeled_data$cal$dates, na.rm = TRUE)
        p <- p + ggplot2::geom_vline(xintercept = boundary, linetype = "dashed", color = "gray", linewidth = 1)
      }

      # Use the actual taxon from the data, not the parameter
      actual_taxon <- unique(df$species)[1]
      if (is.na(actual_taxon) || actual_taxon == "") {
        actual_taxon <- unique(df$taxon)[1]
      }
      if (is.na(actual_taxon) || actual_taxon == "") {
        actual_taxon <- taxon  # fallback to parameter
      }
      
      # Create a more informative title that shows both the actual taxon and model ID
      title_text <- paste("Hindcast vs Truth —", plot_id, "(", site_id, ")")
      if (actual_taxon != taxon && actual_taxon != model_id) {
        subtitle_text <- paste("Taxon:", actual_taxon, "| Model ID:", model_id, "| Functional Group:", taxon)
      } else {
        subtitle_text <- paste("Taxon:", actual_taxon, "| Model:", model_id)
      }
      
      p <- p +
        ggplot2::labs(
          title = title_text,
          subtitle = subtitle_text,
          x = "Date", y = "Abundance"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      # Create filename that includes both actual taxon and model ID for clarity
      if (actual_taxon != taxon && actual_taxon != model_id) {
        filename <- paste0("hindcast_", plot_id, "_", actual_taxon, "_", gsub("env_cycl_", "", model_id), ".png")
      } else {
        filename <- paste0("hindcast_", plot_id, "_", actual_taxon, ".png")
      }
      outfile <- file.path(out_dir, filename)
      # Force a headless-safe device
      ggplot2::ggsave(filename = outfile, plot = p, width = 12, height = 6, dpi = 300, device = "png")
    }
  }

  invisible(out_dir)
}

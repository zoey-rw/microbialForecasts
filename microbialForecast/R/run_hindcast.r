# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Debug logging - set HINDCAST_DEBUG=true to enable verbose output
.hindcast_debug <- identical(Sys.getenv("HINDCAST_DEBUG"), "true")
.dcat <- function(...) if (.hindcast_debug) cat(...)

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
  eta_prev <- numeric(Nmc)  # AR(1) on eta scale (not log(mu))
  eta <- numeric(Nmc)
  shape1 <- numeric(Nmc)
  shape2 <- numeric(Nmc)
  
  # CRITICAL: Initialize eta_prev from initial mu using forward link
  # Model structure: AR(1) is on eta (linear predictor), so we need eta[0]
  ic_values <- pmin(0.999, pmax(0.001, initial_conditions))
  if (use_cloglog) {
    # cloglog forward link: eta = log(-log(1 - mu))
    eta_prev[] <- log(-log(1 - ic_values))
  } else {
    # logit forward link: eta = log(mu / (1 - mu))
    eta_prev[] <- log(ic_values / (1 - ic_values))
  }
  
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

      # legacy term for this time index (scalar)
      legacy_t <- if (!is.null(legacy) && length(legacy) >= time_idx && is.finite(legacy[time_idx]))
                    legacy[time_idx] else 0

      # CRITICAL FIX: AR(1) is on eta (linear predictor scale), matching model structure
      # Model: eta[t] ~ dnorm(rho * eta[t-1] + intercept + ...)
      # So we use eta_prev directly (already initialized from IC or previous step)
      eta[] <- params$rho * eta_prev + linear_predictor +
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
        1 - exp(-exp(eta_cap))              # cloglog^-1: mu <- 1 - exp(-exp(eta))
      } else {
        # logit^-1: mu <- exp(eta) / (1 + exp(eta)) = 1 / (1 + exp(-eta))
        # More numerically stable: use 1 / (1 + exp(-eta)) instead of exp(eta) / (1 + exp(eta))
        exp_neg_eta <- exp(-eta_cap)
        pmin(0.999, pmax(0.001, 1 / (1 + exp_neg_eta)))
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
      
      # CRITICAL: Store current eta as eta_prev for next iteration (AR(1) on eta scale)
      eta_prev[] <- eta[]
    }
  
  all_predictions
}

#' @export
#' Build a column index cache for param_samples to avoid repeated grep() calls.
#' Call once per taxon, reuse across all fcast_logit_beta calls.
#' @param param_samples Parameter samples matrix
#' @param model_name Model type name
#' @return Named list of column indices
build_col_index <- function(param_samples, model_name) {
  cn <- colnames(param_samples)
  rho_idx   <- grep("^rho(\\[|$)", cn)
  beta_idx  <- grep("^beta\\[", cn)
  int_idx   <- grep("^intercept(\\[|$)", cn)
  legacy_idx <- grep("^legacy_effect(\\[|$)", cn)

  # Precision: try precision, phi, kappa in order
  prec_col <- grep("^precision(\\[|$)", cn, value = TRUE)
  if (length(prec_col) == 0) prec_col <- grep("^phi(\\[|$)", cn, value = TRUE)
  if (length(prec_col) == 0) prec_col <- grep("^kappa(\\[|$)", cn, value = TRUE)
  if (length(prec_col) == 0) {
    stop("Could not find precision/phi/kappa parameter in samples. Available columns: ",
         paste(head(cn, 20), collapse=", "))
  }
  if (length(prec_col) > 1) {
    prec_col <- if (any(prec_col == "precision")) "precision" else prec_col[1]
  }

  # Select beta columns based on model type
  n_beta <- switch(model_name,
    cycl_only = 2L,
    env_cov   = 6L,
    env_cycl  = 8L,
    length(beta_idx)
  )

  list(
    rho = rho_idx,
    beta = beta_idx[seq_len(min(n_beta, length(beta_idx)))],
    intercept = int_idx,
    precision = prec_col,
    legacy = legacy_idx
  )
}

#' Pre-extract and cache parameters with proper dimension handling
#' @param col_idx Optional pre-computed column index from build_col_index()
extract_all_parameters <- function(param_samples, row_samples, model_name, col_idx = NULL) {
  if (is.null(col_idx)) col_idx <- build_col_index(param_samples, model_name)

  rho <- param_samples[row_samples, col_idx$rho]
  if (is.matrix(rho)) rho <- as.numeric(rho)

  betas <- param_samples[row_samples, col_idx$beta, drop = FALSE]

  intercept <- param_samples[row_samples, col_idx$intercept]
  if (is.matrix(intercept)) intercept <- as.numeric(intercept)

  precision <- param_samples[row_samples, col_idx$precision, drop = TRUE]
  if (is.matrix(precision)) precision <- as.numeric(precision)

  site_effects <- rep(0, length(row_samples))

  legacy_effect <- if (length(col_idx$legacy) > 0) {
    x <- param_samples[row_samples, col_idx$legacy[1], drop = TRUE]
    if (is.matrix(x)) as.numeric(x) else x
  } else {
    rep(0, length(row_samples))
  }

  list(
    rho = rho,
    betas = betas,
    intercept = intercept,
    precision = pmax(precision, 1e-6),
    site_effects = site_effects,
    legacy_effect = legacy_effect
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
  # CRITICAL FIX: Ensure plotID and siteID are single values to avoid "condition has length > 1" errors
  # This is a minimal defensive check - if they're vectors, take the first element
  if (length(plotID) > 1) {
    warning("plotID has length > 1, using first element: ", plotID[1])
    plotID <- plotID[1]
  }
  if (length(siteID) > 1) {
    warning("siteID has length > 1, using first element: ", siteID[1])
    siteID <- siteID[1]
  }
  
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
    # CRITICAL FIX: Ensure plotID is single value before condition checks
    plotID_check <- if (length(plotID) == 1) plotID else plotID[1]
    if (!is.null(model.inputs$pH) && !is.null(plotID_check) && is.character(plotID_check) && length(plotID_check %in% rownames(model.inputs$pH)) == 1 && plotID_check %in% rownames(model.inputs$pH)) {
      # pH is constant over time - extend to full time series length
      pH_value <- model.inputs$pH[plotID_check, 1]
      pH_data[] <- pH_value
    }
    if (!is.null(model.inputs$pC) && !is.null(plotID_check) && is.character(plotID_check) && length(plotID_check %in% rownames(model.inputs$pC)) == 1 && plotID_check %in% rownames(model.inputs$pC)) {
      # pC is constant over time - extend to full time series length
      pC_value <- model.inputs$pC[plotID_check, 1]
      pC_data[] <- pC_value
    }
    if (1 <= relEM_available && !is.null(plotID_check) && is.character(plotID_check) && length(plotID_check %in% rownames(model.inputs$relEM)) == 1 && plotID_check %in% rownames(model.inputs$relEM)) {
      relEM_end <- min(NT, relEM_available)
      relEM_data[1:relEM_end] <- model.inputs$relEM[plotID_check, 1:relEM_end]
    }
    if (!is.null(model.inputs$LAI) && !is.null(plotID_check)) {
      # LAI data has site-level rownames, so extract site ID from plot ID
      siteID_for_LAI <- gsub("_.*", "", plotID_check)  # Extract site ID from plot ID (e.g., "BART_001" -> "BART")
      if (length(siteID_for_LAI) == 1 && siteID_for_LAI %in% rownames(model.inputs$LAI)) {
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
      pH_sd = {
        plotID_check_sd2 <- if (length(plotID) == 1) plotID else plotID[1]
        if (!is.null(model.inputs$pH_sd) && !is.null(plotID_check_sd2) && length(plotID_check_sd2 %in% rownames(model.inputs$pH_sd)) == 1 && plotID_check_sd2 %in% rownames(model.inputs$pH_sd)) {
          # pH_sd is constant over time - extend to full time series length
          pH_sd_value <- model.inputs$pH_sd[plotID_check_sd2, 1]
          rep(pH_sd_value, NT)
        } else rep(0.1, NT)
      },
      pC_sd = {
        plotID_check_sd3 <- if (length(plotID) == 1) plotID else plotID[1]
        if (!is.null(model.inputs$pC_sd) && !is.null(plotID_check_sd3) && length(plotID_check_sd3 %in% rownames(model.inputs$pC_sd)) == 1 && plotID_check_sd3 %in% rownames(model.inputs$pC_sd)) {
          # pC_sd is constant over time - extend to full time series length
          pC_sd_value <- model.inputs$pC_sd[plotID_check_sd3, 1]
          rep(pC_sd_value, NT)
        } else rep(0.1, NT)
      },
      LAI_sd = {
        plotID_check_sd4 <- if (length(plotID) == 1) plotID else plotID[1]
        if (!is.null(model.inputs$LAI_sd) && !is.null(plotID_check_sd4)) {
          # LAI_sd has site-level rownames, so extract site ID from plot ID
          siteID_for_LAI_sd <- gsub("_.*", "", plotID_check_sd4)
          if (length(siteID_for_LAI_sd) == 1 && length(siteID_for_LAI_sd %in% rownames(model.inputs$LAI_sd)) == 1 && siteID_for_LAI_sd %in% rownames(model.inputs$LAI_sd)) {
          LAI_sd_value <- model.inputs$LAI_sd[siteID_for_LAI_sd, 1]
          rep(LAI_sd_value, NT)
          } else {
            rep(0.1, NT)  # Default value if site not found
          }
        } else rep(0.1, NT)
      }
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
  
  # CRITICAL FIX: Covariate order MUST match the fitted model structure
  # Model order (from 01_fitModels_driverUncertainty.R):
  #   env_cycl: beta[1]=sin_mo, beta[2]=cos_mo, beta[3]=temp, beta[4]=mois, beta[5]=pH, beta[6]=pC, beta[7]=relEM, beta[8]=LAI
  #   env_cov:  beta[1]=temp, beta[2]=mois, beta[3]=pH, beta[4]=pC, beta[5]=relEM, beta[6]=LAI
  #   cycl_only: beta[1]=sin_mo, beta[2]=cos_mo
  # Order for 8-covariate array: sin_mo, cos_mo, temp, mois, pH, pC, relEM, LAI (indices 1-8)
  N.beta <- 8
  covariate_names <- c("sin_mo", "cos_mo", "temp", "mois", "pH", "pC", "relEM", "LAI")
  
  # CRITICAL FIX: Create 3D array for covariate samples with comprehensive validation
  # Each Monte Carlo sample gets its own covariate realization
  if (length(covariate_samples) > 0) {
    # Create 3D array: [Nmc, N.beta, NT]
    covariate_samples_array <- array(NA, dim = c(Nmc, N.beta, NT))
    
    # CRITICAL FIX: Validate input dimensions before processing
    if (Nmc <= 0 || N.beta <= 0 || NT <= 0) {
      stop("Invalid dimensions: Nmc=", Nmc, ", N.beta=", N.beta, ", NT=", NT)
    }
    
    
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
          
          # Check for invalid values before rnorm
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
          
          # Covariates with driver uncertainty get per-draw noise to propagate
          # measurement uncertainty, matching the model's dnorm(*_est | *, *_sd) structure.
          # Deterministic covariates (sin_mo, cos_mo, relEM, LAI) use mean values only.
          driver_uncertain_covs <- c("temp", "mois", "pH", "pC")
          tryCatch({
            if (cov_name %in% driver_uncertain_covs) {
              covariate_samples_array[i, j, ] <- rnorm(NT, mean_val, sd_val)
            } else {
              covariate_samples_array[i, j, ] <- mean_val
            }
          }, error = function(e) {
            warning("Error assigning ", cov_name, " for MC sample ", i, ": ", e$message)
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
													col_idx = NULL,
													...) {

  # CRITICAL FIX: Ensure plotID is a single value to avoid "condition has length > 1" errors
  # This is a minimal defensive check - if it's a vector, take the first element
  .dcat("🔍 DEBUG fcast_logit_beta: plotID length:", length(plotID), "value:", paste(plotID, collapse=", "), "\n")

  if (length(plotID) > 1) {
    warning("plotID has length > 1 in fcast_logit_beta, using first element: ", plotID[1])
    plotID <- plotID[1]
    .dcat("🔍 DEBUG fcast_logit_beta: Fixed plotID to single value:", plotID, "\n")
  
  }
  
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
	

	siteID <- substr(plotID, 1, 4)
  # CRITICAL FIX: Ensure siteID is a single value (substr can return vector if plotID was vector)
  if (length(siteID) > 1) {
    siteID <- siteID[1]
  }

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
  
  # Prepare plot_mu samples (samples2) for both observed and unobserved flows
  samples2_matrix <- NULL
  if (!is.null(samples2_from_file)) {
    if (inherits(samples2_from_file, "mcmc.list")) {
      samples2_matrix <- as.matrix(do.call(rbind, samples2_from_file))
    } else if (is.matrix(samples2_from_file)) {
      samples2_matrix <- samples2_from_file
    }
  }
  if (is.null(samples2_matrix) && !is.null(plot_summary) && is.list(plot_summary)) {
    if ("plot_mu" %in% names(plot_summary) && is.matrix(plot_summary$plot_mu)) {
      samples2_matrix <- plot_summary$plot_mu
    }
  }

  # Sample parameters
  # Nmc_large is the maximum pool size to sample from (cap at 3000 for memory efficiency)
  Nmc_large <- min(3000, nrow(param_samples))
  row_samples <- if (Nmc > Nmc_large) {
    sample.int(Nmc_large, Nmc, replace = TRUE)
	} else {
    sample.int(Nmc_large, Nmc, replace = FALSE)
	}
  # OPTIMIZATION 1: Pre-extract all parameters
  params <- extract_all_parameters(param_samples, row_samples, model_name, col_idx = col_idx)
  
  # Handle site effects for new vs observed sites
	if (is_new_site) {
    siteID_check_pred <- if (length(siteID) == 1) siteID else siteID[1]
    if (!is.null(predict_site_effects) && is.data.frame(predict_site_effects) && length(siteID_check_pred %in% predict_site_effects$siteID) == 1 && siteID_check_pred %in% predict_site_effects$siteID) {
      # Use predicted site effects with uncertainty
      # CRITICAL FIX: Extract both fit (mean) and se_fit (standard error) for uncertainty propagation
      site_row <- predict_site_effects[predict_site_effects$siteID == siteID, , drop = FALSE]
      if (nrow(site_row) == 0) {
        stop("Site effect prediction provided for siteID ", siteID, " but no matching row found. This should never happen.")
      } else if (nrow(site_row) > 1) {
        # Multiple matches - this shouldn't happen if filtering is correct upstream
        # But if it does, take the first one and warn
        warning("Multiple site effect matches for siteID ", siteID, " (", nrow(site_row), " rows) - using first row")
        site_row <- site_row[1, , drop = FALSE]
      }
      
      site_effect_mean <- site_row$fit[1]
      site_effect_se <- if ("se_fit" %in% names(site_row) && !is.na(site_row$se_fit[1]) && is.finite(site_row$se_fit[1])) {
        site_row$se_fit[1]
      } else {
        # Fallback: if se_fit is missing or invalid, use a small default uncertainty
        # This should not happen if predictions are generated correctly, but provides robustness
        warning("se_fit missing or invalid for siteID ", siteID, " - using default uncertainty of 0.1")
        0.1
      }
      
      # Validate the predicted values
      if (is.na(site_effect_mean) || !is.finite(site_effect_mean)) {
        stop("Invalid site effect mean for siteID ", siteID, ": ", site_effect_mean, ". This should never happen when predict_site_effects is provided.")
      }
      if (site_effect_se <= 0 || !is.finite(site_effect_se)) {
        warning("Invalid site effect se_fit for siteID ", siteID, ": ", site_effect_se, " - using default uncertainty of 0.1")
        site_effect_se <- 0.1
      }
      
      # Bayesian shrinkage: combine PLSR prediction (likelihood) with prior N(0, sigma_prior)
      # where sigma_prior = site_effect_sd from the MCMC posterior.
      # This naturally handles extrapolation: when se_fit >> sigma_prior (high leverage,
      # poor prediction), the posterior shrinks toward zero (like random effects).
      # When se_fit << sigma_prior (good prediction), the posterior keeps the PLSR fit.
      #
      # posterior_mean = fit * sigma_prior^2 / (sigma_prior^2 + se_fit^2)
      # posterior_sd   = sqrt(sigma_prior^2 * se_fit^2 / (sigma_prior^2 + se_fit^2))

      site_effect_sd_cols <- grep("site_effect_sd", colnames(param_samples))
      if (length(site_effect_sd_cols) > 0) {
        site_effect_sd_per_draw <- param_samples[row_samples, site_effect_sd_cols[1], drop = TRUE]
        if (is.matrix(site_effect_sd_per_draw)) {
          site_effect_sd_per_draw <- as.numeric(site_effect_sd_per_draw)
        }

        # Per-draw shrinkage: each MC draw uses its own posterior site_effect_sd
        sigma_prior2 <- site_effect_sd_per_draw^2
        se_fit2 <- site_effect_se^2

        posterior_mean <- site_effect_mean * sigma_prior2 / (sigma_prior2 + se_fit2)
        posterior_sd <- sqrt(sigma_prior2 * se_fit2 / (sigma_prior2 + se_fit2))

        params$site_effects[] <- rnorm(Nmc, mean = posterior_mean, sd = posterior_sd)
      } else {
        # No site_effect_sd available — fall back to raw PLSR prediction
        params$site_effects[] <- rnorm(Nmc, mean = site_effect_mean, sd = site_effect_se)
      }
		} else {
			# Sample from site effect variance (random effects for new sites)
			site_effect_sd_cols <- grep("site_effect_sd", colnames(param_samples))
      if (length(site_effect_sd_cols) > 0) {
        site_effect_sd_per_draw <- param_samples[row_samples, site_effect_sd_cols[1], drop = TRUE]
        if (is.matrix(site_effect_sd_per_draw)) {
          site_effect_sd_per_draw <- as.numeric(site_effect_sd_per_draw)
        }

        # Make random site effects unique per site by incorporating siteID
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
          get(".Random.seed", envir = .GlobalEnv)
        } else {
          NULL
        }
        site_hash <- sum(utf8ToInt(siteID)) %% 10000
        set.seed(site_hash)
        params$site_effects[] <- rnorm(Nmc, 0, site_effect_sd_per_draw)
        if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
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
  covar <- create_covariate_samples_fixed(model.inputs, plotID = plotID, siteID = siteID,
                                         Nmc_large, Nmc, model_name)
  
  # CRITICAL FIX: Select covariates based on model type, matching the fitted model order
  # Array order is now: sin_mo, cos_mo, temp, mois, pH, pC, relEM, LAI (indices 1-8)
	if (model_name == "cycl_only") {
    # cycl_only: beta[1]=sin_mo, beta[2]=cos_mo
    covar <- covar[, 1:2, , drop = FALSE]  # sin_mo, cos_mo
	} else if (model_name == "env_cov") {
    # env_cov: beta[1]=temp, beta[2]=mois, beta[3]=pH, beta[4]=pC, beta[5]=relEM, beta[6]=LAI
    # CRITICAL FIX: Force garbage collection after subset to free original array memory
    covar_subset <- covar[, 3:8, , drop = FALSE]   # temp, mois, pH, pC, relEM, LAI
    rm(covar)  # Explicitly remove original array
    gc(verbose = FALSE)  # Force garbage collection
    covar <- covar_subset
	} else if (model_name == "env_cycl") {
    # env_cycl: beta[1]=sin_mo, beta[2]=cos_mo, beta[3]=temp, beta[4]=mois, beta[5]=pH, beta[6]=pC, beta[7]=relEM, beta[8]=LAI
    covar <- covar[, 1:8, , drop = FALSE]   # all 8 in order: sin_mo, cos_mo, temp, mois, pH, pC, relEM, LAI
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
  # CRITICAL FIX: Ensure plotID is a single value to avoid "condition has length > 1" error
  .dcat("🔍 DEBUG: Before plot_start check - plotID length:", length(plotID), "value:", paste(plotID, collapse=", "), "\n")

  plotID_single <- if (length(plotID) == 1) plotID else plotID[1]
  condition_result <- plotID_single %in% names(model.inputs$plot_start)
  .dcat("🔍 DEBUG: plot_start condition_result length:", length(condition_result), "value:", paste(condition_result, collapse=", "), "\n")

  if (length(condition_result) != 1) {
    stop("plot_start condition has length ", length(condition_result), " not 1")
  }
  if (condition_result) {
    plot_start_date <- model.inputs$plot_start[plotID_single]
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
      # CRITICAL FIX: Check if truth column exists before filtering
      # If truth column doesn't exist, keep defaults (ic_mean = 0.5, ic_sd = 0.1)
      if ("truth" %in% names(truth_data)) {
        observed_data <- truth_data %>% 
          filter(species == taxon_name, !is.na(truth)) %>% 
          pull(truth)
        if (length(observed_data) > 0) {
          ic_mean <- mean(observed_data)
          if (length(observed_data) > 1) {
            ic_sd <- max(min(sd(observed_data), 0.3), 0.05)
          } else {
            ic_sd <- 0.1
          }
        }
        # If no observed_data found, defaults remain (ic_mean = 0.5, ic_sd = 0.1)
      }
      # If truth column doesn't exist, defaults remain (ic_mean = 0.5, ic_sd = 0.1)
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
        # Get plot index for this plotID
        plot_indices <- unique(truth.plot.long$plot_num[truth.plot.long$plotID == plotID])
        if (length(plot_indices) == 0) {
          stop("Cannot extract initial condition: plotID ", plotID, " not found in truth.plot.long")
        }
        plot_idx <- plot_indices[1]

        # Get calibration timepoints (model timepoint indices)
        if (!is.null(metadata) && is.list(metadata) && "model_data" %in% names(metadata)) {
          calibration_timepoints <- unique(metadata$model_data$truth.plot.long$timepoint)
        } else {
          calibration_timepoints <- NULL
        }

        # Branch: MCMC format [n_draws, n_params] with colnames like plot_mu[plot, time]
        has_plot_mu_cols <- !is.null(colnames(samples2_matrix)) &&
          any(grepl("^plot_mu\\[|^mu\\[", colnames(samples2_matrix)))
        if (has_plot_mu_cols) {
          # samples2 format: rows = MCMC draws, columns = plot_mu[plot_idx, timepoint]
          # Use posterior mean (and sd) at last calibration timepoint so hindcast median starts there
          cal_mu_cols <- grep(paste0("^plot_mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
          if (length(cal_mu_cols) == 0) {
            cal_mu_cols <- grep(paste0("^mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
          }
          if (length(cal_mu_cols) == 0) {
            stop("Cannot extract initial condition: no plot_mu columns for plot_idx ", plot_idx, " in samples2_matrix")
          }
          timepoint_nums <- as.numeric(gsub(".*,\\s*([^]]+)\\]", "\\1", cal_mu_cols))
          if (!is.null(calibration_timepoints)) {
            in_cal <- timepoint_nums %in% calibration_timepoints
          } else {
            in_cal <- rep(TRUE, length(timepoint_nums))
          }
          cal_cols <- cal_mu_cols[in_cal]
          cal_tps <- timepoint_nums[in_cal]
          if (length(cal_cols) == 0) {
            stop("Cannot extract initial condition: no calibration timepoint columns for plotID ", plotID)
          }
          last_tp <- max(cal_tps)
          last_col <- cal_cols[cal_tps == last_tp][1]
          last_draws <- samples2_matrix[, last_col]
          last_draws <- last_draws[!is.na(last_draws) & is.finite(last_draws)]
          if (length(last_draws) == 0) {
            stop("Cannot extract initial condition: all NA/non-finite at last calibration timepoint for plotID ", plotID)
          }
          # Use posterior median so hindcast median at first timepoint equals calibration median at end
          ic_mean <- median(last_draws)
          ic_sd <- sd(last_draws)
          if (is.na(ic_sd) || ic_sd <= 0) ic_sd <- 0.1
        } else {
          # Summary matrix format [n_plots, n_timepoints]: one value per cell
          if (is.null(calibration_timepoints)) {
            calibration_timepoints <- seq_len(ncol(samples2_matrix))
          }
          valid_timepoints <- calibration_timepoints[calibration_timepoints <= ncol(samples2_matrix)]
          plot_mu_values <- samples2_matrix[plot_idx, valid_timepoints]
          last_valid_idx <- which(!is.na(plot_mu_values))
          if (length(last_valid_idx) == 0) {
            stop("Cannot extract initial condition: no valid plot_mu values found for plotID ", plotID, " in samples2_matrix")
          }
          last_idx <- max(last_valid_idx)
          ic_mean <- plot_mu_values[last_idx]
          recent_values <- plot_mu_values[max(1, last_idx - 2):last_idx]
          recent_values <- recent_values[!is.na(recent_values)]
          ic_sd <- if (length(recent_values) > 1) max(sd(recent_values), 0.05) else 0.1
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
              # CRITICAL: For observed sites, we must have valid estimates
              stop("Cannot extract initial condition: no valid timepoints with estimates found for plotID ", plotID, " in plot_summary")
            }
          } else {
            # CRITICAL: For observed sites, we must have plot summary data
            stop("Cannot extract initial condition: no valid quantile column found for plotID ", plotID, " in plot_summary")
          }
        } else {
          # CRITICAL: For observed sites, we must have plot estimates
          stop("Cannot extract initial condition: no plot estimates found for plotID ", plotID, " in plot_summary")
        }
      } else {
        # CRITICAL: For observed sites, we must have either samples2 or plot_summary
        stop("Cannot extract initial condition: neither samples2_matrix nor plot_summary available for observed site plotID ", plotID)
      }
    }, error = function(e) {
      # CRITICAL: For observed sites, initial condition extraction must succeed
      stop("Failed to extract initial condition from final posterior estimate for plotID ", plotID, ": ", e$message)
    })
  }
  
  ic <- truncnorm::rtruncnorm(Nmc, mean = ic_mean, sd = ic_sd, a = 0, b = 1)
  
  # Define timepoints to forecast (full time series - calibration will be replaced with fitted posterior)
  valid_timepoints <- plot_start_date:model.inputs$N.date
  
  # CRITICAL FIX: For observed sites, the IC should be used for the FIRST HINDCAST timepoint,
  # not the first timepoint overall. We need to determine where the hindcast period starts.
  # Determine calibration end timepoint for IC positioning
  cal_end_dateID <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
    as.numeric(metadata$cal_end_dateID)
  } else {
    201801  # legacy fallback
  }
  
  # Find the timepoint corresponding to cal_end_dateID
  cal_end_timepoint <- which(date_key$dateID == cal_end_dateID)
  if (length(cal_end_timepoint) == 0) {
    # Fallback: assume cal_end is around timepoint 56 (Jan 2018)
    cal_end_timepoint <- 56
  } else {
    cal_end_timepoint <- cal_end_timepoint[1]
  }
  
  # For observed sites, the first hindcast timepoint is AFTER the last calibration observation
  if (is_new_site) {
    first_hindcast_timepoint <- plot_start_date
  } else {
    first_hindcast_timepoint <- cal_end_timepoint + 1
    # Ensure we don't go beyond available timepoints
    first_hindcast_timepoint <- min(first_hindcast_timepoint, model.inputs$N.date)
  }
  
  # Transformation: detect from metadata first (most reliable), then model_id, then default to cloglog
  use_cloglog <- if (!is.null(metadata) && is.list(metadata) && "link_function" %in% names(metadata)) {
    # Use metadata if available (passed from calling script)
    identical(metadata$link_function, "cloglog")
  } else if (grepl("cloglog", model_id, ignore.case=TRUE)) {
    # Check model_id for cloglog
    TRUE
  } else if (grepl("logit", model_id, ignore.case=TRUE)) {
    # Check model_id for logit
    FALSE
  } else {
    # Fallback: default to cloglog for backward compatibility
    TRUE
  }
  
  # Debug: log which transformation is being used (helpful for troubleshooting scale mismatches)
  if (!is.null(metadata) && is.list(metadata) && "link_function" %in% names(metadata)) {
    .dcat("Using link function from metadata:", metadata$link_function, "(use_cloglog=", use_cloglog, ")\n")
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
  # CRITICAL FIX: For observed sites, we need to generate predictions starting from the first hindcast timepoint
  # so the IC (from last calibration observation) is used correctly
  # Determine first hindcast timepoint for observed sites
  if (!is_new_site) {
    cal_end_dateID_for_hindcast <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
      as.numeric(metadata$cal_end_dateID)
    } else {
      201801  # legacy fallback
    }
    cal_end_idx_for_hindcast <- which(date_key$dateID == cal_end_dateID_for_hindcast)
    if (length(cal_end_idx_for_hindcast) == 0) cal_end_idx_for_hindcast <- 56  # fallback
    hindcast_start_timepoint <- cal_end_idx_for_hindcast[1] + 1
    hindcast_start_timepoint <- min(hindcast_start_timepoint, model.inputs$N.date)
  }
  
  if (is_new_site) {
    # For unobserved sites, use full time series
    forecast_timepoints <- valid_timepoints
  } else {
    # For observed sites, only forecast hindcast period (IC will be used for first hindcast timepoint)
    forecast_timepoints <- valid_timepoints[valid_timepoints >= hindcast_start_timepoint]
  }
  
  all_predictions <- vectorized_forecast(
		params = params,
		covar = covar,
		initial_conditions = ic,
		timepoints = forecast_timepoints,
		Nmc = Nmc,
		use_cloglog = use_cloglog,
		legacy = legacy_vec
	)
  
  # For observed sites, we need to pad the predictions to match valid_timepoints
  # (calibration period will be filled with fitted posterior later)
  if (!is_new_site && length(forecast_timepoints) < length(valid_timepoints)) {
    # Create full prediction matrix
    full_predictions <- matrix(NA_real_, nrow = Nmc, ncol = length(valid_timepoints))
    # Map hindcast predictions to correct positions
    hindcast_indices <- which(valid_timepoints >= hindcast_start_timepoint)
    full_predictions[, hindcast_indices] <- all_predictions
    all_predictions <- full_predictions
  }

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
  # plot_summary can be:
  # 1. A data frame (pred.quantiles from model_summary[[3]]) with plot_num, timepoint, and quantile columns
  # 2. A summary.mcmc object with $quantiles and $statistics matrices
  # 3. A list with [[1]]=statistics and [[2]]=quantiles matrices
  cal_plot_data <- NULL
  
  # DEBUG: Log what plot_summary structure we received
  .dcat("DEBUG[calibration extraction]: Starting - plotID:", plotID, "is_new_site:", is_new_site, "\n")
  if (!is.null(plot_summary)) {
    .dcat("DEBUG[calibration extraction]: plot_summary structure - class:", paste(class(plot_summary), collapse=", "), 
        if (is.data.frame(plot_summary)) paste("(nrow:", nrow(plot_summary), "ncol:", ncol(plot_summary), ")") else
        if (is.list(plot_summary)) paste("(length:", length(plot_summary), ")") else "",
        "\n")
    if (is.list(plot_summary)) {
      .dcat("DEBUG[calibration extraction]:   List names:", paste(names(plot_summary), collapse=", "), "\n")
      if ("quantiles" %in% names(plot_summary)) {
        .dcat("DEBUG[calibration extraction]:   Has 'quantiles' element\n")
      }
      if ("statistics" %in% names(plot_summary)) {
        .dcat("DEBUG[calibration extraction]:   Has 'statistics' element\n")
      }
      if ("plot_mu" %in% names(plot_summary)) {
        .dcat("DEBUG[calibration extraction]:   Has 'plot_mu' element (matrix)\n")
      }
      if (length(plot_summary) > 0) {
        .dcat("DEBUG[calibration extraction]:   First element class:", paste(class(plot_summary[[1]]), collapse=", "), "\n")
        if (is.data.frame(plot_summary[[1]])) {
          .dcat("DEBUG[calibration extraction]:   First element dim:", paste(dim(plot_summary[[1]]), collapse=" x "), "\n")
          .dcat("DEBUG[calibration extraction]:   First element colnames:", paste(head(colnames(plot_summary[[1]]), 10), collapse=", "), "\n")
        } else if (is.matrix(plot_summary[[1]])) {
          .dcat("DEBUG[calibration extraction]:   First element (matrix) dim:", paste(dim(plot_summary[[1]]), collapse=" x "), "\n")
        }
      }
    }
    if (inherits(plot_summary, "summary.mcmc")) {
      .dcat("DEBUG[calibration extraction]:   summary.mcmc detected\n")
    }
  } else {
    .dcat("DEBUG[calibration extraction]: plot_summary is NULL\n")
  }
  
  if (!is.null(plot_summary)) {
    # Case 1: Already processed data frame (from model_summary[[3]] = pred.quantiles)
    # This is what the production script passes (cal_quants = model_summary[[3]])
    if (is.data.frame(plot_summary) && "plot_num" %in% names(plot_summary) && "timepoint" %in% names(plot_summary)) {
      # This is pred.quantiles - already has plot_num, timepoint, and quantile columns
      cal_plot_data <- plot_summary
      cat("DEBUG: Using Case 1 (data frame) - cal_plot_data has", nrow(cal_plot_data), "rows\n")
      cat("DEBUG:   Columns:", paste(head(colnames(cal_plot_data), 15), collapse=", "), "\n")
      # Check for quantile columns
      quantile_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
      has_quantiles <- all(quantile_cols %in% colnames(cal_plot_data))
      cat("DEBUG:   Has quantile columns:", has_quantiles, "\n")
      if (!has_quantiles) {
        # Try alternative column names (e.g., from different summary formats)
        alt_cols <- c("lo", "lo_25", "med", "hi_75", "hi", "q2.5", "q25", "q50", "q75", "q97.5")
        found_alt <- intersect(alt_cols, colnames(cal_plot_data))
        if (length(found_alt) > 0) {
          cat("DEBUG:   Found alternative quantile columns:", paste(found_alt, collapse=", "), "\n")
        }
      }
    # Case 1b: List with plot_mu matrix (from test script) - check if it's actually a summary.mcmc with plot_mu added
    } else if (is.list(plot_summary) && "plot_mu" %in% names(plot_summary) && inherits(plot_summary, "summary.mcmc")) {
      # This is a summary.mcmc object with plot_mu added by test script
      # Use the summary.mcmc structure (quantiles and statistics)
      cat("DEBUG: Detected summary.mcmc with plot_mu added - using summary.mcmc structure\n")
      plot_summary_for_parsing <- plot_summary
      class(plot_summary_for_parsing) <- "summary.mcmc"
      # Remove plot_mu temporarily to process as summary.mcmc
      plot_summary_for_parsing$plot_mu <- NULL
      # Fall through to Case 2 handling
      plot_summary <- plot_summary_for_parsing
      if (nrow(cal_plot_data) > 0) {
        cat("DEBUG:   Columns:", paste(head(colnames(cal_plot_data), 10), collapse=", "), "\n")
        cat("DEBUG:   plot_num range:", min(cal_plot_data$plot_num, na.rm=TRUE), "to", max(cal_plot_data$plot_num, na.rm=TRUE), "\n")
        cat("DEBUG:   timepoint range:", min(cal_plot_data$timepoint, na.rm=TRUE), "to", max(cal_plot_data$timepoint, na.rm=TRUE), "\n")
        if ("50%" %in% colnames(cal_plot_data)) {
          na_count <- sum(is.na(cal_plot_data$"50%"))
          non_zero_count <- sum(cal_plot_data$"50%" != 0, na.rm=TRUE)
          cat("DEBUG:   50% column - NA:", na_count, "non-zero:", non_zero_count, "out of", nrow(cal_plot_data), "\n")
        }
      }
    } else if (inherits(plot_summary, "summary.mcmc") || (is.list(plot_summary) && "quantiles" %in% names(plot_summary) && "statistics" %in% names(plot_summary))) {
      # Case 2: summary.mcmc object - need to parse rownames
      # Also handle list that has summary.mcmc structure (quantiles and statistics)
      cat("DEBUG: Using Case 2 (summary.mcmc) - parsing rownames\n")
      plot_quantiles_df <- plot_summary$quantiles
      plot_means_df <- plot_summary$statistics
      cat("DEBUG:   quantiles dim:", paste(dim(plot_quantiles_df), collapse=" x "), "\n")
      cat("DEBUG:   statistics dim:", paste(dim(plot_means_df), collapse=" x "), "\n")
      if (!is.null(rownames(plot_quantiles_df))) {
        cat("DEBUG:   First 3 rownames:", paste(head(rownames(plot_quantiles_df), 3), collapse=", "), "\n")
      }
      # Normalize rownames: mu[...] -> plot_mu[...]
      if (!is.null(rownames(plot_quantiles_df)) && any(grepl("^mu\\[", rownames(plot_quantiles_df)))) {
        rownames(plot_quantiles_df) <- sub("^mu\\[", "plot_mu[", rownames(plot_quantiles_df))
        cat("DEBUG:   Normalized rownames (mu -> plot_mu)\n")
      }
      # Load parse_plot_mu_vars if not available
      if (!exists("parse_plot_mu_vars", mode = "function")) {
        # Try to source helperFunctions.r
        helper_paths <- c("microbialForecast/R/helperFunctions.r", "R/helperFunctions.r", 
                         file.path(dirname(getwd()), "microbialForecast/R/helperFunctions.r"))
        for (hp in helper_paths) {
          if (file.exists(hp)) {
            source(hp, local = TRUE)
            break
          }
        }
      }
      # Use parse_plot_mu_vars() if available, otherwise replicate logic
      if (exists("parse_plot_mu_vars", mode = "function")) {
        cat("DEBUG:   Using parse_plot_mu_vars() function\n")
        cal_plot_data <- parse_plot_mu_vars(as.data.frame(plot_quantiles_df))
        cat("DEBUG:   Parsed cal_plot_data has", nrow(cal_plot_data), "rows\n")
        # Add quantile columns
        quantile_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
        if (all(quantile_cols %in% colnames(plot_quantiles_df))) {
          cal_plot_data[, quantile_cols] <- plot_quantiles_df[, quantile_cols]
          cat("DEBUG:   Added quantile columns to cal_plot_data\n")
        } else {
          cat("DEBUG:   WARNING: Missing quantile columns:", paste(setdiff(quantile_cols, colnames(plot_quantiles_df)), collapse=", "), "\n")
        }
        # Add Mean and SD from statistics
        if (!is.null(plot_means_df) && "Mean" %in% colnames(plot_means_df) && "SD" %in% colnames(plot_means_df)) {
          plot_means_parsed <- parse_plot_mu_vars(as.data.frame(plot_means_df))
          cal_plot_data$Mean <- plot_means_parsed$Mean[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                             paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
          cal_plot_data$SD <- plot_means_parsed$SD[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                         paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
        }
      } else if (requireNamespace("stringr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)) {
        # Fallback: replicate parse_plot_mu_vars logic
        cat("DEBUG:   Using fallback parsing (tidyr::separate)\n")
        plot_quantiles_df_conv <- as.data.frame(plot_quantiles_df)
        plot_quantiles_parsed <- plot_quantiles_df_conv %>% 
          mutate(rowname = rownames(.)) %>%
          tidyr::separate(rowname, sep=", ", into=c("plot_num_str","timepoint_str")) %>%
          mutate(
            plot_num = as.integer(gsub("plot_mu\\[|mu\\[", "", plot_num_str)),
            timepoint = as.integer(gsub("\\]", "", timepoint_str))
          )
        cal_plot_data <- plot_quantiles_parsed
        cat("DEBUG:   Parsed cal_plot_data has", nrow(cal_plot_data), "rows\n")
        # Add Mean and SD from statistics
        if (!is.null(plot_means_df) && "Mean" %in% colnames(plot_means_df) && "SD" %in% colnames(plot_means_df)) {
          plot_means_df_conv <- as.data.frame(plot_means_df)
          plot_means_parsed <- plot_means_df_conv %>% 
            mutate(rowname = rownames(.)) %>%
            tidyr::separate(rowname, sep=", ", into=c("plot_num_str","timepoint_str")) %>%
            mutate(
              plot_num = as.integer(gsub("plot_mu\\[|mu\\[", "", plot_num_str)),
              timepoint = as.integer(gsub("\\]", "", timepoint_str))
            )
          cal_plot_data$Mean <- plot_means_parsed$Mean[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                             paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
          cal_plot_data$SD <- plot_means_parsed$SD[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                         paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
        }
      }
    } else if (is.list(plot_summary) && length(plot_summary) >= 2) {
      # Case 3: List with [[1]]=statistics and [[2]]=quantiles
      plot_quantiles_df <- plot_summary[[2]]
      plot_means_df <- plot_summary[[1]]
      # Normalize rownames if needed
      if (!is.null(rownames(plot_quantiles_df)) && any(grepl("^mu\\[", rownames(plot_quantiles_df)))) {
        rownames(plot_quantiles_df) <- sub("^mu\\[", "plot_mu[", rownames(plot_quantiles_df))
      }
      # Use parse_plot_mu_vars() if available
      if (exists("parse_plot_mu_vars", mode = "function")) {
        cal_plot_data <- parse_plot_mu_vars(as.data.frame(plot_quantiles_df))
        quantile_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
        if (all(quantile_cols %in% colnames(plot_quantiles_df))) {
          cal_plot_data[, quantile_cols] <- plot_quantiles_df[, quantile_cols]
        }
        if (!is.null(plot_means_df) && "Mean" %in% colnames(plot_means_df) && "SD" %in% colnames(plot_means_df)) {
          plot_means_parsed <- parse_plot_mu_vars(as.data.frame(plot_means_df))
          cal_plot_data$Mean <- plot_means_parsed$Mean[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                             paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
          cal_plot_data$SD <- plot_means_parsed$SD[match(paste(cal_plot_data$plot_num, cal_plot_data$timepoint),
                                                         paste(plot_means_parsed$plot_num, plot_means_parsed$timepoint))]
        }
      }
    }
  }
  
  # For observed sites, we need to add calibration period from fitted posterior
  # For unobserved sites, ci already has all timepoints
  if (!is_new_site && !is.null(cal_plot_data) && nrow(cal_plot_data) > 0) {
    cat("DEBUG: Starting calibration period extraction for plotID:", plotID, "\n")
    # Determine calibration boundary
    cal_end_dateID_for_cal <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
      as.numeric(metadata$cal_end_dateID)
    } else {
      201801  # legacy fallback
    }
    cat("DEBUG:   cal_end_dateID_for_cal:", cal_end_dateID_for_cal, "\n")
    cat("DEBUG:   plot_start_date:", plot_start_date, "\n")
    
    # Find calibration timepoints (from plot_start_date to cal_end)
    # CRITICAL: Only include timepoints >= plot_start_date (values before plot_start are NA)
    cal_end_idx <- which(date_key$dateID == cal_end_dateID_for_cal)
    if (length(cal_end_idx) == 0) cal_end_idx <- 56  # fallback
    cal_timepoints <- valid_timepoints[valid_timepoints >= plot_start_date & valid_timepoints <= cal_end_idx[1]]
    cat("DEBUG:   cal_timepoints range:", min(cal_timepoints), "to", max(cal_timepoints), "(length:", length(cal_timepoints), ")\n")
    
    if (length(cal_timepoints) > 0) {
      # CRITICAL FIX: Extract fitted posterior for THIS SPECIFIC PLOT only
      # Get plot index for this plotID (same logic as used for IC)
      plot_indices <- unique(truth.plot.long$plot_num[truth.plot.long$plotID == plotID])
      if (length(plot_indices) > 0) {
        plot_idx <- plot_indices[1]  # Use first plot index if multiple
        cat("DEBUG:   plot_idx:", plot_idx, "\n")
        
        # Filter cal_plot_data to this plot and calibration timepoints
        # CRITICAL: Use plotID to filter (not plot_num) because plot_num can differ between cal_quants and truth.plot.long
        # plot_num is assigned based on order of appearance, which can vary between data preparation steps.
        # plotID and plot_num are NOT interchangeable unless a plot_key mapping exists from training model_dat.
        # CRITICAL: Only include timepoints >= plot_start_date to avoid NA values from before plot started
        if ("plotID" %in% names(cal_plot_data)) {
          # Validate that plotID exists in cal_plot_data before filtering
          available_plotIDs <- unique(cal_plot_data$plotID[!is.na(cal_plot_data$plotID)])
          if (!plotID %in% available_plotIDs) {
            cat("DEBUG:   WARNING: plotID", plotID, "not found in cal_plot_data\n")
            cat("DEBUG:     Available plotIDs:", length(available_plotIDs), "unique\n")
            cat("DEBUG:     Sample available:", paste(head(available_plotIDs, 5), collapse=", "), "\n")
            cat("DEBUG:     This plot was not in the training data - calibration extraction will be empty\n")
            cal_plot_filtered <- cal_plot_data[0, ]  # Empty data frame with same structure
          } else {
            # Use plotID for filtering (more reliable than plot_num)
            cal_plot_filtered <- cal_plot_data %>%
              filter(plotID == !!plotID & timepoint %in% cal_timepoints & timepoint >= plot_start_date) %>%
              arrange(timepoint)
          }
        } else {
          # Fallback: use plot_num if plotID not available (should not happen with proper model_summary structure)
          # WARNING: This is unsafe because plot_num can match the wrong plot
          cal_plot_filtered <- cal_plot_data %>%
            filter(plot_num == plot_idx & timepoint %in% cal_timepoints & timepoint >= plot_start_date) %>%
            arrange(timepoint)
          cat("DEBUG:   WARNING: Using plot_num for filtering (plotID not available in cal_plot_data)\n")
          cat("DEBUG:     This is unsafe - plot_num may match the wrong plot if data order differs\n")
          # Validate we got the right plot by checking if plotID exists in filtered result
          if ("plotID" %in% names(cal_plot_filtered) && nrow(cal_plot_filtered) > 0) {
            filtered_plotID <- unique(cal_plot_filtered$plotID[!is.na(cal_plot_filtered$plotID)])
            if (length(filtered_plotID) > 0 && !all(filtered_plotID == plotID)) {
              cat("DEBUG:     ERROR: plot_num filter returned wrong plotID:", paste(filtered_plotID, collapse=", "), 
                  "expected:", plotID, "\n")
              stop("CRITICAL: plot_num filtering returned wrong plot. cal_plot_data must have plotID column.")
            }
          }
        }
        
        cat("DEBUG:   After filtering - cal_plot_filtered has", nrow(cal_plot_filtered), "rows\n")
        if (nrow(cal_plot_filtered) > 0) {
          cat("DEBUG:     timepoint range:", min(cal_plot_filtered$timepoint, na.rm=TRUE), "to", max(cal_plot_filtered$timepoint, na.rm=TRUE), "\n")
          # Extract quantiles and means
          # CRITICAL: Preserve NA values - do not convert to 0
          # NOTE: pred.quantiles (model_summary[[3]]) has quantile columns, pred.means (model_summary[[2]]) has Mean/SD
          # For unobserved sites using pred.means, quantiles won't be available, but calibration extraction
          # is skipped anyway (is_new_site=TRUE). For observed sites using pred.quantiles, quantiles are available.
          quantile_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
          if (all(quantile_cols %in% names(cal_plot_filtered))) {
            cat("DEBUG:     Extracting quantiles from columns:", paste(quantile_cols, collapse=", "), "\n")
            cal_quantiles <- as.matrix(cal_plot_filtered[, quantile_cols])
            cal_timepoints_actual <- cal_plot_filtered$timepoint
            
            # Check for NA values in quantiles
            na_count_50 <- sum(is.na(cal_quantiles[, "50%"]))
            non_zero_count_50 <- sum(cal_quantiles[, "50%"] != 0, na.rm=TRUE)
            cat("DEBUG:     50% column - NA:", na_count_50, "non-zero:", non_zero_count_50, "out of", nrow(cal_quantiles), "\n")
            if (nrow(cal_quantiles) > 0 && non_zero_count_50 > 0) {
              cat("DEBUG:     Sample 50% values (first 5 non-NA):", 
                  paste(head(cal_quantiles[!is.na(cal_quantiles[, "50%"]) & cal_quantiles[, "50%"] != 0, "50%"], 5), collapse=", "), "\n")
            }
            
            # Get means and SDs (preserve NA values)
            cal_mean <- if ("Mean" %in% names(cal_plot_filtered)) cal_plot_filtered$Mean else rep(NA_real_, nrow(cal_plot_filtered))
            cal_sd <- if ("SD" %in% names(cal_plot_filtered)) cal_plot_filtered$SD else rep(NA_real_, nrow(cal_plot_filtered))
            
            # Create calibration period dataframe
            # CRITICAL: Preserve NA values - do not convert to 0
            # CRITICAL: Ensure all CI columns are numeric (not lists) to prevent plotting issues
            cal_ci <- data.frame(
              lo = as.numeric(cal_quantiles[, "2.5%"]),
              lo_25 = as.numeric(cal_quantiles[, "25%"]),
              med = as.numeric(cal_quantiles[, "50%"]),
              hi_75 = as.numeric(cal_quantiles[, "75%"]),
              hi = as.numeric(cal_quantiles[, "97.5%"]),
              mean = as.numeric(cal_mean),
              sd = as.numeric(cal_sd),
              plotID = plotID,
              siteID = siteID,
              species = taxon_name,
              new_site = is_new_site,
              date_num = cal_timepoints_actual,
              stringsAsFactors = FALSE
            )
            
            cat("DEBUG:     Created cal_ci with", nrow(cal_ci), "rows\n")
            cat("DEBUG:     cal_ci median range:", min(cal_ci$med, na.rm=TRUE), "to", max(cal_ci$med, na.rm=TRUE), "\n")
            cat("DEBUG:     cal_ci CI columns - lo (non-NA):", sum(!is.na(cal_ci$lo)), "hi (non-NA):", sum(!is.na(cal_ci$hi)), "\n")
            
            # CRITICAL FIX: Ensure ci has date_num before rbind
            # For observed sites, ci was padded to match valid_timepoints, but only hindcast period has predictions
            # We need to extract only the hindcast rows from ci and assign correct date_num
            if (!"date_num" %in% names(ci)) {
              if (!is_new_site && exists("forecast_timepoints") && exists("hindcast_start_timepoint")) {
                # ci has been padded to length(valid_timepoints), but only hindcast rows have predictions
                # Extract only hindcast rows (non-NA medians) and assign forecast_timepoints
                hindcast_mask <- !is.na(ci$med) & is.finite(ci$med)
                if (sum(hindcast_mask) == length(forecast_timepoints)) {
                  ci$date_num <- NA_real_
                  ci$date_num[hindcast_mask] <- forecast_timepoints
                  # Keep only hindcast rows
                  ci <- ci[hindcast_mask, , drop=FALSE]
                } else {
                  # Fallback: use forecast_timepoints directly (assuming ci matches)
                  ci$date_num <- forecast_timepoints[1:nrow(ci)]
                }
              } else if (exists("valid_timepoints")) {
                ci$date_num <- valid_timepoints[1:nrow(ci)]
              } else {
                # Fallback: use row numbers
                ci$date_num <- seq_len(nrow(ci))
              }
            } else {
              # If date_num already exists, ensure we only keep hindcast rows
              if (!is_new_site && exists("hindcast_start_timepoint")) {
                # Keep only rows with date_num >= hindcast_start_timepoint
                ci <- ci[ci$date_num >= hindcast_start_timepoint, , drop=FALSE]
              }
            }
            
            cat("DEBUG:     ci (hindcast) has", nrow(ci), "rows before rbind\n")
            
            # Combine calibration and hindcast periods
            ci <- rbind(cal_ci, ci)
            # Sort by date_num to ensure correct order
            ci <- ci[order(ci$date_num), ]
            cat("DEBUG:     Combined ci has", nrow(ci), "rows (calibration:", nrow(cal_ci), "+ hindcast:", nrow(ci)-nrow(cal_ci), ")\n")
          } else {
            cat("DEBUG:     WARNING: Missing quantile columns:", paste(setdiff(quantile_cols, names(cal_plot_filtered)), collapse=", "), "\n")
          }
        } else {
          cat("DEBUG:   WARNING: cal_plot_filtered is empty after filtering\n")
          cat("DEBUG:     This may indicate a plot_num mismatch between cal_plot_data and truth.plot.long\n")
          cat("DEBUG:     cal_plot_data plot_num range:", 
              if (nrow(cal_plot_data) > 0) paste(min(cal_plot_data$plot_num, na.rm=TRUE), "to", max(cal_plot_data$plot_num, na.rm=TRUE)) else "N/A",
              "\n")
          cat("DEBUG:     Expected plot_idx:", plot_idx, "for plotID:", plotID, "\n")
        }
      } else {
        cat("DEBUG:   WARNING: No plot_indices found for plotID:", plotID, "\n")
      }
    } else {
      cat("DEBUG:   WARNING: No calibration timepoints found\n")
    }
  } else {
    if (is_new_site) {
      cat("DEBUG: Skipping calibration extraction - is_new_site=TRUE\n")
    } else if (is.null(cal_plot_data)) {
      cat("DEBUG: Skipping calibration extraction - cal_plot_data is NULL\n")
    } else if (nrow(cal_plot_data) == 0) {
      cat("DEBUG: Skipping calibration extraction - cal_plot_data has 0 rows\n")
    }
  }
  
  if (!is_new_site && (is.null(cal_plot_data) || nrow(cal_plot_data) == 0)) {
    # Fallback: If cal_plot_data not available, try old approach with raw samples2
    # This requires plot_idx and cal_timepoints to be defined
    if (!is.null(samples2_matrix) && is.matrix(samples2_matrix) && nrow(ci) > 0) {
      # Determine calibration boundary for fallback
      cal_end_dateID_for_cal <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
        as.numeric(metadata$cal_end_dateID)
      } else {
        201801  # legacy fallback
      }
      cal_end_idx <- which(date_key$dateID == cal_end_dateID_for_cal)
      if (length(cal_end_idx) == 0) cal_end_idx <- 56  # fallback
      cal_timepoints <- valid_timepoints[valid_timepoints >= plot_start_date & valid_timepoints <= cal_end_idx[1]]
      
      # Get plot index
      plot_indices <- unique(truth.plot.long$plot_num[truth.plot.long$plotID == plotID])
      if (length(plot_indices) > 0 && length(cal_timepoints) > 0) {
        plot_idx <- plot_indices[1]
        
        if (!is.null(colnames(samples2_matrix)) && any(grepl("^plot_mu\\[|^mu\\[", colnames(samples2_matrix)))) {
          # samples2 format: samples (rows) x parameters (columns like plot_mu[plot, time])
          # Extract columns matching plot_mu[plot_idx, timepoint] for each timepoint
          cal_mu_cols <- grep(paste0("^plot_mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
          if (length(cal_mu_cols) == 0) {
            cal_mu_cols <- grep(paste0("^mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
          }
          
          # Filter to calibration timepoints (only >= plot_start_date to avoid NA values)
          if (length(cal_mu_cols) > 0) {
            # Extract timepoint numbers from column names (handle space after comma: "plot_mu[1, 1]")
            timepoint_nums <- as.numeric(gsub(".*,\\s*([^]]+)\\]", "\\1", cal_mu_cols))
            
            # Filter to timepoints >= plot_start_date and in cal_timepoints
            valid_timepoint_mask <- timepoint_nums >= plot_start_date & timepoint_nums %in% cal_timepoints
            valid_cols <- cal_mu_cols[valid_timepoint_mask]
            
            if (length(valid_cols) > 0) {
              # Extract samples for these columns (rows are samples, columns are timepoints)
              cal_mu_samples <- samples2_matrix[, valid_cols, drop = FALSE]
              
              # Calculate quantiles across samples (rows) for each timepoint (column)
              # CRITICAL: Preserve NA values - na.rm=TRUE only removes NA from quantile calculation, not from output
              cal_quantiles <- t(apply(cal_mu_samples, 2, quantile, 
                                      c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
              cal_mean <- apply(cal_mu_samples, 2, mean, na.rm = TRUE)
              cal_sd <- apply(cal_mu_samples, 2, sd, na.rm = TRUE)
              
              # Map back to cal_timepoints
              cal_timepoints_actual <- timepoint_nums[valid_timepoint_mask]
              
              # Create calibration period dataframe
              # CRITICAL: Preserve NA values - do not convert to 0
              cal_ci <- data.frame(
                lo = cal_quantiles[, 1],
                lo_25 = cal_quantiles[, 2],
                med = cal_quantiles[, 3],
                hi_75 = cal_quantiles[, 4],
                hi = cal_quantiles[, 5],
                mean = cal_mean,
                sd = cal_sd,
                plotID = plotID,
                siteID = siteID,
                species = taxon_name,
                new_site = is_new_site,
                date_num = cal_timepoints_actual,
                stringsAsFactors = FALSE
              )
              
              # Combine calibration and hindcast periods
              ci <- rbind(cal_ci, ci)
              # Sort by date_num to ensure correct order
              ci <- ci[order(ci$date_num), ]
            }
          }
        }
      }
    }
  }
  
  # Add date information
  # Use the actual timepoints that were forecast (forecast_timepoints for observed sites, valid_timepoints for unobserved)
  # CRITICAL FIX: If date_num already exists (from calibration period rbind), don't overwrite it
  if (!"date_num" %in% names(ci)) {
    if (is_new_site) {
      ci$date_num <- valid_timepoints
    } else {
      # For observed sites, we need to reconstruct the full timepoint sequence
      # (calibration period will be added from fitted posterior)
      ci$date_num <- forecast_timepoints
    }
  } else {
    # date_num already exists (from calibration period rbind) - ensure it's complete
    # For observed sites, if some rows are missing date_num, fill them from forecast_timepoints
    if (!is_new_site && exists("forecast_timepoints") && any(is.na(ci$date_num))) {
      na_mask <- is.na(ci$date_num)
      if (sum(na_mask) == length(forecast_timepoints)) {
        ci$date_num[na_mask] <- forecast_timepoints
      }
    }
  }
  
  ci <- left_join(ci, date_key, by = "date_num") %>%
    mutate(
      dates = as.Date(paste0(dateID, "01"), format = "%Y%m%d"),
      dateID = as.numeric(as.character(dateID))  # Ensure consistent numeric type
    )
  
  # Add truth values - ensure dateID format matches
  # CRITICAL FIX: Check if truth column exists in truth_data before using it
  if ("truth" %in% names(truth_data)) {
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
  } else {
    ci$truth <- NA
  }
  
  # Note: Raw observation data (sampleID level) is available in model.inputs$sample_values
  # but is not included in the main output to avoid breaking downstream rbind operations
  
  # Add metadata and start date information
  # CRITICAL: Extract cal_end BEFORE mutate to ensure it's evaluated correctly
  cal_end <- if (!is.null(metadata) && "cal_end_dateID" %in% names(metadata)) {
    as.numeric(metadata$cal_end_dateID)
  } else {
    # legacy fallback if not provided
    201801
  }
  
  # CRITICAL: Assign fcast_period BEFORE mutate to ensure correct comparison
  # Convert dateID to numeric explicitly to ensure comparison works
  ci$dateID <- as.numeric(as.character(ci$dateID))
  ci$fcast_period <- if (is_new_site) {
    # For unobserved sites, entire period is hindcast (2013-2020)
    rep("hindcast", nrow(ci))
  } else {
    # For observed sites: use cal_end extracted above
    # Direct vectorized assignment for reliability
    ifelse(!is.na(ci$dateID) & ci$dateID <= cal_end, "calibration", "hindcast")
  }
  
  # CRITICAL FIX: Ensure siteID is a single value before using in mutate
  siteID_single <- if (length(siteID) == 1) siteID else siteID[1]
  site_start_val <- if (is_new_site) {
    1
  } else {
    if (siteID_single %in% names(model.inputs$site_start)) {
      model.inputs$site_start[siteID_single]
    } else {
      NA
    }
  }
  
  ci <- ci %>% mutate(
    model_id = model_id,
    model_name = model_name,
    plot_start = plot_start_date,  # Already defined earlier in function
    site_start = site_start_val
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

  # Pick one plot per site (original behavior)
  sel_plots <- hindcast_data |>
    dplyr::group_by(siteID) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(plotID)

  # If nothing to plot, still create a placeholder figure
  if (!length(sel_plots)) {
    placeholder <- data.frame(x = Sys.Date() + 0:1, y = c(NA_real_, NA_real_))
    p <- ggplot2::ggplot(placeholder, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_blank() +
      ggplot2::labs(title = "No plots available", subtitle = paste("Model:", model_id, "Taxon:", taxon)) +
      ggplot2::theme_minimal()
    outfile <- file.path(out_dir, "hindcast_no_plots.png")
    ggplot2::ggsave(filename = outfile, plot = p, width = 8, height = 4, dpi = 200, device = "png")
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
    
    # Check if this is a new site (both types present and no calibration data)
    is_new_site <- has_site_effects && nrow(modeled_data$cal) == 0 && nrow(random_data$cal) == 0 && 
                   nrow(modeled_data$hin) > 0 && nrow(random_data$hin) > 0
    
    # Use the actual taxon from the data, not the parameter
    actual_taxon <- unique(df$species)[1]
    if (is.na(actual_taxon) || actual_taxon == "") {
      actual_taxon <- unique(df$taxon)[1]
    }
    if (is.na(actual_taxon) || actual_taxon == "") {
      actual_taxon <- taxon  # fallback to parameter
    }
    
    # For new sites with both site effect types, create a side-by-side comparison
    # DISABLED for unobserved sites to speed up processing
    # if (is_new_site && has_site_effects) {
    if (FALSE && is_new_site && has_site_effects) {
      # Create combined data frame with site effect type labels
      combined_df <- dplyr::bind_rows(
        dplyr::mutate(modeled_df, site_effect_type = "Modeled (PLSR)"),
        dplyr::mutate(random_df, site_effect_type = "Random")
      )
      
      # Add truth validation
      if (has_truth && any(is.finite(combined_df$truth))) {
        truth_values <- combined_df$truth
        if (is.numeric(truth_values) && any(truth_values > 10000, na.rm = TRUE)) {
          cat("WARNING: Truth values appear to be corrupted with dateID values. Setting truth values to NA.\n")
          combined_df$truth <- NA_real_
          has_truth <- FALSE
        } else if (is.numeric(truth_values) && any(truth_values < 0 | truth_values > 1, na.rm = TRUE)) {
          combined_df$truth[combined_df$truth < 0 | combined_df$truth > 1] <- NA_real_
        }
      }
      
      # Create side-by-side comparison plot
      p_compare <- ggplot2::ggplot(combined_df, ggplot2::aes(x = dates)) +
        ggplot2::facet_wrap(~ site_effect_type, ncol = 2, scales = "free_y") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = site_effect_type), alpha = 0.3) +
        ggplot2::geom_line(ggplot2::aes(y = med, color = site_effect_type), linewidth = 1.2) +
        ggplot2::scale_fill_manual(values = c("Modeled (PLSR)" = "steelblue", "Random" = "coral"),
                                   name = "Site Effect Type") +
        ggplot2::scale_color_manual(values = c("Modeled (PLSR)" = "darkblue", "Random" = "darkred"),
                                   name = "Site Effect Type") +
        ggplot2::labs(
          title = paste("Site Effect Comparison —", plot_id, "(", site_id, ")"),
          subtitle = paste("Taxon:", actual_taxon, "| Model:", model_id),
          x = "Date", y = "Abundance"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
        )
      
      # Add truth if available
      if (has_truth && any(!is.na(combined_df$truth))) {
        p_compare <- p_compare +
          ggplot2::geom_point(ggplot2::aes(y = truth), color = "black", size = 2, alpha = 0.8) +
          ggplot2::geom_line(ggplot2::aes(y = truth), color = "black", linewidth = 0.8, alpha = 0.6, linetype = "dotted")
      }
      
      # Save comparison plot
      filename_compare <- paste0("hindcast_comparison_", plot_id, "_", actual_taxon, ".png")
      outfile_compare <- file.path(out_dir, filename_compare)
      ggplot2::ggsave(filename = outfile_compare, plot = p_compare, width = 14, height = 6, dpi = 300, device = "png")
      cat("  ✅ Created side-by-side comparison plot for new site:", filename_compare, "\n")
    }
    
    # Create overlay plot (original behavior) for comparison
    {
      p <- ggplot2::ggplot(modeled_df, ggplot2::aes(x = dates))
      
      # Add modeled effects (blue/green)
      if (nrow(modeled_data$cal) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = modeled_data$cal, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightblue") +
          ggplot2::geom_line(  data = modeled_data$cal, ggplot2::aes(y = med), linewidth = 1, color = "blue")
      }
      if (nrow(modeled_data$hin) > 0) {
        p <- p + 
          ggplot2::geom_ribbon(data = modeled_data$hin, ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightgreen") +
          ggplot2::geom_line(  data = modeled_data$hin, ggplot2::aes(y = med), linewidth = 1, color = "darkgreen", linetype = "solid")
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
      if (has_truth && any(is.finite(modeled_df$truth))) {
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
            ggplot2::geom_point(data = modeled_df, ggplot2::aes(y = truth), color = "black", size = 2, alpha = 0.8) +
            ggplot2::geom_line( data = modeled_df, ggplot2::aes(y = truth), color = "black", linewidth = 0.8, alpha = 0.6, linetype = "dotted")
        }
      }

      # boundary line
      if (nrow(modeled_data$cal) > 0) {
        boundary <- max(modeled_data$cal$dates, na.rm = TRUE)
        if (is.finite(boundary)) {
          p <- p + ggplot2::geom_vline(xintercept = boundary, linetype = "dashed", color = "gray", linewidth = 1)
        }
      }
      
      # Create a more informative title that shows both the actual taxon and model ID
      title_text <- paste("Hindcast vs Truth —", plot_id, "(", site_id, ")")
      if (has_site_effects && is_new_site) {
        title_text <- paste(title_text, "— Overlay View")
      }
      if (actual_taxon != taxon && actual_taxon != model_id) {
        subtitle_text <- paste("Taxon:", actual_taxon, "| Model ID:", model_id, "| Functional Group:", taxon)
      } else {
        subtitle_text <- paste("Taxon:", actual_taxon, "| Model:", model_id)
      }
      
      # Add legend for site effect types if both are present
      p <- p +
        ggplot2::labs(
          title = title_text,
          subtitle = subtitle_text,
          x = "Date", y = "Abundance",
          color = "Site Effect",
          linetype = "Site Effect"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      if (has_site_effects && is_new_site) {
        # Add manual legend for overlay plot
        p <- p + ggplot2::scale_linetype_manual(
          values = c("Modeled" = "solid", "Random" = "dashed"),
          guide = "legend"
        )
      }
      
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

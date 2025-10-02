# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Vectorized Monte Carlo simulation for improved performance
vectorized_forecast <- function(params, covar, initial_conditions, timepoints, Nmc) {
  n_timepoints <- length(timepoints)
  all_predictions <- matrix(NA, nrow = Nmc, ncol = n_timepoints)
  
  # Set initial conditions
  all_predictions[, 1] <- initial_conditions
  
  # Pre-allocate vectors for efficiency
  linear_predictor <- numeric(Nmc)
  log_x_prev <- numeric(Nmc)
  log_Ex_mean <- numeric(Nmc)
  shape1 <- numeric(Nmc)
  shape2 <- numeric(Nmc)
  
  # Vectorized computation across all timepoints
  for (t in 2:n_timepoints) {
    time_idx <- timepoints[t]
    
    if (time_idx <= dim(covar)[3]) {
      # Handle both 2D (old) and 3D (new) covariate arrays
      if (length(dim(covar)) == 3) {
        # New 3D array: [Nmc, N.beta, NT] - each MC sample has its own covariate realization
        Z <- covar[, , time_idx]  # Shape: [Nmc, N.beta]
      } else {
        # Old 2D array: [N.beta, NT] - same covariate realization for all MC samples
        Z <- matrix(covar[, time_idx], nrow = Nmc, ncol = nrow(covar), byrow = TRUE)
      }
      
      # OPTIMIZATION 1: Vectorized linear predictor calculation
      # Handle dimension compatibility gracefully
      n_beta_use <- min(ncol(Z), ncol(params$betas))
      linear_predictor[] <- rowSums(Z[, 1:n_beta_use] * params$betas[, 1:n_beta_use])
      
      # OPTIMIZATION 2: Vectorized AR(1) process
      log_x_prev[] <- log(pmax(0.01, pmin(0.99, all_predictions[, t-1])))
      log_Ex_mean[] <- params$rho * log_x_prev + linear_predictor + params$site_effects + params$intercept
      
      # OPTIMIZATION 3: Vectorized beta distribution sampling
      mu <- 1 / (1 + exp(-log_Ex_mean))
      shape1[] <- mu * params$precision
      shape2[] <- (1 - mu) * params$precision
      
      # Ensure positive parameters
      shape1[] <- pmax(shape1, 1e-6)
      shape2[] <- pmax(shape2, 1e-6)
      
      # Sample from beta distribution
      all_predictions[, t] <- rbeta(Nmc, shape1, shape2)
    }
  }
  
  return(all_predictions)
}

#' Pre-extract and cache parameters with proper dimension handling
extract_all_parameters <- function(param_samples, row_samples, model_name) {
  # Extract parameters with proper dimension handling
  rho <- param_samples[row_samples, grep("rho", colnames(param_samples))]
  if (is.matrix(rho)) rho <- as.numeric(rho)
  
  # CRITICAL FIX: Extract betas based on model type to maintain proper matrix structure
  beta_cols <- grep("beta\\[", colnames(param_samples))
  if (model_name == "cycl_only") {
    # cycl_only uses only 2 covariates (sin_mo, cos_mo)
    betas <- param_samples[row_samples, beta_cols[1:2], drop = FALSE]
  } else if (model_name == "env_cov") {
    # env_cov uses 6 environmental covariates
    betas <- param_samples[row_samples, beta_cols[1:6], drop = FALSE]
  } else if (model_name == "env_cycl") {
    # env_cycl uses all 8 covariates
    betas <- param_samples[row_samples, beta_cols[1:8], drop = FALSE]
  } else {
    # Fallback: use all available beta columns
    betas <- param_samples[row_samples, beta_cols, drop = FALSE]
  }
  
  intercept <- param_samples[row_samples, grep("intercept", colnames(param_samples))]
  if (is.matrix(intercept)) intercept <- as.numeric(intercept)
  
  precision <- param_samples[row_samples, grep("precision|sigma", colnames(param_samples))]
  if (is.matrix(precision)) precision <- as.numeric(precision)
  
  # Handle site effects - check if they exist
  site_effect_cols <- grep("site_effect", colnames(param_samples))
  if (length(site_effect_cols) > 0) {
    site_effects <- param_samples[row_samples, site_effect_cols]
    if (is.matrix(site_effects)) site_effects <- as.numeric(site_effects)
  } else {
    site_effects <- rep(0, length(row_samples))
  }
  
  # Handle legacy effects - check if they exist
  legacy_cols <- grep("legacy_effect", colnames(param_samples))
  if (length(legacy_cols) > 0) {
    legacy_effect <- param_samples[row_samples, legacy_cols]
    if (is.matrix(legacy_effect)) legacy_effect <- as.numeric(legacy_effect)
  } else {
    legacy_effect <- rep(0, length(row_samples))
  }
  
  return(list(
    rho = rho,
    betas = betas,
    intercept = intercept,
    precision = precision,
    site_effects = site_effects,
    legacy_effect = legacy_effect
  ))
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
    if (!is.null(model.inputs$pH) && !is.null(plotID) && plotID %in% rownames(model.inputs$pH)) {
      # pH is constant over time - extend to full time series length
      pH_value <- model.inputs$pH[plotID, 1]
      pH_data[] <- pH_value
    }
    if (!is.null(model.inputs$pC) && !is.null(plotID) && plotID %in% rownames(model.inputs$pC)) {
      # pC is constant over time - extend to full time series length
      pC_value <- model.inputs$pC[plotID, 1]
      pC_data[] <- pC_value
    }
    if (1 <= relEM_available && !is.null(plotID) && plotID %in% rownames(model.inputs$relEM)) {
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
  
  # CRITICAL FIX: Create 3D array for covariate samples instead of overwriting
  # Each Monte Carlo sample gets its own covariate realization
  if (length(covariate_samples) > 0) {
    # Create 3D array: [Nmc, N.beta, NT]
    covariate_samples_array <- array(NA, dim = c(Nmc, N.beta, NT))
    
    for (i in 1:Nmc) {
      for (j in 1:N.beta) {
        cov_name <- covariate_names[j]
        if (cov_name %in% names(covariate_samples)) {
          # Get mean and sd for this covariate
          mean_val <- covariate_samples[[cov_name]]
          sd_val <- covariate_samples[[paste0(cov_name, "_sd")]]
          
          # Ensure we're working with vectors, not lists
          if (is.list(mean_val)) {
            mean_val <- unlist(mean_val)
          }
          if (is.list(sd_val)) {
            sd_val <- unlist(sd_val)
          }
          
          if (is.null(sd_val)) {
            sd_val <- rep(0.1, length(mean_val))
          }
          
          # Sample from normal distribution - handle NA values
          if (any(is.na(mean_val)) || any(is.na(sd_val))) {
            # Use default values for NA entries
            mean_val[is.na(mean_val)] <- 0
            sd_val[is.na(sd_val)] <- 0.1
          }
          
          # Debug: check for invalid values before rnorm
          if (any(is.infinite(mean_val)) || any(is.infinite(sd_val))) {
            message("WARNING: Infinite values in ", cov_name, " - replacing with defaults")
            mean_val[is.infinite(mean_val)] <- 0
            sd_val[is.infinite(sd_val)] <- 0.1
          }
          
          if (any(sd_val <= 0)) {
            message("WARNING: Non-positive sd values in ", cov_name, " - replacing with 0.1")
            sd_val[sd_val <= 0] <- 0.1
          }
          
          # Store samples in the 3D array for this specific MC sample
          covariate_samples_array[i, j, ] <- rnorm(NT, mean_val, sd_val)
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

  # Convert param_samples if needed
	if (inherits(param_samples, "mcmc.list")) {
		param_samples <- as.matrix(param_samples)
	}

	siteID <- substr(plotID, 1, 4)

  # Determine if site is observed or unobserved
  unobserved_sites <- c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", 
                        "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", 
                        "WREF", "YELL")
  is_new_site <- siteID %in% unobserved_sites || 
                 !(siteID %in% names(model.inputs$site_start))
  
  # Process truth data
  truth_data <- if (is.data.frame(truth.plot.long)) {
    truth.plot.long
  } else if (is.list(truth.plot.long)) {
    if ("data" %in% names(truth.plot.long)) truth.plot.long$data else truth.plot.long[[1]]
		} else {
    stop("Invalid truth.plot.long format")
		}
	truth_data <- as.data.frame(truth_data)
	
  # Extract taxon name
  taxon_name <- if (!is.null(metadata) && "species" %in% names(metadata)) {
    metadata$species
  } else if (length(unique(truth_data$species)) > 0) {
    unique(truth_data$species)[1]
	} else {
    "unknown_taxon"
  }
  
  # CRITICAL FIX: Apply rank name mapping fix for truth data extraction
  if (!is.null(rank.name) && !is.null(taxon_name) && taxon_name != "unknown_taxon") {
    # Load all_ranks data to find the correct rank name
    all_ranks <- c(readRDS("./data/clean/groupAbundances_16S_2023.rds"),
                   readRDS("./data/clean/groupAbundances_ITS_2023.rds"))
    
    # Find which data source contains this taxon
    found_sources <- c()
    for (source_name in names(all_ranks)) {
      if (taxon_name %in% colnames(all_ranks[[source_name]])) {
        found_sources <- c(found_sources, source_name)
      }
    }
    
    if (length(found_sources) > 0) {
      # Use the first source found
      correct_source <- found_sources[1]
      
      # For taxonomic taxa, we need to map to the full rank name
      if (grepl("_bac$", correct_source)) {
        correct_rank <- correct_source  # Keep the full name like "genus_bac"
      } else if (grepl("_fun$", correct_source)) {
        correct_rank <- correct_source  # Keep the full name like "phylum_fun"
	} else {
        correct_rank <- correct_source  # For functional groups, use as-is
      }
      
      # Update the rank.name parameter with the correct value
      rank.name <- correct_rank
      message("DEBUG: Corrected rank.name to ", rank.name, " for taxon ", taxon_name)
    }
  }
  
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
  Nmc_large <- min(3000, nrow(param_samples))
  row_samples <- if (Nmc > Nmc_large) {
    sample.int(Nmc_large, Nmc, replace = TRUE)
	} else {
    sample.int(Nmc_large, Nmc, replace = FALSE)
	}

  # OPTIMIZATION 1: Pre-extract all parameters
  params <- extract_all_parameters(param_samples, row_samples, model_name)
  
  # Handle site effects for new vs observed sites
	if (is_new_site) {
    if (!is.null(predict_site_effects) && siteID %in% predict_site_effects$siteID) {
      # Use predicted site effects
      site_effect_val <- predict_site_effects$fit[predict_site_effects$siteID == siteID]
      params$site_effects <- rep(site_effect_val, Nmc)
		} else {
			# Sample from site effect variance
			site_effect_sd_cols <- grep("site_effect_sd", colnames(param_samples))
      if (length(site_effect_sd_cols) > 0) {
        site_effect_sd <- mean(param_samples[row_samples, site_effect_sd_cols], na.rm = TRUE)
        params$site_effects <- rnorm(Nmc, 0, site_effect_sd)
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
      params$site_effects <- as.numeric(param_samples[row_samples, site_param])
    }
  }
  
  # OPTIMIZATION 2: Get cached covariates
  covar <- get_site_covariates_cached(siteID, model.inputs, model_name, plotID, Nmc_large, Nmc)
  
  # Select covariates based on model type
	if (model_name == "cycl_only") {
    covar <- covar[, 7:8, ]  # sin_mo, cos_mo
	} else if (model_name == "env_cov") {
    covar <- covar[, 1:6, ]   # temp, mois, pH, pC, relEM, LAI
	} else if (model_name == "env_cycl") {
    covar <- covar[, 1:8, ]   # all 8
  }
  
  # Set initial conditions and start date using model.inputs
  plot_start_date <- model.inputs$plot_start[plotID]
  if (is.na(plot_start_date)) {
    # Fallback if plot_start not found
    plot_start_date <- if (is_new_site) 1 else 56
  }
  
  if (is_new_site) {
    # Use taxon-specific initial condition
    observed_data <- truth_data %>% 
      filter(species == taxon_name, !is.na(truth)) %>% 
      pull(truth)
    ic_mean <- if (length(observed_data) > 0) mean(observed_data) else 0.5
    ic_sd <- if (length(observed_data) > 0) sd(observed_data) else 0.2
    ic_sd <- max(min(ic_sd, 0.3), 0.05)
  } else {
    # For observed sites: EXTRACT FROM PLOT_SUMMARY
    if (!is.null(plot_summary)) {
      # Convert plot_summary to expected format if needed
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
        if (!is.null(metadata) && "model_data" %in% names(metadata)) {
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
            
            # Start from calibration end for hindcast
            # plot_start_date already set from model.inputs$plot_start[plotID]
          } else {
            # Fallback if no valid estimates
            ic_mean <- 0.5
            ic_sd <- 0.1
            # plot_start_date already set from model.inputs$plot_start[plotID]
          }
        } else {
          # No plot summary available, use fallback
          ic_mean <- 0.5
          ic_sd <- 0.1
          # plot_start_date already set from model.inputs$plot_start[plotID]
        }
      } else {
        # No plot summary, use fallback
        ic_mean <- 0.5
        ic_sd <- 0.1
        # plot_start_date already set from model.inputs$plot_start[plotID]
      }
    } else {
      # No plot_summary provided at all
      ic_mean <- 0.5
      ic_sd <- 0.1
      # plot_start_date already set from model.inputs$plot_start[plotID]
    }
  }
  
  ic <- truncnorm::rtruncnorm(Nmc, mean = ic_mean, sd = ic_sd, a = 0, b = 1)
  
  # Define timepoints to forecast
  valid_timepoints <- plot_start_date:model.inputs$N.date
  
  # OPTIMIZATION 3: Vectorized Monte Carlo simulation
  all_predictions <- vectorized_forecast(
		params = params,
		covar = covar,
		initial_conditions = ic,
		timepoints = valid_timepoints,
		Nmc = Nmc
	)

	# Create output dataframe
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
  
  colnames(ci)[1:5] <- c("lo", "lo_25", "med", "hi_75", "hi")
  
  # Add date information
  ci$date_num <- valid_timepoints
  
  date_sequence <- seq.Date(from = as.Date("2013-06-01"), by = "month", 
                            length.out = model.inputs$N.date)
  date_key <- data.frame(
		date_num = 1:model.inputs$N.date,
    dateID = as.numeric(format(date_sequence, "%Y%m"))
  )
  
  ci <- left_join(ci, date_key, by = "date_num") %>%
    mutate(dates = as.Date(paste0(dateID, "01"), format = "%Y%m%d"))
  
  # Add truth values - ensure dateID format matches
  plot_obs <- truth_data %>% filter(plotID == !!plotID)
  if (nrow(plot_obs) > 0) {
    # Debug output
    cat("    DEBUG: Found", nrow(plot_obs), "truth observations for", plotID, "\n")
    cat("    DEBUG: Truth dateID range:", range(plot_obs$dateID), "\n")
    cat("    DEBUG: Prediction dateID range:", range(ci$dateID), "\n")
    
    # Ensure dateID is numeric in both datasets
    truth_subset <- plot_obs %>%
      mutate(dateID = as.numeric(as.character(dateID))) %>%
      group_by(dateID) %>%
      summarise(truth = mean(truth, na.rm = TRUE), .groups = "drop")
    
    cat("    DEBUG: Aggregated truth data:", nrow(truth_subset), "rows\n")
    cat("    DEBUG: Aggregated truth dateIDs:", paste(truth_subset$dateID, collapse = ", "), "\n")
    
    # Ensure ci$dateID is also numeric before joining
    ci <- ci %>%
      mutate(dateID = as.numeric(as.character(dateID))) %>%
      left_join(truth_subset, by = "dateID")
    
    # Debug join results
    truth_matches <- sum(!is.na(ci$truth))
    cat("    DEBUG: Truth values successfully joined:", truth_matches, "out of", nrow(ci), "predictions\n")
  } else {
    cat("    DEBUG: No truth observations found for", plotID, "\n")
    ci$truth <- NA
  }
  
  # Note: Raw observation data (sampleID level) is available in model.inputs$sample_values
  # but is not included in the main output to avoid breaking downstream rbind operations
  # Raw data can be accessed separately if needed for analysis
  
  # Add metadata and start date information
  ci <- ci %>% mutate(
    model_id = model_id,
    model_name = model_name,
    fcast_period = if (is_new_site) {
      # For unobserved sites, entire period is hindcast (2013-2020)
      "hindcast"
    } else {
      # For observed sites, split by date
      ifelse(dates <= as.Date("2018-01-01"), "calibration", "hindcast")
    },
    plot_start = plot_start_date,  # Already defined earlier in function
    site_start = ifelse(is_new_site, 1, 
                       ifelse(siteID %in% names(model.inputs$site_start), 
                              model.inputs$site_start[siteID], NA))
  )
  
  return(ci)
}

# =============================================================================
# DIAGNOSTIC FUNCTIONS
# =============================================================================

#' Generate diagnostic visualizations for hindcast data
#' @param hindcast_data Data frame with hindcast predictions
#' @param model_id Model identifier
#' @param taxon Taxon name
#' @export
generate_hindcast_diagnostics <- function(hindcast_data, model_id, taxon) {
  cat("Creating diagnostic visualizations for", taxon, "hindcasts...\n")
  
  # Get unique sites and select one plot from each site
  unique_sites <- unique(hindcast_data$siteID)
  cat("Found", length(unique_sites), "unique sites for diagnostics\n")
  
  selected_plots <- hindcast_data %>%
    group_by(siteID) %>%
    slice_head(n = 1) %>%
    pull(plotID)
  
  cat("Selected", length(selected_plots), "plots (one per site) for diagnostics\n")
  
  # Create visualization directory based on site type
  # Check if this is for observed or unobserved sites based on the data
  unobserved_sites <- c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", 
                        "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", 
                        "WREF", "YELL")
  
  # Determine if this is for observed or unobserved sites
  unique_sites_in_data <- unique(hindcast_data$siteID)
  is_unobserved_data <- any(unique_sites_in_data %in% unobserved_sites)
  
  site_type <- if (is_unobserved_data) "unobserved_sites" else "observed_sites"
  viz_dir <- here("figures", "hindcast_diagnostics", site_type, model_id)
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Generate plots for each selected plot
  for (i in seq_along(selected_plots)) {
    plot_id <- selected_plots[i]
    site_id <- unique_sites[i]
    
    # Filter data for this plot - keep ALL data, not just non-NA predictions
    plot_data <- hindcast_data %>%
      filter(plotID == plot_id) %>%
      mutate(plot_date = as.Date(dates))
    
    if (nrow(plot_data) == 0) {
      cat("  Skipping", plot_id, "- no data\n")
      next
    }
    
    # Separate calibration and hindcast data for plotting
    calibration_data <- plot_data %>%
      filter(fcast_period == "calibration" & !is.na(med) & !is.na(lo) & !is.na(hi))
    
    hindcast_period_data <- plot_data %>%
      filter(fcast_period == "hindcast" & !is.na(med) & !is.na(lo) & !is.na(hi))
    
    if (nrow(calibration_data) == 0 && nrow(hindcast_period_data) == 0) {
      cat("  Skipping", plot_id, "- no valid predictions in either period\n")
      next
    }
    
    # For unobserved sites, limit plot data to actual hindcast period (2018+) to reduce plot density
    if (nrow(calibration_data) == 0 && nrow(hindcast_period_data) > 0) {
      # This is an unobserved site - limit to hindcast period only
      plot_data <- plot_data %>%
        filter(plot_date >= as.Date("2018-01-01"))
      
      # Re-filter hindcast_period_data from the updated plot_data to avoid duplicates
      hindcast_period_data <- plot_data %>%
        filter(fcast_period == "hindcast" & !is.na(med) & !is.na(lo) & !is.na(hi))
    }
    
    cat("  Plotting", plot_id, "with", nrow(calibration_data), "calibration predictions and", nrow(hindcast_period_data), "hindcast predictions out of", nrow(plot_data), "total timepoints\n")
    
    # Create time series plot with both calibration and hindcast periods
    p1 <- ggplot(plot_data, aes(x = plot_date)) +
      # Add calibration period confidence intervals and predictions
      {if(nrow(calibration_data) > 0) {
        list(
          geom_ribbon(data = calibration_data, aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightblue"),
          geom_line(data = calibration_data, aes(y = med), color = "blue", linewidth = 1),
          geom_point(data = calibration_data, aes(y = med), color = "blue", size = 1, alpha = 0.7)
        )
      }} +
      # Add hindcast period confidence intervals and predictions
      {if(nrow(hindcast_period_data) > 0) {
        list(
          geom_ribbon(data = hindcast_period_data, aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "lightgreen"),
          geom_line(data = hindcast_period_data, aes(y = med), color = "green", linewidth = 1),
          geom_point(data = hindcast_period_data, aes(y = med), color = "green", size = 1, alpha = 0.7)
        )
      }}
    
    # Add truth values if available (for both periods)
    all_valid_data <- bind_rows(calibration_data, hindcast_period_data)
    if(sum(!is.na(all_valid_data$truth)) > 0) {
      p1 <- p1 + 
        geom_point(data = all_valid_data, aes(y = truth), color = "red", size = 2, alpha = 0.8) +
        geom_line(data = all_valid_data, aes(y = truth), color = "red", linewidth = 1, alpha = 0.8)
    }
    
    # Add vertical line to separate calibration and hindcast periods
    if(nrow(calibration_data) > 0 && nrow(hindcast_period_data) > 0) {
      # Find the boundary date between calibration and hindcast
      boundary_date <- max(calibration_data$dateID)
      p1 <- p1 + geom_vline(xintercept = boundary_date, linetype = "dashed", color = "gray", linewidth = 1)
    }
    
    p1 <- p1 +
      labs(title = paste("Hindcast Predictions vs Truth for", plot_id, "(", site_id, ")"),
           subtitle = paste("Taxon:", taxon, "| Model:", model_id, "|", 
                           nrow(calibration_data), "calibration predictions |", 
                           nrow(hindcast_period_data), "hindcast predictions |", 
                           sum(!is.na(all_valid_data$truth)), "truth values"),
           x = "Date", y = "Abundance",
           caption = "Blue: Calibration period | Green: Hindcast period | Red: Observed truth") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Create prediction distribution plot (combine both periods)
    p2 <- ggplot(all_valid_data, aes(x = med, fill = fcast_period)) +
      geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("calibration" = "lightblue", "hindcast" = "lightgreen")) +
      labs(title = paste("Distribution of Predictions for", plot_id),
           x = "Predicted Abundance", y = "Count", fill = "Period") +
      theme_minimal()
    
    # Create uncertainty plot (combine both periods)
    all_valid_data$uncertainty <- all_valid_data$hi - all_valid_data$lo
    p3 <- ggplot(all_valid_data, aes(x = dateID, y = uncertainty, color = fcast_period)) +
      geom_line(linewidth = 1) +
      geom_point(size = 1) +
      scale_color_manual(values = c("calibration" = "blue", "hindcast" = "green")) +
      labs(title = paste("Prediction Uncertainty for", plot_id),
           x = "Date", y = "Uncertainty (95% CI width)", color = "Period") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Add vertical line to separate periods in uncertainty plot
    if(nrow(calibration_data) > 0 && nrow(hindcast_period_data) > 0) {
      boundary_date <- max(calibration_data$dateID)
      p3 <- p3 + geom_vline(xintercept = boundary_date, linetype = "dashed", color = "gray", linewidth = 1)
    }
    
    # Create validation plot (predictions vs truth) if truth values are available
    if(sum(!is.na(all_valid_data$truth)) > 0) {
      p4 <- ggplot(all_valid_data, aes(x = truth, y = med, color = fcast_period)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "dashed") +
        scale_color_manual(values = c("calibration" = "blue", "hindcast" = "green")) +
        labs(title = paste("Predictions vs Truth for", plot_id),
             x = "Observed Truth", y = "Predicted Abundance", color = "Period") +
        theme_minimal()
      
      # Combine plots with validation plot
      combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1, 1, 1))
		} else {
      # Combine plots without validation plot
      combined_plot <- gridExtra::grid.arrange(p1, p2, p3, ncol = 1, heights = c(2, 1, 1))
    }
    
    # Save plot
    plot_file <- file.path(viz_dir, paste0("hindcast_", plot_id, "_", taxon, ".png"))
    ggsave(plot_file, combined_plot, width = 12, height = 10, dpi = 300)
    
    cat("  Saved visualization for", plot_id, "to", plot_file, "\n")
  }
  
  cat("Diagnostic visualizations saved to:", viz_dir, "\n")
}

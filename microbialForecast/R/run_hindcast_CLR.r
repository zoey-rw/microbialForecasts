# Function to forecast functional groups at all NEON sites using CLR models
# Returns a dataframe with CIs and observed truth values (plot means)

#' @title fcast_clr
#' @description Forecast functional groups at NEON plots and sites using CLR model output
#' @export
fcast_clr <- function(plotID,
                      model.inputs,
                      param_samples,
                      truth.plot.long,
                      plot_summary,
                      Nmc = 1000,
                      drop_other = T,
                      predict_site_effects = NULL,
                      rank.name=NULL,
                      model_id,
                      ...) {

  siteID <- substr(plotID, 1, 4)

  # Debug: Check if site exists in model inputs
  if (!siteID %in% names(model.inputs$site_start)) {
    message("Warning: Site ", siteID, " not found in model.inputs$site_start")
    message("Available sites: ", paste(names(model.inputs$site_start), collapse = ", "))
    return(NULL)
  }

  # Prep MCMC sampling IDs
  # Initial condition uncertainty
  # For CLR models, we work with log-ratio values, so initial conditions should be around 0
  if (!exists("ic")) {
    # This happens for sites not in original model data
    # For CLR models, use small values around 0 since we're working with log-ratios
    ic <- rnorm(Nmc, mean = 0, sd = 0.1)
  }
  
  # Validate initial condition
  if (any(is.na(ic)) || any(is.infinite(ic))) {
    message("Warning: Invalid initial condition for plot ", plotID, " - skipping")
    return(NULL)
  }

  NT = model.inputs$N.date
  Nmc_large <- max(nrow(param_samples)) # Larger sample number for covariate/IC set of values
  
  # Handle case where we want more samples than available
  if (Nmc > Nmc_large) {
    # If we want more samples than available, sample with replacement
    row_samples <- sample.int(Nmc_large, Nmc, replace = TRUE)
    cat("Warning: Requested", Nmc, "samples but only", Nmc_large, "available. Sampling with replacement.\n")
  } else {
    # If we have enough samples, sample without replacement
    row_samples <- sample.int(Nmc_large, Nmc, replace = FALSE)
  }

  date_key <- model.inputs$truth.plot.long %>%
    select(dateID, date_num) %>% distinct()

  # Sample covariate data - use the existing function since CLR models use the same covariates
  tryCatch({
    covar <- create_covariate_samples(model.inputs, plotID, siteID,
                                     Nmc_large, Nmc)
  }, error = function(e) {
    if (grepl("missing value where TRUE/FALSE needed", e$message)) {
      message("Error in forecast for plot ", plotID, ": missing value where TRUE/FALSE needed")
      message("This usually means site ", siteID, " is missing from model.inputs$site_start")
      return(NULL)
    } else {
      stop(e)  # Re-throw other errors
    }
  })
  
  # Check if covar is NULL (error occurred)
  if (is.null(covar)) {
    return(NULL)
  }
  
  plot_obs <- model.inputs$truth.plot.long %>%
    filter(plotID==!!plotID) %>%
    select(-c(plot_num,site_num))

  taxon_names <- model.inputs$truth.plot.long %>% select(species) %>% distinct() %>% unlist()
  taxon_name <- taxon_names[1]

  # Check whether there's already an estimated site effect. If not, we'll sample!
  is_new_site <- ifelse(siteID %in% truth.plot.long$siteID, FALSE, TRUE)

  if (is_new_site) {
    # For new sites, we need to predict site effects
    if (is.null(predict_site_effects)) {
      message("Warning: No predicted site effects provided for new site ", siteID, " - skipping")
      return(NULL)
    }
    
    # Get predicted site effect for this site
    site_effect_pred <- predict_site_effects %>%
      filter(siteID == !!siteID, model_id == !!model_id)
    
    if (nrow(site_effect_pred) == 0) {
      message("Warning: No predicted site effect found for site ", siteID, " and model ", model_id, " - skipping")
      return(NULL)
    }
    
    # Sample from predicted site effect distribution
    site_effect <- rnorm(Nmc, mean = site_effect_pred$fit, sd = site_effect_pred$se.fit)
  } else {
    # For existing sites, get the estimated site effect
    site_effect <- param_samples[row_samples,] %>% select(grep("site_effect", colnames(.))) %>% unlist()
    
    # Filter to the specific site effect for this plot
    plot_site_num <- model.inputs$plot_site_num[which(model.inputs$plotID == plotID)]
    if (length(plot_site_num) == 0) {
      message("Warning: No plot_site_num found for plot ", plotID, " - skipping")
      return(NULL)
    }
    
    # Get the site effect for this specific site
    site_effect_col <- paste0("site_effect[", plot_site_num, "]")
    site_effect <- param_samples[row_samples, site_effect_col]
    
    if (length(site_effect) == 0) {
      message("Warning: No site effect found!")
      return(NULL)  # Skip this plot entirely
    }
  }
  
  # Check for valid site effect values
  if (any(is.na(site_effect))) {
    message("Warning: NA values in site effects for plot ", plotID, " - skipping")
    return(NULL)
  }

  ### Get other parameter estimates for CLR models
  ### Betas (for CLR models, these are the coefficients for predictors)
  betas <- param_samples[row_samples,] %>% select(grep("beta", colnames(.)))
  ### Intercept
  intercept <- param_samples[row_samples,] %>% select(grep("intercept", colnames(.))) %>% unlist()
  ### Process error (sigma for CLR models)
  sigma <- param_samples[row_samples,] %>% select(grep("sigma", colnames(.))) %>% unlist()

  # Validate that we have all required parameters
  if (length(betas) == 0 || length(intercept) == 0 || length(sigma) == 0) {
    message("Warning: Missing required parameters for plot ", plotID, " - skipping")
    return(NULL)
  }
  
  # Check for valid parameter values
  if (any(is.na(betas)) || any(is.na(intercept)) || any(is.na(sigma))) {
    message("Warning: NA values in parameters for plot ", plotID, " - skipping")
    return(NULL)
  }

  # For CLR models, we need to handle different model types
  # Extract model type from model_id
  if (grepl("^cycl_only", model_id)) {
    # Seasonal only model - use sin/cos predictors
    covar <- covar[, c(7, 8), ]  # sin_mo, cos_mo
  } else if (grepl("^env_cov", model_id)) {
    # Environmental only model - use environmental predictors
    covar <- covar[, c(1:6), ]   # temperature, moisture, pH, carbon, relative EM, LAI
  } else if (grepl("^env_cycl", model_id)) {
    # Environmental + seasonal model - use all predictors
    # covar remains unchanged (all 8 predictors)
  }

  all_tax_abs <- matrix(NA, Nmc, NT)
  # Initialize first timepoint with initial condition
  plot_start_date <- model.inputs$plot_start[which(model.inputs$plotID == plotID)]
  if (length(plot_start_date) == 0) {
    plot_start_date <- 1
  }
  
  all_tax_abs[, plot_start_date] <- ic
  ## simulate
  x <- ic

  # Only forecast for valid timepoints (where we have data or can reasonably predict)
  valid_timepoints <- plot_start_date:NT
  
  # Ensure we don't try to forecast beyond available data
  if (max(valid_timepoints) > NT) {
    valid_timepoints <- valid_timepoints[valid_timepoints <= NT]
  }

  for (time in valid_timepoints) {
    if (time == plot_start_date) {
      # Skip the initial condition timepoint
      next
    }
    
    Z <- covar[, , time]
    
    # For CLR models, we use a different prediction equation
    # CLR models predict log-ratios directly, not through logit transformation
    if (grepl("^cycl_only", model_id)) {
      # Seasonal only: Ex = intercept + beta[1]*sin_mo + beta[2]*cos_mo + site_effect
      Ex <- intercept + apply(Z * betas, 1, sum) + site_effect
    } else if (grepl("^env_cov", model_id)) {
      # Environmental only: Ex = intercept + sum(beta[i]*env_predictor[i]) + site_effect
      Ex <- intercept + apply(Z * betas, 1, sum) + site_effect
    } else if (grepl("^env_cycl", model_id)) {
      # Environmental + seasonal: Ex = intercept + sum(beta[i]*predictor[i]) + site_effect
      Ex <- intercept + apply(Z * betas, 1, sum) + site_effect
    } else {
      # Default: assume all predictors
      Ex <- intercept + apply(Z * betas, 1, sum) + site_effect
    }

    # For CLR models, we add process error directly to the predicted log-ratio
    # The process error represents uncertainty in the log-ratio prediction
    x <- Ex + rnorm(Nmc, mean = 0, sd = sigma)

    # Save to array
    all_tax_abs[, time] <- x
  }

  # Create output dataframe
  predict <- all_tax_abs
  ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025, .25, 0.5, .75, 0.975), na.rm=T))) %>%
    mutate(mean = apply(predict, 2, mean, na.rm=T),
           sd = apply(predict, 2, sd, na.rm=T),
           date_num = as.numeric(1:NT),
           plotID = plotID,
           siteID = siteID,
           taxon_name = taxon_name,
           species = taxon_name,
           taxon = taxon_name,
           new_site = ifelse(is_new_site, T, F))
  colnames(ci)[1:5] <- c("lo","lo_25","med","hi_75","hi")

  ci <- left_join(ci, date_key, by=c("date_num"))
  ci$dates <- fixDate(ci$dateID)
  ci <- ci %>% filter(!is.na(med))

  coalesce_by_column <- function(df) {
    return(coalesce(df[1], df[2]))
  }

  # Add on truth values and/or model estimates
  if (!is_new_site){
    plot_obs_simple = plot_obs %>%
      select(siteID, plotID, dateID, species,
             truth, lo, lo_25, med, hi_75, hi,
             mean = Mean, sd = SD) %>% mutate(truth=as.numeric(truth))
    ci <- full_join(ci, plot_obs_simple, by = intersect(colnames(ci), colnames(plot_obs_simple)))
    ci = ci %>%
      group_by(dateID) %>%
      summarise_all(coalesce_by_column)
  } else {
    plot_obs_simple = plot_obs %>%
      select(siteID, plotID, dateID, species,
             truth) %>% mutate(truth=as.numeric(truth))
    ci <- full_join(ci, plot_obs_simple, by = intersect(colnames(ci), colnames(plot_obs_simple)))
    ci = ci %>%
      group_by(dateID) %>%
      summarise_all(coalesce_by_column)
  }

  ci$date_num <- coalesce(ci$date_num, date_key$date_num)

  obs_date_num = ci[which(!is.na(ci$truth) & !is.na(ci$date_num)),]$date_num
  obs_date_num = obs_date_num[which(obs_date_num > plot_start_date)]
  
  if (length(obs_date_num) > 0){
    crps_df = data.frame()
    for (i in obs_date_num){
      obs = ci[which(ci$date_num==i),]$truth
      # Ensure obs is numeric for CRPS calculation
      obs = as.numeric(obs)
      X = predict[,i]
      X = X[complete.cases(X)]

      # Skip CRPS calculation if no valid forecast samples
      if (length(X) == 0) {
        message("Warning: No valid forecast samples for timepoint ", i, " - skipping CRPS calculation")
        next
      }

      crps_val = cbind(crps = crps_sample(obs, X), date_num = i)
      crps_df = rbind(crps_df, crps_val)
    }

    # Only join CRPS data if there are results to join
    if (nrow(crps_df) > 0) {
      ci <- left_join(ci, crps_df, by="date_num")
    }

    return(ci)
  } # End of if (length(obs_date_num) > 0) block
  
  # If no valid observations, return the CI without CRPS
  return(ci)
}

# Helper function to create covariate samples for CLR models
# This function should be adapted from the existing create_covariate_samples function
# to handle CLR model inputs properly

#' @title create_covariate_samples_clr
#' @description Create covariate samples for CLR model forecasting
#' @export
create_covariate_samples_clr <- function(model.inputs, plotID, siteID, Nmc_large, Nmc) {
  # This is a placeholder - the actual implementation would need to be
  # adapted from the existing create_covariate_samples function
  # to handle CLR model inputs properly
  
  # For now, return NULL to indicate this function needs to be implemented
  message("Warning: create_covariate_samples_clr function not yet implemented")
  return(NULL)
}

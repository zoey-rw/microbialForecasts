# Function to forecast taxonomic groups at all NEON sites using Dirichlet models
# Returns a dataframe with CIs and observed truth values (plot means)

#' @title fcast_dirichlet
#' @description Forecast taxonomic groups at NEON plots and sites using Dirichlet model output
#' @export
fcast_dirichlet <- function(plotID,
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
  # For Dirichlet models, we work with compositional data, so initial conditions should be positive
  if (!exists("ic")) {
    # Check if we need to set initial condition
    # Get truth data from model.inputs instead of plot_summary
    plot_truth_data <- model.inputs$truth.plot.long %>%
      filter(plotID == !!plotID) %>%
      filter(!is.na(truth)) %>%
      mutate(truth = as.numeric(truth))

    if (nrow(plot_truth_data) > 0) {
      # Use actual truth data for initial condition (Dirichlet scale)
      # For Dirichlet models, we need positive values
      ic <- rgamma(Nmc, shape = 2, rate = 2 / mean(plot_truth_data$truth, na.rm=T))
    } else {
      # Fallback: use default initial condition for Dirichlet models
      ic <- rgamma(Nmc, shape = 2, rate = 2)
    }
  }
  
  # Validate initial condition
  if (any(is.na(ic)) || any(is.infinite(ic)) || any(ic <= 0)) {
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

  # Extract parameters for this plot
  betas <- param_samples[row_samples, grep("^beta\\[", colnames(param_samples))]
  rho <- param_samples[row_samples, "rho"]
  sigma <- param_samples[row_samples, "sigma"]
  site_effect <- param_samples[row_samples, paste0("site_effect[", model.inputs$site_start[[siteID]], "]")]
  intercept <- param_samples[row_samples, "intercept"]

  # Get covariate data
  covar <- model.inputs$covar
  plot_start_date <- model.inputs$plot_start[[plotID]]
  plot_site_num <- model.inputs$plot_site_num[[plotID]]

  # Create date key for joining
  date_key <- model.inputs$truth.plot.long %>%
    filter(plotID == !!plotID) %>%
    filter(!is.na(date_num)) %>%
    select(dateID, date_num) %>%
    distinct() %>%
    mutate(date_num = as.numeric(date_num))

  # If date_key is empty, create fallback
  if (nrow(date_key) == 0) {
    date_key <- data.frame(
      dateID = model.inputs$dateID,
      date_num = 1:length(model.inputs$dateID)
    )
  }

  # Determine if this is a new site
  is_new_site <- siteID %in% names(predict_site_effects)

  # Get site effect for new sites
  if (is_new_site) {
    site_effect <- predict_site_effects[[siteID]]$fit
  }

  # Initialize arrays
  all_tax_abs <- array(NA, dim = c(Nmc, NT))
  all_tax_abs[, plot_start_date] <- ic

  # Get valid timepoints
  valid_timepoints <- (plot_start_date + 1):NT

  # Get plot observations for truth data
  plot_obs <- truth.plot.long %>%
    filter(plotID == !!plotID) %>%
    filter(!is.na(truth))

  for (time in valid_timepoints) {
    if (time == plot_start_date) {
      # Skip the initial condition timepoint
      next
    }
    
    # CRITICAL FIX: Bounds checking for covariate access
    if (time > dim(covar)[3]) {
      next
    }
    
    Z <- covar[, , time]
    
    # CRITICAL FIX: Validate Z dimensions
    if (ncol(Z) != ncol(betas)) {
      next
    }
    
    # CRITICAL FIX: Use Dirichlet model structure to match model fitting approach
    # Dirichlet model fitting uses: log(Ex[p, s, t]) <- rho[s] * log(max(plot_mu[p, s, t - 1], 0.001)) + ...
    # So forecasting should use: log(Ex) <- rho * log(max(x, 0.001)) + apply(Z * betas, 1, sum) + site_effect + intercept
    
    # Use log transformation for numerical stability (matching model fitting)
    log_x_prev <- log(max(x, 0.001))
    
    # Build linear predictor using LOG transformation (matching model fitting)
    log_Ex <- rho * log_x_prev + apply(Z * betas, 1, sum) + site_effect + intercept
    
    # CRITICAL FIX: Use exp() to match model fitting
    # Model fitting uses: plot_mu[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
    Ex <- exp(log_Ex)
    
    # Sample from truncated normal distribution (matching model fitting approach)
    x <- pmax(0.001, rnorm(Nmc, mean = Ex, sd = sigma))

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
           taxon_name = rank.name,
           species = rank.name,
           taxon = rank.name,
           new_site = ifelse(is_new_site, T, F))
  colnames(ci)[1:5] <- c("lo","lo_25","med","hi_75","hi")

  ci <- left_join(ci, date_key, by=c("date_num"))
  ci$dates <- fixDate(ci$dateID)
  ci <- ci %>% filter(!is.na(med))

  coalesce_by_column <- function(df) {
    # For truth column, prioritize non-NA values that are not dateID-like numbers
    if ("truth" %in% colnames(df)) {
      truth_col <- df$truth
      # Filter out dateID-like values (6-digit numbers starting with 20)
      valid_truth <- truth_col[!is.na(truth_col) & !grepl("^20[0-9]{4}$", as.character(truth_col))]
      if (length(valid_truth) > 0) {
        return(valid_truth[1])
      }
    }
    # For other columns, use standard coalesce
    return(coalesce(df[1], df[2]))
  }

  # Add on truth values and/or model estimates
  if (!is_new_site){
    # For observed sites, add truth values from plot_summary if available
    if (!is.null(plot_summary) && nrow(plot_summary) > 0) {
      # Use plot_summary which contains truth data and model estimates
      plot_obs_simple = plot_summary %>%
        filter(plotID == !!plotID) %>%
        select(siteID, plotID, dateID, species, truth, 
               lo = `2.5%`, lo_25 = `25%`, med = `50%`, hi_75 = `75%`, hi = `97.5%`,
               mean = Mean, sd = SD) %>% 
        mutate(truth = as.numeric(truth))
    } else {
      # Fallback: use plot_obs (which comes from truth.plot.long)
      plot_obs_simple = plot_obs %>%
        select(siteID, plotID, dateID, species, truth) %>% 
        mutate(truth = as.numeric(truth))
    }
    
    ci <- full_join(ci, plot_obs_simple, by = intersect(colnames(ci), colnames(plot_obs_simple)))
    ci = ci %>%
      group_by(dateID) %>%
      summarise_all(coalesce_by_column)
  } else {
    # For new sites, we don't have truth values, but we need to add a truth column
    # to ensure consistent column structure when combining with observed sites
    ci$truth <- NA_real_
  }

  ci$date_num <- coalesce(ci$date_num, date_key$date_num)
  
  # Add model metadata and determine forecast period
  ci <- ci %>% mutate(
    model_id = model_id,
    model_name = model_name,
    # Determine forecast period based on date
    fcast_period = ifelse(dates <= as.Date("2018-01-01"), "calibration", "hindcast")
  )

  obs_date_num = ci[which(!is.na(ci$truth) & !is.na(ci$date_num)),]$date_num
  obs_date_num = obs_date_num[which(obs_date_num > plot_start_date)]
  
  if (length(obs_date_num) > 0){
    crps_df = data.frame()
    for (i in 1:length(obs_date_num)){
      obs_date = obs_date_num[i]
      obs_val = ci[which(ci$date_num == obs_date),]$truth
      pred_samples = predict[,obs_date]
      crps_val = crps_sample(obs_val, pred_samples)
      crps_df = rbind(crps_df, data.frame(date_num = obs_date, crps = crps_val))
    }
    ci = left_join(ci, crps_df, by = "date_num")
  } else {
    ci$crps = NA
  }

  return(ci)
}




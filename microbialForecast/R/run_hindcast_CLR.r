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

  # Check if site exists in model inputs (for existing sites only)
  # New sites (unobserved sites) won't be in model.inputs$site_start, which is expected
  is_observed_site <- siteID %in% names(model.inputs$site_start)
  if (!is_observed_site) {
    message("Site ", siteID, " not found in model.inputs$site_start - treating as new site")
  }

  # Define plot_obs outside conditional blocks to avoid scope issues
  # CRITICAL FIX: Convert to data frame before using dplyr functions
  truth_data <- as.data.frame(model.inputs$truth.plot.long)
  plot_obs <- truth_data %>%
    filter(plotID==!!plotID) %>%
    select(-c(plot_num,site_num))

  # Prep MCMC sampling IDs
  # Initial condition will be set later in the function based on site type
  # Skip validation here - it will be done after ic is properly set

  # Create date_key from model.inputs since truth.plot.long data may be incomplete
  # Use the actual dates from model.inputs which has complete time series data
  if (is.null(truth_data) || nrow(truth_data) == 0 || 
      all(is.na(truth_data$date_num))) {
    # Fallback: create date_key from model.inputs date information
    NT <- model.inputs$N.date
    date_key <- data.frame(
      date_num = 1:NT,
      dateID = unique(truth_data$dateID)[1:NT]
    ) %>%
    filter(!is.na(dateID)) %>%
    mutate(date_num = as.numeric(date_num))
  } else {
    # Use date_key from truth_data if available
    date_key <- truth_data %>%
      select(dateID, date_num) %>% 
      distinct() %>%
      filter(!is.na(date_num)) %>%
      mutate(date_num = as.numeric(date_num))
    NT = length(unique(date_key$date_num))
  }
  Nmc_large <- max(nrow(param_samples)) # Larger sample number for covariate/IC set of values
  
  # Handle case where we want more samples than available
  if (Nmc > Nmc_large) {
    # If we want more samples than available, sample with replacement
    row_samples <- sample.int(Nmc_large, Nmc, replace = TRUE)
  } else {
    # If we have enough samples, sample without replacement
    row_samples <- sample.int(Nmc_large, Nmc, replace = FALSE)
  }

  # Parse model_id to determine model type EARLY (moved from later in function)
  model_info <- tryCatch({
    parse_model_id(model_id)
  }, error = function(e) {
    message("Warning: Error parsing model_id '", model_id, "': ", e$message)
    # Fallback parsing
    info <- strsplit(model_id, "_")[[1]]
    if (length(info) >= 2) {
      if (info[1] == "cycl" && info[2] == "only") {
        model_name <- "cycl_only"
      } else if (info[1] == "env" && info[2] == "cov") {
        model_name <- "env_cov"
      } else if (info[1] == "env" && info[2] == "cycl") {
        model_name <- "env_cycl"
      } else {
        model_name <- info[1]
      }
    } else {
      model_name <- "unknown"
    }
    list(model_name = model_name)
  })
  
  # Extract model name safely
  if (is.list(model_info) && length(model_info) >= 6) {
    # parse_model_id succeeded and returned a list with 9 elements (including driver uncertainty flag)
    model_name <- model_info[[6]]
    has_driver_uncertainty <- if (length(model_info) >= 9) model_info[[9]] else FALSE
  } else if (is.list(model_info) && "model_name" %in% names(model_info)) {
    # Fallback parsing was used and returned a list with model_name
    has_driver_uncertainty <- FALSE
    model_name <- model_info$model_name
  } else {
    # Final fallback: try to parse from model_id directly
    info <- strsplit(model_id, "_")[[1]]
    if (length(info) >= 2) {
      if (info[1] == "cycl" && info[2] == "only") {
        model_name <- "cycl_only"
      } else if (info[1] == "env" && info[2] == "cov") {
        model_name <- "env_cov"
      } else if (info[1] == "env" && info[2] == "cycl") {
        model_name <- "env_cycl"
      } else {
        model_name <- info[1]
      }
    } else {
      model_name <- "unknown"
    }
  }
  
  # Debug: ensure model_name is set
  if (!exists("model_name") || is.null(model_name)) {
    model_name <- "unknown"
  }
  
  # Store model_name in a more persistent way to avoid scope issues
  current_model_name <- model_name
  
  #Sample covariate data - FAIL FAST: This will stop execution if data is missing
  tryCatch({
    covar <- create_covariate_samples(model.inputs, plotID, siteID,
                                      Nmc_large, Nmc, model_type = current_model_name)
  }, error = function(e) {
    message("ERROR in create_covariate_samples for plot ", plotID, " site ", siteID, ":")
    message("  ", e$message)
    message("  This indicates missing or invalid data that needs investigation.")
    message("  Check:")
    message("    1. Site mapping between calibration and validation periods")
    message("    2. Array dimensions and rownames")
    message("    3. Data quality (NA, infinite values)")
    message("    4. Time period compatibility")
    stop("Hindcasting stopped due to data issues - investigate and fix root cause")
  })
  
  # Validation - covar should always be valid now since function fails fast
  if (!is.array(covar) || length(dim(covar)) != 3) {
    stop("CRITICAL ERROR: create_covariate_samples returned invalid array - this should never happen!")
  }
  
  # Validate we have enough covariates for the model type
  expected_covars <- if (current_model_name == "cycl_only") 2 else if (current_model_name == "env_cov") 6 else 8
  if (dim(covar)[2] < expected_covars) {
    stop("CRITICAL ERROR: covar has ", dim(covar)[2], " covariates but model expects ", expected_covars, " for plot ", plotID)
  }

  # Use the already converted truth_data
  taxon_names <- truth_data %>% select(species) %>% distinct() %>% unlist()
  taxon_name <- taxon_names[1]  # Use first taxon name for output

  # Check whether there's already an estimated site effect (following beta regression pattern)
  # A site is "new" if it doesn't have an estimated site effect from the calibration data
  # This is determined by whether the site is in the original model.inputs$site_start
  is_new_site <- !is_observed_site
  # Site type determined

  if (is_new_site) {
    # For new sites, we need to predict site effects
    if (!is.null(predict_site_effects)) {
      # Use predicted site effects if available
      site_effect <- predict_site_effects$fit
      # Ensure we have the right number of samples
      if (length(site_effect) == 1) {
        site_effect <- rep(site_effect, Nmc)
      } else if (length(site_effect) != Nmc) {
        # Sample from the available predictions
        site_effect <- sample(site_effect, Nmc, replace = TRUE)
      }
    } else {
      # Sample from site effect variance
      site_effect_sd <- param_samples[row_samples,] %>% select(grep("site_effect_sd", colnames(.))) %>% unlist()
      site_effect_sd <- mean(site_effect_sd, na.rm = TRUE)
      site_effect <- rnorm(Nmc, 0, site_effect_sd)
    }

    # For new sites, use the proper start date calculated by prepCLRData
    # prepCLRData with full_timeseries=TRUE should have calculated plot_start for all sites
    plot_start_date <- model.inputs$plot_start[plotID]
    if (is.na(plot_start_date)) {
      message("Warning: No plot_start found for plot ", plotID, " - this should not happen with full_timeseries=TRUE")
      plot_start_date <- 1  # Fallback only if prepCLRData failed
    }
    
    # For new sites, use taxon-specific initial condition based on observed sites
    # Calculate mean and sd of this taxon across all observed sites
    observed_taxon_data <- truth.plot.long %>%
      filter(species == taxon_name) %>%
      filter(!is.na(truth)) %>%
      pull(truth)
    
    if (length(observed_taxon_data) > 0) {
      # Convert to CLR scale for initial condition
      truth_clr <- log(observed_taxon_data / (1 - observed_taxon_data))
      taxon_mean <- mean(truth_clr, na.rm = TRUE)
      taxon_sd <- sd(truth_clr, na.rm = TRUE)
      # Ensure sd is not too small or too large
      taxon_sd <- max(min(taxon_sd, 0.5), 0.1)
    } else {
      # Fallback if no observed data available
      taxon_mean <- 0
      taxon_sd <- 0.2
    }
    
    ic <- rnorm(Nmc, mean = taxon_mean, sd = taxon_sd)
    message("Using taxon-specific initial condition for new site plot ", plotID, " (mean=", round(taxon_mean, 3), ", sd=", round(taxon_sd, 3), ") at timepoint ", plot_start_date)

  } else {
    site_num <- unique(truth.plot.long[truth.plot.long$siteID==siteID,]$site_num)
    site_param <- paste0("site_effect[", site_num, "]")
    site_effect <- param_samples[row_samples,] %>% select(!!site_param) %>% unlist()

    # WORKAROUND: Use final observation from truth data instead of plot estimates
    # Get the final observation for this plot from the truth data
    plot_truth_data <- truth_data %>%
      filter(plotID == !!plotID) %>%
      filter(!is.na(truth)) %>%
      mutate(truth = as.numeric(truth)) %>%
      arrange(date_num)
    
    if (nrow(plot_truth_data) > 0) {
      # Use the final observation as initial condition
      final_obs <- plot_truth_data[nrow(plot_truth_data), ]
      
      # Use timepoint instead of date_num if date_num is NA
      if (!is.na(final_obs$date_num)) {
        plot_start_date <- final_obs$date_num
      } else if (!is.na(final_obs$timepoint)) {
        plot_start_date <- final_obs$timepoint
      } else {
        # Fallback to model.inputs plot_start or use row number
        plot_start_date <- model.inputs$plot_start[plotID]
        if (is.na(plot_start_date)) {
          # Use the row number as a proxy for timepoint
          plot_start_date <- nrow(plot_truth_data)
        }
      }
      
      # Convert to CLR scale for initial condition
      truth_clr <- log(final_obs$truth / (1 - final_obs$truth))
      
      # Calculate standard deviation for initial condition uncertainty
      truth_sd <- sd(plot_truth_data$truth, na.rm=T)
      if (is.na(truth_sd) || truth_sd == 0) {
        truth_sd <- 0.1  # Default small standard deviation
      }
      
      ic <- rnorm(Nmc, truth_clr, truth_sd)
      
      # Initial condition generated from final observation
      message("Using final observation from truth data for plot ", plotID, ": CLR value = ", round(truth_clr, 3), " at timepoint ", plot_start_date)
    } else {
      # Fallback: try to get from plot_summary if available
      if (!is.null(plot_summary) && nrow(plot_summary) > 0) {
        plot_est <- plot_summary %>%
          filter(plotID == !!plotID) %>%
          select(-c(plot_num, site_num, dateID, dates, truth, rank))
        
        if (nrow(plot_est) > 0) {
          plot_obs <- left_join(plot_obs, plot_est, by = intersect(colnames(plot_obs), colnames(plot_est)))
          last_obs <- plot_est %>% filter(timepoint==max(timepoint))
          plot_start_date <- last_obs$timepoint
          # Convert to CLR scale
          truth_clr <- log(last_obs$`50%` / (1 - last_obs$`50%`))
          ic <- rnorm(Nmc, truth_clr, 0.1)
          message("Using plot estimates for plot ", plotID, ": CLR value = ", round(truth_clr, 3), " at timepoint ", plot_start_date)
        } else {
          # Final fallback: use default initial condition
          plot_start_date <- model.inputs$plot_start[plotID]
          if (is.na(plot_start_date)) {
            plot_start_date <- 1
          }
          ic <- rnorm(Nmc, mean = 0, sd = 0.1)
          message("Using default initial condition for plot ", plotID, " at timepoint ", plot_start_date)
        }
      } else {
        # Final fallback: use default initial condition
        plot_start_date <- model.inputs$plot_start[plotID]
        if (is.na(plot_start_date)) {
          plot_start_date <- 1
        }
        ic <- rnorm(Nmc, mean = 0, sd = 0.1)
        message("Using default initial condition for plot ", plotID, " at timepoint ", plot_start_date)
      }
    }
  }

  # Validate initial condition after it's been set
  if (any(is.na(ic)) || any(is.infinite(ic)) || length(ic) == 0) {
    message("Warning: Invalid initial condition for plot ", plotID, " - skipping")
    return(NULL)
  }

  # CRITICAL FIX: Validate plot_start_date before using it
  if (is.na(plot_start_date) || is.null(plot_start_date) || plot_start_date < 1) {
    message("Warning: Invalid plot_start_date ", plot_start_date, " for plot ", plotID, " - using default")
    plot_start_date <- 1
  }

  # NT is already set correctly from the samples file above

  # Ensure plot_start_date is within bounds
  if (plot_start_date > NT) {
    message("Warning: plot_start_date ", plot_start_date, " exceeds N.date ", NT, " for plot ", plotID, " - adjusting")
    plot_start_date <- 1
  }

  if (length(site_effect)==0) {
    message("No site effect found!")
    return(NULL)  # Skip this plot entirely
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
    # Environmental covariates only
    covar <- covar[, c(1:6), ]  # temp, mois, pH, pC, relEM, LAI
  } else if (grepl("^env_cycl", model_id)) {
    # Environmental + cyclical: all 8 covariates
    covar <- covar[, c(1:8), ]  # temp, mois, pH, pC, relEM, LAI, sin_mo, cos_mo
  } else {
    # Default fallback: assume all covariates
    message("Warning: Unknown model type '", current_model_name, "', using all covariates")
  }

  all_tax_abs <- matrix(NA, Nmc, NT)
  # Initialize first timepoint with initial condition
  all_tax_abs[,plot_start_date] <- ic
  ## simulate
  x <- ic

  # Only forecast for valid timepoints (where we have data or can reasonably predict)
  valid_timepoints <- plot_start_date:NT
  
  # CRITICAL FIX: Ensure we don't try to forecast beyond available data
  if (max(valid_timepoints) > NT) {
    valid_timepoints <- valid_timepoints[valid_timepoints <= NT]
  }
  
  # plot_start_date bounds checking already done earlier in the function

  for (time in valid_timepoints) {
    if (time == plot_start_date) {
      # Skip the initial condition timepoint
      next
    }
    
    # CRITICAL FIX: Bounds checking for covariate access
    if (time > dim(covar)[3]) {
      message("Warning: Time ", time, " exceeds covariate array bounds for plot ", plotID, " - skipping")
      next
    }
    
    Z  <- covar[, ,time]
    
    # CRITICAL FIX: Validate Z dimensions
    if (ncol(Z) != ncol(betas)) {
      message("Warning: Covariate dimension mismatch at time ", time, " for plot ", plotID, " - skipping")
      next
    }
    
    # CRITICAL FIX: Use CLR transformation to match model fitting approach
    # For CLR models: Ex <- rho * x + apply(Z * betas, 1, sum) + site_effect + intercept
    
    # Build linear predictor using CLR transformation (matching model fitting)
    Ex <- x + apply(Z * betas, 1, sum) + site_effect + intercept
    
    # Add process error (matching model fitting approach)
    x <- Ex + rnorm(Nmc, mean = 0, sd = sigma)

    # Save to array
    all_tax_abs[,time] <- x
  }

  # Create output dataframe
  predict <- all_tax_abs
  ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,.25,0.5,.75,0.975), na.rm=T))) %>%
    mutate(mean = apply(predict, 2, mean, na.rm=T),
           sd = apply(predict, 2, sd, na.rm=T),
           plotID = plotID,
           siteID = siteID,
           taxon_name = taxon_name,
           species = taxon_name,
           taxon = taxon_name,
           taxon_name = taxon_name,
           new_site = ifelse(is_new_site, T, F))
  colnames(ci)[1:5] <- c("lo","lo_25","med","hi_75","hi")

  # Assign sequential date_num values and join with date_key
  ci$date_num <- as.numeric(1:NT)
  
  # Ensure date_key$date_num is also numeric
  date_key$date_num <- as.numeric(date_key$date_num)
  
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

  # Add model metadata and determine forecast period
  ci <- ci %>% mutate(
    model_id = model_id,
    model_name = current_model_name,
    # Determine forecast period based on date
    fcast_period = ifelse(dates <= as.Date("2018-01-01"), "calibration", "hindcast"),
    plot_start = plot_start_date,
    site_start = ifelse(is_new_site, NA, model.inputs$site_start[siteID])
  )

  return(ci)
}

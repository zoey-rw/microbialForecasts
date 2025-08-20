# Helper functions for Dirichlet models
# Adapted from the main microbialForecast package functions

#' Create initial values for Dirichlet MCMC runs
#' @param constants List of model constants
#' @param type Type of model ("tax" for taxonomic)
#' @return List of initial values
initsFun_dirichlet <- function(constants, type = "tax") {
  # For Dirichlet models, we need different initialization
  if (type == "tax") {
    # Initialize compositional data parameters
    y_init <- matrix(rep(rep(1/constants$N.spp, constants$N.spp), constants$N.core),
                     ncol = constants$N.spp, nrow = constants$N.core)
    
    # Initialize plot means for each taxon and time
    plot_mu_init <- array(rnorm(constants$N.plot * constants$N.spp * constants$N.date, 0.5, 0.1),
                          dim = c(constants$N.plot, constants$N.spp, constants$N.date))
    
    # Ensure positive values for gamma distribution
    plot_mu_init <- pmax(plot_mu_init, 0.001)
    
    # Initialize relative abundances
    plot_rel_init <- array(rep(rep(rep(1/constants$N.spp, constants$N.spp), constants$N.plot), constants$N.date),
                           dim = c(constants$N.plot, constants$N.spp, constants$N.date))
    
    # Initialize beta parameters (effect sizes)
    beta_init <- matrix(rnorm(constants$N.beta * constants$N.spp, 0, 0.1), 
                       nrow = constants$N.spp, ncol = constants$N.beta)
    
    # Initialize rho (temporal persistence)
    rho_init <- rep(0.5, constants$N.spp)  # Start at 0.5 for stability
    
    # Initialize sigma (process error)
    sigma_init <- rep(0.1, constants$N.spp)
    
    # Initialize intercept
    intercept_init <- rep(0, constants$N.spp)
    
    # Initialize site effects
    site_effect_init <- matrix(rnorm(constants$N.site * constants$N.spp, 0, 0.1),
                              nrow = constants$N.site, ncol = constants$N.spp)
    
    # Initialize sig (site effect scale)
    sig_init <- 0.1
    
    # Initialize expected values
    Ex_init <- plot_mu_init
    
    return(list(
      y = y_init,
      plot_mu = plot_mu_init,
      plot_rel = plot_rel_init,
      beta = beta_init,
      rho = rho_init,
      sigma = sigma_init,
      intercept = intercept_init,
      site_effect = site_effect_init,
      sig = sig_init,
      Ex = Ex_init
    ))
  } else {
    stop("Only 'tax' type supported for Dirichlet models")
  }
}

#' Check if samples file already exists
#' @param file_path Path to check
#' @return Logical indicating if file exists
samples_exists <- function(file_path) {
  # Remove chain number and check if combined file exists
  base_path <- gsub("_chain[1234567]", "", file_path)
  return(file.exists(base_path))
}

#' Fast summary of MCMC chains
#' @param chains List of MCMC chains
#' @return Summary statistics
fast.summary.mcmc <- function(chains) {
  if (length(chains) == 0) {
    return(list(statistics = data.frame(), quantiles = data.frame()))
  }
  
  # Combine all chains
  all_samples <- do.call(rbind, chains)
  
  # Calculate statistics
  means <- colMeans(all_samples, na.rm = TRUE)
  sds <- apply(all_samples, 2, sd, na.rm = TRUE)
  
  # Calculate quantiles
  quantiles <- apply(all_samples, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
  
  # Create statistics dataframe
  stats <- data.frame(
    Mean = means,
    SD = sds,
    row.names = names(means)
  )
  
  # Create quantiles dataframe
  quants <- data.frame(
    t(quantiles),
    row.names = names(means)
  )
  colnames(quants) <- c("2.5%", "25%", "50%", "75%", "97.5%")
  
  return(list(
    statistics = stats,
    quantiles = quants
  ))
}

#' Parse model ID to extract components
#' @param model_id Model identifier string
#' @return List of parsed components
parse_model_id <- function(model_id) {
  # Split by underscores
  parts <- strsplit(model_id, "_")[[1]]
  
  # Extract components based on expected format
  # Format: model_name_taxon_time_period
  if (length(parts) >= 3) {
    model_name <- parts[1]
    if (parts[2] == "env" && parts[3] == "cycl") {
      model_name <- "env_cycl"
      taxon <- parts[4]
      time_period <- paste(parts[5:6], collapse = "_")
    } else if (parts[2] == "env" && parts[3] == "cov") {
      model_name <- "env_cov"
      taxon <- parts[4]
      time_period <- paste(parts[5:6], collapse = "_")
    } else {
      taxon <- parts[2]
      time_period <- paste(parts[3:4], collapse = "_")
    }
  } else {
    model_name <- "unknown"
    taxon <- "unknown"
    time_period <- "unknown"
  }
  
  # Determine rank and group
  if (grepl("phylum_fun", taxon)) {
    rank <- "phylum_fun"
    group <- "ITS"
    summary_type <- "taxonomic"
  } else if (grepl("phylum_bac", taxon)) {
    rank <- "phylum_bac"
    group <- "16S"
    summary_type <- "taxonomic"
  } else {
    rank <- "unknown"
    group <- "unknown"
    summary_type <- "unknown"
  }
  
  return(list(
    rank.name = rank,
    time_period = time_period,
    rank_only = rank,
    species = taxon,
    group = group,
    model_name = model_name,
    fcast_type = summary_type
  ))
}

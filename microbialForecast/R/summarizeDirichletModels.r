#' @title summarize_dirichlet_model
#' @description Summarize NIMBLE Dirichlet regression models for microbial taxa
#' Assumes input RDS files contain a list of:
#' MCMC samples, parameter summaries, latent state summaries, and model-fitting metadata
#' @param file_path Path to the samples file
#' @param save_summary Whether to save the summary file
#' @param overwrite Whether to overwrite existing summary files
#' @param drop_other Whether to drop "other" category from results
#' @export
summarize_dirichlet_model <- function(file_path, save_summary = NULL, overwrite = NULL, drop_other = TRUE) {
  require(stringr)
  
  if (summary_exists(file_path)) {
    if (is.null(overwrite)) {
      return("Summary file already exists")
    }
  }
  
  # Read in file
  read_in <- readRDS(file_path)
  
  # Read in samples
  samples <- read_in$samples
  param_summary <- read_in$param_summary
  plot_summary <- read_in$plot_summary
  truth.plot.long <- read_in$metadata$model_data
  
  # Extract run information
  info <- basename(file_path) %>% str_split("_") %>% unlist()
  model_id <- basename(file_path) %>% str_replace("samples_dirichlet_", "") %>% str_replace(".rds", "")
  
  parsed_id = parse_model_id(model_id)
  rank.name.eval <- parsed_id[[1]]
  model_name <- parsed_id[[6]]
  summary_type <- parsed_id[[3]]
  group <- parsed_id[[5]]
  time_period <- parsed_id[[2]]
  species <- parsed_id[[4]]
  rank_only <- parsed_id[[3]]
  
  # Add columns based on type
  if (summary_type == "functional") {
    rank.name <- rank.name.eval
    rank_only <- summary_type
    fg_cat <- assign_fg_categories(species)
    group <- assign_fg_kingdoms(fg_cat)
  } else {
    taxa_key = stack(microbialForecast:::rank_spec_names) %>%
      select(species = values, rank.name = ind)
    
    rank.name <- taxa_key[match(species, taxa_key$species), ]$rank.name
    rank_only <- rank.name %>% str_split("_") %>% unlist() %>% head(1)
  }
  taxon.name = species
  
  message("Summarizing Dirichlet model: ", species, ", ", rank.name, ", ", time_period, ", ", model_name)
  
  # Get covariate key
  cov_key <- switch(model_name,
                    "all_covariates" = microbialForecast:::all_covariates_key,
                    "env_cov" = microbialForecast:::all_covariates_key,
                    "env_cycl" = microbialForecast:::all_covariates_key,
                    "cycl_only" = microbialForecast:::cycl_only_key)
  
  # Get taxon and site keys - handle empty data case
  if (nrow(truth.plot.long) == 0) {
    # Create minimal keys for empty data
    taxon_key <- c(species)
    names(taxon_key) <- "1"
    site_key <- c("SITE1")
    names(site_key) <- "1"
  } else {
    taxon_key <- unique(truth.plot.long$species)
    names(taxon_key) <- seq(1, length(taxon_key))
    
    sites <- truth.plot.long %>% select(site_num, siteID) %>% unique()
    site_key <- sites[["siteID"]]
    names(site_key) <- sites[["site_num"]]
  }
  
  # Add metadata to observational data - handle empty data case
  if (nrow(truth.plot.long) == 0) {
    # Create minimal metadata structure
    truth.plot.long <- data.frame(
      dates = as.Date("2015-11-01"),
      truth = 0,
      model_name = model_name,
      taxon = species,
      rank = rank.name,
      group = group,
      rank_only = rank_only,
      time_period = time_period,
      fcast_type = summary_type,
      pretty_group = ifelse(group %in% c("16S", "bac"), "Bacteria", "Fungi"),
      model_id = model_id,
      stringsAsFactors = FALSE
    )
  } else {
    truth.plot.long <- truth.plot.long %>%
      mutate(dates = fixDate(dateID),
             truth = as.numeric(truth),
             model_name = !!model_name,
             taxon = species,
             rank = rank.name,
             group = !!group,
             rank_only = !!rank_only,
             time_period = !!time_period,
             fcast_type = !!summary_type,
             pretty_group = ifelse(group %in% c("16S", "bac"), "Bacteria", "Fungi"),
             model_id = !!model_id) %>%
      mutate(time_period = recode(as.character(!!time_period), !!!microbialForecast:::date_recode))
  }
  
  if (summary_type == "functional") {
    truth.plot.long <- truth.plot.long %>%
      mutate(fg_cat = !!fg_cat,
             fcast_type = "Functional")
  } else {
    truth.plot.long <- truth.plot.long %>%
      mutate(fg_cat = NA,
             fcast_type = "Taxonomic")
  }
  
  if (drop_other) {
    truth.plot.long <- truth.plot.long %>% filter(species != "other")
  }
  
  if (nchar(species) > 0) {
    taxon_key[[1]] = species
  }
  
  # Calculate plot median and quantiles - handle empty plot summary
  if (nrow(plot_summary[[1]]) > 0 && nrow(plot_summary[[2]]) > 0) {
    pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
      merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)
    
    # For scoring the predictions, need mean and SD
    pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
      merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)
    
    pred.quantiles$Mean <- pred.means$Mean
    pred.quantiles$SD <- pred.means$SD
  } else {
    # Create empty plot summaries
    pred.quantiles <- data.frame()
    pred.means <- data.frame()
  }
  
  # Get mean values for parameters
  means <- param_summary[[1]]
  quantiles <- param_summary[[2]]
  
  # Extract Dirichlet-specific parameters
  # Dirichlet models have beta[i, j] where i=covariate, j=taxon
  eff_list <- lapply(c("sigma", "sig$", "intercept", "rho", "core_sd"),
                     function(x) extract_summary_row(means, var = x)) %>%
    plyr::rbind.fill() %>%
    mutate(taxon = !!species)
  
  eff_list2 <- lapply(c("sigma", "sig$", "intercept", "rho", "core_sd"),
                      function(x) extract_summary_row(quantiles, var = x)) %>%
    plyr::rbind.fill() %>%
    mutate(taxon = !!species)
  
  if (nrow(eff_list2) > 0) {
    eff_list$Median = eff_list2[, "50%"]
  } else {
    eff_list$Median = NA
  }
  
  # Get site effect sizes - handle empty case
  site_eff_out <- tryCatch({
    site_data <- extract_summary_row(means, var = "site") %>%
      extract_bracketed_vals(varname1 = "site_num")
    
    if (nrow(site_data) > 0 && length(site_key) > 0 && !any(is.null(site_key))) {
      site_data <- site_data %>%
        mutate(taxon = !!species,
               siteID = recode(site_num, !!!site_key))
    } else {
      site_data <- data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
    }
    site_data
  }, error = function(e) {
    data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
  })
  
  site_eff_out2 <- tryCatch({
    site_data <- extract_summary_row(quantiles, var = "site") %>%
      extract_bracketed_vals(varname1 = "site_num")
    
    if (nrow(site_data) > 0 && length(site_key) > 0 && !any(is.null(site_key))) {
      site_data <- site_data %>%
        mutate(taxon = !!species,
               siteID = recode(site_num, !!!site_key))
    } else {
      site_data <- data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
    }
    site_data
  }, error = function(e) {
    data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
  })
  
  if (nrow(site_eff_out2) > 0 && "50%" %in% colnames(site_eff_out2)) {
    site_eff_out$Median = site_eff_out2[, "50%"]
  } else {
    site_eff_out$Median = NA
  }
  
  # Get beta parameters - Dirichlet format is beta[covariate, taxon]
  beta_out <- tryCatch({
    beta_data <- extract_summary_row(means, var = "beta") %>%
      extract_bracketed_vals(varname1 = "covariate_num", varname2 = "taxon_num")
    
    if (nrow(beta_data) > 0 && length(cov_key) > 0 && length(taxon_key) > 0) {
      beta_data <- beta_data %>%
        mutate(beta = recode(covariate_num, !!!cov_key),
               taxon = recode(taxon_num, !!!taxon_key))
    } else {
      beta_data <- data.frame(taxon = species, beta = "UNKNOWN", stringsAsFactors = FALSE)
    }
    beta_data
  }, error = function(e) {
    data.frame(taxon = species, beta = "UNKNOWN", stringsAsFactors = FALSE)
  })
  
  # Use quantiles to assign significance to beta parameters
  beta_ci <- tryCatch({
    beta_data <- extract_summary_row(param_summary[[2]], var = "beta") %>%
      extract_bracketed_vals(varname1 = "covariate_num", varname2 = "taxon_num")
    
    if (nrow(beta_data) > 0 && length(cov_key) > 0 && length(taxon_key) > 0) {
      beta_data <- beta_data %>%
        mutate(beta = recode(covariate_num, !!!cov_key),
               taxon = recode(taxon_num, !!!taxon_key))
    } else {
      beta_data <- data.frame(taxon = species, beta = "UNKNOWN", `2.5%` = NA, `97.5%` = NA, stringsAsFactors = FALSE)
    }
    beta_data
  }, error = function(e) {
    data.frame(taxon = species, beta = "UNKNOWN", `2.5%` = NA, `97.5%` = NA, stringsAsFactors = FALSE)
  })
  
  # Only calculate significance if we have valid data
  if (nrow(beta_out) > 0 && nrow(beta_ci) > 0 && all(!is.na(beta_ci$`2.5%`)) && all(!is.na(beta_ci$`97.5%`))) {
    beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
    beta_out$effSize <- abs(beta_out$Mean)
  } else {
    beta_out$significant <- NA
    beta_out$effSize <- NA
  }
  
  # Combine parameter estimates into summary
  summary_df <- plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>% 
    mutate(rank = rank.name) %>%
    left_join(truth.plot.long[, c("model_name", "rank", "group", "rank_only", "time_period",
                                 "fcast_type", "pretty_group", "model_id")] %>% distinct())
  
  # Calculate gelman diagnostics to assess convergence
  gd <- add_gelman(read_in, rank.name) %>% 
    mutate(rank = rank.name, taxon = !!species) %>%
    left_join(truth.plot.long[, colnames(truth.plot.long) %in% c("model_name", "rank", "group", "rank_only", "time_period",
                                                                "fcast_type", "pretty_group", "model_id")] %>% distinct())
  
  if (drop_other && "taxon" %in% colnames(summary_df)) {
    summary_df <- summary_df %>% filter(taxon != "other")
  }
  
  if (drop_other && "taxon" %in% colnames(pred.quantiles) && nrow(pred.quantiles) > 0) {
    pred.quantiles <- pred.quantiles %>% filter(taxon != "other")
  }
  
  out <- list(summary_df, pred.quantiles, gd)
  if (!is.null(save_summary)) {
    savePath <- gsub("samples", "summary", file_path)
    saveRDS(out, savePath)
    message("Saved Dirichlet summary to ", savePath)
    return(TRUE)
  } else {
    return(out)
  }
}

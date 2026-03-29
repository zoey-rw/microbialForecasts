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
  model_data <- read_in$metadata$model_data

  # Dirichlet model_data is a list with y, keep_taxa, truth.plot.long, plot_num, timepoint.
  # Build long-format truth: one row per (plot_num, timepoint, taxon) with truth = observed relative abundance.
  required <- c("y", "keep_taxa", "truth.plot.long", "plot_num", "timepoint")
  if (!is.list(model_data) || is.data.frame(model_data) || !all(required %in% names(model_data))) {
    stop("Dirichlet metadata$model_data must be a list with: ", paste(required, collapse = ", "))
  }
  base <- as.data.frame(model_data$truth.plot.long)
  base$plot_num <- model_data$plot_num
  base$timepoint <- model_data$timepoint
  base$siteID <- model_data$siteID
  base$site_num <- model_data$plot_site_num[base$plot_num]
  keep_taxa <- model_data$keep_taxa
  y <- as.matrix(model_data$y)
  n_obs <- nrow(y)
  n_taxa <- length(keep_taxa)
  truth.plot.long <- base[rep(seq_len(n_obs), each = n_taxa), , drop = FALSE]
  truth.plot.long$species <- rep(keep_taxa, times = n_obs)
  truth.plot.long$truth <- as.numeric(t(y))

  # Extract run information DIRECTLY from saved metadata to avoid string parsing errors
  meta <- read_in$metadata
  model_id <- meta$model_id
  rank.name <- meta$rank.name
  model_name <- meta$model_name
  time_period <- paste0(meta$min.date, "_", meta$max.date)
  species <- meta$species
  has_driver_uncertainty <- isTRUE(meta$has_driver_uncertainty)

  # Derive categorical helpers safely
  summary_type <- ifelse(is.null(meta$fcast_type), "Taxonomic", meta$fcast_type)
  rank_only <- str_split(rank.name, "_")[[1]][1]

  # Group assignment fallback
  if (!is.null(meta$group)) {
    group <- meta$group
  } else {
    group <- ifelse(grepl("fun", rank.name), "ITS", "16S")
  }

  if (summary_type == "functional") {
    fg_cat <- assign_fg_categories(species)
    group <- assign_fg_kingdoms(fg_cat)
  }

  taxon.name = species

  message("Summarizing Dirichlet model: ", species, ", ", rank.name, ", ", time_period, ", ", model_name)
  
  # Get covariate key
  cov_key <- switch(model_name,
                    "all_covariates" = microbialForecast:::env_cycl_covariates_key,
                    "env_cov" = microbialForecast:::env_cov_covariates_key,
                    "env_cycl" = microbialForecast:::env_cycl_covariates_key,
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
             truth = if ("truth" %in% names(.)) as.numeric(truth) else NA_real_,
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
  
  # taxon_key is already correctly built from keep_taxa above;
  # do NOT overwrite taxon 1 with the model-level species name.

  # Calculate plot median and quantiles - handle empty plot summary
  if (nrow(plot_summary[[1]]) > 0 && nrow(plot_summary[[2]]) > 0) {

    # Custom parser for 3D plot_rel[plot, species, time]
    parse_plot_rel <- function(df) {
      df$parameter <- rownames(df)
      extracted <- str_match(df$parameter, "plot_rel\\[(\\d+),\\s*(\\d+),\\s*(\\d+)\\]")
      df$plot_num <- as.numeric(extracted[, 2])
      df$taxon <- as.character(taxon_key[as.character(extracted[, 3])])
      df$timepoint <- as.numeric(extracted[, 4])
      return(df)
    }

    pred.quantiles <- parse_plot_rel(plot_summary[[2]]) %>%
      merge(truth.plot.long, by = c("plot_num", "timepoint", "taxon"), all = TRUE)

    pred.means <- parse_plot_rel(plot_summary[[1]]) %>%
      merge(truth.plot.long, by = c("plot_num", "timepoint", "taxon"), all = TRUE)

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
  
  # Get beta parameters - NIMBLE model outputs beta[taxon, covariate]
  beta_out <- tryCatch({
    beta_data <- extract_summary_row(means, var = "beta") %>%
      extract_bracketed_vals(varname1 = "taxon_num", varname2 = "covariate_num")

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
      extract_bracketed_vals(varname1 = "taxon_num", varname2 = "covariate_num")
    
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
  
  # Only calculate significance if we have valid data and matching row counts
  if (nrow(beta_out) > 0 && nrow(beta_ci) == nrow(beta_out) && all(!is.na(beta_ci$`2.5%`)) && all(!is.na(beta_ci$`97.5%`))) {
    sig <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
    if (length(sig) == nrow(beta_out)) {
      beta_out$significant <- sig
      beta_out$effSize <- abs(beta_out$Mean)
    } else {
      beta_out$significant <- NA
      beta_out$effSize <- NA
    }
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

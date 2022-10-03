
#'  @title summary_exists
#'  @description Check if summary file already exists
#'
#' @export
summary_exists <- function(file_path) {
	save_path <- gsub("samples","summary",file_path)
	if (file.exists(save_path)){
	sample_info <- file.info(file_path) %>% select(mtime) %>% unlist()
	summary_info <- file.info(save_path) %>% select(mtime) %>% unlist()
	if (summary_info > sample_info) {
		return(TRUE)
	} else return(FALSE)
	} else return(FALSE)
}


#'  @title add_gelman
#'  @description Calculate gelman diagnostics to assess convergence
#'
#' @export
add_gelman <- function(read_in, rank.name) {
if ("gelman" %in% names(read_in)) {
	gd <- read_in$gelman %>% cbind.data.frame(effSize = effectiveSize(read_in$samples))
} else {
	gd <- cbind.data.frame(gelman.diag(read_in$samples, multivariate = FALSE)[[1]],
												 effSize = effectiveSize(read_in$samples))
}
gd <- gd %>% mutate(rank.name = !!rank.name,
										#taxon.name = !!taxon.name,
										niter = read_in$metadata[[2]],
										nburnin = read_in$metadata[[3]])
return(gd)
}



#'  @title summarize_fg_div_model
#'  @description summarize_fg_div_model
#'
#' @export
summarize_fg_div_model <- function(file_path, drop_other = TRUE){
	require(stringr)

	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	#list2env(read_in,globalenv())

	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	truth.plot.long <- read_in$metadata$model_data

	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_name <- basename(dirname(file_path))
	rank.name <- info %>% head(-2) %>% tail(-1) %>% paste0(collapse = "_")
	time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")

	message("\nSummarizing ", rank.name, ", ", time_period, ", ", model_name)

	cov_key <- switch(model_name,
										"all_covariates" = all_covariates_key,
										"cycl_only" = cycl_only_key)


	# Add some info to observational data for merging
	truth.plot.long <- truth.plot.long %>%
		mutate(dates = fixDate(dateID),
					 truth = as.numeric(truth),
					 model_name = !!model_name,
					 taxon = rank.name,
					 time_period = !!time_period)

	if (grepl("functional_groups", file_path)) {
		truth.plot.long <- truth.plot.long %>% mutate(
			fcast_type = "Functional group",
			rank = "functional_group",
			fg_cat = assign_fg_categories(species),
			group = assign_fg_kingdoms(fg_cat),
			pretty_group = ifelse(group == "16S", "Bacteria", "Fungi"))
	} else if (grepl("div", file_path)) {
		truth.plot.long <- truth.plot.long %>% mutate(
			fcast_type = "Diversity",
			rank = "diversity",
			group = ifelse(grepl("16S", rank.name), "16S", "ITS"),
			pretty_group = ifelse(group == "16S", "Bacteria", "Fungi"))
	} else {
		message("File path does not contain 'div' or 'functional_groups'")
	}

	if (drop_other) {
		truth.plot.long <- truth.plot.long %>% filter(species != "other")
	}

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)

	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)

	# Get mean values for parameters
	means <- param_summary[[1]]
	eff_list <- lapply(c("sigma", "sig$", "core", "intercept"),
										 function(x) extract_summary_row(means, var = x)) %>%
		plyr::rbind.fill()

	# Get site effect sizes per rank
	site_eff_out <- extract_summary_row(means, var = "site")  %>%
		extract_bracketed_vals("site_num") %>%
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num), ]$siteID)

	# Get beta sizes per rank
	beta_out <- extract_summary_row(means, var = "beta|rho") %>%
		extract_bracketed_vals("beta_num") %>%
		mutate(beta = recode(beta_num, !!!cov_key))
	beta_out[grep("rho", beta_out$rowname), ]$beta = "rho"
	beta_out[grep("rho", beta_out$rowname), ]$beta_num = "0"


	# Use quantiles to assign significance to beta parameters.
	beta_ci <- extract_summary_row(param_summary[[2]], var = "beta|rho") %>%
		extract_bracketed_vals("beta_num")
	beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	summary_df <-
		plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>%
		cbind.data.frame(truth.plot.long[1, 11:18])

	## Calculate gelman diagnostics to assess convergence
	gd <- add_gelman(read_in, rank.name)


	# if (drop_other) {
	# 	summary_df <- summary_df %>% filter(taxon != "other")
	# 	pred.means <- pred.means %>% filter(taxon != "other")
	# 	pred.quantiles <- pred.quantiles %>% filter(taxon != "other")
	# }
	out <- list(summary_df, pred.means, pred.quantiles, gd)

	return(out)
}


#' @title summarize_tax_model
#' @description summarize_tax_model
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
summarize_tax_model <- function(file_path, save_summary = NULL, drop_other = TRUE){
	require(stringr)
 if(summary_exists(file_path)) return("Summary file already exists")

	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	truth.plot.long <- read_in$metadata$model_data

	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_name <- basename(dirname(file_path))
	if (tail(info,1) == "summary.rds") {
		if (info[[1]] == "refit") {
			time_period <- "2013-06_2020-01"
			info <- info[c(2,3,4,5,6)]

		} else if (info[[1]] == "calibration") {
			time_period <- "2013-06_2017-01"
			info <- info[c(2,3,4,5,6)]
		} else 	{
			info <- info %>% head(-1)
			time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
		}
	} else {
		time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
	}
	rank.name <- info %>% head(-2) %>% tail(-1) %>% paste0(collapse = "_")
	species <- info %>% head(-2) %>% tail(-3) %>% paste0(collapse = "_")
	rank_only <- info[[2]]
	group <- ifelse(info[[3]] == "bac", "16S", "ITS")

	message("\nSummarizing ",group,", ", rank.name, ", ", time_period, ", ", model_name)

	if (model_name == "convergence_testing") {
		model_name = "cycl_only"
		taxon.name = info[[4]]
	}


	cov_key <- switch(model_name,
										"all_covariates" = microbialForecast:::all_covariates_key,
										"cycl_only" = microbialForecast:::cycl_only_key)

	taxon_key <- unique(truth.plot.long$species)
	names(taxon_key) <- seq(1, length(taxon_key))

	sites <- truth.plot.long %>% select(site_num, siteID) %>% unique()
	site_key <- sites[["siteID"]]
	names(site_key) <- sites[["site_num"]]

	# Add some info to observational data for merging
	truth.plot.long <- truth.plot.long %>%
		mutate(dates = fixDate(dateID),
					 truth = as.numeric(truth),
					 model_name = !!model_name,
					 taxon = species,
					 rank = rank.name,
					 group = !!group,
					 rank_only = !!rank_only,
					 time_period = !!time_period,
					 fcast_type = "Taxonomic",
					 pretty_group = ifelse(group == "16S", "Bacteria", "Fungi"))

	if (drop_other) {
		truth.plot.long <- truth.plot.long %>% filter(species != "other")
	}

	if (nchar(species) > 0) {
		# truth.plot.long <- truth.plot.long %>%
		# 	mutate(taxon = !!species,
		# 				 species = !!species)
		taxon_key[[1]] = species
	}

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = recode(species_num, !!!taxon_key)) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)


	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = recode(species_num, !!!taxon_key)) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)

	pred.quantiles$Mean <- pred.means$Mean
	pred.quantiles$SD <- pred.means$SD

	# Get mean values for parameters
	means <- param_summary[[1]]
	eff_list <- lapply(c("sigma", "sig$", "intercept", "rho"),
										 function(x) extract_summary_row(means, var = x)) %>%
		plyr::rbind.fill() %>% extract_bracketed_vals("taxon_num") %>%
		mutate(taxon = recode(taxon_num, !!!taxon_key))

	# Get site effect sizes
	site_eff_out <- extract_summary_row(means, var = "site") %>%
		extract_bracketed_vals(varname1 = "site_num", varname2 = "taxon_num") %>%
		mutate(taxon = recode(taxon_num, !!!taxon_key),
					 siteID = recode(site_num, !!!site_key))

	# Get beta sizes per rank
	beta_out <- extract_summary_row(means, var = "beta") %>%
		extract_bracketed_vals(varname1 = "taxon_num", varname2 = "beta_num") %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = recode(taxon_num, !!!taxon_key))

	# Use quantiles to assign significance to beta parameters.
	beta_ci <- extract_summary_row(param_summary[[2]], var = "beta") %>%
		extract_bracketed_vals(varname1 = "taxon_num", varname2 = "beta_num") %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = recode(taxon_num, !!!taxon_key))
	beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	if (length(unique(beta_out$taxon)) > 2){
		summary_df <-
			plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>%
			left_join(truth.plot.long[, c("model_name", "taxon", "rank", "group", "rank_only", "time_period",
																		"fcast_type", "pretty_group")] %>% distinct())
	} else {
	summary_df <-
		plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>%
		left_join(truth.plot.long[1, 11:18], by="taxon")
	}


	## Calculate gelman diagnostics to assess convergence
	gd <- add_gelman(read_in, rank.name)

	if (drop_other) {
		summary_df <- summary_df %>% filter(taxon != "other")
		#pred.means <- pred.means %>% filter(taxon != "other")
		pred.quantiles <- pred.quantiles %>% filter(taxon != "other")
	}

	out <- list(summary_df, pred.quantiles, gd)
	if (!is.null(save_summary)) {
		savePath <- gsub("samples","summary",file_path)
		saveRDS(out, savePath)
		message("Saved summary to ", savePath)
		return(TRUE)
	} else {
	return(out)
	}
}





#'  @title summarize_fg_beta_model
#'  @description Summarize beta regression models for functional groups
#'
#' @export
summarize_fg_beta_model <- function(file_path, save_summary = NULL, drop_other = TRUE){
	require(stringr)
	if(summary_exists(file_path)) return("Summary file already exists")

	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	truth.plot.long <- read_in$metadata$model_data

	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_name <- basename(dirname(file_path))
	if (tail(info,1) == "summary.rds") {
		if (info[[1]] == "refit") {
			time_period <- "2013-06_2020-01"
			info <- info[c(2,3,4,5,6)]

		} else if (info[[1]] == "calibration") {
			time_period <- "2013-06_2017-01"
			info <- info[c(2,3,4,5,6)]
		} else 	{
			info <- info %>% head(-1)
			time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
		}
	} else {
		time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
	}
	rank.name <- info %>% head(-2) %>% tail(-2) %>% paste0(collapse = "_")
	species <- info %>% head(-2) %>% tail(-2) %>% paste0(collapse = "_")
	rank_only <- info[[3]]
	fg_cat <- assign_fg_categories(species)
	group <- assign_fg_kingdoms(fg_cat)

	message("\nSummarizing ",group,", ", rank.name, ", ", time_period, ", ", model_name)

	if (model_name == "convergence_testing") {
		model_name = "cycl_only"
		taxon.name = info[[4]]
	}


	cov_key <- switch(model_name,
										"all_covariates" = microbialForecast:::all_covariates_key,
										"cycl_only" = microbialForecast:::cycl_only_key)

	taxon_key <- unique(truth.plot.long$species)
	names(taxon_key) <- seq(1, length(taxon_key))

	sites <- truth.plot.long %>% select(site_num, siteID) %>% unique()
	site_key <- sites[["siteID"]]
	names(site_key) <- sites[["site_num"]]

	# Add some info to observational data for merging
	truth.plot.long <- truth.plot.long %>%
		mutate(dates = fixDate(dateID),
					 truth = as.numeric(truth),
					 model_name = !!model_name,
					 taxon = species,
					 rank = rank.name,
					 group = !!group,
					 rank_only = !!rank_only,
					 time_period = !!time_period,
					 fg_cat = !!fg_cat,
					 fcast_type = "Functional",
					 pretty_group = ifelse(group == "16S", "Bacteria", "Fungi"))

	if (drop_other) {
		truth.plot.long <- truth.plot.long %>% filter(species != "other")
	}

	if (nchar(species) > 0) {
		# truth.plot.long <- truth.plot.long %>%
		# 	mutate(taxon = !!species,
		# 				 species = !!species)
		taxon_key[[1]] = species
	}

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = recode(species_num, !!!taxon_key)) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)


	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = recode(species_num, !!!taxon_key)) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)

	pred.quantiles$Mean <- pred.means$Mean
	pred.quantiles$SD <- pred.means$SD

	# Get mean values for parameters
	means <- param_summary[[1]]
	eff_list <- lapply(c("sigma", "sig$", "intercept", "rho"),
										 function(x) extract_summary_row(means, var = x)) %>%
		plyr::rbind.fill() %>% extract_bracketed_vals("taxon_num") %>%
		mutate(taxon = recode(taxon_num, !!!taxon_key))

	# Get site effect sizes
	site_eff_out <- extract_summary_row(means, var = "site") %>%
		extract_bracketed_vals(varname1 = "site_num", varname2 = "taxon_num") %>%
		mutate(taxon = recode(taxon_num, !!!taxon_key),
					 siteID = recode(site_num, !!!site_key))

	# Get beta sizes per rank
	beta_out <- extract_summary_row(means, var = "beta") %>%
		extract_bracketed_vals(varname1 = "taxon_num", varname2 = "beta_num") %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = recode(taxon_num, !!!taxon_key))

	# Use quantiles to assign significance to beta parameters.
	beta_ci <- extract_summary_row(param_summary[[2]], var = "beta") %>%
		extract_bracketed_vals(varname1 = "taxon_num", varname2 = "beta_num") %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = recode(taxon_num, !!!taxon_key))
	beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	if (length(unique(beta_out$taxon)) > 2){
		summary_df <-
			plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>%
			left_join(truth.plot.long[, c("model_name", "taxon", "rank", "group", "rank_only", "time_period",
																		"fcast_type", "pretty_group")] %>% distinct())
	} else {
		summary_df <-
			plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>%
			left_join(truth.plot.long[1, 11:18], by="taxon")
	}


	## Calculate gelman diagnostics to assess convergence
	gd <- add_gelman(read_in, rank.name)

	if (drop_other) {
		summary_df <- summary_df %>% filter(taxon != "other")
		#pred.means <- pred.means %>% filter(taxon != "other")
		pred.quantiles <- pred.quantiles %>% filter(taxon != "other")
	}

	out <- list(summary_df, pred.quantiles, gd)
	if (!is.null(save_summary)) {
		savePath <- gsub("samples","summary",file_path)
		saveRDS(out, savePath)
		message("Saved summary to ", savePath)
		return(TRUE)
	} else {
		return(out)
	}
}




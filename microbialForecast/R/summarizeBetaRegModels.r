

#'  @title summarize_beta_regression
#'  @description Summarize NIMBLE state-space beta regression models for microbial taxa and functional groups
#'  Assumes input RDS files contain a list of:
#'  MCMC samples, parameter summaries, latent state summaries, and model-fitting metadata
#'	@param overwrite want to save new summary files even if there's an existing, recent summary file
#' @export
summarize_beta_model <- function(file_path, save_summary = NULL, overwrite=NULL, drop_other = TRUE){
	require(stringr)
	if(summary_exists(file_path)) { # checks that a summary is needed (samples files are new)
		if (is.null(overwrite)) {
			return("Summary file already exists")
		}
	}
	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	
	# Handle different data structures
	# Newer models have truth.plot.long nested under model_data
	# Older models have the data directly in model_data
	if("truth.plot.long" %in% names(read_in$metadata$model_data)) {
		truth.plot.long <- read_in$metadata$model_data$truth.plot.long
	} else {
		truth.plot.long <- read_in$metadata$model_data
	}


	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_id <- basename(file_path) %>%  str_replace("samples_", "") %>%  str_replace(".rds", "")

	parsed_id = parse_model_id(model_id)
	rank.name.eval <- parsed_id[[1]]
	model_name <- parsed_id[[6]]
	summary_type <- parsed_id[[3]]
	group  <- parsed_id[[5]]
	time_period <- parsed_id[[2]]

	if (tail(info,1) == "summary.rds") {
			info <- info %>% head(-1)
	}

	species <- parsed_id[[4]]
	rank_only <- parsed_id[[3]]

# Add columns based on
	if (summary_type=="functional") {
		rank.name <- rank.name.eval
		rank_only <- summary_type
	fg_cat <- assign_fg_categories(species)
	group <- assign_fg_kingdoms(fg_cat)
	} else {

		taxa_key = stack(microbialForecast:::rank_spec_names) %>%
			select(species = values, rank.name = ind)

		rank.name <- taxa_key[match(species, taxa_key$species),]$rank.name
		rank_only <-  rank.name  %>% str_split("_") %>% unlist() %>% head(1)
		}
	taxon.name = species

	message("Summarizing ", species, ", ", rank.name, ", ", time_period, ", ", model_name)


	cov_key <- switch(model_name,
										"all_covariates" = microbialForecast:::all_covariates_key,
										"env_cov" = microbialForecast:::all_covariates_key,
										"env_cycl" = microbialForecast:::all_covariates_key,
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
					 fcast_type = !!summary_type,
					 pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"),
					 model_id = !!model_id) %>%
		mutate(time_period =
					 	recode(as.character(!!time_period), !!!microbialForecast:::date_recode))

	if (summary_type=="functional") {
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

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)


	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)

	pred.quantiles$Mean <- pred.means$Mean
	pred.quantiles$SD <- pred.means$SD

	# Get mean values for parameters
	means <- param_summary[[1]]
	quantiles <- param_summary[[2]]
	eff_list <- lapply(c("sigma", "sig$", "intercept", "rho", "core_sd"),
										 function(x) extract_summary_row(means, var = x)) %>%
		plyr::rbind.fill() %>%
		mutate(taxon = !!species)
	eff_list2 <- lapply(c("sigma", "sig$", "intercept", "rho", "core_sd"),
											function(x) extract_summary_row(quantiles, var = x)) %>%
		plyr::rbind.fill() %>%
		mutate(taxon = !!species)
	eff_list$Median = eff_list2[,"50%"]

	# Get site effect sizes
	site_eff_out <- extract_summary_row(means, var = "site") %>%
		extract_bracketed_vals(varname1 = "site_num")  %>%
		mutate(taxon = !!species,
					 siteID = recode(site_num, !!!site_key))
	site_eff_out2 <- extract_summary_row(quantiles, var = "site") %>%
		extract_bracketed_vals(varname1 = "site_num")  %>%
		mutate(taxon = !!species,
					 siteID = recode(site_num, !!!site_key))
	site_eff_out$Median = site_eff_out2[,"50%"]

	# Get beta sizes per rank
	beta_out <- extract_summary_row(means, var = "beta") %>%
		extract_bracketed_vals(varname1 = "beta_num")  %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = !!species)
					# taxon = recode(taxon_num, !!!taxon_key))

	# Use quantiles to assign significance to beta parameters.
	beta_ci <- extract_summary_row(param_summary[[2]], var = "beta") %>%
		extract_bracketed_vals(varname1 = "beta_num")  %>%
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = !!species)
					# taxon = recode(taxon_num, !!!taxon_key))
	beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	#if (length(unique(beta_out$taxon)) > 2){
	summary_df <-
		plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>% mutate(rank = rank.name) %>%
		left_join(truth.plot.long[, c("model_name", #"taxon",
																	"rank", "group", "rank_only", "time_period",
																	"fcast_type","pretty_group","model_id")] %>% distinct())
	# } else {
	# 	summary_df <-
	# 		plyr::rbind.fill(beta_out, eff_list, site_eff_out)  %>% mutate(rank = rank.name) %>%
	# 		left_join(truth.plot.long[1, 11:18])
	# }


	## Calculate gelman diagnostics to assess convergence
	gd <- add_gelman(read_in, rank.name) %>% mutate(rank = rank.name,  taxon = !!species) %>%
		left_join(truth.plot.long[, colnames(truth.plot.long) %in% c("model_name", #"taxon",
																	"rank", "group", "rank_only", "time_period",
																	"fcast_type", "pretty_group","model_id")]  %>% distinct())

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

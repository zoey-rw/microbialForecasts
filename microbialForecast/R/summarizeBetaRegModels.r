

#'  @title summarize_beta_regression
#'  @description Summarize NIMBLE state-space beta regression models for microbial taxa and functional groups
#'  Assumes input RDS files contain a list of:
#'  MCMC samples, parameter summaries, latent state summaries, and model-fitting metadata
#'
#' @export
summarize_beta_model <- function(file_path, save_summary = NULL, overwrite=NULL, drop_other = TRUE){
	require(stringr)
	if(summary_exists(file_path)) {
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

	summary_type <- "taxon"
	# TODO: insert a real check here
	if (summary_type=="functional") {
	rank.name <- info %>% head(2) %>% tail(2) %>% paste0(collapse = "_")
	species <- info %>% tail(2) %>% head(2) %>% paste0(collapse = "_")
	rank_only <- info[[3]]
	fg_cat <- assign_fg_categories(species)
	group <- assign_fg_kingdoms(fg_cat)
	} else {
		rank.name <- info %>% head(3) %>% tail(2) %>% paste0(collapse = "_")
		rank_only <-  info %>% head(2) %>% tail(1)
		species <- info %>% tail(3) %>% head(1) %>% paste0(collapse = "_")
		group <- info %>% head(3) %>% tail(1)
		}

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
					 pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

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
		# truth.plot.long <- truth.plot.long %>%
		# 	mutate(taxon = !!species,
		# 				 species = !!species)
		taxon_key[[1]] = species
	}

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = !!species) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)


	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		mutate(taxon = !!species) %>%
		#mutate(taxon = recode(species_num, !!!taxon_key)) %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint","taxon"), all = T)

	pred.quantiles$Mean <- pred.means$Mean
	pred.quantiles$SD <- pred.means$SD

	# Get mean values for parameters
	means <- param_summary[[1]]
	eff_list <- lapply(c("sigma", "sig$", "intercept", "rho"),
										 function(x) extract_summary_row(means, var = x)) %>%
		plyr::rbind.fill() %>%
		mutate(taxon = !!species) #%>%
		#mutate(taxon = recode(taxon_num, !!!taxon_key))

	# Get site effect sizes
	site_eff_out <- extract_summary_row(means, var = "site") %>%
		extract_bracketed_vals(varname1 = "site_num")  %>%
		mutate(taxon = !!species, #) %>%

		#mutate(taxon = recode(taxon_num, !!!taxon_key),
					 siteID = recode(site_num, !!!site_key))

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
																	"fcast_type","pretty_group")] %>% distinct())
	# } else {
	# 	summary_df <-
	# 		plyr::rbind.fill(beta_out, eff_list, site_eff_out)  %>% mutate(rank = rank.name) %>%
	# 		left_join(truth.plot.long[1, 11:18])
	# }


	## Calculate gelman diagnostics to assess convergence
	gd <- add_gelman(read_in, rank.name) %>% mutate(rank = rank.name,  taxon = !!species) %>%
		left_join(truth.plot.long[, colnames(truth.plot.long) %in% c("model_name", #"taxon",
																	"rank", "group", "rank_only", "time_period",
																	"fcast_type", "pretty_group")])

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

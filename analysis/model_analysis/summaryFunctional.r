library(scales)
library(coda)
library(tidyverse)

#source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepFunctionalData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Create lists for outputs.
plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()


rank.name <- "cellulolytic"
scenario <- "full_uncertainty"
# scenario <- "no_uncertainty"
# scenario <- "temporal_uncertainty"
# scenario <- "spatial_uncertainty"
out_scenarios <- list() 
missing_scenarios <- matrix(ncol = 2, dimnames = list(NULL, c("scenario", "rank.name")))

# # Loop through all ranks/scenarios
# for (scenario in c("no_uncertainty","temporal_uncertainty","spatial_uncertainty","full_uncertainty")){
# 
 	for (rank.name in keep_fg_names) {
		cat(paste0("\nSummarizing ", scenario, "...")); print(rank.name)
	# Subset to one rank.
	#rank.df <- d[[rank.name]]  
	group.model.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario,".rds")
	
	# Checking to see if any models didn't save successfully
# 	if (!file.exists(group.model.path)) {
# 		missing_scenarios <- rbind(missing_scenarios, c(rank.name, scenario))
# 		cat("missing: ",  c(rank.name, scenario))
# 		next()
# 	} else {next()}}}
	
	read_in <- readRDS(group.model.path)
	
#		if (length(read_in) < 4) {
			# Calculate summary and save output.
		# 	param_summary <- fast.summary.mcmc(read_in$samples$samples)
		# 	plot_summary <- fast.summary.mcmc(read_in$samples$samples2)
		# 	cat(paste0("\nSummarized output for ", rank.name))
		# 	out <- list(samples = read_in$samples$samples,
		# 							param_summary = param_summary,
		# 							metadata = read_in$metadata,
		# 							plot_summary = plot_summary)
		# 	saveRDS(out, group.model.path, compress = F)
		# 	cat(paste0("\nSaved summarized output for ", rank.name,", ", scenario))
		# # #} else {
		#	next() }}
	# }
	# Read in samples for visualization
	# full_summary <- read_in[[2]]
	# means <- full_summary[[2]]
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	# TO DO: remove once output is fixed
	truth.plot.long <- read_in$metadata$model_data
		
		# cbind.data.frame(plotID = read_in$metadata$model_data.plotID,
		# 																	plot_num = read_in$metadata$model_data.plot_num,
		# 																	timepoint = read_in$metadata$model_data.timepoint,
		# 																	truth = read_in$metadata$model_data.truth,
		# 																	dateID = read_in$metadata$model_data.dateID,
		# 																	site_num = read_in$metadata$model_data.site_num,
		# 																	siteID = read_in$metadata$model_data.siteID)
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d"),
					 truth = as.numeric(truth))
	
	
	# Calculate plot means and actual values per rank
	plot_est <- plot_summary[[2]]
	pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)),
					 rank = "functional_group",
					 scenario = scenario,
					 taxon = rank.name)
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	plot_est_df_all[[rank.name]] <- allplots
	
	
	# For scoring the predictions (need mean and SD)
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)),
					 scenario = scenario,
					 taxon = rank.name)
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores.list[[rank.name]] <- allplots.scores
	
	
	# Get mean values for parameters
	means <- param_summary[[1]]
	beta_out <- means[grep("beta", rownames(means)),]
	rho_out <- means[grep("rho", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% mutate(beta = "rho")
	sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sig_out <- means[grep("sig$", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	core_out <- means[grep("core", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	int_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	
	print(int_out)
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID,
					 taxon = rank.name, rank = "functional_group")

	# Get beta sizes per rank
	beta_out <-  means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(beta_num,
												 "1" = "Temperature",
												 "2" = "Moisture",
												 "3" = "pH",
												 "4" = "pC",
												 "5" = "Plant species richness",
												 "6" = "Ectomycorrhizal trees",
												 "7" = "% invasive species"),
					 taxon = rank.name, rank = "functional_group")
	
	# Use quantiles to assign significance to beta parameters.
	beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) 
	beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
																 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
																 1, 0)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	summary_df <- plyr::rbind.fill(beta_out, sigma_out, rho_out, sig_out, core_out, int_out, site_eff_out) %>%
		mutate(scenario = scenario,
					 taxon = rank.name)
	summary_df_all[[rank.name]] <- summary_df
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[rank.name]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
	
}
	
	plot_est_df_all_ranks <- plyr::rbind.fill(plot_est_df_all) %>% 
		mutate(fcast_type = "Functional group", fcast_period = "calibration")

	scores.list_all_ranks <- plyr::rbind.fill(allplots.scores.list) %>% 
		mutate(fcast_type = "Functional group", fcast_period = "calibration")

	summary_df_all_ranks <- do.call(rbind, summary_df_all) %>% 
		mutate(fcast_type = "Functional group", fcast_period = "calibration", group = "16S")

summary_df_all_ranks$fg_cat <- assign_fg_categories(summary_df_all_ranks$taxon)
summary_df_all_ranks[grepl("Troph", summary_df_all_ranks$fg_cat),]$group = "ITS"

out_scenarios[[scenario]] <- list(plot_est = plot_est_df_all_ranks,
						 summary_df = summary_df_all_ranks,
						 gelman_list = gelman_list,
						 scores.list = scores.list_all_ranks)
}

saveRDS(out_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")

# saveRDS(out_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries_fixed_subset.rds")

# missing_scenarios
# write.csv(missing_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/missing_scenarios.csv")

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




file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups/", recursive = T,
												#pattern = "^refit_samples",
												full.names = T)
f <- file.list[[18]]
f <- file.list[[98]]

 	for (f in file.list) {
	# Subset to one rank.
	#rank.df <- d[[rank.name]]  
	group.model.path <- f
	
	# Checking to see if any models didn't save successfully
# 	if (!file.exists(group.model.path)) {
# 		missing_scenarios <- rbind(missing_scenarios, c(rank.name, scenario))
# 		cat("missing: ",  c(rank.name, scenario))
# 		next()
# 	} else {next()}}}
	
	read_in <- readRDS(group.model.path)
	
	info <- basename(f) %>% str_split("_")
	model_name <- basename(dirname(f))
	time_period <- info[[1]][1]
	last <- length(info[[1]])
	scenario <- paste0(info[[1]][(last-1):last], collapse = "_") %>% str_replace(".rds", "")
	rank.name <- paste0(info[[1]][3:(last-2)], collapse = "_") 

	
	cat(paste0("\nSummarizing ", model_name, ", ", time_period, ", ", rank.name))
	
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
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	allplots$scenario <- scenario
	allplots$time_period <- time_period
	allplots$model_name <- model_name
	allplots$taxon <- rank.name
	allplots$rank <- "functional_group"
	
	plot_est_df_all[[f]] <- allplots
	
	
	# For scoring the predictions (need mean and SD)
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores$scenario <- scenario
	allplots.scores$time_period <- time_period
	allplots.scores$model_name <- model_name
	allplots.scores$taxon <- rank.name
	allplots.scores.list[[f]] <- allplots.scores
	
	
	
	
	all_covariates_key <- c("1" = "Temperature",
													"2" = "Moisture",
													"3" = "pH",
													"4" = "pC",
													"5" = "Ectomycorrhizal trees",
													"6" = "LAI",
													"7" = "sin",
													"8" = "cos",
													"NA" = "NA")
	
	cycl_only_key <- list("1" = "sin",
												"2" = "cos")
	
	if(model_name == "all_covariates") cov_key <- all_covariates_key
	if(model_name == "cycl_only") cov_key <- cycl_only_key
	
	
	# Get mean values for parameters
	means <- param_summary[[1]]
	sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sig_out <- means[grep("sig$", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	core_out <- means[grep("core", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	int_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID)

	# Get beta sizes per rank
	beta_out <-  means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta|rho", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(beta_num, !!!cov_key),
					 taxon = rank.name, rank = "functional_group")
	beta_out[grep("rho", beta_out$rowname),]$beta = "rho"
	beta_out[grep("rho", beta_out$rowname),]$beta_num = "0"
	
	
	# Use quantiles to assign significance to beta parameters.
	beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta|rho", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) 
	beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
																 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
																 1, 0)
	beta_out$effSize <- abs(beta_out$Mean)

	# Combine parameter estimates into summary
	summary_df <- plyr::rbind.fill(beta_out, sigma_out, #rho_out, 
																 sig_out, core_out, int_out, site_eff_out) 
	summary_df$scenario <- scenario
	summary_df$time_period <- time_period
	summary_df$model_name <- model_name
	summary_df$taxon <- rank.name
	summary_df_all[[f]] <- summary_df
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[f]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
	
}
	
	plot_est_df_all_ranks <- plyr::rbind.fill(plot_est_df_all) %>% 
		mutate(fcast_type = "Functional group")

	scores.list_all_ranks <- plyr::rbind.fill(allplots.scores.list) %>% 
		mutate(fcast_type = "Functional group")

	summary_df_all_ranks <- do.call(rbind, summary_df_all) %>% 
		mutate(fcast_type = "Functional group", group = "16S",
					 pretty_name = "Functional group",
					 only_rank = "functional")

summary_df_all_ranks$fg_cat <- assign_fg_categories(summary_df_all_ranks$taxon)
summary_df_all_ranks[grepl("Troph", summary_df_all_ranks$fg_cat),]$group = "ITS"
summary_df_all_ranks$pretty_group <- ifelse(summary_df_all_ranks$group=="16S", "Bacteria", "Fungi")


rownames(summary_df_all_ranks) <- NULL

out_scenarios[[scenario]] <- list(plot_est = plot_est_df_all_ranks,
						 summary_df = summary_df_all_ranks,
						 gelman_list = gelman_list,
						 scores.list = scores.list_all_ranks)
#}

saveRDS(out_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")

# saveRDS(out_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries_fixed_subset.rds")

# missing_scenarios
# write.csv(missing_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/missing_scenarios.csv")


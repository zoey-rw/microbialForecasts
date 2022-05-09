
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

	file.list <- list.files("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/convergence_testing", pattern = "20151101_20180101", full.names = T)
	f <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/convergence_testing/samples_class_bac_acidimicrobiia_20151101_20180101.rds"
	
	
	testing <- summarize_tax_model(f)
	file_summaries <- lapply(file.list[1:2], summarize_tax_model)
	
	cl <- makeCluster(28, outfile="")
	registerDoParallel(cl)
	#Run for multiple chains, in parallel (via PSOCK)
	file_summaries = foreach(f=file.list,
													 .errorhandling = 'pass') %dopar% {
													 	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
													 	
													 	require(stringr)
													 	out <- summarize_tax_model(f)
													 	return(out)
													 }
	stopCluster(cl)
	
	
	summary_df <- map(file_summaries, 1) %>% plyr::rbind.fill() %>% 
		mutate(time_period = recode(as.character(time_period), !!!date_recode))
	plot_est <- map(file_summaries, 3) %>% plyr::rbind.fill() %>% 
		mutate(time_period = recode(as.character(time_period), !!!date_recode))
	gelman_list <- map(file_summaries, 4) 
	names(gelman_list) <- paste(basename(dirname(file.list)),basename(file.list), sep = "_")
	scores.list <- map(file_summaries, 2) %>% plyr::rbind.fill() %>% 
		mutate(time_period = recode(as.character(time_period), !!!date_recode))
	
	
	saveRDS(list(summary_df = summary_df,
							 plot_est = plot_est,
							 scores.list = scores.list,
							 gelman.summary = gelman_list), 
					"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxon_convergence_summaries_201511_201801.rds")
	
	
	
	
	
	
	
	
	
	# old
	
	plot_est_df_all_ranks <- plyr::rbind.fill(plot_est_df_all)
	scores.list_all_ranks <- plyr::rbind.fill(allplots.scores.list)
	summary_df_all_ranks <- do.call(rbind, summary_df_all)
	
	
	summary_df_all_ranks <- summary_df_all_ranks %>% mutate(fcast_type = "Taxonomic_convergence",
																	fcast_period = "calibration",
																	fg_cat = "Bacterial taxa",
																	group = "16S") 
	summary_df_all_ranks[grepl("fun$", summary_df_all_ranks$rank.name),]$fg_cat = "Fungal taxa"
	summary_df_all_ranks[grepl("fun$", summary_df_all_ranks$rank.name),]$group = "ITS"

	gelman.summary <- do.call(rbind, gelman_list)
	rownames(gelman.summary) <- NULL
	
	
saveRDS(list(summary_df = summary_df_all_ranks,
						 gelman.summary = gelman.summary), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/allrank_convergence_summaries.rds")



gelman_list2 <- list()
for (f in file.list){
	
	# Prep model inputs/outputs.
	model.path <- f 
	
	split <- str_split(basename(f), "_|\\.")
	rank_only <- split[[1]][2]
	kingdom <- split[[1]][3]
	rank.name <- paste0(rank_only, "_", kingdom)
	taxon.name <- split[[1]][4]
	
	cat(paste("Summarizing convergence output of ", rank.name, ", group: ", taxon.name))
	
	# Read in samples
	read_in <- readRDS(model.path)
	samples <- read_in$samples
	
	## Calculate gelman diagnostics to assess convergence
	gd.orig <- gelman.diag(samples)
	gd <- cbind.data.frame(gd.orig[[1]], effSize = effectiveSize(samples))
	gd$rank.name <- rank.name
	gd$taxon.name <- taxon.name
	gd$niter <- read_in$metadata[[2]]
	gd$nburnin <- read_in$metadata[[3]]
	
	gd2 <- cbind.data.frame(gd.orig[[2]])
	gd2$rank.name <- rank.name
	gd2$taxon.name <- taxon.name
	gd2$niter <- read_in$metadata[[2]]
	gd2$nburnin <- read_in$metadata[[3]]
	gelman_list2[[f]] <- gd2
}
	
gelman.summary2 <- do.call(rbind, gelman_list2)



gelman.summary_remove <- gelman.summary[gelman.summary$`Upper C.I.` > 20,]



df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/allrank_convergence_summaries.rds")

gelman.summary <- df$gelman.summary
summary.df <- df$summary_df

keep_taxa <- list()


by.rank <- gelman.summary %>% 
	group_by(rank.name, taxon.name) %>% 
	summarize(median_gbr = median(`Point est.`))

remove <- by.rank[by.rank$median_gbr > 2.7,]$taxon.name
keep <- by.rank[by.rank$median_gbr <= 2.2,]$taxon.name
keep <- by.rank[by.rank$median_gbr <= 2.2,]


keep_list <- split(keep, keep$rank.name)
#by.rank[by.rank$taxon.name %in% keep,]

saveRDS(keep_list, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")

gelman.summary_remove <- gelman.summary[gelman.summary$taxon.name %in% remove,]


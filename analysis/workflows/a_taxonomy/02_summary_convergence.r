# Summarize MCMC output from calibration, cyclical-covariate single-taxon models
# To determine which should be included in dirichlet models

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

file.list = intersect(list.files(here("data/model_outputs/single_taxon/cycl_only"),recursive = T,
																	pattern = "20151101_20180101", full.names = T),
											 list.files(here("data/model_outputs/single_taxon/cycl_only"), recursive = T,
											 					 pattern = "samples", full.names = T))
#	out <- summarize_tax_model(file.list[[340]], save_summary=T, drop_other = T)

	cl <- makeCluster(28, outfile="")
	registerDoParallel(cl)
	#Run for multiple chains, in parallel (via PSOCK)
	file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
		source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
													 	out <- summarize_tax_model(f, save_summary=T, drop_other = F)
													 	return(out)
													 }
	stopCluster(cl)



	file.list2 = intersect(list.files(here("data/model_outputs/single_taxon/cycl_only"),
																		pattern = "20151101_20180101", full.names = T),
												 list.files(here("data/model_outputs/single_taxon/cycl_only"),
												 					 pattern = "summary", full.names = T))

	file_summaries <- purrr::map(file.list2, readRDS)

	#in1 <- readRDS(file.list2[[1]])


	summary_df <- map_df(file_summaries, 1) %>%
		mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
	plot_est <- map_df(file_summaries, 2) %>%
		mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
	gelman_list <- map_df(file_summaries, 3)


	saveRDS(list(summary_df = summary_df,
							 plot_est = plot_est,
							 gelman.summary = gelman_list),
					"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxon_convergence_summaries_201511_201801.rds")



	convergence_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxon_convergence_summaries_201511_201801.rds")

	scores.list <- convergence_in$scores.list


gelman.summary <- convergence_in$gelman.summary

gelman.summary <- gelman.summary %>%
	separate(rank.name, sep = "_", into = c("only_rank", "group", NA), remove = F) %>%
	separate(rank.name, sep = "_bac_|_fun_", into = c(NA, "taxon.name"), remove = F) %>%
	mutate(rank_group = paste(only_rank, group, sep = "_"))

by.group <- gelman.summary %>%
	group_by(group, rank_group, only_rank, rank.name, taxon.name) %>%
	dplyr::summarize(median_gbr = median(`Point est.`))

ggplot(by.group) + geom_jitter(aes(x = only_rank, y = median_gbr, color = group)) + ylim(c(0,10))



keep_taxa <- list()

by.rank <- gelman.summary %>%
	group_by(rank_group, taxon.name) %>%
	summarize(median_gbr = median(`Point est.`))

keep <- by.group[by.group$median_gbr <= 2,]
keep_list <- split(keep, keep$rank_group)

rerun <- by.group[by.group$median_gbr > 1.5,]
rerun_list <- split(rerun, rerun$rank_group)
gelman.summary_remove <- gelman.summary[gelman.summary$taxon.name %in% rerun$taxon.name,]

#by.rank[by.rank$taxon.name %in% keep,]

# This was actually run with cycl_only models
saveRDS(keep_list, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")


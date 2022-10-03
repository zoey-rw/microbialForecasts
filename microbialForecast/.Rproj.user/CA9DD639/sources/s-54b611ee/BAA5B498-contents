# Summarize MCMC output from all single-taxon models

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

file.list = intersect(list.files(here("data/model_outputs/single_taxon/"),recursive = T,
																 pattern = "20151101_20180101", full.names = T),
											list.files(here("data/model_outputs/single_taxon/"), recursive = T,
																 pattern = "samples", full.names = T))


file.list = list.files(here("data/model_outputs/single_taxon/"), recursive = T,
																 pattern = "samples", full.names = T)
#	out <- summarize_tax_model(file.list[[340]], save_summary=T, drop_other = F)

cl <- makeCluster(28, outfile="")
registerDoParallel(cl)
#Run summary function for multiple chains, in parallel (via PSOCK)
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	out <- summarize_tax_model(f, save_summary=T, drop_other = F)
	return(out)
}
stopCluster(cl)


# Combine summary files for all taxa calibration models
file.list2 = intersect(list.files(here("data/model_outputs/single_taxon/"),recursive = T,
																	pattern = "20151101_20180101", full.names = T),
											 list.files(here("data/model_outputs/single_taxon/"), recursive = T,
											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(file.list2, readRDS)
summary_df <- map_df(file_summaries, 1) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% distinct() %>%
	filter(taxon != "other")
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% distinct()
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

plot_est$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", plot_est$rank)
summary_df$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", summary_df$rank)
gelman_list$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", gelman_list$rank)

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/single_taxon_summaries_201511_201801.rds"))






# Combine summary files for all taxa refit models
file.list3 = intersect(list.files(here("data/model_outputs/single_taxon/"),recursive = T,
																	pattern = "20151101_20200101", full.names = T),
											 list.files(here("data/model_outputs/single_taxon/"), recursive = T,
											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(file.list3, readRDS)
summary_df <- map_df(file_summaries, 1) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct() %>%
	filter(taxon != "other")
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct() %>%
	filter(taxon != "other")
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

plot_est$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", plot_est$rank)
summary_df$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", summary_df$rank)
gelman_list$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", gelman_list$rank)

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/single_taxon_summaries_201511_202001.rds"))


ggplot(gelman_list) + geom_jitter(aes(x = `rank.name`,y = `Point est.`), alpha=.3) + facet_wrap(~rank_name, scales="free")







scores.list <- convergence_in$scores.list
gelman.summary <- gelman_list
gelman.summary <- gelman.summary %>%
	separate(rank.name, sep = "_", into = c("only_rank", "group", NA), remove = F) %>%
	separate(rank.name, sep = "_bac_|_fun_", into = c(NA, "taxon.name"), remove = F) %>%
	mutate(rank_group = paste(only_rank, group, sep = "_"))

by.group <- gelman.summary %>%
	group_by(group, rank_group, only_rank, rank.name, taxon.name, niter) %>%
	dplyr::summarize(median_gbr = median(`Point est.`))

ggplot(by.group) + geom_jitter(aes(x = only_rank, y = median_gbr, color = group)) + ylim(c(0,10))



keep_taxa <- list()

by.rank <- gelman.summary %>%
	group_by(rank_group, taxon.name) %>%
	summarize(median_gbr = median(`Point est.`))

keep <- by.group[by.group$median_gbr <= 2,]
keep_list <- split(keep, keep$rank_group)

rerun <- by.group[by.group$median_gbr > 1.5 & by.group$niter==1000000,]
rerun_list <- split(rerun, rerun$rank_group)
gelman.summary_remove <- gelman.summary[gelman.summary$taxon.name %in% rerun$taxon.name,]

# Whereas this is all_cov models
saveRDS(rerun_list, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")

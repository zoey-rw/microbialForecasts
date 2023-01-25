# Summarize MCMC output from all single-taxon models

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

file.list = intersect(list.files(here("data/model_outputs/single_taxon/"),recursive = T,
																 pattern = "20151101_20180101|20151101_20200101", full.names = T),
											list.files(here("data/model_outputs/single_taxon/"), recursive = T,
																 pattern = "samples", full.names = T))
file.list = file.list[!grepl("chain", file.list)]
file.list = file.list[grepl("2018", file.list)]


cl <- makeCluster(27, outfile="")
registerDoParallel(cl)
#Run summary function for multiple chains, in parallel (via PSOCK)
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	out <- summarize_tax_model(f, save_summary=T, drop_other = F, overwrite = TRUE)
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
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T) %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T)  %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

plot_est$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", plot_est$rank)
summary_df$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", summary_df$rank)
gelman_list$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", gelman_list$rank)

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/single_taxon_summaries_201511_201801.rds"))


current_summary = readRDS(here("data/summary/single_taxon_summaries_201511_201801.rds"))
mean(current_summary$gelman.summary[current_summary$gelman.summary$rank_name=="family_fun",]$`Point est.`, na.rm = T)
mean(gelman_list[gelman_list$rank_name=="family_fun",]$`Point est.`, na.rm = T)



# Combine summary files for all taxa refit models
file.list3 = intersect(list.files(here("data/model_outputs/single_taxon/"),recursive = T,
																	pattern = "20151101_20200101", full.names = T),
											 list.files(here("data/model_outputs/single_taxon/"), recursive = T,
											 					 pattern = "summary", full.names = T))
file.list3 = file.list3[!grepl("samples", file.list3)]

file_summaries <- purrr::map(file.list3, readRDS)

summary_df <- map_df(file_summaries, 1) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T) %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T) %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

plot_est$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", plot_est$rank)
summary_df$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", summary_df$rank)
gelman_list$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", gelman_list$rank)

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/single_taxon_summaries_201511_202001.rds"))


ggplot(gelman_list) + geom_jitter(aes(x = `rank.name`,y = `Point est.`), alpha=.3) + facet_wrap(~rank_name, scales="free")





gelman.summary <- gelman_list
# or
convergence_in = readRDS(here("data/summary/single_taxon_summaries_201511_201801.rds"))
gelman.summary <- convergence_in$gelman.summary

gelman.summary <- gelman.summary %>%
	separate(rank, sep = "_bac_|_fun_", into = c(NA, "taxon.name"), remove = F) %>%
	mutate(taxon_model_rank = paste(taxon.name, model_name, rank_name))


by.group <- gelman.summary %>%
	group_by(group, rank_name, rank_only, model_name, taxon.name) %>%
	dplyr::summarize(median_gbr = median(`Point est.`))

ggplot(by.group %>% filter(model_name=="all_covariates")) + geom_jitter(aes(x = rank_only, y = median_gbr, color = group), alpha=.2, height = 0) + ylim(c(0,10)) + geom_hline(yintercept = 1)


keep_taxa <- list()

by.rank <- gelman.summary %>%
	group_by(rank_name, taxon.name) %>%
	summarize(median_gbr = median(`Point est.`))

beta_only <- gelman.summary %>% filter(grepl("beta|rho|intercept", parameter)) %>%
	group_by(group,rank_name, rank_only,taxon.name, model_name, taxon_model_rank) %>%
	summarize(median_gbr = median(`Point est.`))

keep <- beta_only[beta_only$median_gbr <= 3,]
keep_list <- split(keep, keep$rank_name)

rerun <- beta_only[beta_only$median_gbr > 2,]
rerun_list <- split(rerun, rerun$rank_name)
gelman.summary_remove <- gelman.summary[gelman.summary$taxon.name %in% rerun$taxon.name,]

ggplot(beta_only %>% filter(model_name=="all_covariates")) + geom_jitter(aes(x = rank_only, y = median_gbr, color = group), alpha=.5, height = 0, size=3) + ylim(c(0,10)) + geom_hline(yintercept = 1)

# Whereas this is all_cov models
saveRDS(rerun, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
saveRDS(keep, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")

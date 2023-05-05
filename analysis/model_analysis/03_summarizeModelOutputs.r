# Summarize MCMC output from all single-taxon models
# Assumes input files have already had MCMC chains combined

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")


file.list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
																 pattern = "20130601_20151101|20151101_20180101|20151101_20200101", full.names = T),
											list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
																 pattern = "samples", full.names = T))

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files
info <- file.info(file.list)
newer <- rownames(info[which(info$mtime > "2023-04-30 00:00:00 EDT"), ])
#newer <- rownames(info[which(info$mtime > "2023-03-09 00:00:00 EDT"), ])
file.list <- file.list[file.list %in% newer]

f=file.list[[1]]
cl <- makeCluster(28, outfile="")
registerDoParallel(cl)

#Run summary function for multiple groups, in parallel
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/summarizeBetaRegModels.r")
	out <- summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
	return(out)
}

stopCluster(cl)


summary_file_list = list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
															 pattern = "summary", full.names = T)

# Subset to newest output files
info <- file.info(summary_file_list)
newer <- rownames(info[which(info$mtime > "2023-03-09 00:00:00 EDT"), ])
summary_file_list <- summary_file_list[summary_file_list %in% newer]

# Combine summary files for all models/time-periods
# summary_file_list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
# 																	pattern = "20151101_20180101", full.names = T),
# 											 list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
# 											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(summary_file_list, readRDS)
summary_df <- map_df(file_summaries, 1)
plot_est <- map_df(file_summaries, 2)
gelman_list <- map_df(file_summaries, 3)

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/logit_beta_regression_summaries.rds"))




#summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))


gelman.summary <- gelman_list %>%
	filter(model_name != "all_covariates") %>%
	mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd", parameter), T, F))
#mutate(model_id2 = paste(model_name, rank, taxon, time_period, sep = "_"))

by_rank <- gelman.summary %>%
	group_by(model_id, is_major_param) %>%
	dplyr::mutate(median_gbr = median(`Point est.`,na.rm=T),
								quant_95 = quantile(`Point est.`, c(.95),na.rm=T),
								min_es = min(es, na.rm=T),
								median_es = min(es, na.rm=T),
								mean_gbr = mean(`Point est.`,na.rm=T)) %>%
	distinct(.keep_all = T)

model_median = by_rank %>% select(c("rank.name", "is_major_param","niteration", "rank", "taxon", "model_name", "group", "rank_only",
																		"time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
	distinct(.keep_all = T)
model_median = model_median %>% pivot_wider(values_from = c("mean_gbr","median_gbr","quant_95","min_es","median_es"),names_from = is_major_param)


ggplot(model_median) + geom_jitter(aes(x = niteration, y = median_gbr_TRUE, color = group)) + ylim(c(0,5)) + geom_hline(yintercept = 1)



keep_models <- model_median %>%
	group_by(model_id) %>%
	filter(median_gbr_TRUE <= 1.1) %>%
	filter(mean_gbr_TRUE <= 1.2) %>%
	filter(mean_gbr_FALSE <= 1.5) %>%
	filter(min_es_TRUE > 75)
	#	filter(median_gbr_FALSE < 1.1)
keep_list <- unique(keep_models$model_id)


keep_models_weak <- model_median %>%
group_by(model_id) %>%
	filter(median_gbr_TRUE <= 1.15) %>%
	filter(mean_gbr_TRUE <= 1.5) %>%
	filter(mean_gbr_FALSE <= 2) %>%
	filter(min_es_TRUE > 15)
keep_list_weak <- unique(keep_models_weak$model_id)

rerun <- model_median %>% filter(!model_id %in% keep_list_weak)
rerun_list <- unique(rerun$model_id)

#rerun %>% ungroup %>% select(c(12:21)) %>% as.matrix %>% pairs

saveRDS(keep_list, here("data/summary/converged_taxa_list.rds"))
saveRDS(keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))
saveRDS(rerun_list, here("data/summary/unconverged_taxa_list.rds"))









#
#
#
# gelman_list[gelman_list$`Point est.` > 1.1,]$taxon
#
# fg_cycl_refit = readRDS(here("data/summary/beta_fg_summaries_20151101_20200101.rds"))
# fg_cycl_refit_plots = fg_cycl_refit$plot_est
# fg_cycl_refit_sum = fg_cycl_refit$summary_df
#
# summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
# hindcast_dat = summaries$plot_est
# #hindcast_dat = plot_est
# select_plots <- c("HARV_033","HARV_004","KONZ_001","KONZ_002")
# select_plots <- c("HARV_034","HARV_035","BLAN_033","BLAN_032","BART_071","BART_024")
#
#
# hindcast_dat <- hindcast_dat[hindcast_dat$date_num != 1,]
#
# ggplot(hindcast_dat %>% filter(plotID %in% select_plots &
# 															 	taxon == "assim_nitrate_reduction"),
# 			 aes(fill=species, x = dates, group=plotID)) +
# 	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
# 	facet_grid(siteID~ model_id, scales="free") #+ ylim(c(0,0.4))
#
# ggplot(hindcast_dat %>% filter(plotID %in% select_plots &
# 															 	taxon == "ascomycota"),
# 			 aes(fill=species, x = dates, y =Mean, group=plotID)) +
# 	geom_line(aes(x = dates, y = `Mean`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +facet_grid(model_id~siteID, scales="free")
#
#
#
#
# ggplot(hindcast_in2 %>% filter(plotID %in% select_plots &
# 															 	taxon == "ascomycota"),
# 			 aes(fill=species, x = dates, y =Mean, group=plotID)) +
# 	geom_line(aes(x = dates, y = `Mean`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +facet_grid(model_name~siteID)
#
# new_cellulo = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/model_outputs/logit_beta_regression/cycl_only/samples_cellulolytic_cellulolytic_20151101_20180101_chain4.rds")
#
# old_cellulo = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/model_outputs/logit_beta_regression/cycl_only/samples_cellulolytic_cellulolytic_20151101_20180101.rds")
#
#
#
# # For scoring the predictions, need mean and SD
# pred.means <- new_cellulo$samples2 %>% fast.summary.mcmc
#
#
# old_means = old_cellulo$plot_summary[[1]] %>%  parse_plot_mu_vars()
# old_quantiles <- old_cellulo$plot_summary[[2]] %>% parse_plot_mu_vars()
# old_quantiles$Mean <- old_means$Mean
# old_quantiles$SD <- old_means$SD
#
# # Calculate plot median and quantiles
# new_means = pred.means[[1]] %>%  parse_plot_mu_vars()
# new_quantiles <- pred.means[[2]] %>% parse_plot_mu_vars()
# new_quantiles$Mean <- new_means$Mean
# new_quantiles$SD <- new_means$SD
#
#
# ggplot(new_quantiles %>% filter(plot_num %in% c(1:5)),
# 			 aes(x = timepoint, group=plot_num)) +
# 	geom_line(aes(x = timepoint, y = `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = timepoint, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = timepoint, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	facet_grid(~plot_num, scales="free") #+ ylim(c(0,0.4))
#
#
#
#
#
#
#
# ggplot(plot_est %>% filter(plotID %in% c("HARV_004"),
# 													 #,"BART_004") &
# 													 taxon %in% c("ectomycorrhizal","chitin_complex") &
# 													 	time_period %in% c("20130601_20151101","2015-11_2018-01")),
# 			 aes(fill=species, x = dates, y = `50%`, group=plotID)) +
# 	geom_line(aes(x = dates, y =  `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
# 	facet_grid(model_name+time_period~taxon)
#
#
#
#
#
#
#
# ggplot(plot_est %>% filter(plotID %in% c("HARV_004"),
# 													 #,"BART_004") &
# 													 taxon %in% c("ectomycorrhizal") &
# 													 	model_name == "all_covariates" &
# 													 	time_period %in% c("20130601_20151101","2015-11_2018-01")),
# 			 aes(fill=species, x = dates, y = `50%`, group=plotID)) +
# 	geom_line(aes(x = dates, y =  `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
# 	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +
#
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
# 	facet_grid(time_period~taxon)

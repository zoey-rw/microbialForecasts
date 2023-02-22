# Summarize MCMC output from all single-taxon models
# Assumes input files have already had MCMC chains combined

source("./source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast/R/summarizeBetaRegModels.r")


file.list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
																 pattern = "20151101_20180101|20151101_20200101", full.names = T),
											list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
																 pattern = "samples", full.names = T))

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# For now, just run on calibration models
file.list = file.list[grepl("2018", file.list)]

f = file.list[[1]]
cl <- makeCluster(4, outfile="")
registerDoParallel(cl)
#Run summary function for multiple groups, in parallel
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast/R/summarizeBetaRegModels.r")
	out <- summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
	return(out)
}

stopCluster(cl)



# Combine summary files for all taxa calibration models
summary_file_list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
																	pattern = "20151101_20180101", full.names = T),
											 list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(summary_file_list, readRDS)
summary_df <- map_df(file_summaries, 1) %>%
	mutate(time_period =
				 	recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T)
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period =
				 	recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>%
	distinct(.keep_all=T)
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				here("data/summary/logit_beta_regression_summaries.rds"))


fg_cycl_refit = readRDS(here("data/summary/beta_fg_summaries_20151101_20200101.rds"))
fg_cycl_refit_plots = fg_cycl_refit$plot_est
fg_cycl_refit_sum = fg_cycl_refit$summary_df

summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
hindcast_in = summaries$plot_est %>% mutate(model_type = "logit_beta")


summaries2 = readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_201801.rds")
hindcast_in2 = summaries2$plot_est %>% mutate(model_type = "dirichlet")

hindcast_dat = rbindlist(list(hindcast_in, hindcast_in2), fill=T)
select_plots <- c("HARV_033","HARV_004","KONZ_001","KONZ_002")


ggplot(hindcast_dat %>% filter(plotID %in% select_plots &
															 	taxon == "acidobacteriota"),
			 aes(fill=species, x = dates, y =Mean, group=plotID)) +
	geom_line(aes(x = dates, y = `Mean`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +

	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +facet_grid(siteID+ model_name ~ model_type, scales="free") + ylim(c(0,0.4))

ggplot(hindcast_dat %>% filter(plotID %in% select_plots &
																taxon == "ascomycota"),
			 aes(fill=species, x = dates, y =Mean, group=plotID)) +
	geom_line(aes(x = dates, y = `Mean`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +

		theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +facet_grid(model_type~model_name:siteID, scales="free")




ggplot(hindcast_in2 %>% filter(plotID %in% select_plots &
															 	taxon == "ascomycota"),
			 aes(fill=species, x = dates, y =Mean, group=plotID)) +
	geom_line(aes(x = dates, y = `Mean`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.2) +
	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`),fill="red", alpha=0.4) +

	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +facet_grid(model_name~siteID)

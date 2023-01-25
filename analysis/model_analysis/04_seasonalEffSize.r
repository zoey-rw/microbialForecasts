# Uses model parameter sizes to estimate seasonal amplitude
# Produces seasonalAmplitude.rds

# Visualize effect size estimates (beta covariates) from all model

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

sum.all <- readRDS(here("data/summary/all_fcast_effects.rds"))

sum.all_params <- sum.all %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
seas_vals <- sum.all_params %>% pivot_wider(id_cols = c("taxon","model_name","fcast_type","time_period",
																												"pretty_name","pretty_group","rank","only_rank"),
																						names_from = beta,
																						values_from = "Mean")

# Couldn't figure out how to vectorize.
out <- list()
for (i in 1:nrow(seas_vals)) {
	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
}
out <- rbindlist(out)
seas_vals <- cbind.data.frame(seas_vals, out)


keep = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
keep = keep[keep$median_gbr < 1.3,]
seas_vals_converged = seas_vals %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	filter(taxon_model_rank %in% keep$taxon_model_rank)

cycl_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
																		model_name == "cycl_only")
all_cov_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
																		 	model_name == "all_covariates")

ggplot(data = seas_vals_converged %>% filter(model_name=="cycl_only"),
			 aes(x = pretty_group,y = amplitude,
														 color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(alpha=.3, size=3, show.legend = F) +
	ylab("Seasonal amplitude") +
	theme_bw() +
	stat_compare_means(label = "p.signif", show.legend = F) + # Add significance levels
	stat_compare_means(label.y = .3, show.legend = F) +
	xlab(NULL) +
	theme_bw(base_size = 18)

e <- ggplot(data = cycl_vals, aes(x = pretty_group,y = amplitude,
														 color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(alpha=.3, size=3, show.legend = F) +
	ylab("Seasonal amplitude") +
	theme_bw() +
	stat_compare_means(label = "p.signif", show.legend = F) + # Add significance levels
	stat_compare_means(label.y = .3, show.legend = F) +
	xlab(NULL) +
	theme_bw(base_size = 18)

f <- ggplot(data = all_cov_vals, aes(x = pretty_group,y = amplitude,
																color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(alpha=.3, size=3, show.legend = F) +
	ylab("Residual seasonal amplitude") +
	theme_bw() +
	stat_compare_means(label = "p.signif", show.legend = F) + # Add significance levels
	stat_compare_means(label.x = 1.5, label.y = .65, show.legend = F) +
	xlab(NULL) +
	theme_bw(base_size = 18)


a <- ggplot(cycl_vals) +
	geom_jitter(aes(x = only_rank,y = amplitude,
									color = pretty_group),
							width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	ggtitle("Seasonal amplitude") +
	xlab(NULL) + ylab("Seasonal amplitude") +
	facet_grid(cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() + theme(text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		strip.text.y = element_text(size=12,face="bold")) +
	geom_smooth(aes(x = as.numeric(only_rank), y = amplitude), show.legend = F)

b <- ggplot(all_cov_vals) +
	geom_jitter(aes(x = only_rank,y = amplitude,
									color = pretty_group),
							width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	ggtitle("Residual seasonal amplitude (with other covariates)") +
	xlab("Rank") + ylab("Seasonal amplitude") +
	facet_grid(cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	geom_smooth(aes(x = as.numeric(only_rank), y = amplitude), show.legend = F)
ggarrange(a,b, nrow=2)

# Functional vs taxonomic effect sizes
c <- ggplot(all_cov_vals, aes(x = fcast_type,y = amplitude,
															color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F)+
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	#ggtitle("Residual seasonal amplitude") +
	xlab(NULL) +
	ylab("Residual seasonal amplitude") +
	facet_grid(cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold"))
c


# Test for differences in residual seasonal amplitude sizes between fcast types
amplitude_tukey_group = list()
for (group in c("Fungi", "Bacteria")){
	df_group=all_cov_vals %>% filter(pretty_group==!!group & fcast_type != "Diversity")
	out = df_group %>%
		#group_by(beta) %>%
		aov(amplitude~fcast_type,.) %>%
		emmeans::emmeans(object = ., pairwise ~ "fcast_type", adjust = "tukey") %>% .$emmeans %>%
		multcomp::cld(object = ., Letters = letters) %>% as.data.frame() %>%
		rename(Letters_Tukey = `.group`) %>%
		#rownames_to_column("beta") %>%
		mutate(pretty_group = !!group)
	amplitude_tukey_group[[group]] = out
}
tukey_amplitude = do.call(rbind, amplitude_tukey_group)
tukey_amplitude$tot = tukey_amplitude$upper.CL + .3
tukey_amplitude

amplitude_plot <- c +
	geom_text(data = tukey_amplitude,

						aes(x = fcast_type, y = tot, label = Letters_Tukey), show.legend = F, color = 1, size =6)
amplitude_plot

# Functional vs taxonomic effect sizes
d <- ggplot(cycl_vals, aes(x = fcast_type,y = amplitude,
															color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F)+
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	#ggtitle("Residual seasonal amplitude") +
	xlab(NULL) +
	ylab("Seasonal amplitude") +
	facet_grid(cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold"))
d
ggarrange(c,d, nrow=1)

saveRDS(list(all_cov_vals, cycl_vals), here("data/summary/seasonalAmplitude.rds"))


seas_in = readRDS(here("data/summary/seasonalAmplitude.rds"))
all_cov_vals = seas_in[[1]]
cycl_vals = seas_in[[2]]

site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))


hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds")) %>% mutate(only_rank = rank_only)
hindcast_data$month <- lubridate::month(hindcast_data$dates)
hindcast_data$month_label <- lubridate::month(hindcast_data$dates, label = T)
hindcast_data_site <- merge(hindcast_data, site_descr)  %>% filter(time_period ==  "2015-11_2020-01")

hindcast_data_site$latitude_bin = cut(hindcast_data_site$latitude, breaks = 10)

# Merge seasonality with hindcast data
# All cov model
hindcast_all_cov <- merge(hindcast_data_site %>% filter(model_name == "all_covariates"),
												all_cov_vals, by = c("taxon", "model_name", #"fcast_type",
																						 "pretty_name", "pretty_group",
																						 "only_rank"))
phylum_all_cov <-  hindcast_all_cov %>%
	filter(pretty_name == "Phylum") %>%
	distinct()

calibration_all_cov = hindcast_all_cov %>% filter(fcast_period == "calibration")


# Cycl model
hindcast_cycl_only <- merge(hindcast_data_site %>% filter(model_name == "cycl_only") %>% select(-rank_only),
														cycl_vals, by = c("taxon", "model_name", #"fcast_type",
																							 "pretty_name", "pretty_group",
																							 "only_rank"))
phylum_cycl_only <-  hindcast_cycl_only %>%
	filter(pretty_name == "Phylum") %>%
	distinct()

phylum_calibration_all_cov = phylum_all_cov %>% filter(fcast_period == "calibration")
phylum_calibration_cycl_only = phylum_cycl_only %>% filter(fcast_period == "calibration")

ggplot(phylum_all_cov %>% filter(taxon %in% c("ascomycota","basidiomycota",#"mortierellomycota",
																							"glomeromycota"#,"chytridiomycota","rozellomycota"
																							))) +
	geom_smooth(aes(x = month, y = med,
									#group = interaction(ecoregion, nlcd),

									group = siteID,
									color = ecoregion),
							size=1, alpha = .2, na.rm = T, se = F) +
	geom_point(aes(x = month, y = truth,
								 #group = ecoregion,
								 color = siteID),
						 size=1, alpha = .2) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)") +
	theme_bw(base_size = 18)

by_ecoregion <- phylum_calibration_all_cov %>% filter(taxon %in% c("ascomycota","basidiomycota",#"mortierellomycota",
																										"glomeromycota"#,"chytridiomycota","rozellomycota"
)) %>% filter(!is.na(truth)) %>% group_by(ecoregion) %>% tally() %>% arrange(n)

for_plot = phylum_calibration_all_cov %>% filter(ecoregion %in% c("D01","D03","D14"))
for_plot = calibration_all_cov %>% filter(ecoregion %in% c("D01","D03","D14"))

ggplot(for_plot %>% filter(
	taxon %in% c("acidobacteriota","actinobacteriota","firmicutes",
							 "bacteroidota","chloroflexi","cyanobacteria"))) +
	geom_smooth(aes(x = month, y = med,
									#group = interaction(ecoregion, nlcd),

									#group = ecoregion,
									color = ecoregion),
							size=1, alpha = .2, na.rm = T, se = F) +
	geom_jitter(aes(x = month, y = truth,
								 #group = ecoregion,
								 color = ecoregion),
						 size=1, alpha = .2) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)") +
	theme_bw(base_size = 18)


refit_summary <- readRDS(here("data/summary/single_taxon_summaries_201511_202001.rds"))
refit_plot_est <- refit_summary$plot_est %>% mutate(month = lubridate::month(dates),
																	month_label = lubridate::month(dates, label = T)) %>%
	filter(taxon != "other") %>% mutate(truth = as.numeric(truth))
refit_plot_est <- merge(refit_plot_est, site_descr)

southeast <- refit_plot_est %>% filter(ecoregion %in% c("D03","D02","D08"))
southeast_phyla <- southeast %>% filter(rank_name == "phylum_fun")

southeast_fg <- sum.all.fg_refit$plot_est %>% mutate(month = lubridate::month(dates),
												month_label = lubridate::month(dates, label = T)) %>%
	filter(taxon != "other") %>% mutate(truth = as.numeric(truth)) %>% filter(pretty_group == "Fungi")
southeast_fg <- merge(southeast_fg, site_descr) %>% filter(ecoregion %in% c("D03","D02","D08"))


southeast_div<- sum.div.all$plot_est %>% mutate(month = lubridate::month(dates),
																										 month_label = lubridate::month(dates, label = T)) %>%
	filter(taxon != "other") %>% mutate(truth = as.numeric(truth)) %>% filter(pretty_group == "Fungi")
southeast_div <- merge(southeast_div, site_descr) %>% filter(ecoregion %in% c("D03","D02","D08"))

calibration_asco <- phylum_all_cov %>% filter(taxon %in% c("ascomycota","basidiomycota"))
refit_asco <- refit_plot_est %>% filter(taxon %in% c("ascomycota","basidiomycota"))

ggplot() +
	geom_smooth(data = refit_asco,
							aes(x = month, y = `50%`,
									group = interaction(ecoregion, nlcd_clean),
									color = nlcd_clean),
							size=1, na.rm = T, se = F) +
	geom_jitter(data = refit_asco,
							aes(x = month, y = `50%`,
									group = interaction(ecoregion, nlcd_clean),
									color = nlcd_clean),
							size=1, alpha = .4, width = .5, height = 0) +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_grid(rows=vars(taxon), scales="free") +
	theme_bw(base_size = 18) + scale_y_sqrt() + guides(color = guide_legend("NLCD Class"))


label_df <- southeast_phyla %>%
	select(siteID, taxon, month, `50%`,nlcd_clean) %>%
	group_by(siteID, taxon) %>%
	mutate(label = if_else(!is.na(`50%`), as.character(siteID), NA_character_))  %>%
	filter(!is.na(`50%`)) %>%
	#mutate(label = if_else(month == max(month) & !is.na(`50%`), as.character(siteID), NA_character_))  %>%
	arrange(desc(month)) %>%
	filter(row_number()==1)
ggplot(data = southeast_phyla,
			 aes(x = month, y = `50%`,
			 		group = interaction(siteID, nlcd_clean),
			 		color = nlcd_clean)) +
	geom_smooth(size=1, na.rm = T, se = F) +
	geom_jitter(size=1, alpha = .4, width = .5, height = 0) +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_wrap(~taxon, scales="free") +
	theme_bw(base_size = 18) + scale_y_sqrt() + guides(color = guide_legend("NLCD Class")) +
geom_label_repel(data = label_df,
	aes(label = label),
								 nudge_x = 1,
								 na.rm = TRUE)


refit_asco <- refit_plot_est %>% group_by(siteID, taxon,model_name,group) %>%
	mutate(Mean_scale = scale(Mean)) %>%
	as.data.frame() %>% filter(rank_name == "phylum_fun")

ggplot() +
	geom_smooth(data = refit_asco,
							aes(x = month, y = Mean_scale,
									group = interaction(ecoregion, nlcd),
									color = nlcd_clean),
							size=1, na.rm = T, se = F) +
	geom_jitter(data = refit_asco,
							aes(x = month, y = Mean_scale,
									group = interaction(ecoregion, nlcd),
									color = nlcd_clean),
							size=1, alpha = .4, width = .5, height = 0) +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_grid(rows=vars(taxon), scales="free") +
	theme_bw(base_size = 18) + scale_y_sqrt()




ggplot(phylum_all_cov %>% filter(siteID %in% c("HARV","BART") &
																 	taxon %in% c("acidobacteriota","actinobacteriota","firmicutes",
																 							 "bacteroidota","chloroflexi","cyanobacteria"))) +
	geom_jitter(aes(x = month, y = med,
									color = siteID),
							size=1, alpha = .2, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)





#### -----------------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------------



vals_to_merge <- vals %>% select(-c(pretty_name))
vals_to_merge$taxon <- recode(vals_to_merge$taxon, div_ITS = "diversity_ITS", div_16S = "diversity_16S")
all_cov_vals_to_merge <- all_cov_vals %>% select(-c(pretty_name))
all_cov_vals_to_merge$taxon <- recode(all_cov_vals_to_merge$taxon, div_ITS = "diversity_ITS", div_16S = "diversity_16S")

cycl_hindcast <- hindcast_data %>% filter(time_period ==  "2015-11_2018-01") %>%
	filter(model_name == "cycl_only" & taxon != "other" & fcast_type != "Diversity")
cycl_hindcast_vals <- cycl_hindcast %>% merge(vals_to_merge)
cycl_hindcast_vals <- left_join(cycl_hindcast, vals_to_merge)
cycl_hindcast_vals <- cycl_hindcast_vals %>% mutate(med = ifelse(is.na(`50%`), med, `50%`),
																lo = ifelse(is.na(`2.5%`), lo,`2.5%`),
																hi = ifelse(is.na(`97.5%`), hi, `97.5%`))
cycl_hindcast_vals$peak_season  <- as.numeric(Hmisc::cut2(cycl_hindcast_vals$max, g=4))

all_cov_hindcast <- hindcast_data %>% filter(time_period ==  "2015-11_2018-01") %>%
	filter(model_name == "all_covariates" & taxon != "other" & fcast_type != "Diversity")
#all_cov_hindcast_vals <- all_cov_hindcast %>% merge(all_cov_vals_to_merge)
all_cov_hindcast_vals <- left_join(all_cov_hindcast, all_cov_vals_to_merge)



phylum <- cycl_hindcast_vals[cycl_hindcast_vals$pretty_name == "Phylum",] %>%
	filter(!is.na(med)) %>%
	distinct()
phylum_all_cov <-  all_cov_hindcast_vals %>% filter(time_period ==  "2015-11_2018-01") %>%
	filter(pretty_name ==
				 	"Phylum") %>%
	filter(!is.na(med)) %>%
	distinct()






ggplot(phylum) +
	geom_line(aes(x = dates, y = med,
									color = max,
								group = plotID),
						size=1, alpha = .4, show.legend = F, na.rm = T)

ggplot(phylum) +
	geom_line(aes(x = dates, y = med,
								color = peak_season,
								group = interaction(plotID, taxon)),
						size=1, alpha = .1, na.rm = T) +
	#scale_x_date(breaks = "1 month") +
	theme_bw() +
	binned_scale(aesthetics = "color",
							 scale_name = "stepsn",
							 palette = function(x) c("red", "green", "yellow", "red"),
							 breaks = c(11, 2, 5, 8),
							 show.limits = TRUE,
							 guide = "colorsteps"
	)




ggplot(phylum %>% filter(siteID=="HARV" & taxon %in% c("acidobacteriota","actinobacteriota","ascomycota","basidiomycota"))) +
	geom_line(aes(x = month, y = med,
								color = peak_season,
								group = interaction(plotID, taxon)),
						size=1, alpha = .1, na.rm = T) +
	#scale_x_date(breaks = "1 month") +
	theme_bw() +
	binned_scale(aesthetics = "color",
							 scale_name = "stepsn",
							 palette = function(x) c("red", "green", "yellow", "red"),
							 breaks = c(11, 2, 5, 8),
							 show.limits = TRUE,
							 guide = "colorsteps"
	)



ggplot(phylum %>% filter(#siteID=="HARV" &
												 	taxon %in% c("ascomycota","basidiomycota","mortierellomycota",
												 							 "glomeromycota","chytridiomycota","rozellomycota"))) +
	geom_jitter(aes(x = month, y = med,
								color = siteID),
						size=1, alpha = .1, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)


ggplot(phylum %>% filter(fcast_period=="calibration" &
	taxon %in% c("ascomycota","basidiomycota","mortierellomycota",
							 "glomeromycota","chytridiomycota","rozellomycota"))) +
	geom_jitter(aes(x = month, y = med,
									color = siteID),
							size=1, alpha = .1, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends", subtitle = "calibration estimates only") +
	theme_bw(base_size = 18)


ggplot(phylum_all_cov %>% filter(
																 	taxon %in% c("ascomycota","basidiomycota","mortierellomycota",
																 							 "glomeromycota","chytridiomycota","rozellomycota"))) +
	geom_jitter(aes(x = month, y = med,
									color = siteID),
							size=1, alpha = .2, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)





ggplot(phylum %>% filter(model_name == "cycl_only" &
												 	taxon %in% c("acidobacteriota","actinobacteriota","firmicutes",
												 							 "bacteroidota","chloroflexi","cyanobacteria"))) +
	geom_jitter(aes(x = month, y = med,
									color = siteID),
							size=1, alpha = .1, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)


ggplot(phylum_all_cov %>% filter(siteID %in% c("HARV","BART") &
												 	taxon %in% c("acidobacteriota","actinobacteriota","firmicutes",
												 							 "bacteroidota","chloroflexi","cyanobacteria"))) +
	geom_jitter(aes(x = month, y = med,
									color = siteID),
							size=1, alpha = .2, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)




scoring_metrics_in <- readRDS("./data/summary/scoring_metrics_cv.rds")
cv_metric_scaled_in <- scoring_metrics_in$cv_metric_scaled
cv_metric_scaled <- cv_metric_scaled_in %>% filter(cv_type == "mean_per_plot_cv")

scoring_metrics_cv_site <- scoring_metrics_in$scoring_metrics_cv_site
#cv_metric_scaled <- cv_metric_scaled_in %>% filter(cv_type == "mean_per_site_cv")

# hindcast_metric_merged <- merge(phylum_all_cov,
# 																cv_metric_scaled,
# 																by = c("pretty_name", "taxon",
# 																			 "pretty_group", "model_name", "fcast_type"))


hindcast_metric_merged_mean <- left_join(phylum_all_cov,
																		cv_metric_scaled)
hindcast_metric_merged_site <- left_join(phylum_all_cov,
																		scoring_metrics_cv_site)

ggplot(hindcast_metric_merged_site %>%
			 	filter(pretty_group == "Bacteria" & #metric == "CRPS" &
			 				 	!is.na(per_site_cv))) +
	geom_jitter(aes(x = month, y = med,
									#color = siteID,
									color = per_site_cv),
							size=1, alpha = .2, na.rm = T) +
	#scale_x_date(breaks = "2 month") +
	ylab("Abundance estimates") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)", subtitle = "calibration & hindcast estimates") +
	theme_bw(base_size = 18)





unique_df <- hindcast_metric_merged_mean %>% filter(!is.na(CV_scale)) %>%
	select("siteID","fcast_type","rank_name", "species", "rank", "taxon", "fg_cat", "category",
				 "taxon_name",
				 "rank.name", "only_rank", "newsite", "pretty_name", "MAT", "MAP",
				 "nlcd", "ecoregion", "sin", "cos", "max", "amplitude", "metric",
				 "score", "cv_type", "cv", "CV_scale", "metric_scale") %>%
	distinct(.keep_all = T)



unique_df <- hindcast_metric_merged_site %>% #filter(!is.na(CV_scale)) %>%
	select("siteID","fcast_type","rank_name", "species", "rank", "taxon", "fg_cat", "category",
				 "taxon_name",
				 "rank.name", "only_rank", "newsite", "pretty_name", "MAT", "MAP",
				 "nlcd", "ecoregion", "sin", "cos", "max", "amplitude",
				 "per_site_cv", "CRPS","RMSE") %>%
	distinct(.keep_all = T)

ggplot(unique_df) + # %>%
			 	# filter(#pretty_group == "Bacteria" & #
			 	# 	metric == "RMSE")) +
	geom_point(aes(x = RMSE, y = per_site_cv,
								 #color = siteID,
								 color = ecoregion),
						 size=1, alpha = .2, na.rm = T)

ggplot(unique_df %>%
			 	filter(#pretty_group == "Bacteria" & #
			 		metric == "RMSE")) +
	geom_point(aes(x = CV_scale, y = metric_scale,
									#color = siteID,
									color = ecoregion),
							size=1, alpha = .2, na.rm = T) #+ facet_wrap(~taxon, scales = "free")

summary(lm(RMSE ~ ecoregion + per_site_cv, df_merged))
summary(lm(metric ~ ecoregion + CV_scale, hindcast_metric_merged_mean[hindcast_metric_merged_mean$metric=="RMSE",], na.action = "omit"))


summary(lm(amplitude ~ ecoregion + per_site_cv, unique_df, na.action = na.omit))

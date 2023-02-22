# Compare functional group forecasts by category/kingdom

pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Read in hindcasts
hindcast_201511_201801 <- readRDS("./data/summary/beta_hindcast_fg_2015-11_2018-01.rds")
#hindcast_201511_201801 <- readRDS("./data/summary/hindcast_fg.rds")

# Read in hindcast scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))

# Subset to functional groups
fcast_info_simple <- scores_list$scoring_metrics_site_lon %>% ungroup %>%
	select(fcast_type, pretty_group, model_name, pretty_name, taxon) %>%
	distinct()
fg_rsq <- scores_list$scoring_metrics_long %>%
	filter(metric %in% c("RSQ.1","RSQ")) %>%
	mutate(score = ifelse(score < 0, 0, score)) %>%
	filter(model_name == "all_covariates" & fcast_type == "Functional group") %>%
	distinct()  %>%
	merge(fcast_info_simple, all.x=T, all.y=F)
fg_rsq$fg_category <- microbialForecast:::assign_fg_categories(fg_rsq$taxon)
fg_rsq$fg_source <- assign_fg_sources(fg_rsq$taxon)



# Looking at assignment method
ggplot(fg_rsq %>%
			 	filter(metric %in% "RSQ.1" & site_prediction == "New time (observed site)"),
			 aes(x = fg_source, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.5) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw(base_size = 18) + facet_grid(metric~site_prediction) +
	theme(
		axis.text.x=element_text(
			angle = 320, vjust=1, hjust = -0.05)	)

# Looking at functional group category
ggplot(fg_rsq %>%  filter(metric %in% "RSQ.1" & site_prediction == "New time (observed site)"), aes(x = fg_category, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter( size=3, width = .1, height=0, alpha=.1) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw()  + theme_bw(base_size = 18) +
	facet_grid(metric~site_prediction) +
	theme(
		axis.text.x=element_text(
			angle = 320, vjust=1, hjust = -0.05)	)


# Not restricted by site-pred type or metric
ggplot(fg_rsq, aes(x = fg_category, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.1) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw() + facet_grid(metric~site_prediction)

ggplot(fg_rsq, aes(x = fg_source, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.1) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw() + facet_grid(metric~site_prediction)



# truth_data <- hindcast_fg_201306_201701 %>%
# 	select(species, dateID, plotID, siteID, truth) %>%
# 	filter(!is.na(truth) & species != "other")
# hindcast_201511_201801$truth <- NULL
# hindcast_201511_201801 <- left_join(hindcast_201511_201801, truth_data)
# saveRDS(hindcast_201511_201801, "./data/summary/hindcast_fg_2015-11_2018-01.rds")


both_hindcast_periods <- plyr::rbind.fill(hindcast_fg_201306_201701, hindcast_201511_201801)
both_hindcast_periods <- plyr::rbind.fill(hindcast_201511_201801)

fungi_allcov <- both_hindcast_periods  %>%
	filter(species %in% c("animal_pathogen", "ectomycorrhizal", "endophyte", "lichenized",
												"plant_pathogen", "saprotroph"),
				 model_name == "all_covariates")

fungi_cycl <- both_hindcast_periods  %>%
	filter(species %in% c("animal_pathogen", "ectomycorrhizal", "endophyte", "lichenized",
												"plant_pathogen", "saprotroph"),
				 model_name == "cycl_only")


# Get best plots for examples (most calibration AND validation points)
not_na_hindcast <- both_hindcast_periods %>% filter(fcast_period == "hindcast" & !is.na(truth))
not_na_calibration <- both_hindcast_periods %>% filter(fcast_period == "calibration" & !is.na(truth))
plot_hindcast_counts <- sort(table(not_na_hindcast$plotID))
plot_calibration_counts <- sort(table(not_na_calibration$plotID))
top_hindcast_plots <- names(tail(plot_hindcast_counts, 30))
top_calibration_plots <- names(tail(plot_calibration_counts, 30))
top_plots <- intersect(top_hindcast_plots, top_calibration_plots)
top_plots



ggplot(fungi_allcov %>% filter(plotID=="HARV_013")) +
	facet_grid(rows=vars(species),
	cols = vars(time_period), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')




top_plots
fungi_allcov_top_plots <- fungi_allcov %>%
	filter(plotID %in% top_plots & time_period == "2015-11_2018-01")
pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_fungal_functional.pdf", height = 10)
for(i in 1:length(top_plots)){
	print(ggplot(fungi_allcov_top_plots) +
					# facet_grid(rows=vars(taxon),
					# 					 cols = vars(time_period), drop=T, scales="free") +
					geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
					geom_line(aes(x = dates, y = `50%`), show.legend = F) +
					geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
					geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
					theme_bw()+
					scale_fill_brewer(palette = "Paired") +
					theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
								legend.position = "bottom",legend.title = element_text(NULL),
								plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
					geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
					facet_grid_paginate(
						species~plotID,
						drop=T, scales="free",
						ncol = 1, nrow = 6, page = i))
}
dev.off()



library(scoringRules)
both_hindcast_periods$newsite <- ifelse(both_hindcast_periods$new_site, "New site", "Observed site")

scored_hindcasts <- both_hindcast_periods %>%
	filter(!is.na(truth) & fcast_period == "hindcast") %>%
	mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd))


scored_hindcasts_mean <- scored_hindcasts %>% group_by(fcast_type, category,group,time_period, model_name, newsite,species) %>%
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))

ggplot(scored_hindcasts_mean %>% filter(newsite=="Observed site")) +
	facet_grid(~model_name) +
	geom_violin(aes(x = time_period, y = crps_mean), draw_quantiles = c(.5)) +
geom_jitter(aes(x = time_period, y = crps_mean, color = category))


ggplot(scored_hindcasts_mean  %>%  filter(model_name=="all_covariates" &
																						time_period == "2015-11_2018-01" &
																						newsite=="Observed site")) +
	geom_violin(aes(x = category, y = crps_mean), draw_quantiles = c(.5)) +
	geom_jitter(aes(x = category, y = crps_mean, color = category))

# Only 16S groups since ITS is own category

ggplot(scored_hindcasts_mean  %>%  filter(model_name=="all_covariates" &
                                            time_period == "2015-11_2018-01" &
                                            newsite=="Observed site" &
                                            group == "16S")) +
  geom_violin(aes(x = reorder(category, crps_mean), y = crps_mean), draw_quantiles = c(.5)) +
  geom_jitter(aes(x = reorder(category, crps_mean), y = crps_mean, color = category), size=3, alpha = .8, show.legend = F) +
  theme_bw(base_size = 18) +xlab("Functional category") +ylab("CRPS score") +
  ggtitle("Predictability of bacterial functional categories")


# Skill score from Scavia 2021
skill_score <- scored_hindcasts_mean %>%
	pivot_wider(id_cols = c("fcast_type","category","group","model_name","time_period","species"), values_from = "crps_mean", names_from = "newsite") %>%
	mutate(skill_score = (1 - (`New site`/`Observed site`)))


ggplot(skill_score %>% filter(model_name=="all_covariates" & time_period == "2015-11_2018-01")) +
	#facet_grid(~model_name) +
	geom_violin(aes(x = category, y = skill_score), draw_quantiles = c(.5)) +
	geom_jitter(aes(x = category, y = skill_score, color = category))






ggplot(hindcast_201511_201801 %>% filter(plotID=="BART_002" & taxon == "oligotroph")) +
	facet_grid(#rows=vars(taxon),
		cols = vars(model_name), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')


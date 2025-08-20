# Compare functional group forecasts by category/kingdom

pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce, ggrepel)
source("source.R")

# Read in hindcasts
hindcast_data <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))

# Read in hindcast scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr2.rds")))

# Subset to functional groups
fcast_info_simple <- scores_list$scoring_metrics_site_lon %>% ungroup %>%
	select(fcast_type, pretty_group, model_name, pretty_name, taxon) %>%
	distinct()

fg_rsq <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% scores_list$converged_strict_list) %>%
	filter(metric %in% c("RSQ.1","RSQ","RMSE.norm")) %>%
	mutate(score = ifelse(score < 0, 0, score)) %>%
	filter(#model_name == "env_cycl" &
				 	pretty_name == "functional") %>%
	distinct()  %>%
	merge(fcast_info_simple, all.x=T, all.y=F)
fg_rsq$fg_category <- microbialForecast:::assign_fg_categories(fg_rsq$taxon) %>% make.names
fg_rsq$fg_source <- assign_fg_sources(fg_rsq$taxon) #%>% make.names



pretty_names <- list("cellulolytic" = "Cellulose degraders",
										 "assim_nitrite_reduction" = "Assimilatory nitrite reducers",
										 "dissim_nitrite_reduction" = "Dissimilatory nitrite reducers",
										 "assim_nitrate_reduction" = "Assimilatory nitrate reducers",
										 "n_fixation" = "Nitrogen fixers",
										 "dissim_nitrate_reduction" = "Dissimilatory nitrate reducers",
										 "nitrification" = "Nitrifiers",
										 "denitrification" = "Denitrifiers",
										 "chitinolytic" = "Chitin degraders",
										 "lignolytic" = "Lignin degraders",
										 "methanotroph" = "Methanotrophs",
										 "copiotroph" = "Copiotrophs",
										 "oligotroph" = "Oligotrophs",
										 "benomyl_antibiotic" = "Benomyl-resistant",
										 "glucose_simple" = "Glucose-enriched",
										 "pyruvate_simple" = "Pyruvate-enriched",
										 "streptomycin_antibiotic" = "Streptomycin-resistant",
										 "sucrose_complex"  = "Sucrose-enriched",
										 "acetogen_anaerobic" = "Acetogen anaerobic",
										 "chloramphenicol_antibiotic"  = "Chloramphenicol-resistant",
										 "erythromycin_antibiotic"  = "Erythromycin-resistant",
										 "gentamycin_antibiotic"  = "Gentamycin-resistant",
										 "glycerol_simple"   = "Glycerol-enriched",
										 "acetate_simple"  = "Acetate-enriched",
										 "acidic_stress"   = "Acidic stress-tolerant",
										 "cellobiose_complex"   = "Cellobiose-enriched",
										 "cellulose_complex"   = "Cellulose-enriched",
										 "chitin_complex"   = "Chitin-enriched",
										 "galactose_simple"   = "Galactose-enriched",
										 "xylose_simple"   = "Xylose-enriched",
										 "salt_stress" = "Salt stress-tolerant",
										 "herbicide_stress" = "Herbicide stress-tolerant",
										 "osmotic_stress" = "Osmotic stress-tolerant",
										 "heat_stress" = "Heat stress-tolerant",
										 "light_stress" = "Light stress-tolerant",
										 "arbuscular" = "Arbuscular mycorrhizae",
										 "endophyte" = "Endophyte",
										 "litter_saprotroph" = "Litter saprotrophs",
										 "lichenized" = "Lichenized fungi",
										 "animal_pathogen" = "Animal pathogens",
										 "plant_pathogen" = "Plant pathogens",
										 "saprotroph" = "Saprotrophs",
										 "wood_saprotroph" = "Wood saprotrophs",
										 "ectomycorrhizal" = "Ectomycorrhizae"
)
fg_rsq$pretty_fg_names <- recode(fg_rsq$taxon, !!!pretty_names)# %>% make.names
#fg_rsq$pretty_fg_names <- recode(fg_rsq$taxon, !!!microbialForecast:::pretty_names)


stat_pvalue_fg_source <- fg_rsq %>%
	filter(metric %in% "RMSE.norm" &
				 	site_prediction == "New time (observed site)" &
				 	model_name == "cycl_only") %>%
	rstatix::tukey_hsd(score ~ fg_source) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .02) #%>%
#	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

stat_pvalue_fg_source <- fg_rsq %>%
	group_by(model_name) %>%
	filter(metric %in% "RMSE.norm" &
				 	site_prediction == "New time (observed site)") %>%
	rstatix::tukey_hsd(score ~ fg_source) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(y.trans=log10) #%>%
#	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

stat_pvalue_fg_category <- fg_rsq %>%
	group_by(model_name) %>%
	filter(metric %in% "RMSE.norm" &
				 	site_prediction == "New time (observed site)") %>%
	rstatix::tukey_hsd(score ~ fg_category) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .02) #%>%
	#mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

pos <- position_jitter(width = 0.1, height=0, seed = 1)

# Looking at assignment method
fg_source = ggplot(fg_rsq %>%
			 	filter(metric %in% "RMSE.norm" &
			 				 	site_prediction == "New time (observed site)" &
			 				 	model_name == "cycl_only"),
			 aes(x = fg_source, y = as.numeric(score)))  +
	#geom_boxplot(alpha=.5) +
	geom_jitter(aes(color=pretty_group),
							size=4, #width = .3, height=0,
							alpha=.5,
							position=pos) +
	xlab(NULL) + labs(fill='') +
	theme_classic(base_size = 22) +
	theme(axis.text.x=element_text(size=20,
			angle = 320, vjust=1, hjust = -0.05)	) +
	ylab("Forecast error (nRMSE)") + labs(color=NULL) +
	geom_text_repel(aes(label=pretty_fg_names), size=6,
									 max.overlaps = 10,
									 position=pos) +
	ggpubr::stat_pvalue_manual(stat_pvalue_fg_source %>%
														 	filter(model_name == "cycl_only"),
														 label = "p.adj.signif",
														 #bracket.nudge.y = -1.1,
size=7, hide.ns = T) +
	scale_y_continuous(trans = "log10")


png(here("figures","functional_group_source.png"), width = 800, height=1000)
print(fg_source)
dev.off()

# Looking at functional group category
ggplot(fg_rsq %>%  filter(metric %in% "RMSE.norm" &
														site_prediction == "New time (observed site)",
													model_name == "cycl_only"),
			 aes(x = fg_category, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter( size=3, width = .1, height=0, alpha=.1) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw()  + theme_bw(base_size = 18) +
	facet_grid(metric~site_prediction) +
	theme(
		axis.text.x=element_text(
			angle = 320, vjust=1, hjust = -0.05)	) + ylab("RSQ") +
	#geom_text(aes(label=pretty_fg_names), hjust=0, vjust=0) +
	ggpubr::stat_pvalue_manual(stat_pvalue_fg_category %>%
														 	filter(model_name == "cycl_only"),
														 label = "p.adj.signif", #, bracket.nudge.y = -.4,
														 size=7, hide.ns = T)


# Not restricted by site-pred type or metric
ggplot(fg_rsq, aes(x = fg_category, y = as.numeric(score)))  +
	geom_boxplot() +
	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.1) +
	xlab(NULL) + labs(fill='') +
	stat_compare_means() + theme_bw() + facet_wrap(metric~site_prediction)

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


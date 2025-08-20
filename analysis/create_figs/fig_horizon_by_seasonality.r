# read in seasonal values
library(lubridate)
library(ggrepel)
source("source.R")
source("microbialForecast/R/assignPhenology.r")

sum.in <- readRDS(here("data", "summary/logit_beta_regression_summaries.rds"))

# Read in and split up plot-level estimates
plot_estimates = sum.in$plot_est %>% filter(model_name != "all_covariates") %>% filter(model_id %in% sum.in$keep_models$model_id)
plot_estimates$month = lubridate::month(plot_estimates$dates)
plot_estimates$year = lubridate::year(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))

model_info_key <- plot_estimates %>% select(c("model_id","fcast_type","time_period","pretty_group","rank_only","model_name","taxon")) %>% distinct()

read_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
full_phenophase_abundance = read_in[[1]] %>% merge(model_info_key, all.x=T)
pheno_info_key <- full_phenophase_abundance[,c("siteID","year","month","site_cat")] %>% distinct()
# This factor should be ordered (so that phenophases are sequential)
pheno_info_key$sampling_season =
	factor(pheno_info_key$site_cat, ordered = T, levels = c("dormancy","greenup","peak","greendown"))


# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))
site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups


site_descr$latitude_category = ifelse(site_descr$latitude > 44, "High-latitude",
																			ifelse(site_descr$latitude < 31, "Low-latitude",
																						 "Mid-latitude"))
# This factor should be ordered
site_descr$latitude_category =
	factor(site_descr$latitude_category, ordered = T, levels = c("Low-latitude",  "Mid-latitude","High-latitude"))
pheno_info_key <- merge(pheno_info_key, site_descr, all=T)


plot_phenophase_abundance <- merge(plot_estimates, pheno_info_key, all.x=T)


plot_phenophase_abundance_cal <- plot_phenophase_abundance %>%  #full_phenophase_abundance %>%
	filter(time_period=="2015-11_2018-01")

for_stats <- plot_phenophase_abundance_cal  %>%
	filter(!is.na(`50%`) & !is.na(site_cat))
tukey_median_pheno = for_stats %>%
	group_by(fcast_type,pretty_group,model_name,taxon,model_id) %>%
	summarize(tukey(x = sampling_season, y = `50%`, y.offset = 0)) %>%
	rename(sampling_season = x)
tukey_median_pheno_sig <- tukey_median_pheno %>%
	group_by(fcast_type,pretty_group,model_name,taxon,model_id) %>%
	mutate(significant_diff = ifelse(n_distinct(Letters_Tukey) > 1, T, F)) %>% select(model_id, significant_diff) %>% distinct()

# Get percent of groups that differ across phenophases: 91%
tukey_median_pheno_sig %>% #filter(model_name=="env_cycl") %>%
	ungroup() %>%
	group_by(model_name) %>%
	add_count(name = "total") %>%
	group_by(model_name,significant_diff) %>%
	add_count(name="group") %>% select(model_name,group,total) %>%
	distinct %>% mutate(pct = group/total) %>% filter(significant_diff) %>% print()

# Get percentage of site-times assigned to each category
for_stats %>% select(siteID, plotID, dateID, sampling_season) %>% distinct() %>%
	select(sampling_season) %>% table()


horizon_in = readRDS(here("data", paste0("summary/fcast_horizon_clean.rds")))
fcast_horizon_x_RSQ = horizon_in[[2]]
fcast_horizon_x_CRPS = horizon_in[[4]]

fcast_horizon_long_seas <- merge(fcast_horizon_x_CRPS, tukey_median_pheno_sig, all.y = T, by=c("taxon","pretty_group", "model_name", "model_id")) %>%
	filter(!is.na(pretty_group)) %>% mutate(significant_diff = recode(as.character(significant_diff),
																																		"TRUE" = "Microbes that vary \nwith plant phenophase",
																																		"FALSE" = "Microbes that do not vary\n with plant phenophase"))



stat_pvalue <- fcast_horizon_long_seas %>%
	group_by(model_name) %>%
	rstatix::t_test(months_since_obs ~ significant_diff,detailed = T) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .2) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

fig3g <- ggplot(fcast_horizon_long_seas %>% filter(model_name %in% c("cycl_only")),
																				#										"env_cycl")),# %>% filter(parameter_type=="horizon" & horizon_parameter=="rsq_fcast_horizon"),
			 aes(x = significant_diff,
			 		y = months_since_obs)) +

	geom_violin(draw_quantiles=c(.5), alpha=.5, show.legend = F, outlier.shape = NA) +
	#facet_grid(~model_name) +
	geom_point(aes(color = pretty_group),
		size=3, alpha=.4, #show.legend = F,
		position=position_jitter(height=.2, width=.2),
		show.legend = T
		#position=position_jitterdodge(jitter.width = .2, jitter.height = 0, dodge.width = 1)
		) + #facet_grid(#pretty_group
																																			#					~model_name) +
	coord_flip() + theme_bw(base_size = 16) +
	ylab("Forecast horizon (months of skilled predictions)") + xlab(NULL) +
	#ggtitle("Forecast horizon") +
	ggpubr::stat_pvalue_manual(stat_pvalue,
																													 label = "p = {p}",
																													 bracket.nudge.y = -.1,
																													 size=4, hide.ns = T) +
	labs(color="Kingdom") + #	scale_color_discrete(name = "Kingdom") +
	theme(legend.box="horizontal",legend.position = "top")
#																										legend.position=c(.7,.15)

#fig3g <- tag_facet(fig3g, tag_pool = "G", y = 5)

fig3g
fig3g <- ggarrange(fig3g, labels="G")
png(here("figures","fcast_horizon_by_seasonality.png"), width = 1200, height=400)
print(fig3g)
dev.off()

long_horizon_subset1 = fcast_horizon_x_RSQ[which(fcast_horizon_x_RSQ$months_since_obs > 10),]
long_horizon_subset2 = fcast_horizon_x_CRPS[which(fcast_horizon_x_CRPS$months_since_obs > 10),]
long_horizon_subset <- inner_join(long_horizon_subset1, long_horizon_subset2)
intersect(long_horizon_subset1$model_id, long_horizon_subset2$model_id)
print(long_horizon_subset)


fit_combined = horizon_in[[1]]


fg_list =  c(#"ectomycorrhizal",
						 "plant_pathogen")
						 #"saprotroph","endophyte",#"oligotroph",
					#	 "copiotroph")
stats_for_plot = tukey_median_pheno %>% filter(taxon %in% fg_list & model_name=="cycl_only")
stats_for_plot$tot = ifelse(stats_for_plot$tot > .55, .55, stats_for_plot$tot)
plot_median_abun <- ggplot(for_stats %>%
													 	filter(taxon %in% fg_list & model_name=="cycl_only"),
													 aes(y = `50%`,color = sampling_season,
													 		x = sampling_season)) +

	geom_violin(alpha=.5, show.legend = F, outlier.shape = NA) +
	# geom_point(#aes(color = siteID),
	# 	size=1, alpha=.1, #show.legend = F,
	# 	position=position_jitter(height=0), show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	#facet_grid(rows=vars(taxon), scales="free") +
	facet_wrap(~taxon, scales="free") +
	ylab("Median abundance across all sites") +
	xlab("Plant phenophase")  +
	geom_text(data = stats_for_plot,
						aes(x = sampling_season, y = .1, label = Letters_Tukey),
						show.legend = F, color = 1, size =6) + ylim(c(0,.15))
plot_median_abun



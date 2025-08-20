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

phenophase_fg_abundance_fungi <- plot_phenophase_abundance %>%  #full_phenophase_abundance %>%
	filter(pretty_group=="Fungi" & fcast_type == "Functional" &
				 	time_period=="2015-11_2018-01")


for_stats <- plot_phenophase_abundance_cal  %>%
	filter(!is.na(`50%`) & !is.na(site_cat))
tukey_median_pheno = for_stats %>%
	group_by(fcast_type,model_name,taxon,model_id) %>%
	summarize(tukey(x = sampling_season, y = `50%`)) %>%
	rename(sampling_season = x)
tukey_median_pheno_sig <- tukey_median_pheno %>%
	group_by(fcast_type,model_name,taxon,model_id) %>%
	mutate(significant_diff = ifelse(n_distinct(Letters_Tukey) > 1, T, F)) %>% select(model_id, significant_diff) %>% distinct()

# core_cal_hindcast_seas <- merge(core_cal_hindcast,
# 																tukey_median_pheno_sig, all = T, by=c("taxon", "model_name", "model_id", "fcast_type"))


fcast_horizon = readRDS(here("data/summary/fcast_horizon_df_core.rds"))

core_horizon_seas <- merge(fcast_horizon[[1]],
																tukey_median_pheno_sig, all.y = T, by=c("taxon", "model_name", "model_id"))
ggplot(core_horizon_seas,
			 aes(x = significant_diff,#color = latitude_category,
			 		y = rsq_fcast_horizon)) +

	geom_boxplot(alpha=.5, show.legend = F) +
	geom_point(#aes(color = siteID),
		size=1, alpha=.3, #show.legend = F,
		position=position_jitter(height=0), show.legend = F) + facet_grid(pretty_group~model_name) +coord_flip()


for_stats <- phenophase_fg_abundance_fungi  %>%
	filter(model_name=="cycl_only" & !is.na(`50%`) & !is.na(site_cat) & !is.na(pretty_group))
for_stats <- merge(for_stats, site_descr, all.x=T)


tukey_median_pheno_evergreen = for_stats %>%
	group_by(fcast_type,model_name,taxon,latitude_category) %>%
	summarize(tukey(x = sampling_season, y = `50%`)) %>%
	rename(sampling_season = x)
plot_median_abun <- ggplot(for_stats %>% filter(taxon %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
													 aes(y = `50%`,#color = latitude_category,
													 		x = sampling_season)) +

	geom_boxplot(alpha=.5, show.legend = F) +
	geom_point(#aes(color = siteID),
						 size=1, alpha=.3, #show.legend = F,
						 position=position_jitter(height=0), show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~taxon, scales="free") +
	ylab("Median abundance across all sites") +
	xlab("Plant phenophase")  +
	geom_text(data = tukey_median_pheno_evergreen %>% filter(taxon %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
						aes(x = sampling_season, y = tot-.1, label = Letters_Tukey),
						show.legend = F, color = 1, size =4)
plot_median_abun




site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups




core_scores_in = readRDS(here("data", paste0("summary/scoring_metrics_core_level.rds")))
core_cal_hindcast = core_scores_in[[4]]


# Wide format taxonomic
truth_vals <- core_cal_hindcast %>%
	filter(!is.na(core_truth)) %>%
	filter(model_name == "env_cycl")

truth_vals <- merge(truth_vals, site_descr, all.x=T, by="siteID")

pheno_info_key <- full_phenophase_abundance[,c("siteID","year","month","site_cat")] %>% distinct()

truth_vals$month = month(truth_vals$dates)
truth_vals$year = year(truth_vals$dates)
truth_vals <- merge(truth_vals, pheno_info_key, all.x=T, by=c("siteID","year","month"))
# This factor should be ordered (so that phenophases are sequential)
truth_vals$sampling_season =
	factor(truth_vals$site_cat, ordered = T, levels = c("dormancy","greenup","peak","greendown"))

saveRDS(truth_vals, here("data", paste0("summary/truth_vals_phenophase.rds")))


ggplot(truth_vals %>%
			 	filter(!is.na(sampling_season)) %>%
			 	filter(taxon_name %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
													 aes(y = core_truth,color = latitude_category,
													 		x = sampling_season)) +

	geom_boxplot(alpha=.5, show.legend = F) +
	geom_point(#aes(color = siteID),
		size=1, alpha=.3, #show.legend = F,
		position=position_jitter(height=0), show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(taxon~latitude_category, scales="free") +
	ylab("Median abundance across all sites") +
	xlab("Plant phenophase")
	# geom_text(data = tukey_median_pheno_evergreen %>% filter(taxon %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
	# 					aes(x = sampling_season, y = tot-.1, label = Letters_Tukey),
	# 					show.legend = F, color = 1, size =4)




scores_list = readRDS(here("data", "summary/scoring_metrics_cv.rds"))


site_cv = scores_list$scoring_metrics_cv_site
site_cv_simple <- merge(site_cv, site_descr, all.x=T, by="siteID") %>%
	select(-c(metric,score,mean_crps_sample)) %>% distinct()



ggplot(site_cv_simple,
			 aes(y = per_site_cv,color = pretty_group,
			 		x = latitude)) +

	geom_boxplot(alpha=.5, show.legend = F) +
	geom_point(#aes(color = siteID),
		size=1, alpha=.3, #show.legend = F,
		position=position_jitter(height=0), show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) #+
	#facet_grid(taxon~latitude_category, scales="free")

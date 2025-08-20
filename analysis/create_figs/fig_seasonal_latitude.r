library(ggallin)
library(ggpmisc) # for polynomial plotting function
library(broom)
source("source.R")


# Scoring metrics aggregated by taxon/model
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr2.rds")))
converged = scores_list$converged_list
converged = scores_list$converged_strict_list

hindcast_filter <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(model_name != "all_covariates" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm") &
				 	site_prediction == "New time (observed site)")

# Raw hindcast crps values
hindcast_data <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))
hindcast_data_df = hindcast_data %>%
	filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period=="hindcast" & new_site==FALSE)

site_scores_allmetrics = scores_list$scoring_metrics_site_long %>%
	filter(!siteID %in% "MLBS") %>%
	filter(model_id %in% converged &
				 	site_prediction == "New time (observed site)")


# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))

site_scores_allmetrics <- merge(site_scores_allmetrics, site_descr)


# Predictability for latitude - fungi and bacteria
# Using Rsq to compare across groups at site level
# Because RMSE.norm is biased by zero-abundance values at the site level
ggplot(site_scores_allmetrics %>%
			 	filter(model_name=="env_cycl") %>%
			 	filter(metric=="RSQ"), aes(x = latitude,
			 														 y = score)) +
	geom_point(aes(color=pretty_group), alpha=.4, size=3,
						 position=position_jitter(height=0, width = .3), show.legend = F) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") + ylab("RSQ") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")),
					 size=6, label.x.npc = .01) +
	facet_grid(rows=vars(pretty_group), scales = "free_y")

# Model predictability by latitude for individual groups
# Fig S9
model_results = site_scores_allmetrics %>% filter(metric=="RMSE.norm") %>%
	# filter(model_id %in% c("env_cov_aspergillus_20151101_20180101","cycl_only_sordariales_20151101_20180101","env_cycl_aspergillus_20151101_20180101")) %>%
	group_by(pretty_group,model_name,model_id) %>%
	nest() %>%
	mutate(model = map(data, lm(score ~ latitude, data = .)),
				 tidied = map(model, tidy))

# actually swtiching back to RSQ to make the figs consistent
model_results = site_scores_allmetrics %>% filter(metric=="RSQ") %>%
	group_by(pretty_group,model_name,model_id) %>%
	mutate(glance(lm(score ~ latitude))) %>%
	dplyr::select(siteID:taxon,score,latitude,adj.r.squared,p.value)
sig_model_results = model_results %>% filter(p.value < .1)

# Fig S10
ggplot(sig_model_results %>%
			 	filter(model_name=="env_cycl"),
			 # "salt_stress","oligotroph",
			 # "heat_stress")),
			 aes(x = latitude, color=pretty_group,
			 		y = score)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm", se=F) +
	theme_minimal(base_size = 20) + labs(color = "Kingdom") +
	xlab("Latitude") + ylab("Forecast accuracy (RSQ)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")),
					 size=5, label.x.npc = .45, color=1, p.accuracy = .001) +
	facet_wrap(pretty_name~taxon, scales = "free") + theme(legend.position = c(.9,.1))







##### OLD #####


seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))

seas_amplitude_long <- seas_in[[1]] %>% filter(model_id %in% converged)
seas_amplitude_long$is_seasonal <- ifelse(seas_amplitude_long$significant_sin==1 |
																						seas_amplitude_long$significant_cos == 1, T, F)

# seas_amplitude_simple <- seas_amplitude_long %>% distinct(model_id, is_seasonal)
high_seas_taxa = seas_amplitude_long %>%
	#filter( model_name=="cycl_only") %>%
	filter( model_name=="cycl_only" & amplitude > .05) %>%
	filter(significant_sin==1 | significant_cos == 1)

seas_taxa = unique(high_seas_taxa$taxon)
ggplot(high_seas_taxa) + geom_point(aes(x = rank_only, y = amplitude, color = pretty_group))

site_scores_allmetrics$is_seasonal <- ifelse(site_scores_allmetrics$taxon %in% seas_taxa, 1, 0)
ggplot(site_scores_allmetrics %>%
			 	filter(model_name=="env_cycl") %>%
			 	filter(metric=="RSQ"), aes(x = latitude,
			 														 y = score)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_wrap(~is_seasonal, scales = "free")


core_scores_in = readRDS(here("data", paste0("summary/scoring_metrics_core_level.rds")))

core_site_scores_allmetrics = core_scores_in[[6]] %>% filter(model_id %in% converged)
core_site_scores_allmetrics$is_seasonal <- ifelse(core_site_scores_allmetrics$taxon %in% seas_taxa, 1, 0)
core_site_scores_allmetrics <- merge(core_site_scores_allmetrics, site_descr, all.x=T)

ggplot(core_site_scores_allmetrics %>%
			 	filter(model_name=="env_cycl"), aes(x = latitude,
			 														 y = RMSE.norm)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_grid(site_prediction~is_seasonal, scales = "free")


## Fig for ESA 2023
ggplot(core_site_scores_allmetrics %>%
			 	filter(model_name=="env_cycl" & fcast_type=="Functional") %>%
			 	filter(taxon %in% c("ectomycorrhizal","plant_pathogen",
			 											"endophyte",
			 											"lignolytic")),
			 											# "salt_stress","oligotroph",
			 											# "heat_stress")),
			 aes(x = latitude,
			 																																 y = RMSE.norm)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm", se=F) +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") + ylab("Forecast error (nRMSE)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_wrap(~taxon, scales = "free", ncol=1)



plot_site_scores_allmetrics = core_scores_in[[6]] %>% filter(model_id %in% converged)
plot_site_scores_allmetrics$is_seasonal <- ifelse(plot_site_scores_allmetrics$taxon %in% seas_taxa, 1, 0)
plot_site_scores_allmetrics <- merge(plot_site_scores_allmetrics, site_descr, all.x=T)

ggplot(plot_site_scores_allmetrics %>%
			 	filter(model_name=="env_cycl"), aes(x = latitude,
			 																			y = mean_crps_sample)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_wrap(pretty_group~is_seasonal, scales = "free")


site_descr$only_deciduous_evergreen = ifelse(site_descr$deciduous=="Deciduous" & site_descr$evergreen == "",
																				 "Deciduous",
																				 ifelse(site_descr$evergreen=="Evergreen" & site_descr$deciduous == "", "Evergreen", ""))



site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups
site_descr$latitude_category = ifelse(site_descr$latitude > 44, "High-latitude",
																			ifelse(site_descr$latitude < 31, "Low-latitude",
																						 "Mid-latitude"))
# This factor should be ordered
site_descr$latitude_category =
	factor(site_descr$latitude_category, ordered = T, levels = c("Low-latitude",  "Mid-latitude","High-latitude"))

#####

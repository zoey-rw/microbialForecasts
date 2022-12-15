# Read in model estimates for seasonal cycles, with and without other covariates (temperature, moisture, etc)
seas_in = readRDS(here("data/summary/seasonalAmplitude.rds"))
# Combine seasonality estimates from two model versions
seasonality_estimates = rbind(seas_in[[1]],seas_in[[2]])

# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))
site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups

# Read in plot-level model states from longer calibration period ("refit")
# All the data up till Jan 2020 has been assimilated
refit_summary <- readRDS(here("data/summary/single_taxon_summaries_201511_202001.rds"))
refit_fg_summary <- readRDS(here("data/summary/fg_summaries_20151101_20200101.rds"))

# Add better date labels
refit_plot_est <- rbindlist(list(refit_summary$plot_est, refit_fg_summary$plot_est %>% mutate(rank_name = "Functional")), fill = T)


refit_plot_est <- refit_plot_est %>%
	mutate(month = lubridate::month(dates),
				 month_label = lubridate::month(dates, label = T)) %>%
	filter(!grepl("other", taxon)) %>% select(-c("site_num","plot_num","species_num","date_num")) %>%
	distinct(.keep_all = T)




# Read in hindcast data
hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds")) %>% mutate(only_rank = rank_only) %>%
	mutate(month = lubridate::month(dates),
				 month_label = lubridate::month(dates, label = T))
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




# Subset
phylum_all_cov <-  hindcast_all_cov %>%
	filter(pretty_name == "Phylum") %>%
	distinct()
phylum_cycl_only <-  hindcast_cycl_only %>%
	filter(pretty_name == "Phylum") %>%
	distinct()

for_plot = phylum_calibration_all_cov %>% filter(ecoregion %in% c("D01","D03","D14"))



# Add to estimated values for each plot/date
# Takes a minute to merge
refit_estimates <- merge(refit_plot_est, site_descr, by="siteID") %>%
	filter(time_period ==  "2015-11_2020-01") %>%
	merge(seasonality_estimates %>% select(-fcast_type), by = c("taxon", "model_name", "pretty_group"))


scoring_metrics_in <- readRDS("./data/summary/scoring_metrics_cv.rds")
#taxon_site_cv = scoring_metrics_in$scoring_metrics_cv_site %>% select(taxon, siteID, per_site_cv) %>% distinct()


# Get variation per site, then average per taxon

truth_vals <- refit_estimates %>%
	filter(model_name == "all_covariates" & !is.na(truth))
cv_tax_per_site <- truth_vals %>%
	group_by(fcast_type, pretty_group, rank_name, taxon, siteID, plotID) %>%
	summarize(per_plot_cv = calc_cv(truth),
						median_plot_abundance = median(truth, na.rm=T)) %>%
	ungroup %>%
	group_by(fcast_type, pretty_group, rank_name, taxon, siteID) %>%
	dplyr::summarize(per_site_cv = mean(per_plot_cv, na.rm=T),
									 median_site_abundance = median(median_plot_abundance, na.rm=T))


taxon_site_cv = merge(cv_tax_per_site, site_descr, all.x=T)  %>%
	group_by(fcast_type, pretty_group, rank_name, taxon, siteID)
taxon_site_cv

ecoregion = "D01"
for (ecoregion in unique(taxon_site_cv$ecoregion)){

	ecoregion_cv = taxon_site_cv %>%
		filter(ecoregion == !!ecoregion) %>%
		group_by(fcast_type, pretty_group, rank_name, taxon) %>%
		summarize(mean_cv = mean(per_site_cv, na.rm=T),
							median_ecoregion_abundance = median(median_site_abundance, na.rm=T))
	top_cv = ecoregion_cv %>%
		ungroup() %>%
		filter(median_ecoregion_abundance > .1) %>%
		group_by(fcast_type) %>%
		slice_max(order_by = mean_cv, n = 2)

}


seasonality_estimates[seasonality_estimates$model_name=="cycl_only",]  %>%
	group_by(fcast_type) %>%
	slice_max(order_by = amplitude, n = 2)

ggplot(refit_estimates %>% filter(taxon %in% top_cv$taxon &
																		ecoregion %in% !!ecoregion)) +
	geom_smooth(aes(x = month, y = Mean,
									color = siteID),
							size=1, alpha = .2, na.rm = T, se = F, span=.3) +
	geom_jitter(aes(x = month, y = truth,
								 color = siteID),
						 size=3, alpha = .2) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_wrap(pretty_group~taxon, scales="free", drop=T) +
	ggtitle("Seasonal trends (all-covariate model)") +
	theme_bw(base_size = 18)

sum.all <- readRDS(here("data/summary/all_fcast_effects.rds"))
sum.all_params <- sum.all %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
seas_vals <- sum.all_params %>% pivot_wider(id_cols = c("taxon","model_name","fcast_type","time_period",
																																"pretty_name","pretty_group","only_rank"),
																										names_from = beta,
																										values_from = "Mean")

# Couldn't figure out how to vectorize.
out <- list()
for (i in 1:nrow(seas_vals)) {
	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
}
out <- rbindlist(out)
seas_vals <- cbind.data.frame(seas_vals, out)


most_seasonal = seas_vals[seas_vals$model_name=="cycl_only",]  %>%
	filter(fcast_type != "Diversity") %>%
	filter(time_period == "2015-11_2018-01" | fcast_type == "Functional") %>%
	select(-time_period) %>%
	group_by(pretty_group, fcast_type) %>%
	slice_max(order_by = amplitude, n = 10)
most_seasonal = merge(most_seasonal, ecoregion_cv %>% ungroup %>% select(-fcast_type))
most_seasonal[most_seasonal$median_ecoregion_abundance > .01,]
most_seasonal_taxon = c("chytridiomycota","chitin_complex","caulobacterales","endophyte")


northeast_taxon = c("saprotroph",
										"actinobacteriota",
										"cellulose_complex",
										"ascomycota",
										"basidiomycota",
										"n_fixation")

ggplot(refit_estimates %>% filter(taxon %in% northeast_taxon &
																		ecoregion %in% c("D01"))) +
	geom_smooth(aes(x = month, y = `50%`,
									color = pretty_group),
							size=1, na.rm = T, se = F) +
	geom_jitter(aes(x = month, y = truth,
									color = pretty_group),
							size=2, alpha = .3, height = 0) +
	geom_jitter(aes(x = month, y = `50%`,
									color = pretty_group),
							size=1, alpha = .05, height = 0) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_wrap(pretty_group~taxon, scales="free", drop=T) +
	ggtitle("Seasonal trends in Northeast US") +
	theme_bw(base_size = 18)


ggplot(refit_estimates %>% filter(taxon %in% northeast_taxon &
																		ecoregion %in% c("D01","D14"))) +
	geom_smooth(aes(x = month, y = `50%`,
									color = ecoregion),
							size=1, na.rm = T, se = F) +
	geom_jitter(aes(x = month, y = truth,
									color = ecoregion),
							size=2, alpha = .3, height = 0) +
	geom_jitter(aes(x = month, y = `50%`,
									color = ecoregion),
							size=1, alpha = .01, height = 0) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	facet_wrap(pretty_group~taxon, #scales="free",
						 drop=T) +
	ggtitle("Seasonal trends in Southeast and Northeast US") +
	theme_bw(base_size = 18)

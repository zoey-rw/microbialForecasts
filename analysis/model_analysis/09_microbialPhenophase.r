# read in seasonal values
library(lubridate)
library(ggrepel)
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/assignPhenology.r")

sum.in <- readRDS(here("data", "summary/logit_beta_regression_summaries.rds"))
seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))

seas_amplitude_long <- seas_in[[1]] %>% filter(model_id %in% sum.in$keep_models$model_id)

# Extract groups with "max" dates in winter
seas_amplitude_long$max_month = month(seas_amplitude_long$max_y_date)
seas_amplitude_long$simulated_month = month(seas_amplitude_long$dates)
seas_amplitude <- seas_amplitude_long %>%
	select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)
seas_amplitude_cycl_only = seas_amplitude %>%
	filter(model_name == "cycl_only")

# Read in and split up plot-level estimates
plot_estimates = sum.in$plot_est %>% filter(model_name != "all_covariates") %>% filter(model_id %in% sum.in$keep_models$model_id)

plot_estimates$month = lubridate::month(plot_estimates$dates)
plot_estimates$year = lubridate::year(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cov_est = plot_estimates %>% filter(grepl("env_cov",model_name))

# Read in reference dataframe of phenophase categories for all site/dates
pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
pheno_categories_long = pheno_categories_in[[2]] %>%
	filter(ID %in% unique(plot_estimates$siteID))

model_info_key <- seas_amplitude_long %>% select(c("model_id","fcast_type","time_period","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")) %>% distinct()


max_dates = seas_amplitude# %>% mutate(dates = max_y_date)
#max_dates$sampling_season =  apply(max_dates, 1, assign_pheno_category)
#max_dates_cycl_only$sampling_season = factor(max_dates_cycl_only$sampling_season, ordered = T, levels = c("dormancy","dormancy_greenup", "greenup","greenup_peak", "peak", "greendown_peak", "greendown","greendown_dormancy"))

#mutate(site_cat = assign_pheno_date(dates))


##### #####
# Approach 1: use sine/cosine wave from model coefficients. Captures shared annual trend across sites.

max_cycl_abundance <- seas_amplitude %>%
	mutate(dates = ymd(paste0("2016-", seas_amplitude$max_month, "-15")))# %>%
#mutate(site_cat = assign_pheno_date(dates))

# Assign phenophase to those. Should only take a few seconds
phenophase_cycl_abundance <- lapply(max_cycl_abundance$dates, assign_pheno_date) %>%
	as.matrix() %>% unlist() %>%
	cbind.data.frame(sampling_season=., max_cycl_abundance)
phenophase_cycl_abundance$sampling_season = factor(phenophase_cycl_abundance$sampling_season, ordered = T, levels = c("dormancy","dormancy_greenup",
																																																											"greenup","greenup_peak",
																																																											"peak", "greendown_peak",
																																																											"greendown","dormancy_greendown"
))
phenophase_cycl_abundance <- merge(phenophase_cycl_abundance, model_info_key, all.x=T)

##### #####
# Approach 2: Get phenophase of the month of peak modeled abundances, which vary each site/year

# Calculate mean abundances for every site/year
site_estimates = plot_estimates %>%
	group_by(pretty_group,fcast_type,rank_only, model_id,taxon,time_period, siteID, year, month, dates) %>%
	# Use median plot estimates; take site/date mean of these values
	summarize(mean_modeled_abun = mean(`50%`, na.rm=T),
						sd_modeled_abun = mean(`50%`, na.rm=T),
						median_modeled_abun = median(`50%`, na.rm=T)) %>%
	filter(!is.na(median_modeled_abun))

# Subset model estimates to years 2016-2018, which have data for most sites
max_modeled_abun = site_estimates %>%
	# Use median plot estimates; take site/date mean of these values
	filter(time_period == "20130601_20151101" & year %in% c("2014") |
				 	time_period == "2015-11_2020-01" & year %in% c("2016","2017","2018","2019") |
				 	time_period == "2015-11_2018-01" & year %in% c("2016","2017")) %>%
	# Use median plot estimates; take site/date mean of these values
	ungroup %>%
	group_by(model_id, siteID, year) %>%
	# Only retain highest value per site-year
	filter(mean_modeled_abun == max(mean_modeled_abun, na.rm=T))

# Assign the phenophases
# Approx 60k rows; takes a few minutes to assign phenophases to each site/date
phenophase_modeled_abundance <- max_modeled_abun %>%
	mutate(dates = ymd(paste0(year, "-", month, "-15"))) %>%
	filter(siteID %in% pheno_categories_long$ID) %>%
	mutate(site_cat = assign_pheno_site_date(siteID, dates))
# Note: About 6 sites are missing all MODIS phenophase data; these "NA" values will be excluded
# Across all sites, about 7000 additional dates are missing data

# This factor should be ordered (so that phenophases are sequential)
phenophase_modeled_abundance$sampling_season =
	factor(phenophase_modeled_abundance$site_cat, ordered = T, levels = c("dormancy",#"dormancy_greenup",
																																						 "greenup",#"greenup_peak",
																																						 "peak", #"greendown_peak",
																																						 "greendown"#,"dormancy_greendown"
))

# Get most frequently-assigned phenophase per group
seasonality_mode = phenophase_modeled_abundance %>%
	filter(!is.na(sampling_season)) %>%
	group_by(model_id, sampling_season) %>%
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>%
	select(-n) #%>%
#pivot_wider(names_from="sampling_season", values_from = "freq", values_fill = 0)
seasonality_mode2 = seasonality_mode %>% group_by(model_id) %>% filter(freq == max(freq))
seasonality_mode2 <- merge(seasonality_mode2, model_info_key, all.x=T)

seasonality_mode3 = seasonality_mode %>% group_by(model_id) %>% filter(freq > .5)
seasonality_mode3 <- merge(seasonality_mode3,  model_info_key, all.x=T)

max_modeled_abun_to_plot <- merge(max_modeled_abun,  model_info_key, all.x=T)

seasonality_mode_to_plot <- merge(seasonality_mode,  model_info_key, all.x=T)

# saveRDS(list(seasonality_mode2, max_abun,seasonality_mode3, max_abun_to_plot, seasonality_mode_to_plot), here("data/clean/group_peak_phenophases.rds"))
#saveRDS(list(seasonality_mode2, max_abun,seasonality_mode3, max_abun_to_plot, seasonality_mode_to_plot), here("data/clean/pheno_group_peak_phenophases.rds"))



phenophase_mode_filtered = phenophase_modeled_abundance %>%
filter(!is.na(sampling_season)) %>%
	group_by(model_id, sampling_season) %>%
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>%
	select(-n)
seasonality_mode_filtered <- merge(phenophase_mode_filtered, model_info_key, all.x=T)

seasonality_mode_to_plot = seasonality_mode_filtered



freq_phenophase = seasonality_mode_to_plot %>%
	pivot_wider(values_from=freq, names_from=sampling_season, values_fill = 0)
freq_phenophase$greenup_peak = freq_phenophase$greenup + freq_phenophase$peak
freq_phenophase$greenup_peak_greendown = freq_phenophase$greenup + freq_phenophase$peak + freq_phenophase$greendown
freq_phenophase$max_phase_green = ifelse(freq_phenophase$greenup_peak > .5, 1, 0)
freq_phenophase$max_phase_growing = ifelse(freq_phenophase$greenup_peak_greendown > .5, 1, 0)
freq_phenophase_long = freq_phenophase %>%
	pivot_longer(cols=c(greenup_peak, greendown, dormancy),
							 names_to = "sampling_season", values_to="freq")
freq_phenophase_long$sampling_season = ordered(freq_phenophase_long$sampling_season, levels = c("greenup_peak","greendown","dormancy"))


phenophase_proportion = freq_phenophase %>%
	group_by(model_name, pretty_group, rank_only) %>%
	add_tally(name = "total") %>%
	ungroup() %>%
	group_by(model_name, pretty_group, rank_only, total, max_phase_growing) %>%
	tally(name = "count") %>%
	filter(max_phase_growing == 1 ) %>%
	mutate(proportion_growing_season = count/total) %>% select(-c(max_phase_growing, count))


ggplot(freq_phenophase_long %>% filter(#fcast_type == "Functional" &
	model_name=="cycl_only"),
	aes(fill=pretty_group, y = freq,
			x = sampling_season, group=model_id)) +
	geom_line(aes(color = pretty_group), alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~pretty_group, scales="free") +
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase") +
	geom_label_repel(aes(label=taxon), size=4)

# Takes a few seconds (assigns phenophase to every site/date combo - first vectorizes function)
phenophase_key <- plot_estimates %>%
	filter(!is.na(Mean)) %>%
	select(siteID, year, month) %>% distinct() %>%
	mutate(dates = ymd(paste0(year, "-", month, "-15"))) %>%
	filter(siteID %in% pheno_categories_long$ID)
assign_pheno_site_date_v <- Vectorize(assign_pheno_site_date)
phenophase_key <- phenophase_key %>%
	mutate(site_cat = assign_pheno_site_date_v(siteID, dates)) %>% select(-dates)



saveRDS(list(phenophase_key, phenophase_modeled_abundance,seasonality_mode, seasonality_mode2,seasonality_mode3, site_estimates,phenophase_mode_filtered,phenophase_cycl_abundance,phenophase_proportion), here("data/clean/pheno_group_peak_phenophases.rds"))


read_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
#pheno_info_key = read_in[[1]] %>% merge(model_info_key, all.x=T)
phenophase_key = read_in[[1]]
plot_phenophase_abundance <- merge(plot_estimates, phenophase_key, all.x=T)

phenophase_fg_abundance_fungi <- plot_phenophase_abundance %>%  #full_phenophase_abundance %>%
	filter(pretty_group=="Fungi" & fcast_type == "Functional" & time_period=="2015-11_2018-01") %>% 	filter(!is.na(Mean))

# This factor should be ordered (so that phenophases are sequential)
phenophase_fg_abundance_fungi$sampling_season =
	factor(phenophase_fg_abundance_fungi$site_cat, ordered = T, levels = c("dormancy",#"dormancy_greenup",
																																				"greenup",#"greenup_peak",
																																				"peak", #"greendown_peak",
																																				"greendown"#,"dormancy_greendown"
	))


mean_phenophase_abundance = phenophase_fg_abundance_fungi %>%
	filter(!is.na(sampling_season)) %>%
	group_by(model_id, model_name, sampling_season) %>%
	summarize(mean_pheno_abundance = mean(Mean, na.rm=T),
				 sd_pheno_abundance = mean(SD, na.rm=T))
mean_phenophase_abundance <- merge(mean_phenophase_abundance, model_info_key, all.x=T)


ggplot(mean_phenophase_abundance  %>% filter(
	model_name=="cycl_only"),
			 aes(color=taxon, y = mean_pheno_abundance,
			 		x = sampling_season)) +
	geom_point(aes(color = taxon), size=3, alpha=.5, #show.legend = F,
						 position=position_jitter(height=0, width=.1)) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	#facet_grid(taxon~model_name, scales="free") +
	ylab("Mean abundance in phenophase") +
	xlab("Phenophase") #+


x <- ggplot(phenophase_fg_abundance_fungi  %>% filter(
	model_name=="cycl_only"),
	aes(color=siteID, y = mean_modeled_abun,
			x = sampling_season)) +
	geom_point(aes(color = taxon), size=3, alpha=.5, #show.legend = F,
						 position=position_jitter(height=0)) +
	geom_boxplot(aes(color = taxon), size=1, alpha=.5, show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~taxon, scales="free") +
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase") #+



library(rstatix)


tukey_mean_pheno = phenophase_fg_abundance_fungi %>%
	group_by(fcast_type,model_name,taxon) %>%
	summarize(tukey(x = sampling_season, y = mean_modeled_abun)) %>%
	rename(sampling_season = x)



tukey_median_pheno = phenophase_fg_abundance_fungi %>%
	filter(model_name=="cycl_only") %>%
	group_by(fcast_type,model_name,taxon) %>%
	summarize(tukey(x = sampling_season, y = median_modeled_abun)) %>%
	rename(sampling_season = x)

stat.test2 <- phenophase_fg_abundance_fungi %>%
	group_by(fcast_type,model_name,taxon) %>%
	t_test(median_modeled_abun ~ sampling_season) %>%
	#adjust_pvalue(method = "bonferroni") %>%
	add_significance()
stat.test2 <- stat.test2 %>% add_xy_position(x = "sampling_season")

x + stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T)

x +
	geom_text(data = tukey_median_pheno,
						aes(x = sampling_season, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6)



# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))
site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups



phenophase_abundance_site <- merge(phenophase_fg_abundance_fungi, site_descr, all=T)

tukey_median_pheno_ecoregion = phenophase_abundance_site %>%
	filter(model_name=="cycl_only") %>%
	group_by(fcast_type,model_name,taxon,deciduous) %>%
	summarize(tukey(x = sampling_season, y = median_modeled_abun)) %>%
	rename(sampling_season = x)
ggplot(phenophase_fg_abundance_fungi  %>% filter(
	model_name=="cycl_only"),
	aes(color=siteID, y = mean_modeled_abun,
			x = sampling_season)) +
	geom_point(aes(color = taxon), size=3, alpha=.5, #show.legend = F,
						 position=position_jitter(height=0)) +
	geom_boxplot(aes(color = taxon), size=1, alpha=.5, show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(deciduous~taxon, scales="free") +
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase")  +
	geom_text(data = tukey_median_pheno_ecoregion,
						aes(x = sampling_season, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6)

site_descr$latitude_category = ifelse(site_descr$latitude > 44, "High-latitude",
																			ifelse(site_descr$latitude < 31, "Low-latitude",
																						 "Mid-latitude"))

# This factor should be ordered (so that phenophases are sequential)
site_descr$latitude_category =
	factor(site_descr$latitude_category, ordered = T, levels = c("Low-latitude",  "Mid-latitude","High-latitude"))

for_stats <- phenophase_fg_abundance_fungi  %>%
	filter(model_name=="cycl_only" & !is.na(`50%`) & !is.na(site_cat))
for_stats <- merge(for_stats, site_descr, all.x=T)


tukey_median_pheno_evergreen = for_stats %>%
	group_by(fcast_type,model_name,taxon,latitude_category) %>%
	summarize(tukey(x = sampling_season, y = `50%`)) %>%
	rename(sampling_season = x)
plot_median_abun <- ggplot(for_stats %>% filter(taxon %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
	aes(y = `50%`,color = sampling_season,
			x = sampling_season)) +
	# geom_point(#aes(color = siteID),
	# 					 size=1, alpha=.3, #show.legend = F,
	# 					 position=position_jitter(height=0), show.legend = F) +
	geom_boxplot(alpha=.5, show.legend = F) +
	#geom_line(aes(color = taxon), size=3, alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(taxon~latitude_category, scales="free") +
	ylab("Median abundance across all sites") +
	xlab("Plant phenophase")  +
	geom_text(data = tukey_median_pheno_evergreen %>% filter(taxon %in% c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")),
						aes(x = sampling_season, y = tot-.1, label = Letters_Tukey),
						show.legend = F, color = 1, size =4)
plot_median_abun



ecto <- for_stats %>% filter(taxon=="ectomycorrhizal")
ecto %>%
	group_by(fcast_type,model_name,taxon,latitude_bin) %>%
	summarize(tukey(x = sampling_season, y = `50%`)) %>%
	rename(sampling_season = x)

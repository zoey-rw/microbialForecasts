# read in seasonal values
library(lubridate)
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/assignPhenology.r")

sum.in <- readRDS(here("data", "summary/logit_beta_regression_summaries.rds"))
seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))

sum.in_pheno <- readRDS(here("data", "summary/pheno_summaries.rds"))
seas_in_pheno = readRDS(here("data/summary/pheno_seasonal_amplitude.rds"))

seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% sum.in$keep_models_weak$model_id)
seas_vals_wide <- seas_in[[2]] %>% filter(model_id %in% sum.in$keep_models_weak$model_id)

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>%
	select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)
cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")

plot_estimates = sum.in$plot_est %>% filter(model_name != "all_covariates")

plot_estimates$month = lubridate::month(plot_estimates$dates)
plot_estimates$year = lubridate::year(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cov_est = plot_estimates %>% filter(grepl("env_cov",model_name))


pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
pheno_categories_long = pheno_categories_in[[2]] %>% filter(ID %in% unique(cycl_only_est$siteID))


max_dates = seas_vals_short %>% mutate(dates = max_y_date)
#max_dates$sampling_season =  apply(max_dates, 1, assign_pheno_category)
#max_dates_cycl_only$sampling_season = factor(max_dates_cycl_only$sampling_season, ordered = T, levels = c("dormancy","dormancy_greenup", "greenup","greenup_peak", "peak", "greendown_peak", "greendown","greendown_dormancy"))

max_cycl <- max_dates %>%
	mutate(dates = ymd(paste0("2014-", max_dates$max_month, "-15")))# %>%
	#mutate(site_cat = assign_pheno_date(dates))

# Subset model estimates to years 2016-2018, which have data for most sites
new_max_abun = plot_estimates %>%
	filter(year %in% c("2014","2016","2017")) %>%
	#filter(time_period == "2015-11_2018-01") %>%
	group_by(model_id, siteID, year, month, dates) %>%
	# Use median plot estimates; take site/date mean of these values
	summarize(mean_modeled_abun = mean(`50%`, na.rm=T)) %>%
	ungroup %>%
	group_by(model_id, siteID, year) %>%
	# Only retain highest value per site-year
	filter(mean_modeled_abun == max(mean_modeled_abun, na.rm=T))

# Approx 60k rows; takes a few minute to assign phenophases to each site/date
max_abun <- new_max_abun %>%
	mutate(dates = ymd(paste0(year, "-", month, "-15"))) %>%
	filter(siteID %in% pheno_categories_long$ID) %>%
	mutate(site_cat = assign_pheno_site_date(siteID, dates))
# About 6 sites are missing all MODIS phenophase data; these "NA" values will be excluded
# About 7000 additional dates have some missing data
max_abun$sampling_season = factor(max_abun$site_cat, ordered = T, levels = c("dormancy",#"dormancy_greenup",
																																						 "greenup",#"greenup_peak",
																																						 "peak", #"greendown_peak",
																																						 "greendown"#,"dormancy_greendown"
))

# temp <- max_cycl
# temp$rowname = 1:nrow(temp)
# cat_list2 = temp %>%
# 	split( .$rowname) %>%
# 	map(assign_pheno_date, dates)
# site_cats2 = as.matrix(cat_list2)
# max_abun <- cbind.data.frame(new_max_abun, site_cats2)

site_cats2 = lapply(max_cycl$dates, assign_pheno_date) %>% as.matrix() %>% unlist()
max_cycl <- cbind.data.frame(max_cycl, site_cats2)
max_cycl$sampling_season = factor(max_cycl$site_cats2, ordered = T, levels = c("dormancy","dormancy_greenup",
																																						 "greenup","greenup_peak",
																																						 "peak", "greendown_peak",
																																						 "greendown","dormancy_greendown"
))
max_cycl <- merge(max_cycl, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")], all.x=T)


# Get most frequently-assigned phenophase per group
seasonality_mode = max_abun %>%
	filter(!is.na(sampling_season)) %>%
	group_by(model_id, sampling_season) %>%
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>%
	select(-n) #%>%
#pivot_wider(names_from="sampling_season", values_from = "freq", values_fill = 0)
seasonality_mode2 = seasonality_mode %>% group_by(model_id) %>% filter(freq == max(freq))
seasonality_mode2 <- merge(seasonality_mode2, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")], all.x=T)

seasonality_mode3 = seasonality_mode %>% group_by(model_id) %>% filter(freq > .5)
seasonality_mode3 <- merge(seasonality_mode3, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")], all.x=T)

max_abun_to_plot <- merge(max_abun, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")])

seasonality_mode_to_plot <- merge(seasonality_mode, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")])

# saveRDS(list(seasonality_mode2, max_abun,seasonality_mode3, max_abun_to_plot, seasonality_mode_to_plot), here("data/clean/group_peak_phenophases.rds"))
saveRDS(list(seasonality_mode2, max_abun,seasonality_mode3, max_abun_to_plot, seasonality_mode_to_plot), here("data/clean/pheno_group_peak_phenophases.rds"))



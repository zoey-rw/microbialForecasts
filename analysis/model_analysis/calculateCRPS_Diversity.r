# Calculate and visualize CRPS for diversity forecasts
pacman::p_load(scoringRules,reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

hindcast_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_div.rds")
#hindcast_data <- hindcast_data %>% separate(scenario, sep = "\\.", into = c("scenario","plot","datenum"))


# hindcast_data$scenario < factor(hindcast_data$scenario, levels = c("no_uncertainty","temporal_uncertainty","spatial_uncertainty","full_uncertainty"))
hindcast_data <- hindcast_data %>% filter(dates > "2016-12-31"  & !is.na(hindcast_data$truth)) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2))


hindcast_data$uncert <- factor(hindcast_data$uncert, levels = c("no_uncertainty",
																																		"temporal_uncertainty",
																																		"spatial_uncertainty",
																																		"full_uncertainty"))

scored_hindcasts <- hindcast_data %>% mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd),
																																									 fcast_type = "Diversity")
saveRDS(scored_hindcasts, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_div.rds")



# COMPARE CRPS BY ERRORS-IN-VARIABLES UNCERTAINTY
# All values

ggplot(scored_hindcasts, aes(x = scenario, y = crps)) + 
	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~group) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = scenario, y = crps), width=.2, height = 0, alpha = .3) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

not_na <- scored_hindcasts[!is.na(scored_hindcasts$truth),]

# Compare means
scored_hindcasts_mean <- scored_hindcasts %>% group_by(group, uncert, new_site) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))
ggplot(scored_hindcasts_mean, aes(x = uncert, y = crps_mean, color = group)) + 
#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~new_site) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = uncert, y = crps_mean), width=.2, height = 0, alpha = .3, size=4) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)


# Do T-tests
newsites <- scores_out[scores_out$newsites=="New sites",]
newsites_bac <- newsites[newsites$group=="16S",]
newsites_fun <- newsites[newsites$group=="ITS",]
oldsites <- scores_out[scores_out$newsites=="Observed sites",]
oldsites_bac <- oldsites[oldsites$group=="16S",]
oldsites_fun <- oldsites[oldsites$group=="ITS",]
t.test(oldsites_bac$crps, oldsites_fun$crps)
t.test(newsites_bac$crps, newsites_fun$crps)

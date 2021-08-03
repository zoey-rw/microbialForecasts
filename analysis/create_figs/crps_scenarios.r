library(scoringRules)
library(tidyverse)


full_uncert <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data.rds") %>% mutate(truth = as.numeric(truth)) %>% group_by(group)  %>% mutate(scenario="full_uncertainty")

spatial_uncert <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data_spatial.rds") %>% mutate(truth = as.numeric(truth)) %>% group_by(group) %>%  mutate(scenario="spatial_uncertainty")

temporal_uncert <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data_temporal.rds") %>% mutate(truth = as.numeric(truth)) %>% group_by(group) %>%  mutate(scenario="temporal_uncertainty") 

no_uncert <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data_no_uncert.rds")  %>% mutate(truth = as.numeric(truth)) %>% group_by(group) %>%  mutate(scenario="no_uncertainty") 

observed_points <- do.call(plyr::rbind.fill, list(full_uncert, spatial_uncert, temporal_uncert, no_uncert))

observed_points$scenario <- factor(observed_points$scenario, levels = c("no_uncertainty","spatial_uncertainty","temporal_uncertainty","full_uncertainty"), ordered=T)
observed_points_subset <- observed_points[observed_points$plotID %in% c("HARV_001", "DSNY_001", "ORNL_001"),]
ggplot(observed_points_subset[observed_points_subset$pretty_group=="Bacteria",]) + 
	facet_wrap(plotID ~ scenario) + 
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = pretty_group, alpha = fcast_period
),
show.legend = F) +	scale_alpha_discrete(range = c(.7, .3)) +
geom_point(aes(x = dates, y = as.numeric(truth))) 






scored_scenarios <- observed_points %>% group_by(group, scenario) %>% filter(fcast_period=="Hindcast") %>% 
	mutate(crps = crps_norm(truth, mean_EDP, sd_EDP)) 
	
p <- ggplot(scored_scenarios, aes(x = pretty_group, y = crps)) + 
	geom_violin(aes(fill = pretty_group), draw_quantiles = c(0.5), show.legend=F) + 
	facet_grid(~scenario) +   coord_trans(y = "log10") +
	geom_jitter(aes(x = pretty_group, y = crps), width=.2, height = 0, alpha = .3) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)
p


scored_scenarios
library(ComplexUpset)
# Create parameters to pass	
params = data.frame(scenario = c("no_uncertainty", "spatial_uncertainty",
																 "temporal_uncertainty", "full_uncertainty"),
										temporalDriverUncertainty = c(F, F, T, T),
										spatialDriverUncertainty = c(F, T, F, T))
params$`Soil moisture` <- ifelse(params$temporalDriverUncertainty, T, F)
params$`Soil temperature` <- ifelse(params$temporalDriverUncertainty, T, F)
params$`Soil pH` <- ifelse(params$spatialDriverUncertainty, T, F)
params$`Soil pC` <- ifelse(params$spatialDriverUncertainty, T, F)
params$`Plant richness` <- F
params$`% grass cover` <- F

scored_scenarios <- merge(scored_scenarios, params, all=T)

uncertainties <- c("Soil moisture", "Soil temperature", "Soil pH", "Soil pC"
									 #, 
									 #"Plant richness", "% grass cover"
									 )
u1 <- ComplexUpset::upset(scored_scenarios, uncertainties,
													base_annotations=list(
														'CI'=(
															ggplot() 
															+ facet_grid(rows=vars(pretty_group))
															+ geom_violin(aes(color = pretty_group, y = crps), draw_quantiles = c(0.5), show.legend=F) 
															+ geom_point(aes(color = pretty_group, fill=pretty_group, y = crps),fill="black",alpha = .4, position = position_jitterdodge(
																jitter.width = .3,
																jitter.height = 0,
																dodge.width = 1
															), show.legend = F) 
															+ ylab("Continuous ranked probability score (CRPS)") 
															+ xlab(NULL) 
																	+ coord_trans(y = "log10") 
															+ theme(text = element_text(size=20))
														)
													),
													set_sizes=FALSE,
													
													themes=upset_default_themes(text=element_text(size=20)),
													sort_intersections_by=c('ratio'
																									#'cardinality','degree'
																									),
) + xlab('Errors in variables')
u1


upset_test(scored_scenarios, uncertainties)

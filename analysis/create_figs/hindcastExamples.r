# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Read in hindcast data
hindcast_in <- readRDS("./data/summary/all_hindcasts.rds")

# Read in CRPS scores
crps_in <- readRDS("./data/summary/CRPS_hindcasts.rds")
crps_by_group <- crps_in$scored_hindcasts_taxon
best_per_group <- crps_by_group %>% 
	filter(model_name == "all_covariates" & newsite == "Observed site") %>% 
	group_by(fcast_type) %>% 
	filter(crps_mean == min(crps_mean, na.rm=TRUE))
worst_per_group <- crps_by_group %>% 
	filter(model_name == "all_covariates" & newsite == "Observed site") %>% 
	group_by(fcast_type) %>% 
	filter(crps_mean == max(crps_mean, na.rm=TRUE))


bad_hindcasts <- hindcast_data %>% filter(model_name == "all_covariates" &
																						plotID=="HARV_004" & taxon %in% worst_per_group$taxon)
ggplot(bad_hindcasts) + 
	facet_grid(#rows=vars(taxon), 
		rows = vars(fcast_type), drop=T, scales="free") +
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


good_hindcasts <- hindcast_data %>% filter(model_name == "all_covariates" &
																						plotID=="HARV_004" & taxon %in% best_per_group$taxon)
ggplot(good_hindcasts) + 
	facet_grid(rows=vars(taxon), drop=T, scales="free") +
		#rows = vars(fcast_type), drop=T, scales="free") +
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

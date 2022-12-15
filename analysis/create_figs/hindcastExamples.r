# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Read in hindcast data
hindcast_in <- readRDS("./data/summary/all_hindcasts.rds")
#hindcast_in$fcast_period <- ifelse(hindcast_in$dates > "2018-01-01", "hindcast", "calibration")

# Read in CRPS scores
crps_in <- readRDS("./data/summary/CRPS_hindcasts.rds")
crps_in <- readRDS("./data/summary/scoring_metrics_cv.rds")
crps_by_group <- crps_in$scoring_metrics_site
best_per_group <- crps_by_group %>%
	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
	group_by(fcast_type) %>%
	filter(CRPS_truncated == min(CRPS_truncated, na.rm=TRUE))
worst_per_group <- crps_by_group %>%
	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
	group_by(fcast_type) %>%
	filter(CRPS_truncated == max(CRPS_truncated, na.rm=TRUE))

bad_hindcasts <- hindcast_in %>% filter(model_name == "all_covariates" & taxon %in% worst_per_group$taxon)
good_hindcasts <- hindcast_in %>% filter(model_name == "all_covariates" & taxon %in% best_per_group$taxon)

# Get best plots for examples (most calibration AND validation points)
not_na_hindcast <- hindcast_in %>% filter(fcast_period == "hindcast" & !is.na(truth))
not_na_calibration <- hindcast_in %>% filter(fcast_period == "calibration" & !is.na(truth))
plot_hindcast_counts <- sort(table(not_na_hindcast$plotID))
plot_calibration_counts <- sort(table(not_na_calibration$plotID))
top_hindcast_plots <- names(tail(plot_hindcast_counts, 50))
top_calibration_plots <- names(tail(plot_calibration_counts, 30))
top_plots <- intersect(top_hindcast_plots, top_calibration_plots)
top_plots
# Filter to those plots
good_hindcasts_top_plots <- good_hindcasts[good_hindcasts$plotID %in% top_plots,]
bad_hindcasts_top_plots <- bad_hindcasts[bad_hindcasts$plotID %in% top_plots,]

ggplot(hindcast_in %>% filter(model_name == "all_covariates" & siteID=="SJER" & taxon=="nitrification")) +
	facet_grid(rows=vars(plotID), drop=T, scales="free") +
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
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +ggtitle("'Good' forecast",  "lowest CRPS score, for site x taxon combination")

ggplot(hindcast_in %>% filter(model_name == "all_covariates" & siteID=="OSBS" & taxon=="agaricomycetes")) +
	facet_grid(rows=vars(plotID), drop=T, scales="free") +
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
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +ggtitle("'Bad' forecast",  "highest CRPS score, for site x taxon combination")

# one plot for testing
ggplot(good_hindcasts_top_plots %>% filter(plotID=="CPER_004")) +
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

good_hindcasts_top_plots[good_hindcasts_top_plots$taxon_name == "armatimonadota",]

pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_good_example.pdf", height = 6)
for(i in 1:length(top_plots)){
	print(ggplot(good_hindcasts_top_plots) +
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
						taxon~plotID,
						drop=T, scales="free",
						ncol = 1, nrow = 3, page = i))
}
dev.off()


pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_bad_example.pdf", height = 6)
for(i in 1:length(top_plots)){
	print(ggplot(bad_hindcasts_top_plots) +
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
						taxon~plotID,
						drop=T, scales="free",
						ncol = 1, nrow = 3, page = i))
}
dev.off()


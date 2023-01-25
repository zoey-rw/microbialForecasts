# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Read in hindcast data
hindcast_in <- readRDS("./data/summary/all_hindcasts.rds")
# Remove first timepoint from calibration
hindcast = hindcast_in %>% filter(!(fcast_period=="calibration" & is_any_start_date)) %>%
	filter(fcast_type != "Diversity")


# Add duplicated data row so that plotting ribbons are continuous.
max_cal_date <- hindcast %>%
	group_by(taxon, plotID) %>%
	filter(fcast_period=="calibration") %>%
	filter(timepoint == max(timepoint)) %>%
	mutate(fcast_period="hindcast")
hindcast <- rbind.data.frame(hindcast, max_cal_date)



#hindcast_in$fcast_period <- ifelse(hindcast_in$dates > "2018-01-01", "hindcast", "calibration")

converged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
converged <- converged %>% filter(median_gbr < 1.5 & model_name == "all_covariates")


# Read in CRPS scores
scores_list <- readRDS("./data/summary/scoring_metrics_cv.rds")
crps_by_group <- scores_list$scoring_metrics %>% filter(fcast_type != "Diversity")
best_per_group <- crps_by_group %>%
	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
	group_by(fcast_type) %>%
	filter(RSQ.1 > 0) %>%
	filter(CRPS_truncated == min(CRPS_truncated, na.rm=TRUE))
worst_per_group <- crps_by_group %>%
	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
	group_by(fcast_type) %>%
	filter(CRPS_truncated == max(CRPS_truncated, na.rm=TRUE))


select_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% c("cellulolytic","nitrification","chitin_complex","mortierellales","trichoderma"))


bad_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% worst_per_group$taxon)
good_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% best_per_group$taxon)

# Get best plots for examples (most calibration AND validation points)
not_na_hindcast <- hindcast %>% filter(fcast_period == "hindcast" & !is.na(truth))
not_na_calibration <- hindcast %>% filter(fcast_period == "calibration" & !is.na(truth))
plot_hindcast_counts <- sort(table(not_na_hindcast$plotID))
plot_calibration_counts <- sort(table(not_na_calibration$plotID))
top_hindcast_plots <- names(tail(plot_hindcast_counts, 50))
top_calibration_plots <- names(tail(plot_calibration_counts, 30))
top_plots <- intersect(top_hindcast_plots, top_calibration_plots)
top_plots
# Filter to those plots
good_hindcasts_top_plots <- good_hindcasts[good_hindcasts$plotID %in% top_plots,]
bad_hindcasts_top_plots <- bad_hindcasts[bad_hindcasts$plotID %in% top_plots,]

select_plots <- c("HARV_033","OSBS_026","WOOD_044","KONZ_001")
select_plots <- c("HARV_033","WOOD_044")
select_plots <- c("HARV_033","KONZ_001")
select_plots <- c("HARV_033","HARV_004","HARV_001","HARV_034")

select_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% c("chitin_complex","mortierellales"))
select_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% c("actinobacteriota","acidobacteriota"))

select_hindcasts_top_plots <- select_hindcasts[select_hindcasts$plotID %in% top_plots,]
select_hindcasts_select_plots <- select_hindcasts[select_hindcasts$plotID %in% select_plots,]

select_hindcasts_select_plots

# New facet label names for supp variable
tax.labs <- c("Bacteria (chitin-enriched)", "Fungi (Order: Mortierellales)")
names(tax.labs) <- c("chitin_complex", "mortierellales")


ggplot(select_hindcasts_select_plots, aes(fill=plotID, x = dates, y =med)) +
	facet_grid(#rows=vars(plotID),
						 rows=vars(species),
						 drop=T, scales="free",
						 labeller = labeller(species = tax.labs)
	) +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo, ymax = hi), alpha=0.3) +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo_25, ymax = hi_75), alpha=.6) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo_25, ymax = hi_75), alpha=0.3) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo, ymax = hi), alpha=0.1) +
	geom_line(data = ~filter(.x, fcast_period=="calibration"),
						 alpha=0.8) +
	geom_line(data = ~filter(.x, fcast_period=="hindcast"),
						aes(x = dates, y = med), alpha=0.3) +
	geom_point(aes(x = dates, y = as.numeric(truth), color=plotID), position = position_jitter()) +
	xlab(NULL) + labs(fill='') +
	scale_fill_brewer(palette = "Set2") +
	scale_color_brewer(palette = "Set2") +
	theme(panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	ggtitle("Example plot-level hindcasts") + theme_minimal(base_size = 18) +
	scale_y_log10() + theme(legend.position = "none")

# COlor by species instead
ggplot(select_hindcasts_select_plots, aes(fill=species, x = dates, y =med, group=plotID)) +
	facet_grid(#rows=vars(plotID),
		rows=vars(species),
		drop=T, scales="free",
		 labeller = labeller(species = tax.labs)
	) +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo, ymax = hi), alpha=0.3) +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo_25, ymax = hi_75), alpha=.6) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo_25, ymax = hi_75), alpha=0.3) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo, ymax = hi), alpha=0.1) +
	geom_line(data = ~filter(.x, fcast_period=="calibration"),
						alpha=0.8) +
	geom_line(data = ~filter(.x, fcast_period=="hindcast"),
						aes(x = dates, y = med), alpha=0.3) +
	geom_point(aes(x = dates, y = as.numeric(truth), color=species), position = position_jitter()) +
	xlab(NULL) + labs(fill='') +
	scale_fill_brewer(palette = "Set2") +
	scale_color_brewer(palette = "Set2") +
	theme(panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	ggtitle("Example plot-level hindcasts") + theme_minimal(base_size = 18) + scale_y_log10() +
	theme(legend.position = "none")

ggplot(hindcast_in %>% filter(model_name == "all_covariates" & siteID=="SJER" & taxon=="nitrification")) +
	facet_grid(rows=vars(plotID), drop=T, scales="free") +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = med), show.legend = F) +
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


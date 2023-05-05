

# bad_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% worst_per_group$taxon)
# good_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% best_per_group$taxon)
#
# # Get best plots for examples (most calibration AND validation points)
# not_na_hindcast <- hindcast %>% filter(fcast_period == "hindcast" & !is.na(truth))
# not_na_calibration <- hindcast %>% filter(fcast_period == "calibration" & !is.na(truth))
# plot_hindcast_counts <- sort(table(not_na_hindcast$plotID))
# plot_calibration_counts <- sort(table(not_na_calibration$plotID))
# top_hindcast_plots <- names(tail(plot_hindcast_counts, 50))
# top_calibration_plots <- names(tail(plot_calibration_counts, 30))
# top_plots <- intersect(top_hindcast_plots, top_calibration_plots)
# top_plots
# # Filter to those plots
# good_hindcasts_top_plots <- good_hindcasts[good_hindcasts$plotID %in% top_plots,]
# bad_hindcasts_top_plots <- bad_hindcasts[bad_hindcasts$plotID %in% top_plots,]
#
#
#
#
# # Read in CRPS scores
# scores_list <- readRDS("./data/summary/scoring_metrics_cv.rds")
# crps_by_group <- scores_list$scoring_metrics %>% filter(fcast_type != "Diversity")
# best_per_group <- crps_by_group %>%
# 	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
# 	group_by(fcast_type) %>%
# 	filter(RSQ.1 > 0) %>%
# 	filter(CRPS_truncated == min(CRPS_truncated, na.rm=TRUE))
# worst_per_group <- crps_by_group %>%
# 	filter(model_name == "all_covariates" & grepl("observed", site_prediction)) %>%
# 	group_by(fcast_type) %>%
# 	filter(CRPS_truncated == max(CRPS_truncated, na.rm=TRUE))

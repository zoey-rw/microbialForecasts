# Calculate and visualize scoring metrics and coefficient of variation for taxonomic forecasts
# Produces "scoring_metrics_cv.rds"

source("../../source.R")

# Source the stats functions directly since they're not exported in the package
source(here("microbialForecast/R/statsFunctions.r"))

# Source the robust version of add_scoring_metrics
source("robust_add_scoring_metrics.R")

# Read in hindcast data created by 01_tidyHindcasts.r
# Check if Parquet file exists, otherwise use RDS
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(parquet_file)) {
  cat("Using Parquet file for memory efficiency...\n")
  hindcast_data <- arrow::read_parquet(parquet_file)
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  hindcast_data <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}

# Create some subsets
hindcast_only = hindcast_data %>% filter(fcast_period=="hindcast") %>%
	filter(!is.na(truth) & !is.na(mean))
calibration_only = hindcast_data %>% filter(fcast_period=="calibration") %>%
	filter(!is.na(truth) & !is.na(mean))

# For the calibration, the first observed date per plot is always wonky due to model structure, so we leave it out of our scoring metrics
calibration_only_not_first = calibration_only %>% filter(date_num > plot_start_date)

# Add scoring metrics to various subsets/groupings using data.table for better performance
library(data.table)

# Convert to data.table for faster processing
hindcast_dt <- as.data.table(hindcast_only %>%
	filter(!is.na(site_prediction)))

# Calculate scoring metrics for each group with error handling
scoring_metrics_list <- hindcast_dt[, {
			# Check if we have data for this group
		if(.N == 0) {
			# Return NA values for empty groups
			list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
						RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
						residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
		} else {
			# Calculate scoring metrics for this group
			tryCatch({
				result <- robust_add_scoring_metrics(
					observed = truth,
					mean_predicted = mean,
					median_predicted = med,
					sd_predicted = sd
				)
				as.list(result)
			}, error = function(e) {
				cat("Warning: Error in group", paste(.BY, collapse=", "), ":", e$message, "\n")
				list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
						RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
						residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
			})
		}
}, by = .(model_id, fcast_type, pretty_group, model_name, pretty_name, rank_name, taxon, site_prediction)]

# Convert back to data.frame
scoring_metrics <- as.data.frame(scoring_metrics_list)

# Now add the CRPS summary
scoring_metrics <- scoring_metrics %>%
	left_join(
		hindcast_dt %>%
			group_by(model_id, fcast_type, pretty_group, model_name, pretty_name, rank_name, taxon, site_prediction) %>%
			summarize(mean_crps_sample = mean(crps, na.rm = TRUE), .groups = "drop"),
		by = c("model_id", "fcast_type", "pretty_group", "model_name", "pretty_name", "rank_name", "taxon", "site_prediction")
	)

#%>% select(-c("RSQ.1", "RMSE.norm"))

calibration_metrics <- calibration_only_not_first %>%
	#filter(!is.na(site_prediction)) %>%
	as.data.table() %>%
	.[, {
		# Check if we have data for this group
		if(.N == 0) {
			# Return NA values for empty groups
			list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
						RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
						residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
		} else {
			# Calculate scoring metrics for this group
			tryCatch({
				result <- robust_add_scoring_metrics(
					observed = truth,
					median_predicted = med,
					mean_predicted = mean,
					sd_predicted = sd
				)
				as.list(result)
			}, error = function(e) {
				cat("Warning: Error in calibration group", paste(.BY, collapse=", "), ":", e$message, "\n")
				list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
							RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
							residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
			})
		}
	}, by = .(model_id, fcast_type, pretty_group, model_name, rank_name, pretty_name, taxon)] %>%
	as.data.frame()

calibration_metrics_site <- calibration_only_not_first %>%
	#filter(!is.na(site_prediction)) %>%
	as.data.table() %>%
	.[, {
		# Check if we have data for this group
		if(.N == 0) {
			# Return NA values for empty groups
			list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
						RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
						residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
		} else {
			# Calculate scoring metrics for this group
			tryCatch({
				result <- robust_add_scoring_metrics(
					observed = truth,
					median_predicted = med,
					mean_predicted = mean,
					sd_predicted = sd
				)
				as.list(result)
			}, error = function(e) {
				cat("Warning: Error in calibration site group", paste(.BY, collapse=", "), ":", e$message, "\n")
				list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
							RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
							residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
			})
		}
	}, by = .(model_id, fcast_type, pretty_group, model_name, rank_name, pretty_name, taxon, siteID)] %>%
	as.data.frame()

# This one took forever so I'm running it with data.table for better performance
scoring_metrics_site <- hindcast_only %>% #filter(newsite=="Observed site") %>%
	filter(!is.na(site_prediction)) %>%
	as.data.table() %>%
	.[, mean_crps_sample := mean(crps, na.rm=T), by = .(model_id, siteID, site_prediction)] %>%
	.[, {
		# Check if we have data for this group
		if(.N == 0 || .N == 1) {
			# Return NA values for empty or single-row groups (can't calculate statistics)
			list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
						RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
						residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
		} else {
			# Calculate scoring metrics for this group
			tryCatch({
				result <- robust_add_scoring_metrics(
					observed = truth,
					mean_predicted = mean,
					median_predicted = med,
					sd_predicted = sd
				)
				as.list(result)
			}, error = function(e) {
				cat("Warning: Error in scoring site group", paste(.BY, collapse=", "), ":", e$message, "\n")
				list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_, 
							RMSE.iqr = NA_real_, RMSE.norm = NA_real_, CRPS_truncated = NA_real_,
							residual_variance = NA_real_, predictive_variance = NA_real_, total_PL = NA_real_)
			})
		}
	}, by = .(model_id, siteID, site_prediction, mean_crps_sample)] %>%
	as.data.frame()


# Pivot longer for easier plotting
scoring_metrics_long <- scoring_metrics %>% pivot_metrics()
# Remove infinite values (which come from having only one validation point)
scoring_metrics_site_long <- scoring_metrics_site  %>% pivot_metrics() %>% filter(!is.infinite(score) & !is.nan(score))
calibration_metrics_long <- calibration_metrics  %>% pivot_metrics()
calibration_metrics_site_long <- calibration_metrics_site  %>% pivot_metrics()  %>% filter(!is.infinite(score) & !is.nan(score))

# check correlations between metrics

# library(corrplot)
# for_cor <- scoring_metrics_site_long %>% ungroup() %>% pivot_wider(names_from = metric, values_from = score) %>% select(RMSE, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1, RMSE.norm)
# cor_scores <- cor(for_cor, use = "pairwise.complete.obs")
# corrplot::corrplot(cor_scores)

# Add coefficient of variation

# Wide format taxonomic
truth_vals <- calibration_only %>%
	filter(model_name == "env_cycl")

# Get variation per plot, then average per site
cv_tax_per_plot <- truth_vals %>%
	group_by(pretty_group, rank_name,taxon, siteID, plotID) %>%
	dplyr::summarize(per_plot_cv = calc_cv(truth))
cv_tax_per_plot_site <- cv_tax_per_plot %>% ungroup %>%
	group_by(pretty_group, rank_name,taxon, siteID) %>%
	dplyr::summarize(mean_per_plot_site_cv = mean(per_plot_cv, na.rm=T))
# And now per taxon
cv_tax_per_plot_taxon <- cv_tax_per_plot %>%
	ungroup %>%
	group_by(pretty_group, rank_name,taxon) %>%
	summarize(mean_per_plot_cv = mean(per_plot_cv, na.rm=T))

# Get variation per site, then average per taxon
cv_tax_per_site <- truth_vals %>%
	group_by(pretty_group, rank_name,taxon, siteID) %>%
	summarize(per_site_cv = calc_cv(truth))
cv_tax_per_site_taxon <- cv_tax_per_site %>% ungroup %>%
	group_by(pretty_group, rank_name, taxon) %>%
	summarize(mean_per_site_cv = mean(per_site_cv, na.rm=T))

# Get variation per site, then average per taxon
cv_tax_per_site <- truth_vals %>%
	group_by(pretty_group, rank_name,taxon, siteID) %>%
	summarize(per_site_cv = calc_cv(truth))
cv_tax_per_site_taxon <- cv_tax_per_site %>% ungroup %>%
	group_by(pretty_group, rank_name, taxon) %>%
	summarize(mean_per_site_cv = mean(per_site_cv, na.rm=T))

cv_tax_overall <- truth_vals %>%
	group_by(pretty_group, pretty_name,rank_name, taxon) %>%
	summarize(overall_cv = calc_cv(truth))

cv_tax <- cv_tax_overall %>%
	merge(cv_tax_per_site_taxon) %>%
	merge(cv_tax_per_site)

scoring_metrics_cv <- merge(scoring_metrics_long, cv_tax, all=T)
scoring_metrics_cv_site <- merge(scoring_metrics_site_long, cv_tax_per_site)



### Predictability

scoring_metrics_cv_long <- scoring_metrics_cv %>% pivot_longer(cols = c(overall_cv, per_site_cv, mean_per_site_cv),
																							 names_to = "cv_type", values_to = "cv")


# CV scaled by rank
cv_metric_scaled <- scoring_metrics_cv_long %>%
	group_by(pretty_group, pretty_name, rank_name, cv_type, metric) %>%
	mutate(CV_scale = scale(cv)[,1],
				 metric_scale = scale(score)[,1])


# Use later to filter: remove taxa that did not converge
unconverged <- readRDS(here("data/summary/unconverged_taxa_list.rds"))
converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))



# Skill score from Scavia 2021
skill_score_taxon <- scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(metric=="CRPS_truncated") %>%
	pivot_wider(id_cols = c("model_id","fcast_type","pretty_group","model_name","pretty_name","rank_name","taxon"),
							values_from = "mean_crps_sample", names_from = "site_prediction") %>%
	mutate(skill_score = (1 - (`New time x site (modeled effect)`/`New time (observed site)`)),
				 skill_score_random = (1 - (`New time x site (random effect)`/`New time (observed site)`)),
	)


# Skill score from Scavia 2021, except with RMSE
skill_score_taxon_RMSE <- scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(metric=="RMSE.norm") %>%
	pivot_wider(id_cols = c("model_id","fcast_type","pretty_group","model_name","pretty_name","rank_name","taxon"),
							values_from = "score", names_from = "site_prediction") %>%
	mutate(skill_score = (1 - (`New time x site (modeled effect)`/`New time (observed site)`)),
				 skill_score_random = (1 - (`New time x site (random effect)`/`New time (observed site)`)),
	)

skill_score_rank <- skill_score_taxon %>%
	group_by(model_id, pretty_group, pretty_name, rank_name) %>%
	summarize(mean_skill_score = mean(skill_score, na.rm=T))




to_save <- list(cv_metric_scaled = cv_metric_scaled,
								scoring_metrics_cv = scoring_metrics_cv,
								scoring_metrics_cv_site = scoring_metrics_cv_site,
								scoring_metrics_cv_long = scoring_metrics_cv_long,
								calibration_truth_vals = truth_vals,

								scoring_metrics = scoring_metrics,
								scoring_metrics_long = scoring_metrics_long,
								scoring_metrics_site = scoring_metrics_site,
								scoring_metrics_site_long = scoring_metrics_site_long,

								calibration_metrics=calibration_metrics,
								calibration_metrics_long = calibration_metrics_long,
								calibration_metrics_site = calibration_metrics_site,
								calibration_metrics_site_long = calibration_metrics_site_long,

								skill_score_taxon = skill_score_taxon,
								skill_score_taxon_RMSE = skill_score_taxon_RMSE,
								skill_score_rank = skill_score_rank,
								unconverged_list = unconverged,
								converged_list = converged,
								converged_strict_list = converged_strict)
saveRDS(to_save, here("data", paste0("summary/scoring_metrics_plsr2.rds")))

cat("✅ Scoring metrics calculated and saved successfully!\n")
cat("File saved to:", here("data", "summary/scoring_metrics_plsr2.rds"), "\n")

# Comment out plotting code for now to focus on core functionality
# #### Visualizing
# 
# #### View predictability scores by rank
# ggplot(scoring_metrics_long %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
# 			 aes(x = pretty_name, y = score,
# 			 		color = pretty_name)) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
# 	facet_grid(metric ~ pretty_group, drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Metric scores") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	ggtitle("Predictability decreases at broader ranks")  +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))
# 
# 
# #### View CV across tax ranks
# ggplot(cv_long %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
# 			 aes(x = pretty_name, y = cv,
# 			 		color = pretty_name), alpha = .5) +
# 	geom_point(size = 3, position = position_jitterdodge(#jitter.width = 1
# 	), alpha=.5, show.legend = F) +
# 	facet_grid(pretty_group~cv_type, scales="free", drop =T) + theme_bw(base_size=18) +
# 	ggtitle("Variability decreases with broader taxonomic ranks") +
# 	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + xlab(NULL) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))
# 
# 
# #### View CV across tax ranks, for each plot
# # Just per-plot CV (capturing temporal variability)
# ggplot(cv_long %>% filter(model_name == "all_covariates" & cv_type == "mean_per_plot_cv"),
# 			 aes(x = pretty_name, y = cv,
# 			 		color = pretty_name), alpha = .5) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(#jitter.width = 1
# 	), alpha=.5, show.legend = F) +
# 	facet_grid(rows=vars(pretty_group), scales="free", drop =T) + theme_bw(base_size=18) +
# 	ggtitle("Variability decreases with broader taxonomic ranks") +
# 	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + xlab(NULL)
# 
# 
# 
# # Just RMSE
# ggplot(scoring_metrics_cv %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
# 			 aes(x = mean_per_plot_cv, y = RMSE,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
# 	facet_grid(~pretty_group, drop = T, scales="free") + theme_bw()
# 
# 
# 
# 
# # All metrics
# ggplot(cv_metric_scaled %>%
# 			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & cv_type == "mean_per_plot_cv"),
# 			 aes(x = CV_scale, y = metric_scale,
# 			 		shape = pretty_group,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
# 						 alpha=.5, show.legend = T) +
# 	facet_wrap(~metric,drop = T, scales="free") +
# 	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
# 																		shape=guide_legend(title="Domain")) +
# 	xlab("Coefficient of variation (scaled by rank)") +
# 	ylab("Predictability (scaled by rank)") +
# 	ggtitle("Predictability increases with variability, but only for some metrics")
# 
# # NOt scaled
# ggplot(cv_metric_scaled %>%
# 			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & cv_type == "mean_per_plot_cv"),
# 			 aes(x = cv, y = score,
# 			 		shape = pretty_group,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
# 						 alpha=.5, show.legend = T) +
# 	facet_wrap(~metric,drop = T, scales="free") +
# 	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
# 																		shape=guide_legend(title="Domain")) +
# 	xlab("Coefficient of variation (NOT mean-scaled by rank)") +
# 	ylab("Predictability (NOT mean-scaled by rank)") +
# 	ggtitle("Predictability increases with variability, but only for some metrics")
# 
# # Just CRPS
# ggplot(cv_metric_scaled %>%
# 			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
# 			 				 	cv_type == "mean_per_plot_cv" & metric == "CRPS"),
# 			 aes(x = CV_scale, y = metric_scale,
# 			 		shape = pretty_group,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
# 						 alpha=.5, show.legend = T) +
# 	facet_grid(~metric, drop = T, scales="free") +
# 	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
# 																		shape=guide_legend(title="Domain")) +
# 	#	stat_regline_equation(aes(y = CV_scale, x = metric_scale), inherit.aes = F, label.x = 1, label.y = 2) +
# 	stat_smooth(aes(x = CV_scale, y = metric_scale), inherit.aes = F, method = "loess", span = .9)
# 
# # F vs B grid, CV types (split by rank)
# ggplot(cv_metric_scaled %>%
# 			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & metric=="CRPS"),
# 			 aes(x = pretty_name, y = cv,
# 			 		shape = pretty_group,
# 			 		color = cv_type)) +
# 	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
# 						 alpha=.5, show.legend = T) +
# 	facet_grid(~pretty_group, drop = T, scales="free") +
# 	theme_bw(base_size = 18) + guides(color=guide_legend(title="Coeff. of \nVariation type"),
# 																		shape=guide_legend(NULL)) + xlab(NULL)  +
# 	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))
# 
# 
# 
# # View skill scores
# ggplot(skill_score_taxon %>% filter(model_name == "all_covariates"),
# 			 aes(x = pretty_name, y = skill_score,
# 			 		color = pretty_group)) +
# 	#geom_violin(draw_quantiles = c(.5))+
# 	geom_point(aes(x = pretty_name, y = skill_score), position = position_jitterdodge(jitter.height = 0, jitter.width = .1), alpha = .5, size=4) +
# 	ylab("Skill score (% decrease in CRPS)") + xlab(NULL) +
# 	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") +
# 	theme(text = element_text(size = 20),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=22))  +
# 	guides(color=guide_legend(title="Domain"),
# 				 shape=guide_legend(title="Forecast type"))
# 
# 
# # Predictability within-site
# ggplot(scoring_metrics_site %>%
# 			 	filter(model_name == "all_covariates"),
# 			 aes(x = pretty_group, y = CRPS,
# 			 		color = pretty_name)) +
# 	geom_jitter(size = 4, alpha=.5, show.legend = T) +
# 	#facet_wrap(~beta, drop = T, scales="free") +
# 	theme_bw(base_size = 18) +
# 	xlab("Rank") +
# 	ylab("CRPS") + ggtitle("Mean predictability at each observed site") +
# 	guides(color=guide_legend(title="Rank")) +scale_y_sqrt()
# 
# ggplot(scoring_metrics_site %>% filter(model_name == "all_covariates" &
# 																			 	#newsite == "Observed site" &
# 																			 	pretty_name != "Diversity"),
# 			 aes(x = pretty_name, y = CRPS,
# 			 		color = pretty_group)) +
# 	geom_violin(draw_quantiles = c(.5))+
# 	geom_point(aes(x = pretty_name, y = CRPS), position=position_jitterdodge(jitter.height = 0),
# 						 alpha = .2, size=4) +
# 	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) +
# 	theme_minimal(base_size=20) + ggtitle("Mean predictability at each observed site") +
# 	theme(text = element_text(size = 20),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=22))  +
# 	guides(color=guide_legend(title="Domain"),
# 				 shape=guide_legend(title="Forecast type")) +scale_y_sqrt()
# 
# 
# # F vs B "other" abundances at each rank
# ggplot(hindcast_data %>%
# 			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & species=="other"),
# 			 aes(x = pretty_name, y = truth,
# 			 		shape = pretty_group)) +
# 	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
# 						 alpha=.5, show.legend = T) +
# 	facet_grid(~pretty_group, drop = T, scales="free") +
# 	theme_bw(base_size = 18) + guides(color=guide_legend(title="Coeff. of \nVariation type"),
# 																		shape=guide_legend(NULL)) + xlab(NULL)  +
# 	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))
# 
# 
# 
# 
# eval_df <- scoring_metrics_site_long %>% filter(model_name == "all_covariates" &
# 																					 	pretty_name != "Diversity")
# 
# good_fcasts_bias <- eval_df %>% filter(metric == "BIAS" & (value > -.01 | value < .01))
# bad_fcasts_bias <- eval_df %>% filter(metric == "BIAS" & (value < -.015 | value > .015))
# 
# input_df <- hindcast_data
# plot_model(hindcast_data, siteID = "DSNY", taxon = "ectomycorrhizal", site_plots = "facet")
# plot_model(hindcast_data, siteID = "OSBS", taxon = "ectomycorrhizal", site_plots = "facet")
# plot_model(hindcast_data, siteID = "GUAN", taxon = "ectomycorrhizal", site_plots = "facet")
# 
# 
# good_fcasts_RMSE <- eval_df %>% filter(metric == "RMSE" & (value < .01))
# bad_fcasts_RMSE <- eval_df %>% filter(metric == "RMSE" & (value > .05))
# good_fcasts_RMSE[good_fcasts_RMSE$taxon=="lignolytic",]
# plot_model(hindcast_data, siteID = "DSNY", taxon = "lignolytic", site_plots = "facet")
# plot_model(hindcast_data, siteID = "SCBI", taxon = "lignolytic", site_plots = "facet")
# plot_model(hindcast_data, siteID = "UKFS", taxon = "lignolytic", site_plots = "facet")
# 
# good_fcasts_CRPS <- eval_df %>% filter(metric == "CRPS" & (value < .02))
# bad_fcasts_CRPS <- eval_df %>% filter(metric == "CRPS" & (value > .03))
# bad_fcasts_CRPS[bad_fcasts_CRPS$taxon=="saprotroph",]
# plot_model(hindcast_data, siteID = "HARV", taxon = "oligotroph", site_plots = "facet")
# 
# plot_model(hindcast_data, siteID = "OSBS", taxon = "saprotroph", site_plots = "facet")
# 
# good_CRPS_bad_RMSE <- merge(good_fcasts_CRPS, bad_fcasts_RMSE)
# good_RMSE_bad_CRPS <- merge(good_fcasts_RMSE, bad_fcasts_CRPS)
# good_RMSE_bad_bias <- merge(good_fcasts_RMSE, bad_fcasts_bias)
# good_CRPS_bad_bias <- merge(good_fcasts_CRPS, bad_fcasts_bias)
# 
# 
# plot_model(hindcast_data, siteID = "HARV", taxon = "plant_pathogen", site_plots = "facet")
# 
# 
# # Get good functional groups
# table(good_fcasts_CRPS[good_fcasts_CRPS$pretty_name == "Functional group",]$taxon) %>% sort()
# table(good_fcasts_RMSE[good_fcasts_RMSE$pretty_name == "Functional group",]$taxon) %>% sort()
# 
# 
# 
# 
# 
# for_examples <- cv_metric_scaled %>%
# 	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
# 				 	cv_type == "mean_per_plot_cv" & metric == "CRPS")
# cv_tax_per_site[cv_tax_per_site$taxon=="lignolytic",]
# plot_model(hindcast_data, siteID = "ORNL", taxon = "lignolytic", site_plots = "facet")
# plot_model(hindcast_data, siteID = "UNDE", taxon = "lignolytic", site_plots = "facet")
# 
# scoring_metrics_cv_site[scoring_metrics_cv_site$taxon=="lignolytic",] %>% filter(model_name == "all_covariates")
# cv_tax_per_plot_site[cv_tax_per_plot_site$taxon=="lignolytic",]









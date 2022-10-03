# Calculate and visualize scoring metrics for taxonomic forecasts

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


p_load(Metrics, forestmangr, scoringRules, ggpubr)

hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds")) %>% mutate(truth=as.numeric(truth))

hindcast_data_non_tax = hindcast_data %>% filter(fcast_type != "Taxonomic")
hindcast_data_tax = hindcast_data %>% filter(fcast_type == "Taxonomic")
hindcast_only = hindcast_data %>% filter(#fcast_type == "Taxonomic" &
																				 	fcast_period=="hindcast")
calibration_only = hindcast_data %>% filter(#fcast_type == "Taxonomic" &
																							fcast_period=="calibration")


# pulled model dat from refit script with j=1 or j=6, run this section twice

model.dat_bac = quick_get_rank_df(k = 1)
model.dat_fun = quick_get_rank_df(k = 6)


bac_site_start_dates <- model.dat_bac$site_start %>% stack %>% mutate(pretty_group = "Bacteria")
bac_plot_start_dates <- model.dat_bac$plot_start %>% stack %>% mutate(pretty_group = "Bacteria")
colnames(bac_site_start_dates) <- c("site_start_date", "siteID", "pretty_group")
colnames(bac_plot_start_dates) <- c("plot_start_date", "plotID", "pretty_group")
bac_plot_start_dates$plotID <- as.character(bac_plot_start_dates$plotID)

fun_site_start_dates <- model.dat_fun$site_start %>% stack %>% mutate(pretty_group = "Fungi")
fun_plot_start_dates <- model.dat_fun$plot_start %>% stack %>% mutate(pretty_group = "Fungi")
colnames(fun_site_start_dates) <- c("site_start_date", "siteID", "pretty_group")
colnames(fun_plot_start_dates) <- c("plot_start_date", "plotID", "pretty_group")
fun_plot_start_dates$plotID <- as.character(fun_plot_start_dates$plotID)

#hindcast_data_select <- hindcast_data_tax %>% select(-truth)
#hindcast_data_newtax = merge(hindcast_data_select, truth_plot_select, all.x = T)
#
#not_na <- hindcast_data_newtax %>% filter(!is.na(truth) & dateID > 201701)
#
## Subset to hindcast observations
#hindcast_only <- hindcast_data_merge %>% filter(!is.na(truth) &
#																						fcast_period=="hindcast"
#																					#dates > "2016-12-31"
#) %>%
#	mutate(truth = as.numeric(truth))
calibration_bac <- merge(calibration_only %>% filter(pretty_group=="Bacteria"),
												 bac_site_start_dates, by=c("siteID","pretty_group"))
calibration_bac <- merge(calibration_bac,
												 bac_plot_start_dates, by=c("plotID","pretty_group"))

calibration_fun <- merge(calibration_only %>% filter(pretty_group=="Fungi"),
												 fun_site_start_dates, by=c("siteID","pretty_group"))
calibration_fun <- merge(calibration_fun,
												 fun_plot_start_dates, by=c("plotID","pretty_group"))

calibration_only = rbind.data.frame(calibration_bac, calibration_fun)

calibration_only$is_site_start_date = ifelse(calibration_only$date_num==calibration_only$site_start_date, T, F)
calibration_only$is_plot_start_date = ifelse(calibration_only$date_num==calibration_only$plot_start_date, T, F)
calibration_only$is_any_start_date = ifelse(calibration_only$is_plot_start_date|calibration_only$is_site_start_date, T, F)

calibration_only_not_first = calibration_only %>% filter(timepoint > plot_start_date)


fg_means <- readRDS(here("data", paste0("summary/fg_summaries_20151101_20180101.rds")))
fg_means <- fg_means$scores.list


calibration_only_not_fg <- calibration_only_not_first %>% filter(fcast_type != "Functional group")
calibration_only_fg <- calibration_only_not_first %>% filter(fcast_type == "Functional group")


calibration_only_fg <- merge(calibration_only_fg[,-c("Mean","SD")],
														 fg_means[,c("Mean", "SD", "siteID", "plotID", "dateID", "species",
																															"truth", "dates", "model_name", "taxon", "time_period")],
														 by = c("plotID", "siteID", "dates", "model_name", "time_period", "dateID",
														 			 "species", "truth", "taxon"), all=T)

calibration_only_not_first = rbind.data.frame(calibration_only_not_fg, calibration_only_fg)
#
# # Convert prediction SD to variance
# pl <- test %>% filter(!is.na(mean))
# pl$pred_variance = pl$sd^2
# pl$resid = (pl$mean - pl$truth)^2
# P <- sum(pl$pred_variance)/nrow(pl)
# G <- sum(pl$resid)/nrow(pl)
# total_PL = P+G

scoring_metrics <- hindcast_only %>%
	filter(!is.na(truth) & !is.na(mean)) %>%
	filter(newsite=="Observed site") %>%
	group_by(pretty_group,model_name,pretty_name,taxon) %>%
	summarize(add_scoring_metrics(observed = truth,
														 mean_predicted = mean,
														 sd_predicted = sd))

# ggplot(scoring_metrics, aes(x = CRPS, y = total_PL)) +
# 	geom_point(alpha=.4) +
# 	geom_abline(slope = 1)


# test <- hindcast_only %>% filter(!is.na(truth) & newsite=="Observed site") %>%
# 	group_by(pretty_group,model_name,pretty_name,taxon) %>%
# first()
# test1 = test %>% 	summarize(add_scoring_metrics(observed = truth,
# 																								mean_predicted = Mean,
# 																								sd_predicted = SD))
# ggplot(test, aes(x = mean, y = truth)) + geom_point() +
# 	geom_abline(slope = 1) +
# 	ggtitle("Actinobacteria calibration",
# 					subtitle = paste0("RSQ Colin: ", test1$RSQ.1.colin,
# 														"\nRSQ Mike: ", test1$RSQ.1.mike)) +
# 	geom_smooth(method="lm") + xlab("predicted") + ylab("observed")

calibration_metrics <- calibration_only_not_first %>%
	filter(!is.na(truth) & !is.na(Mean)) %>%
	group_by(pretty_group,model_name,pretty_name,taxon) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = Mean,
																sd_predicted = SD))
calibration_metrics_long <- calibration_metrics %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1, RMSE.norm, residual_variance, predictive_variance, total_PL),
																												 names_to = "metric")


#low_rsq1 <- calibration_metrics[calibration_metrics$model_name=="all_covariates" &
#																	calibration_metrics$RSQ.1 < -10,]$taxon %>% unique
#calibration_metrics_test <- calibration_only_not_first %>%
#	filter(!is.na(truth) & !is.na(Mean)) %>%
#	filter(taxon %in% low_rsq1 & model_name=="all_covariates") %>%
#	group_by(pretty_group,model_name,pretty_name,rank_name,taxon) %>%
#	mutate(add_scoring_metrics(observed = truth,
#																mean_predicted = Mean,
#																sd_predicted = SD)) %>%
#	mutate(label = paste(model_name, taxon, "\nRSQ 1:1: ", RSQ.1, "\nRSQ: ", RSQ))
#
#ggplot(calibration_metrics_test, aes(x = Mean, y = truth)) +
#	geom_point(aes(color=is_plot_start_date), alpha=.4) +
#	geom_abline(slope = 1) +
#	geom_smooth(method="lm") + xlab("predicted") + ylab("observed") +
#	facet_wrap(~label, drop=T, scales="free")


calibration_metrics_site <- calibration_only_not_first %>%
	filter(!is.na(truth) & !is.na(Mean)) %>%
	group_by(pretty_group,model_name,pretty_name,taxon, siteID) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = Mean,
																sd_predicted = SD))
scoring_metrics_site <- hindcast_only %>% #filter(newsite=="Observed site") %>%
	filter(!is.na(truth) & !is.na(mean)) %>%
	group_by(fcast_type, pretty_group,model_name,pretty_name, taxon,siteID) %>%
	#filter(!all(is.na(truth))) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = mean,
																sd_predicted = sd))

scoring_metrics_long <- scoring_metrics %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1, RMSE.norm,residual_variance, predictive_variance, total_PL),
																												 names_to = "metric")
scoring_metrics_site_long <- scoring_metrics_site %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1, RMSE.norm,residual_variance, predictive_variance, total_PL),
																												 names_to = "metric")

# check correlations between metrics

#library(corrplot)
#for_cor <- scoring_metrics_site %>% ungroup() %>% select(RMSE, MAE, CRPS, RSQ, RSQ.1, RMSE.norm)
#cor_scores <- cor(for_cor, use = "pairwise.complete.obs")
# corrplot::corrplot(cor_scores)
# ggplot(for_cor, aes(x = CRPS, y = RSQ.1)) +
# 	geom_point(aes(color=siteID), alpha=.4) +
# 	#geom_abline(slope = 1) +
# 	#ggtitle("Actinobacteria calibration",
# 	#				subtitle = paste0("RSQ Colin: ", test1$RSQ.1.colin,
# 	#													"\nRSQ Mike: ", test1$RSQ.1.mike)) +
# 	geom_smooth(method="lm")

# Add coefficient of variation

# Wide format taxonomic
cv_tax_long <- hindcast_data %>%
	filter(model_name == "all_covariates" & taxon != "other" & !is.na(truth)) %>%
	mutate(truth = as.numeric(truth))

cv_tax_per_plot <- cv_tax_long %>%
	group_by(fcast_type, pretty_name, taxon, siteID, plotID) %>%
	dplyr::summarize(per_plot_cv = calc_cv(truth))
cv_tax_per_plot_site <- cv_tax_per_plot %>% ungroup %>%
	group_by(fcast_type, pretty_name, taxon, siteID) %>%
	dplyr::summarize(mean_per_plot_site_cv = mean(per_plot_cv, na.rm=T))

cv_tax_per_site <- cv_tax_long %>%  group_by(fcast_type, pretty_name, taxon, siteID) %>%
	summarize(per_site_cv = calc_cv(truth))

cv_tax_per_site_taxon <- cv_tax_per_site %>% ungroup %>%
	group_by(fcast_type, pretty_name, taxon) %>%
	summarize(mean_per_site_cv = mean(per_site_cv, na.rm=T))

cv_tax_per_plot_taxon <- cv_tax_long %>%  group_by(fcast_type, pretty_name, taxon, siteID, plotID) %>%
	summarize(per_plot_cv = calc_cv(truth)) %>%
	ungroup %>%
	group_by(fcast_type, pretty_name, taxon) %>%
	summarize(mean_per_plot_cv = mean(per_plot_cv, na.rm=T))

cv_tax_overall <- cv_tax_long %>%  group_by(fcast_type, pretty_name, taxon) %>%
	summarize(overall_cv = calc_cv(truth))

cv_tax <- cv_tax_overall %>% merge(cv_tax_per_plot_taxon) %>% merge(cv_tax_per_site_taxon)

scoring_metrics_cv <- merge(scoring_metrics, cv_tax, all=T)
scoring_metrics_cv_site <- merge(scoring_metrics_site, cv_tax_per_site)



### Predictability
scoring_metrics_cv_long <- scoring_metrics_cv %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1, RMSE.norm, residual_variance, predictive_variance, total_PL),
																															 names_to = "metric", values_to = "score") %>%
	pivot_longer(cols = c(overall_cv, mean_per_plot_cv, mean_per_site_cv),
							 names_to = "cv_type", values_to = "cv")

cv_long <- scoring_metrics_cv %>% pivot_longer(cols = c(overall_cv, mean_per_plot_cv, mean_per_site_cv),
																							 names_to = "cv_type", values_to = "cv")



# CV scaled by rank
cv_metric_scaled <- scoring_metrics_cv_long %>%
	group_by(pretty_group, pretty_name, cv_type, metric) %>%
	mutate(CV_scale = scale(cv)[,1],
				 metric_scale = scale(score)[,1])

crps_site_obs <- hindcast_only %>% #filter(newsite=="Observed site") %>%
	filter(!is.na(truth) & !is.na(mean)) %>%
	group_by(pretty_group,fcast_type,model_name,pretty_name, taxon,newsite) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = mean,
																sd_predicted = sd, type = c("CRPS","CRPS_truncated")))

# Skill score from Scavia 2021
skill_score_taxon <- crps_site_obs %>%
	pivot_wider(id_cols = c("fcast_type","pretty_group","model_name","pretty_name","taxon"),
							values_from = "CRPS_truncated", names_from = "newsite") %>%
	mutate(skill_score = (1 - (`New site`/`Observed site`)))
skill_score_rank <- skill_score_taxon %>%
	group_by(pretty_group, pretty_name) %>%
	summarize(mean_skill_score = mean(skill_score, na.rm=T))


to_save <- list(cv_metric_scaled = cv_metric_scaled,
								cv_tax = cv_tax,
								cv_tax_long = cv_tax_long,
								cv_tax_overall = cv_tax_overall,
								cv_tax_per_plot = cv_tax_per_plot,
								cv_tax_per_plot_site = cv_tax_per_plot_site,
								cv_tax_per_plot_taxon = cv_tax_per_plot_taxon,
								cv_tax_per_site = cv_tax_per_site,
								cv_tax_per_site_taxon = cv_tax_per_site_taxon,
								scoring_metrics = scoring_metrics,
								scoring_metrics_cv = scoring_metrics_cv,
								scoring_metrics_cv_long = scoring_metrics_cv_long,
								scoring_metrics_cv_site = scoring_metrics_cv_site,
								scoring_metrics_long = scoring_metrics_long,
								scoring_metrics_site = scoring_metrics_site,
								scoring_metrics_site_lon = scoring_metrics_site_long,
								calibration_metrics=calibration_metrics,
								calibration_metrics_long = calibration_metrics_long,
								calibration_metrics_site = calibration_metrics_site,
								skill_score_taxon = skill_score_taxon,
								skill_score_rank = skill_score_rank)
saveRDS(to_save, here("data", paste0("summary/scoring_metrics_cv.rds")))




#### Visualizing





#### View predictability scores by rank
ggplot(scoring_metrics_long %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
			 aes(x = pretty_name, y = value,
			 		color = pretty_name)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
	facet_grid(metric ~ pretty_group, drop = T, scales="free") +
	# geom_jitter(aes(x = metric, y = value), width=.1,
	# 						height = 0, alpha = .8, size=4) +
	ylab("Metric scores") + xlab(NULL) +
	theme_bw(base_size=18) +
	ggtitle("Predictability decreases at broader ranks")  +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm"))


#### View CV across tax ranks
ggplot(cv_long %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
			 aes(x = pretty_name, y = cv,
			 		color = pretty_name), alpha = .5) +
	geom_point(size = 3, position = position_jitterdodge(#jitter.width = 1
	), alpha=.5, show.legend = F) +
	facet_grid(pretty_group~cv_type, scales="free", drop =T) + theme_bw(base_size=18) +
	ggtitle("Variability decreases with broader taxonomic ranks") +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + xlab(NULL) +
	theme(plot.margin = margin(1,2,1,1, "cm"))


#### View CV across tax ranks, for each plot
# Just per-plot CV (capturing temporal variability)
ggplot(cv_long %>% filter(model_name == "all_covariates" & cv_type == "mean_per_plot_cv"),
			 aes(x = pretty_name, y = cv,
			 		color = pretty_name), alpha = .5) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(#jitter.width = 1
	), alpha=.5, show.legend = F) +
	facet_grid(rows=vars(pretty_group), scales="free", drop =T) + theme_bw(base_size=18) +
	ggtitle("Variability decreases with broader taxonomic ranks") +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + xlab(NULL)



# Just RMSE
ggplot(scoring_metrics_cv %>% filter(model_name == "all_covariates" & pretty_name != "Diversity"),
			 aes(x = mean_per_plot_cv, y = RMSE,
			 		color = pretty_name)) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
	facet_grid(~pretty_group, drop = T, scales="free") + theme_bw()




# All metrics
ggplot(cv_metric_scaled %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & cv_type == "mean_per_plot_cv"),
			 aes(x = CV_scale, y = metric_scale,
			 		shape = pretty_group,
			 		color = pretty_name)) +
	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
						 alpha=.5, show.legend = T) +
	facet_wrap(~metric,drop = T, scales="free") +
	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
																		shape=guide_legend(title="Domain")) +
	xlab("Coefficient of variation (scaled by rank)") +
	ylab("Predictability (scaled by rank)") +
	ggtitle("Predictability increases with variability, but only for some metrics")

# NOt scaled
ggplot(cv_metric_scaled %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & cv_type == "mean_per_plot_cv"),
			 aes(x = cv, y = score,
			 		shape = pretty_group,
			 		color = pretty_name)) +
	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
						 alpha=.5, show.legend = T) +
	facet_wrap(~metric,drop = T, scales="free") +
	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
																		shape=guide_legend(title="Domain")) +
	xlab("Coefficient of variation (NOT mean-scaled by rank)") +
	ylab("Predictability (NOT mean-scaled by rank)") +
	ggtitle("Predictability increases with variability, but only for some metrics")

# Just CRPS
ggplot(cv_metric_scaled %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
			 				 	cv_type == "mean_per_plot_cv" & metric == "CRPS"),
			 aes(x = CV_scale, y = metric_scale,
			 		shape = pretty_group,
			 		color = pretty_name)) +
	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
						 alpha=.5, show.legend = T) +
	facet_grid(~metric, drop = T, scales="free") +
	theme_bw(base_size = 18) + guides(color=guide_legend(title="Rank"),
																		shape=guide_legend(title="Domain")) +
	#	stat_regline_equation(aes(y = CV_scale, x = metric_scale), inherit.aes = F, label.x = 1, label.y = 2) +
	stat_smooth(aes(x = CV_scale, y = metric_scale), inherit.aes = F, method = "loess", span = .9)

# F vs B grid, CV types (split by rank)
ggplot(cv_metric_scaled %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & metric=="CRPS"),
			 aes(x = pretty_name, y = cv,
			 		shape = pretty_group,
			 		color = cv_type)) +
	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
						 alpha=.5, show.legend = T) +
	facet_grid(~pretty_group, drop = T, scales="free") +
	theme_bw(base_size = 18) + guides(color=guide_legend(title="Coeff. of \nVariation type"),
																		shape=guide_legend(NULL)) + xlab(NULL)  +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))



# View skill scores
ggplot(skill_score_taxon %>% filter(model_name == "all_covariates"),
			 aes(x = pretty_name, y = skill_score,
			 		color = pretty_group)) +
	#geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = pretty_name, y = skill_score), position = position_jitterdodge(jitter.height = 0, jitter.width = .1), alpha = .5, size=4) +
	ylab("Skill score (% decrease in CRPS)") + xlab(NULL) +
	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=22))  +
	guides(color=guide_legend(title="Domain"),
				 shape=guide_legend(title="Forecast type"))


# Predictability within-site
ggplot(scoring_metrics_site %>%
			 	filter(model_name == "all_covariates"),
			 aes(x = pretty_group, y = CRPS,
			 		color = pretty_name)) +
	geom_jitter(size = 4, alpha=.5, show.legend = T) +
	#facet_wrap(~beta, drop = T, scales="free") +
	theme_bw(base_size = 18) +
	xlab("Rank") +
	ylab("CRPS") + ggtitle("Mean predictability at each observed site") +
	guides(color=guide_legend(title="Rank")) +scale_y_sqrt()

ggplot(scoring_metrics_site %>% filter(model_name == "all_covariates" &
																			 	#newsite == "Observed site" &
																			 	pretty_name != "Diversity"),
			 aes(x = pretty_name, y = CRPS,
			 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = pretty_name, y = CRPS), position=position_jitterdodge(jitter.height = 0),
						 alpha = .2, size=4) +
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) +
	theme_minimal(base_size=20) + ggtitle("Mean predictability at each observed site") +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=22))  +
	guides(color=guide_legend(title="Domain"),
				 shape=guide_legend(title="Forecast type")) +scale_y_sqrt()


# F vs B "other" abundances at each rank
ggplot(hindcast_data %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" & species=="other"),
			 aes(x = pretty_name, y = truth,
			 		shape = pretty_group)) +
	geom_point(size = 4, #position = position_jitterdodge(jitter.width = .5),
						 alpha=.5, show.legend = T) +
	facet_grid(~pretty_group, drop = T, scales="free") +
	theme_bw(base_size = 18) + guides(color=guide_legend(title="Coeff. of \nVariation type"),
																		shape=guide_legend(NULL)) + xlab(NULL)  +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))





eval_df <- scoring_metrics_site_long %>% filter(model_name == "all_covariates" &
																					 	pretty_name != "Diversity")

good_fcasts_bias <- eval_df %>% filter(metric == "BIAS" & (value > -.01 | value < .01))
bad_fcasts_bias <- eval_df %>% filter(metric == "BIAS" & (value < -.015 | value > .015))

input_df <- hindcast_data
plot_model(hindcast_data, siteID = "DSNY", taxon = "ectomycorrhizal", site_plots = "facet")
plot_model(hindcast_data, siteID = "OSBS", taxon = "ectomycorrhizal", site_plots = "facet")
plot_model(hindcast_data, siteID = "GUAN", taxon = "ectomycorrhizal", site_plots = "facet")


good_fcasts_RMSE <- eval_df %>% filter(metric == "RMSE" & (value < .01))
bad_fcasts_RMSE <- eval_df %>% filter(metric == "RMSE" & (value > .05))
good_fcasts_RMSE[good_fcasts_RMSE$taxon=="lignolytic",]
plot_model(hindcast_data, siteID = "DSNY", taxon = "lignolytic", site_plots = "facet")
plot_model(hindcast_data, siteID = "SCBI", taxon = "lignolytic", site_plots = "facet")
plot_model(hindcast_data, siteID = "UKFS", taxon = "lignolytic", site_plots = "facet")

good_fcasts_CRPS <- eval_df %>% filter(metric == "CRPS" & (value < .02))
bad_fcasts_CRPS <- eval_df %>% filter(metric == "CRPS" & (value > .03))
bad_fcasts_CRPS[bad_fcasts_CRPS$taxon=="saprotroph",]
plot_model(hindcast_data, siteID = "HARV", taxon = "oligotroph", site_plots = "facet")

plot_model(hindcast_data, siteID = "OSBS", taxon = "saprotroph", site_plots = "facet")

good_CRPS_bad_RMSE <- merge(good_fcasts_CRPS, bad_fcasts_RMSE)
good_RMSE_bad_CRPS <- merge(good_fcasts_RMSE, bad_fcasts_CRPS)
good_RMSE_bad_bias <- merge(good_fcasts_RMSE, bad_fcasts_bias)
good_CRPS_bad_bias <- merge(good_fcasts_CRPS, bad_fcasts_bias)


plot_model(hindcast_data, siteID = "HARV", taxon = "plant_pathogen", site_plots = "facet")


# Get good functional groups
table(good_fcasts_CRPS[good_fcasts_CRPS$pretty_name == "Functional group",]$taxon) %>% sort()
table(good_fcasts_RMSE[good_fcasts_RMSE$pretty_name == "Functional group",]$taxon) %>% sort()








for_examples <- cv_metric_scaled %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
				 	cv_type == "mean_per_plot_cv" & metric == "CRPS")
cv_tax_per_site[cv_tax_per_site$taxon=="lignolytic",]
plot_model(hindcast_data, siteID = "ORNL", taxon = "lignolytic", site_plots = "facet")
plot_model(hindcast_data, siteID = "UNDE", taxon = "lignolytic", site_plots = "facet")

scoring_metrics_cv_site[scoring_metrics_cv_site$taxon=="lignolytic",] %>% filter(model_name == "all_covariates")
cv_tax_per_plot_site[cv_tax_per_plot_site$taxon=="lignolytic",]









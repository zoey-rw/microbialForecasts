

# Calculate forecast horizons
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

core_scores_in = readRDS(here("data", paste0("summary/scoring_metrics_core_level.rds")))



# converged = converged_strict
core_hindcast_data = core_scores_in[[3]]
core_cal_hindcast = core_scores_in[[4]]

hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))
hindcast_data = hindcast_in

# Determine each plot's final observation in calibration
last_obs = hindcast_data %>%
	#filter(model_id %in% converged) %>%
	group_by(model_id, plotID) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	filter(timepoint == max(timepoint, na.rm=T)) %>%
	mutate(last_obs = timepoint) %>%
	select(plotID, last_obs) %>% distinct()

hindcast_data[hindcast_data$truth==0,]$truth <- .0001


# Calculate "time since final observation" for forecasts
fcast_horizon_df = core_hindcast_data %>%
	#filter(model_id %in% converged) %>%
	filter(!is.na(core_truth) & fcast_period=="hindcast" & new_site==FALSE) %>%
	merge(last_obs, by=c("plotID","model_id")) %>%
	group_by(model_id,model_name) %>% 
	mutate(months_since_obs = timepoint - last_obs)


# Calculate average scores by model
fcast_horizon_model_mean = fcast_horizon_df %>%
	group_by(taxon,model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(mean_crps = mean(crps, na.rm=T),
						RMSE = rmse(actual = core_truth, predicted = med),
						RSQ.1 = 1 - (RMSE^2)/var(core_truth),
						MAPE = mape(actual = core_truth, predicted = med),
						RSQ = summary(lm(core_truth ~ med))$r.squared,
						abundance = mean(core_truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
fcast_horizon_model_mean$RSQ.1 = ifelse(fcast_horizon_model_mean$RSQ.1 < 0, 0, fcast_horizon_model_mean$RSQ.1)
fcast_horizon_model_mean$RMSE.norm = ifelse(fcast_horizon_model_mean$RMSE.norm > 5, 5, fcast_horizon_model_mean$RMSE.norm)

# Calculate number of observations per months out
fcast_n_obs = fcast_horizon_df %>%
	group_by(model_id, months_since_obs) %>%
	tally(name = "n_obs")
# Remove metrics calculated from only 1-2 observations
fcast_horizon_model_mean <- merge(fcast_horizon_model_mean, fcast_n_obs) 


# Calculate the historical mean abundances (within and across sites), to be the "null"
# Also tried this with historical site median abundances, but site median predictions were worse
site_mean = core_cal_hindcast %>% group_by(taxon,model_id,model_name,siteID) %>%
	filter(fcast_period != "hindcast" & !is.na(core_truth)) %>%
	summarize(site_mean = mean(core_truth, na.rm=T),
						site_sd = sd(core_truth, na.rm=T)) %>% ungroup
overall_mean = core_cal_hindcast %>% group_by(taxon,model_id) %>%
	filter(fcast_period != "hindcast" & !is.na(core_truth)) %>%
	summarize(overall_mean = mean(core_truth, na.rm=T),
						overall_sd = sd(core_truth, na.rm=T)) %>% ungroup
historical_mean <- merge(site_mean, overall_mean)

# Merge historical means into main df, and use to calculate null scores
fcast_horizon_df2 <- merge(fcast_horizon_df, historical_mean, by=c("taxon","model_id","siteID","model_name"), all.x=T)


fcast_horizon_null_site = fcast_horizon_df2 %>%
	filter(site_sd != 0) %>% 
	group_by(taxon,model_name,pretty_group, rank_only, model_id) %>%
	summarize(null_CRPS = mean(crps_norm(core_truth, site_mean,
																			 site_sd)),
						null_CRPS_truncated = mean(crps(core_truth, family = "tnorm",
																						location = site_mean, scale = site_sd,
																						lower = 0, upper = 1)),
						null_RMSE = rmse(actual = core_truth, predicted = site_mean),
						null_RSQ.1 = 1 - (null_RMSE^2)/var(core_truth),
						null_MAPE = mape(actual = core_truth, predicted = site_mean),
						null_RSQ = summary(lm(core_truth ~ site_mean))$r.squared,
						abundance = mean(core_truth, na.rm = T), RMSE.norm = null_RMSE/abundance) %>% distinct()
fcast_horizon_null_site$null_RSQ.1 = ifelse(fcast_horizon_null_site$null_RSQ.1 < 0, 0, fcast_horizon_null_site$null_RSQ.1)
fcast_horizon_null_site$months_since_obs = Inf


last_obs_values = inner_join(core_cal_hindcast, last_obs %>% mutate(timepoint=last_obs))


last_obs_score <- last_obs_values %>%
	mutate(months_since_obs = 0) %>% 
	group_by(taxon,model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(CRPS = mean(crps_norm(core_truth, med,
																	sd)),
						CRPS_truncated = mean(crps(core_truth, family = "tnorm",
																			 location = med, scale = sd,
																			 lower = 0, upper = 1)),
						RMSE = rmse(actual = core_truth, predicted = med),
						RSQ.1 = 1 - (RMSE^2)/var(core_truth),
						MAPE = mape(actual = core_truth, predicted = med),
						RSQ = summary(lm(core_truth ~ med))$r.squared,
						abundance = mean(core_truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
last_obs_score$RSQ.1 = ifelse(last_obs_score$RSQ.1 < 0, 0, last_obs_score$RSQ.1)
last_obs_score$RMSE.norm = ifelse(last_obs_score$RMSE.norm > 5, 5, last_obs_score$RMSE.norm)


# Combine for plotting!
colnames(fcast_horizon_null_site) = gsub("null_","", colnames(fcast_horizon_null_site))
comparison_scores <- rbindlist(list(fcast_horizon_null_site, last_obs_score), fill=T)
to_plot <- rbindlist(list(comparison_scores, fcast_horizon_model_mean %>% filter(n_obs > 30)), fill=T)
to_plot_null <-  fcast_horizon_null_site %>% filter(model_id %in%  to_plot$model_id)


saveRDS(list(to_plot, to_plot, fcast_horizon_null_site, fcast_horizon_model_mean, last_obs_score, historical_mean),
				here("data/summary/fcast_horizon_input_core.rds"))






in_list <- readRDS(here("data/summary/fcast_horizon_input_core.rds"))
to_plot <- in_list[[1]]
to_plot_null <- in_list[[2]]
fcast_horizon_null_site <-  in_list[[3]]
fcast_horizon_model_mean <-  in_list[[4]]
model_id_list <- unique(to_plot$model_id)
model_id=model_id_list[[1]]


out_df_list = list()
out_figure_list = list()

to_plot <- in_list[[1]]%>% filter(model_id %in% converged & RSQ != 0)
to_plot$mean_crps = ifelse(is.na(to_plot$mean_crps), to_plot$CRPS_truncated, to_plot$mean_crps)
to_plot_null$mean_crps = to_plot_null$CRPS_truncated
to_plot$months_since_obs = to_plot$months_since_obs + .001

to_model = single_tax %>% filter(!is.na(months_since_obs) & !is.infinite(months_since_obs))
to_model_all = to_plot %>% filter(!is.na(months_since_obs) & !is.infinite(months_since_obs))


model_results = to_model_all %>% 
	# filter(model_id %in% c("env_cov_aspergillus_20151101_20180101","cycl_only_sordariales_20151101_20180101","env_cycl_aspergillus_20151101_20180101")) %>% 
	group_by(pretty_group,model_name,model_id) %>% 
	nest() %>% 
	mutate(model = map(data, ~glm(RSQ ~ months_since_obs, data = ., family = gaussian(link = 'log'))),
				 tidied = map(model, tidy)) 


	#left_join(pframe, by = "model_id") %>% 
fit_out=list()
for (i in 1:nrow(model_results)) {
	fit_out[[i]] = predict(model_results$model[[i]], newdata=data.frame(months_since_obs=1:30), type='response') %>% 
		stack() %>% 
		mutate(model_id=model_results$model_id[[i]], 
					 model_name=model_results$model_name[[i]],
					 pretty_group=model_results$pretty_group[[i]],
					 p_value=model_results$tidied[[i]][2,]$p.value,
					 RSQ=values,
					 months_since_obs =as.numeric(ind ))
}

fit_combined <- rbindlist(fit_out) 
unconverged <- fit_combined[fit_combined$RSQ > .4 & fit_combined$months_since_obs > 18 ,]$model_id %>% unique

#unconverged <- fit_combined[fit_combined$p_value > .05,]$model_id %>% unique
#unconverged= c("env_cov_aspergillus_20151101_20180101","cycl_only_sordariales_20151101_20180101","env_cycl_aspergillus_20151101_20180101")

fit_combined <- fit_combined %>% filter(!model_id %in% unconverged)

null_to_plot = fcast_horizon_null_site %>% select(c(1:5, RMSE_null = RMSE.norm,
																										RSQ_null = RSQ, CRPS_null = CRPS_truncated))
fit_combined = merge(fit_combined, null_to_plot)
null_to_plot2 = fit_combined[fit_combined$RSQ > RSQ_null,]
null_to_plot3 = null_to_plot2 %>% group_by(model_id) %>% filter(RSQ == min(RSQ))
																									
decay_plot1 = ggplot(fit_combined,
			 aes(x = months_since_obs, y = RSQ, group=model_id, color=model_name)) +
	xlim(c(0,30)) +
	geom_line(alpha=.3) +
	geom_point(data=null_to_plot3, 
						 aes(y = RSQ_null, x = months_since_obs, 
						 		group=model_id), color="black", 
						 na.rm = T, alpha=.3, size=2, shape=4, position=position_jitter(width = .1, height=.01)) + 
	theme_minimal(base_size = 16) +
	facet_grid(rows=vars(model_name), labeller = labeller(model_name = model.labs))

summary_plot1 = ggplot(null_to_plot3,
			 aes(x = model_name, y = months_since_obs, color = model_name)) +
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitterdodge(jitter.width = .2, jitter.height = .2), alpha=.3) + 
	theme_minimal(base_size = 16) 
ggarrange(decay_plot1, summary_plot1, nrow=2, common.legend = T)


decay_plot2 = ggplot(fit_combined,
										aes(x = months_since_obs, y = RSQ, group=model_id, color=pretty_group)) +
	xlim(c(0,30)) +
	geom_line(alpha=.3) +
	geom_point(data=null_to_plot3, 
						 aes(y = RSQ_null, x = months_since_obs, 
						 		group=model_id), color="black", 
						 na.rm = T, alpha=.3, size=2, shape=4, position=position_jitter(width = .2)) + 
	theme_minimal(base_size = 16) +facet_grid(rows=vars(pretty_group))

summary_plot2 = ggplot(null_to_plot3,
											 aes(x = pretty_group, y = months_since_obs, color=pretty_group)) +
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitterdodge(jitter.width = .2, jitter.height = .2), alpha=.3) + 
	theme_minimal(base_size = 16) 

ggarrange(decay_plot2, summary_plot2, nrow=2)

kingdom_stat_pvalue <- null_to_plot3 %>% 
	#group_by(pretty_group) %>% 
	rstatix::tukey_hsd(months_since_obs ~ pretty_group) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position(step.increase = .4) %>% 
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))


	glm(RSQ~months_since_obs, family = gaussian(link = 'log'), data = .)
glm_fit = glm(RSQ~months_since_obs, family = gaussian(link = 'log'), data = to_model)
## generate prediction frame
pframe <- expand.grid(model_id = unique(to_model_all$model_id), months_since_obs=1:20)
## add predicted values (on response scale) to prediction frame
pframe$RSQ <- predict(glm_fit,newdata=pframe,type="response")

to_plot_null_df = fcast_horizon_null_site %>% filter(model_id %in% converged) 
ggplot(to_plot,
			 aes(x = months_since_obs, y = RSQ, group=model_id)) +
	xlim(c(0,20)) +
	geom_point(alpha=.1) +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =to_plot_null_df, aes(yintercept = RSQ,
																				group=model_name), na.rm = T, alpha=.1, linetype=2) + 
	geom_smooth(aes(color = model_name), method = "glm", formula = y~x,	
							method.args = list(family = gaussian(link = 'log')), 
							fullrange=TRUE, se=F, color="black", alpha=.5)


#Run for multiple ranks, in parallel
cl <- makeCluster(28, type="FORK", outfile = '/projectnb/dietzelab/zrwerbin/microbialForecasts/log.txt')
registerDoParallel(cl)

# for testing
out_parallel_list = foreach(model_id = model_id_list, 
														.errorhandling = 'remove') %dopar% {

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
															
print(model_id)

single_tax = to_plot %>% filter(model_id == !!model_id)
single_tax$mean_crps = ifelse(is.na(single_tax$mean_crps), single_tax$CRPS_truncated, single_tax$mean_crps)

single_tax_null <-  fcast_horizon_null_site %>% filter(model_id == !!model_id)
single_tax_null$mean_crps = single_tax_null$CRPS_truncated
single_tax$months_since_obs = single_tax$months_since_obs + .001

fig_rsq <- ggplot(single_tax,
									aes(x = months_since_obs, y = RSQ)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = RSQ,
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)

fig_crps <- ggplot(single_tax,
									 aes(x = months_since_obs, y = mean_crps)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (CRPS)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = mean_crps,
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)

fig_rmse <- ggplot(single_tax,
									 aes(x = months_since_obs, y = RMSE.norm)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (RMSE)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = RMSE.norm,
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)

save_fig = ggarrange(fig_rsq + rremove("xlab"), 
										 fig_crps + rremove("xlab"), 
										 fig_rmse, labels = c("A","B","C"), ncol=1) + ggtitle(model_id)

# Extract information about plot
horizon_rsq_data = ggplot_build(fig_rsq)$data
rsq_loess_line = horizon_rsq_data[[3]]
rsq_null_line = horizon_rsq_data[[2]]$yintercept

# Get intersection: first point at which the fit is lower than the null Rsq
rsq_fcast_horizon <- rsq_loess_line$x[which(rsq_loess_line$y < rsq_null_line)] %>% min(na.rm=T)

# Extract information about plot
horizon_rmse_data = ggplot_build(fig_rmse)$data
rmse_loess_line = horizon_rmse_data[[3]]
rmse_null_line = horizon_rmse_data[[2]]$yintercept

# Get intersection: first point at which the  fit is lower than the null rmse
rmse_fcast_horizon <- rmse_loess_line$x[which(rmse_loess_line$y > rmse_null_line)] %>% min(na.rm=T)

# Extract information about plot
horizon_crps_data = ggplot_build(fig_crps)$data
crps_loess_line = horizon_crps_data[[3]]
crps_null_line = horizon_crps_data[[2]]$yintercept

# Get intersection: first point at which the fit is lower than the null crps
crps_fcast_horizon <- crps_loess_line$x[which(crps_loess_line$y > crps_null_line)] %>% min(na.rm=T)

out_df = cbind.data.frame(single_tax[1,1:5], 
													rsq_fcast_horizon = rsq_fcast_horizon, 
													rsq_null_line = rsq_null_line,
													rmse_fcast_horizon = rmse_fcast_horizon, 
													rmse_null_line = rmse_null_line,
													crps_fcast_horizon = crps_fcast_horizon, 
													crps_null_line = crps_null_line)

out_df_list[[model_id]] = out_df

return(list(out_df,save_fig))
}

out_df_list = lapply(out_parallel_list, "[[", 1)
out_figure_list = lapply(out_parallel_list, "[[", 2)

fcast_horizon_results = data.table::rbindlist(out_df_list) %>% as.data.frame()



fcast_horizon_long <- fcast_horizon_results %>% 
	pivot_longer(cols = c("rsq_fcast_horizon", "rmse_fcast_horizon", "crps_fcast_horizon",
												# "rsq_r.squared", "rmse_r.squared", "crps_r.squared",
												# "rsq_adj.r.squared", "rmse_adj.r.squared", "crps_adj.r.squared",
												"rsq_null_line", "rmse_null_line", "crps_null_line"),
							 #"rsq_p.value","rmse_p.value","crps_p.value"), 
							 names_to = "horizon_parameter") %>% 
	mutate(metric = lapply(horizon_parameter, function(x) strsplit(x, "_") %>% unlist())  %>%  
				 	lapply(., "[[", 1) %>% unlist()) %>% 
	mutate(parameter_type = lapply(horizon_parameter, function(x) strsplit(x, "_") %>% unlist())  %>%  
				 	lapply(., "[[", 2) %>% unlist(),
				 parameter_type = ifelse(parameter_type=="fcast","horizon",parameter_type))


saveRDS(list(fcast_horizon_results, out_figure_list, fcast_horizon_long),
				here("data/summary/fcast_horizon_df_core.rds"))


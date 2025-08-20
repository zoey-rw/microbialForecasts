# Calculate forecast horizons

# Visualize the forecast horizon
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
# install.packages('egg', dependencies = TRUE)
library(egg)
library(ggpmisc) # for polynomial plotting function

# Check if Parquet file exists, otherwise use RDS
parquet_file <- here("data/summary/parquet/all_hindcasts.parquet")
rds_file <- here("data/summary/all_hindcasts.rds")

if (file.exists(parquet_file)) {
  cat("Using Parquet file for memory efficiency...\n")
  hindcast_in <- arrow::read_parquet(parquet_file)
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  hindcast_in <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
# converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))
# 
# converged = converged_strict
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
fcast_horizon_df = hindcast_data %>%
	#filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period=="hindcast" & new_site==FALSE) %>%
	merge(last_obs) %>%
	group_by(model_id,model_name) %>% 
	mutate(months_since_obs = timepoint - last_obs)


# Calculate average scores by model
fcast_horizon_model_mean = fcast_horizon_df %>%
	group_by(taxon,model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(mean_crps = mean(crps, na.rm=T),
						RMSE = rmse(actual = truth, predicted = mean),
						RSQ.1 = 1 - (RMSE^2)/var(truth),
						MAPE = mape(actual = truth, predicted = mean),
						RSQ = summary(lm(truth ~ mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
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
site_mean = hindcast_data %>% group_by(taxon,model_id,model_name,siteID) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	summarize(site_mean = mean(truth, na.rm=T),
						site_sd = sd(truth, na.rm=T)) %>% ungroup
overall_mean = hindcast_data %>% group_by(taxon,model_id) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	summarize(overall_mean = mean(truth, na.rm=T),
						overall_sd = sd(truth, na.rm=T)) %>% ungroup
historical_mean <- merge(site_mean, overall_mean)

# Merge historical means into main df, and use to calculate null scores
fcast_horizon_df2 <- merge(fcast_horizon_df, historical_mean, by=c("taxon","model_id","siteID","model_name"), all.x=T)


fcast_horizon_null_site = fcast_horizon_df2 %>%
	filter(site_sd != 0) %>% 
	group_by(taxon,model_name,pretty_group, rank_only, model_id) %>%
	summarize(null_CRPS = mean(crps_norm(truth, site_mean,
																			 site_sd)),
						null_CRPS_truncated = mean(crps(truth, family = "tnorm",
																						location = site_mean, scale = site_sd,
																						lower = 0, upper = 1)),
						null_RMSE = rmse(actual = truth, predicted = site_mean),
						null_RSQ.1 = 1 - (null_RMSE^2)/var(truth),
						null_MAPE = mape(actual = truth, predicted = site_mean),
						null_RSQ = summary(lm(truth ~ site_mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = null_RMSE/abundance) %>% distinct()
fcast_horizon_null_site$null_RSQ.1 = ifelse(fcast_horizon_null_site$null_RSQ.1 < 0, 0, fcast_horizon_null_site$null_RSQ.1)
fcast_horizon_null_site$months_since_obs = Inf


last_obs_values = inner_join(hindcast_data, last_obs %>% mutate(timepoint=last_obs))


last_obs_score <- last_obs_values %>%
	mutate(months_since_obs = 0) %>% 
	group_by(taxon,model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(CRPS = mean(crps_norm(truth, mean,
																			 sd)),
						CRPS_truncated = mean(crps(truth, family = "tnorm",
																						location = mean, scale = sd,
																						lower = 0, upper = 1)),
						RMSE = rmse(actual = truth, predicted = mean),
						RSQ.1 = 1 - (RMSE^2)/var(truth),
						MAPE = mape(actual = truth, predicted = mean),
						RSQ = summary(lm(truth ~ mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
last_obs_score$RSQ.1 = ifelse(last_obs_score$RSQ.1 < 0, 0, last_obs_score$RSQ.1)
last_obs_score$RMSE.norm = ifelse(last_obs_score$RMSE.norm > 5, 5, last_obs_score$RMSE.norm)


# Combine for plotting!
colnames(fcast_horizon_null_site) = gsub("null_","", colnames(fcast_horizon_null_site))
comparison_scores <- rbindlist(list(fcast_horizon_null_site, last_obs_score), fill=T)
to_plot <- rbindlist(list(comparison_scores, fcast_horizon_model_mean %>% filter(n_obs > 10)), fill=T)
to_plot_null <-  fcast_horizon_null_site %>% filter(model_id %in%  to_plot$model_id)


saveRDS(list(to_plot, to_plot, fcast_horizon_null_site, fcast_horizon_model_mean),
				here("data/summary/fcast_horizon_input.rds"))

in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
to_plot <- in_list[[1]]
to_plot_null <- in_list[[2]]
fcast_horizon_null_site <-  in_list[[3]]
fcast_horizon_model_mean <-  in_list[[4]]

model_id_list <- unique(to_plot$model_id)
model_id=model_id_list[[1]]

out_df_list = list()
out_figure_list = list()



#Run for multiple ranks, in parallel
cl <- makeCluster(28, type="FORK", 
									outfile="")
registerDoParallel(cl)

# for testing
out_parallel_list = foreach(model_id = model_id_list, 
														.errorhandling = 'remove') %dopar% {
															source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
															
#for(model_id in model_id_list){
	print(model_id)
	
single_tax = to_plot %>% filter(model_id == !!model_id)
single_tax$mean_crps = ifelse(is.na(single_tax$mean_crps), single_tax$CRPS_truncated, single_tax$mean_crps)

single_tax_null <-  fcast_horizon_null_site %>% filter(model_id == !!model_id)
single_tax_null$mean_crps = single_tax_null$CRPS_truncated
single_tax$months_since_obs = single_tax$months_since_obs + .001

fig_rsq <- ggplot(single_tax,
									aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = RSQ,
																				color = pretty_group, 
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)
		#  stat_poly_line(formula = y ~x, se = F, method = 'glm',	
		#  							 method.args = list(family = gaussian(link = 'log'))) +
		# stat_poly_eq(formula = y ~x, use_label("R2")) + guides(color="none")


fig_crps <- ggplot(single_tax,
									aes(x = months_since_obs, y = mean_crps, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (CRPS)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = mean_crps,
																				color = pretty_group, 
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)
# +
# 	 stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
# 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none")


fig_rmse <- ggplot(single_tax,
									 aes(x = months_since_obs, y = RMSE.norm, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (RMSE)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =single_tax_null, aes(yintercept = RMSE.norm,
																				color = pretty_group, 
																				group=model_name), na.rm = T, alpha=.8, linetype=2) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')), fullrange=TRUE)
# +
# 	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
# 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none")

save_fig = ggarrange(fig_rsq + rremove("xlab"), 
										 fig_crps + rremove("xlab"), 
										 fig_rmse, labels = c("A","B","C"), top = model_id)
out_figure_list[[model_id]] <- save_fig

# Extract information about plot
horizon_rsq_data = ggplot_build(fig_rsq)$data
rsq_loess_line = horizon_rsq_data[[3]]
rsq_null_line = horizon_rsq_data[[2]]$yintercept

# Get intersection: first point at which the poly fit is lower than the null Rsq
# If it comes super close, it should count 
rsq_fcast_horizon <- rsq_loess_line$x[which(rsq_loess_line$y < rsq_null_line)] %>% min(na.rm=T)
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# rsq_fit_info = horizon_rsq_data[[4]][,c("adj.r.squared","r.squared","p.value")] %>% 
# 	setNames(., paste0("rsq_",names(.)))


# Extract information about plot
horizon_rmse_data = ggplot_build(fig_rmse)$data
rmse_loess_line = horizon_rmse_data[[3]]
rmse_null_line = horizon_rmse_data[[2]]$yintercept

# Get intersection: first point at which the poly fit is lower than the null rmse
# If it comes super close, it should count 
rmse_fcast_horizon <- rmse_loess_line$x[which(rmse_loess_line$y > rmse_null_line)] %>% min(na.rm=T)
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# rmse_fit_info = horizon_rmse_data[[4]][,c("adj.r.squared","r.squared","p.value")]  %>% 
# 	setNames(., paste0("rmse_",names(.)))



# Extract information about plot
horizon_crps_data = ggplot_build(fig_crps)$data
crps_loess_line = horizon_crps_data[[3]]
crps_null_line = horizon_crps_data[[2]]$yintercept

# Get intersection: first point at which the poly fit is lower than the null crps
# If it comes super close, it should count 
crps_fcast_horizon <- crps_loess_line$x[which(crps_loess_line$y > crps_null_line)] %>% min(na.rm=T)
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# crps_fit_info = horizon_crps_data[[4]][,c("adj.r.squared","r.squared","p.value")] %>% 
# 	setNames(., paste0("crps_",names(.)))


out_df = cbind.data.frame(single_tax[1,1:5], 
													rsq_fcast_horizon = rsq_fcast_horizon, 
													rsq_null_line = rsq_null_line,
													#rsq_fit_info,
													
													rmse_fcast_horizon = rmse_fcast_horizon, 
													rmse_null_line = rmse_null_line,
													#rmse_fit_info,
													
													crps_fcast_horizon = crps_fcast_horizon, 
													crps_null_line = crps_null_line)
													#crps_fit_info)
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
				here("data/summary/fcast_horizon_df.rds"))



# Confirm that all Inf results are from a forecast horizon beyond 12 months:
inf_horizon = fcast_horizon_results[fcast_horizon_results$rsq_fcast_horizon==Inf,]
ggplot(to_plot %>% filter(model_id %in% inf_horizon$model_id),
			 aes(x = months_since_obs, y = RSQ, color = pretty_group, group=model_id)) +
	
	facet_wrap( ~model_id, 
						 labeller = labeller(model_name = model.labs), scales="free") +
	geom_point(show.legend = F) +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =to_plot_null %>% filter(model_id %in% inf_horizon$model_id & is.infinite(months_since_obs)), 
						 aes(yintercept = RSQ,
																				color = pretty_group, 
																				group=model_name), na.rm = T, alpha=.8, show.legend = F) + 
	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')))
# +
# 	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE),se=F, show.legend = F) +
# 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2"), show.legend = F)

# Change infinite horizon to max (12 months)
# fcast_horizon_results <- fcast_horizon_results %>% mutate(horizon = ifelse(horizon==Inf, 12, horizon))

# View overall horizon results!
ggplot(fcast_horizon_results,
			 aes(x = rank_only, y = rsq_fcast_horizon, color = pretty_group)) +
	facet_grid(pretty_group ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") + 
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitter(), alpha=.3) 


ggplot(fcast_horizon_long %>% filter(parameter_type=="null"),
			 aes(x = rank_only, y = value, color = pretty_group)) +
	facet_grid(metric ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") + 
	#coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitter(), alpha=.3) 



to_plot



ggplot(to_plot %>%  filter(rank_only=="functional") %>%
			 	filter(model_id %in% converged) ,
			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	facet_grid(pretty_group ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") + 
	#coord_flip() +
#	geom_boxplot() +
	geom_point(position=position_jitter(), alpha=.3) +
	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none") +
	geom_hline(data = to_plot_null, aes(yintercept = RSQ.1), alpha=.1)



ggplot(to_plot %>%  filter(rank_only=="functional") %>%
			 	filter(model_id %in% converged) ,
			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	facet_grid(pretty_group ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") +
	geom_hline(data = to_plot_null, aes(yintercept = RSQ.1), alpha=.1)

	to_plot_null

to_plot_null


to_plot_null = to_plot_null %>% group_by(pretty_group,model_name) %>% mutate(median_null = median(RSQ)) %>% 
	filter(rank_only %in% c("phylum","functional")) %>% filter(model_name %in% c("env_cycl"))
p <- ggplot(to_plot %>% filter(rank_only %in% c("phylum","functional")) %>% filter(model_name %in% c("env_cycl")),
			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	
	facet_grid(rows=vars(pretty_group), 
						 labeller = labeller(model_name = model.labs), scales="free") +
	geom_point(position=position_jitter(), alpha=.1, show.legend = F) +
	geom_smooth(span=2,  method="loess", show.legend = F, se=F) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .15, label.y.npc = .9, show.legend = F
	) +
	theme_bw(base_size = 20) +
	ylab("Forecast accuracy (RSQ)") +
	xlab("Months since last observation") +
	geom_hline(data = to_plot_null, aes(yintercept = RSQ , group=model_name), alpha=.05) +
	geom_hline(data = to_plot_null, aes(yintercept = median_null, group=model_name),linetype=2, alpha=.5) 

arrows <- 
	tibble(
		x1 = c(8, 12),
		x2 = c(7.2, 10),
		y1 = c(.5, .3),
		y2 = c(.4, .12), 
		pretty_group = c("Bacteria", "Fungi")
	)

dat_text <- data.frame(
	pretty_group = c("Bacteria", "Fungi","Bacteria", "Fungi"),
	months_since_obs     = c(9 , 11.5, 2.2, 2.2),
	RSQ     = c(.55,.42, .44, .18),
	label = c("Forecast horizon: \n7 months","Forecast horizon: \n10 months",
						"Null predictability","Null predictability")
)
p + geom_curve(
	data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
	arrow = arrow(length = unit(0.08, "inch")), size = 0.8, curvature = -.3, color = "gray20") + 
	geom_text(data = dat_text, aes(label=label), color=1, size=5, show.legend = F) + guides(pretty_group=NULL)


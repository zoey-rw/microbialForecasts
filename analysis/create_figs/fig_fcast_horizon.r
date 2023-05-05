# Visualize the forecast horizon
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
# install.packages('egg', dependencies = TRUE)
library(egg)

hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds"))

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

converged = converged_strict

# Determine each plot's final observation in calibration
last_obs = hindcast_data %>%
	filter(model_id %in% converged) %>%
	group_by(model_id, plotID) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	filter(timepoint == max(timepoint, na.rm=T)) %>%
	mutate(last_obs = timepoint) %>%
	select(plotID, last_obs) %>% distinct()

hindcast_data[hindcast_data$truth==0,]$truth <- .0001

# Calculate "time since final observation" for forecasts
fcast_horizon_df = hindcast_data %>%
	filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period=="hindcast" & new_site==FALSE) %>%
	merge(last_obs) %>%
	group_by(model_id,model_name) %>% 
	mutate(months_since_obs = timepoint - last_obs)


site_mean = hindcast_data %>% group_by(model_id,model_name,siteID) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	summarize(site_mean = mean(truth, na.rm=T),
				 site_sd = sd(truth, na.rm=T)) %>% ungroup

overall_mean = hindcast_data %>% group_by(model_id) %>%
	filter(fcast_period != "hindcast" & !is.na(truth)) %>%
	summarize(overall_mean = mean(truth, na.rm=T),
						overall_sd = sd(truth, na.rm=T)) %>% ungroup
historical_mean <- merge(site_mean, overall_mean)

fcast_horizon_df2 <- merge(fcast_horizon_df, historical_mean, by=c("model_id","siteID","model_name"), all.x=T)

# Average by model
fcast_horizon_null = fcast_horizon_df2 %>%
	group_by(model_name,pretty_group, rank_only, model_id) %>%
	summarize(null_CRPS = mean(crps_norm(truth, overall_mean,
																	overall_sd)),
						null_CRPS_truncated = mean(crps(truth, family = "tnorm",
																			 location = overall_mean, scale = overall_sd,
																																						 lower = 0, upper = 1)),
						null_RMSE = rmse(actual = truth, predicted = overall_mean),
						null_RSQ.1 = 1 - (null_RMSE^2)/var(truth),
						null_MAPE = mape(actual = truth, predicted = overall_mean),
						null_RSQ = summary(lm(truth ~ overall_mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = null_RMSE/abundance) %>% distinct()
fcast_horizon_null$null_RSQ.1 = ifelse(fcast_horizon_null$null_RSQ.1 < 0, 0, fcast_horizon_null$null_RSQ.1)


# Average by model
fcast_horizon_model_mean = fcast_horizon_df %>%
	group_by(model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(mean_crps = mean(crps, na.rm=T),
																									RMSE = rmse(actual = truth, predicted = mean),
																									RSQ.1 = 1 - (RMSE^2)/var(truth),
																									MAPE = mape(actual = truth, predicted = mean),
																									RSQ = summary(lm(truth ~ mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
fcast_horizon_model_mean$RSQ.1 = ifelse(fcast_horizon_model_mean$RSQ.1 < 0, 0, fcast_horizon_model_mean$RSQ.1)
fcast_horizon_model_mean$RMSE.norm = ifelse(fcast_horizon_model_mean$RMSE.norm > 5, 5, fcast_horizon_model_mean$RMSE.norm)

fcast_n_obs = fcast_horizon_df %>%
	group_by(model_id, months_since_obs) %>%
	tally(name = "n_obs")
# Remove metrics calculated from only 1-2 observations
fcast_horizon_model_mean <- merge(fcast_horizon_model_mean, fcast_n_obs) %>% filter(n_obs > 2)


# Average by model and site
fcast_horizon_site_mean = fcast_horizon_df %>%
	group_by(model_name,pretty_group, rank_only, model_id, months_since_obs,siteID) %>%
	summarize(mean_crps = mean(crps, na.rm=T),
																														 RMSE = rmse(actual = truth, predicted = mean),
																														 RSQ.1 = 1 - (RMSE^2)/var(truth),
																														 MAPE = mape(actual = truth, predicted = mean),
																														 RSQ = summary(lm(truth ~ mean))$r.squared) %>% distinct()
fcast_horizon_site_mean$RSQ.1 = ifelse(fcast_horizon_site_mean$RSQ.1 < 0, 0, fcast_horizon_site_mean$RSQ.1)


# View mean values per group, RSQ.1
b = ggplot(fcast_horizon_model_mean %>% filter(months_since_obs < 14), # %>% filter(rank_only=="phylum"),
			 aes(x = months_since_obs, y = RSQ.1, color = pretty_group)) +
	
	facet_grid(pretty_group ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") +
	geom_smooth(span=2,  method="loess", show.legend = F, se=F) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .15, label.y.npc = .15
	) +
	theme_minimal(base_size = 20) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	#	geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
	geom_hline(data =fcast_horizon_null, aes(yintercept = mean(null_RSQ.1)), na.rm = T) 
b

b <- tag_facet(b, size=7)

png(here("figures","forecast_horizon.png"), width = 1400, height=1000)
print(b)

dev.off()



# Plot CRPS values for all observations
to_plot <- fcast_horizon_df %>% filter(months_since_obs < 16) #%>% filter(rank_only=="genus")
a <- ggplot(to_plot,
						aes(x = months_since_obs, y = crps, color = pretty_group)) +
#	geom_point(alpha=.3, position=position_jitter(height=0, width=.3), size=3) +
	facet_grid(rows=vars(model_name),
						 labeller = labeller(model_name = model.labs)) +
	geom_smooth(method="loess", span=1, se=F) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .65, label.y.npc = .75
	) +
	theme_minimal(base_size = 18) +
	ylab("Forecast error (CRPS)") +
	xlab("Months since last observation") +
	labs(color="Kingdom") + ylim(c(0,.25))  +
	scale_y_log10()
a

# View site-means (since the number of observations varies by site)
ggplot(fcast_horizon_site_mean %>% filter(months_since_obs < 16# & rank_only=="genus"
																					),
			 aes(x = months_since_obs, y = mean_crps, color = pretty_group)) +
	# geom_point(alpha=.003, position=position_jitter(height=0, width=.1), size=3) +
	facet_grid(cols=vars(pretty_group),
						 rows=vars(model_name),
						 labeller = labeller(model_name = model.labs), scales = "free") +
	geom_smooth(span=.9,  method="loess") +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .15#, label.y = .25
		) +
	theme_minimal(base_size = 18) +
	ylab("Forecast error") +
	xlab("Months since last observation") +
	labs(color="Kingdom") +
	#geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
	scale_y_log10()

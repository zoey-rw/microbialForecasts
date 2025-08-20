# Visualize the forecast horizon
source("source.R")


fcast_horizon_results <- readRDS("data/summary/fcast_horizon_df_core.rds")



weak_converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
strict_converged <- readRDS(here("data/summary/converged_taxa_list.rds"))
converged = weak_converged
converged = strict_converged

# converged_2018 <- converged[grepl("2018",converged)]
# converged_all_mods =
# 	lapply(sum.in$keep_models_weak$taxon, function(x) {
# 	if (sum(grepl(x, converged_2018)) == 3) return(x)
# 		}) %>% unlist %>% unique


fcast_horizon_df = fcast_horizon_results[[1]]
fcast_horizon_plotting = fcast_horizon_results[[2]]
fcast_horizon_long = fcast_horizon_results[[3]]
fcast_horizon_long = fcast_horizon_results[[3]]

fcast_horizon_df <- fcast_horizon_df %>% filter(model_id %in% converged) 
fcast_horizon_long <- fcast_horizon_long %>% filter(model_id %in% converged)
#fcast_horizon_df <- fcast_horizon_df[fcast_horizon_df$model_id %in% converged_strict,]
#


# Remove infinite horizons 
# fcast_horizon_long <- fcast_horizon_long %>%
# 	filter(parameter_type == "horizon" &
# 				 	value !=Inf) 

# fcast_horizon_df <- fcast_horizon_df %>%
# 	filter(!is.infinite(crps_fcast_horizon) & 
# 				 	!is.infinite(rmse_fcast_horizon) & 
# 				 	!is.infinite(crps_fcast_horizon)) 


# # Set infinite horizons to the max value (20)
# fcast_horizon_long <- fcast_horizon_long %>%
# 	filter(model_id %in% converged) %>%
# 	mutate(value = ifelse(parameter_type == "horizon" &
# 																															 	value ==Inf,
# 																															 20, value))
# 
# fcast_horizon_df <- fcast_horizon_df %>%
# 	filter(model_id %in% converged) %>%
# 	mutate(crps_fcast_horizon = ifelse(crps_fcast_horizon ==Inf,
# 												20, crps_fcast_horizon),
# 				 rmse_fcast_horizon = ifelse(rmse_fcast_horizon ==Inf,
# 				 														20, rmse_fcast_horizon),
# 				 rsq_fcast_horizon = ifelse(rsq_fcast_horizon ==Inf,
# 				 														20, rsq_fcast_horizon))
ggplot(to_plot,
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


fcast_horizon_df$seasonal_predictors = as.factor(ifelse(fcast_horizon_df$model_name %in% c("env_cycl","cycl_only"), 1, 0))
fcast_horizon_df$environmental_predictors = as.factor(ifelse(fcast_horizon_df$model_name %in% c("env_cycl","env_cov"), 1, 0))
fcast_horizon_long$seasonal_predictors = as.factor(ifelse(fcast_horizon_long$model_name %in% c("env_cycl","cycl_only"), 1, 0))
fcast_horizon_long$environmental_predictors = as.factor(ifelse(fcast_horizon_long$model_name %in% c("env_cycl","env_cov"), 1, 0))

mean_fcast_horizon = fcast_horizon_df %>% 
	group_by(model_id, model_name, pretty_group, rank_only, seasonal_predictors, environmental_predictors) %>% 
	summarize(mean_horizon=mean(crps_fcast_horizon, rmse_fcast_horizon, 
															#rsq_fcast_horizon, 
															na.rm=T)) %>% filter(!is.infinite(mean_horizon))
ggplot(mean_fcast_horizon,# %>% filter(rank_only %in% c("functional","phylum")),
			 aes(x = model_name, y = mean_horizon, color = model_name)) +
	facet_grid(~pretty_group, 
						 labeller = labeller(model_name = model.labs), scales="free") + 
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw()

mean_fcast_horizon %>% 
	group_by(pretty_group) %>%
	summarize(tukey(x = model_name, y = mean_horizon)) %>%
	rename(model_name = x)


for_stats = fcast_horizon_long %>% 
	filter(parameter_type=="horizon" & metric=="rsq") %>% 
	filter(!is.na(value) & !is.infinite(value))
for_stats %>% 
	group_by(pretty_group) %>%
	summarize(tukey(x = model_name, y = value)) %>%
	rename(model_name = x)


model_stat_pvalue <- for_stats %>% 
	#group_by(pretty_group) %>% 
	rstatix::tukey_hsd(value ~ model_name) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position(step.increase = .4) %>% 
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

# View horizon results by model type
horizon_by_model_type <- ggplot(fcast_horizon_long %>% filter(parameter_type=="horizon" & metric=="rsq")  %>% 
			 	filter(!is.na(value) & !is.infinite(value)),
			 aes(x = model_name, y = value, color = model_name)) +
	# facet_grid(~pretty_group, 
	# 					 labeller = labeller(model_name = model.labs), scales="free") + 
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw() +
	ggpubr::stat_pvalue_manual(model_stat_pvalue, label = "p.adj.signif", #bracket.nudge.y = -.4, 
														 size=4, hide.ns = T) +
	coord_flip() 
horizon_by_model_type




kingdom_stat_pvalue <- for_stats %>% 
	#group_by(pretty_group) %>% 
	rstatix::tukey_hsd(value ~ pretty_group) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position(step.increase = .4) %>% 
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

# View horizon results by model type
horizon_by_kingdom <- ggplot(fcast_horizon_long %>% filter(parameter_type=="horizon" & metric=="rsq")  %>% 
																	filter(!is.na(value) & !is.infinite(value)),
																aes(x = pretty_group, y = value, color = pretty_group)) +
	# facet_grid(~pretty_group, 
	# 					 labeller = labeller(model_name = model.labs), scales="free") + 
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw() +
	ggpubr::stat_pvalue_manual(kingdom_stat_pvalue, label = "p.adj.signif", #bracket.nudge.y = -.4, 
														 size=4, hide.ns = T) +
	coord_flip() 
horizon_by_kingdom






fcast_horizon_long %>% 
	filter(rank_only %in% c("functional","phylum")) %>% 
	group_by(pretty_group, metric) %>% 
	rstatix::tukey_hsd(value ~ environmental_predictors)


model_parameters_lm <-lm(rsq_fcast_horizon ~ 
												 	model_name*pretty_group,
												 data = fcast_horizon_df)
sjPlot::plot_model(model_parameters_lm,  main = "Effects on predictability (RSQ)", intercept = F)  



stat_pvalue_env <- fcast_horizon_df %>% 
	filter(rank_only %in% c("functional","phylum","genus")) %>% 
	group_by(pretty_group) %>% 
	rstatix::tukey_hsd(rsq_fcast_horizon ~ environmental_predictors) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position(step.increase = .4) %>% 
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))


stat_pvalue <- fcast_horizon_df %>% 
	filter(rank_only %in% c("functional","phylum","genus")) %>% 
	group_by(pretty_group) %>% 
	rstatix::tukey_hsd(rsq_fcast_horizon ~ rank_only) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position(step.increase = .4) %>% 
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))



ggplot(fcast_horizon_df %>% filter(rank_only %in% c("functional","phylum","genus")),
			 aes(x = rank_only, y = rsq_fcast_horizon, color = model_name)) +
	facet_grid(~pretty_group, 
						 labeller = labeller(model_name = model.labs), scales="free") + 
	geom_boxplot() +
	ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.adj.signif", bracket.nudge.y = -.4, size=4, hide.ns = T) +
coord_flip() +
	theme_minimal(base_size = 20) + 
	geom_point(position=position_jitterdodge(), alpha=.3) + 
	scale_y_continuous(n.breaks=10)


# View horizon results by model type
ggplot(fcast_horizon_long %>% filter(parameter_type=="horizon") %>% 
			 filter(rank_only %in% c("phylum")),
			 aes(x = model_name, y = value, color = pretty_group)) +
	facet_grid(metric~pretty_group,  scales="free") + 
	coord_flip() +
	geom_boxplot() +
	geom_point(position=position_jitter(), alpha=.3) 




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



fcast_horizon_df2 <- merge(fcast_horizon_df, historical_mean, by=c("model_id","siteID","model_name"), all.x=T)

fcast_horizon_null_site = fcast_horizon_df2 %>%
	filter(site_sd != 0) %>% 
	group_by(model_name,pretty_group, rank_only, model_id) %>%
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
fcast_horizon_null_site$null_RSQ.1 = ifelse(fcast_horizon_null$null_RSQ.1 < 0, 0, fcast_horizon_null$null_RSQ.1)



# Plot CRPS values for all observations

# View mean values per group, RSQ, with fcast horizon
c = ggplot(to_plot, # %>% filter(rank_only=="phylum"),
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
	geom_hline(data =fcast_horizon_null_site, aes(yintercept = mean(null_RSQ.1)), na.rm = T) 
c


test_HARV <- fcast_horizon_df2 %>% filter(siteID=="HARV",model_id=="env_cycl_russulaceae_20151101_20180101")
test_mod <- fcast_horizon_df2 %>% filter(model_id=="env_cycl_russulaceae_20151101_20180101")
test_out_HARV <- fcast_horizon_null_site %>% filter(siteID=="HARV",model_id=="env_cycl_russulaceae_20151101_20180101")


last_obs_values = inner_join(hindcast_data, last_obs %>% mutate(timepoint=last_obs))


last_obs_score <- last_obs_values %>%
	mutate(months_since_obs = 0) %>% 
	group_by(model_name,pretty_group, rank_only, model_id, months_since_obs) %>%
	summarize(mean_crps = mean(crps, na.rm=T),
						RMSE = rmse(actual = truth, predicted = mean),
						RSQ.1 = 1 - (RMSE^2)/var(truth),
						MAPE = mape(actual = truth, predicted = mean),
						RSQ = summary(lm(truth ~ mean))$r.squared,
						abundance = mean(truth, na.rm = T), RMSE.norm = RMSE/abundance) %>% distinct()
last_obs_score$RSQ.1 = ifelse(last_obs_score$RSQ.1 < 0, 0, last_obs_score$RSQ.1)
last_obs_score$RMSE.norm = ifelse(last_obs_score$RMSE.norm > 5, 5, last_obs_score$RMSE.norm)


to_plot <- merge(fcast_horizon_null_site, fcast_horizon_model_mean, all=T) %>% filter(months_since_obs < 16)
to_plot <- rbindlist(list(to_plot, last_obs_score), fill=T)

mean_null
ggplot(to_plot %>% filter(rank_only=="phylum"),
			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	
	facet_grid(pretty_group ~model_name, 
						 labeller = labeller(model_name = model.labs), scales="free") +
	geom_point() +
	geom_smooth(span=2,  method="loess", show.legend = F, se=F) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .15, label.y.npc = .15
	) +
	theme_minimal(base_size = 20) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(aes(yintercept = mean(null_RSQ), group=model_name))
	#	geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
	# geom_hline(data =fcast_horizon_null_site, aes(yintercept = mean(null_RSQ),
	# 																							 color = pretty_group, group=model_name), na.rm = T) 

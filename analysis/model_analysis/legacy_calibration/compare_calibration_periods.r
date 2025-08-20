# Compare hindcast accuracy/bias depending on calibration period
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

pacman::p_load(scoringRules, data.table) 

recent_hindcast <- readRDS("./data/summary/hindcast_div.rds")
old_hindcast <- readRDS("./data/summary/hindcast_div_legacy_incl.rds")
old_hindcast[old_hindcast$time_period=="calibration",]$time_period  <- "2013-06_2017-01"
hindcast_data <- rbindlist(list(recent_hindcast, 
																old_hindcast), fill = T)

# Add prettier data values 
hindcast_data$newsite <- ifelse(hindcast_data$new_site, "New site", "Observed site")
hindcast_data$pretty_group <- ifelse(hindcast_data$group=="16S", "Bacteria", "Fungi")
hindcast_data <- hindcast_data %>% 
	mutate(calibration_label = recode(as.character(time_period), !!!calibration_label),
				 truth = as.numeric(truth))

# Subset to hindcast observations
not_na_hindcast <- hindcast_data %>% 
	filter(!is.na(hindcast_data$truth) & fcast_period == "hindcast") %>% 
	mutate(in_95 = ifelse(truth > lo & truth < hi, 1, 0))
not_na_calibration <- hindcast_data %>% 
	filter(!is.na(hindcast_data$truth) & fcast_period == "calibration") %>% 
	mutate(in_95 = ifelse(truth > `2.5%` & truth < `97.5%`, 1, 0))
# Score using continuous ranked probability score with estimated mean and SD
scored_hindcasts <- not_na_hindcast %>% 
	mutate(crps = crps_norm(truth, mean, sd))


# Compare CRPS scores for all_cov model
scored_hindcasts_mean <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, model_name, newsite, calibration_label) %>% 
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))

ggplot(scored_hindcasts %>% filter(newsite == "Observed site" & 
																	 	model_name == "all_covariates"), 
			 aes(x = pretty_group, y = crps, 
														 color = calibration_label)) + 
	coord_trans(y = "log10") +
	geom_violin(aes(color = calibration_label), 
							draw_quantiles = c(0.5), show.legend=F) + 
	geom_jitter(alpha = .1, size=4, 
							position=position_jitterdodge(dodge.width = 1)) + 
	#facet_wrap(~pretty_group, drop = T) +  
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + scale_color_manual(values = c(1,2))


# Skill score from Scavia 2021
skill_score_mean <- scored_hindcasts_mean %>% 
	pivot_wider(id_cols = c("fcast_type","pretty_group","model_name","calibration_label"), 
							values_from = "crps_mean", names_from = "newsite") %>% 
	mutate(skill_score = (1 - (`New site`/`Observed site`)))
ggplot(skill_score_mean %>% filter(model_name == "all_covariates"), 
			 aes(x = pretty_group, y = skill_score, 
			 		color = calibration_label)) + 
	geom_jitter(alpha = 1, size=4, 
							position=position_jitterdodge(dodge.width = 1)) + 
	#facet_wrap(~pretty_group, drop = T) +  
	ylab("Skill score (% change in CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + scale_color_manual(values = c(1,2))





# Calculate overdispersion for calibration & hindcasts (by site)
observed <- rbind(not_na_hindcast,not_na_calibration) %>% 
	filter(newsite == "Observed site" & model_name == "all_covariates")
overdispersion <-  observed %>%
	group_by(fcast_period, fcast_type, pretty_group, model_name, newsite, calibration_label, siteID) %>% 
	mutate(countT= n(),
				in_95_count = sum(in_95, na.rm=T),
				in_95_pct = in_95_count/countT) %>% 
	select(fcast_period, fcast_type, pretty_group, model_name, newsite, time_period,
				 in_95_pct, calibration_label, siteID) %>%
	distinct() %>% 
	mutate(outlier = ifelse(in_95_pct < .5, siteID, NA)) 
overdispersion$outside_confidence_interval <- 1-overdispersion$in_95_pct

# Plot overdispersion for calibration & hindcast periods by site
ggplot(overdispersion, aes(x = pretty_group, y = outside_confidence_interval, 
													 color = calibration_label)) + 
	geom_violin(draw_quantiles = c(0.5), show.legend=F) + 
	geom_point(aes(color = calibration_label), 
						 alpha = .3, size=4,
						 position=position_jitterdodge(dodge.width = 1, jitter.width = .2)) + 
	facet_grid(rows=vars(fcast_period), drop=T, scales = "free") +  
	ylab("% observations outside 95% CI") + xlab(NULL) + 
	theme_minimal(base_size=20) + ggtitle("Overdispersion in model estimates") + 
	theme(text = element_text(size = 26),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1)) + 
	guides(color=guide_legend(title=NULL)) + 
	geom_hline(yintercept = .05, linetype=2) + 
	annotate("text", x = 1.7, y = .025, label = "Target value (5%)", size = 7) + scale_color_manual(values = c(1,2))
	# geom_text(aes(label = outlier), na.rm = TRUE, 
	# 					position=position_jitterdodge(dodge.width = 1)) + scale_y_sqrt()


# View actual models & hindcasts for a single site
ggplot(hindcast_data %>% filter(#plotID=="HARV_001" & 
	group=="ITS" & model_name == "all_covariates" & siteID == "CPER")) + 
	facet_grid(rows=vars(plotID), 
						 cols = vars(calibration_label), drop=T, scales="free") +
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
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ggtitle("Fungal diversity hindcasts at NEON site: CPER")

ggplot(hindcast_data %>% filter(
	group=="16S" & model_name == "all_covariates" & siteID == "CPER")) + 
	facet_grid(rows=vars(plotID), 
						 cols = vars(calibration_label), drop=T, scales="free") +
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
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ggtitle("Bacterial diversity hindcasts at NEON site: WOOD")


p_load(ganttrify)

data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")

model_date_summaries <- data_in$plot_est %>% 
	group_by(time_period) %>% 
	filter(!is.na(truth)) %>% 
	summarise(start_date = min(dates, na.rm = T),
						end_date = max(dates, na.rm = T)) %>% 
	mutate(wp = recode(as.character(time_period), !!!calibration_label),
				 activity = recode(as.character(time_period), !!!calibration_label))

ganttrify(project = model_date_summaries,
month_number_label = F,
by_date = T,
month_breaks = 6,
					font_family = "Roboto Condensed")


p_load(Metrics, forestmangr)

scoring_metrics <- not_na_hindcast %>% filter(newsite=="Observed site") %>% 
	group_by(calibration_label, pretty_group,model_name, siteID) %>% 
	mutate(crps_point = crps_norm(truth, mean, sd)) %>% 
	summarise(
		RMSE = rmse(actual = truth, predicted = mean),
		BIAS = bias(actual = truth, predicted = mean),
		MAE = mae(actual = truth, predicted = mean),
		#CRPS = crps_norm(truth, mean, sd))
		#BIAS = percent_bias(actual = truth, predicted = mean),
		CRPS = mean(crps_point))
scoring_metrics_long <- scoring_metrics %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS), names_to = "metric")

# df <-  not_na_hindcast %>% filter(newsite=="Observed site", 
# 																	group == "16S", 
# 																	model_name == "all_covariates", 
# 																	calibration_label=="Including legacy data") 

ggplot(scoring_metrics_long %>% filter(model_name == "all_covariates"), 
			 aes(x = calibration_label, y = value, 
			 		color = calibration_label)) + 
	geom_violin(draw_quantiles = c(0.5), show.legend=F) + 
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
	facet_grid(metric~pretty_group, drop = T, scales="free") +  
	# geom_jitter(aes(x = metric, y = value), width=.1, 
	# 						height = 0, alpha = .8, size=4) + 
	ylab("Metric scores") + xlab(NULL) + 
	theme_bw(base_size=18) + 
	ggtitle("Hindcast scores by calibration period")  + 
	scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) + 
	#guides(color=guide_legend(title=NULL)) +	
	geom_hline(yintercept = 0, linetype=2) 




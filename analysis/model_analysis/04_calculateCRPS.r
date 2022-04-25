# Calculate and visualize CRPS for diversity forecasts
pacman::p_load(scoringRules, reshape2, parallel, lubridate, nimble, coda) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Combine forecast types
div_hindcast_data <- readRDS("./data/summary/hindcast_div.rds") %>% mutate(fcast_type = "Diversity")
fg_hindcast_data <- readRDS("./data/summary/hindcast_fg.rds") %>% mutate(fcast_type = "Functional group")
tax_hindcast_data <- readRDS("./data/summary/hindcast_tax_test.rds") %>% mutate(fcast_type = "Taxonomic") # replace
hindcast_data <- rbindlist(list(div_hindcast_data, fg_hindcast_data, tax_hindcast_data), fill = T)

hindcast_data <- hindcast_data %>% filter(!is.na(hindcast_data$truth) & dates > "2016-12-31")
# hindcast_data <- hindcast_data %>% filter(dates > "2016-12-31")
# hindcast_data <- hindcast_data %>% filter(dates < "2016-12-31")

# Add prettier data values 
hindcast_data$newsite <- ifelse(hindcast_data$new_site, "New site", "Observed site")
# 
# cal <- hindcast_data %>% filter(dates < "2016-12-31")
# val <- hindcast_data %>% filter(dates > "2016-12-31")

harv <- fg_hindcast_data %>% filter(plotID=="HARV_004" & taxon == "copiotroph" & model_name == "cycl_only")

# Score using continuous ranked probability score with estimated mean and SD
scored_hindcasts <- hindcast_data %>% mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd))
scored_hindcasts_mean <- scored_hindcasts %>% group_by(fcast_type, group, model_name, newsite) %>% 
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))

# Skill score from Scavia 2021
skill_score <- scored_hindcasts_mean %>% 
	pivot_wider(id_cols = c("fcast_type","group","model_name"), values_from = "crps_mean", names_from = "newsite") %>% 
	mutate(skill_score = (1 - (`New site`/`Observed site`)))

# Save output before visualizing
saveRDS(scored_hindcasts, "./data/summary/CRPS_hindcasts.rds")





# COMPARE CRPS BY ERRORS-IN-VARIABLES UNCERTAINTY
# All values
ggplot(scored_hindcasts) + 
	geom_violin(aes(x = model_name, y = crps, color = group), 
							draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~newsite) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = model_name, y = crps, color = group), #width=.2, #height = 0, 
							alpha = .3, size=4, 
							position=position_jitterdodge(dodge.width = 1, jitter.width = .1)) + 
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

not_na <- scored_hindcasts[!is.na(scored_hindcasts$truth),]

# Compare means
ggplot(scored_hindcasts_mean, aes(x = model_name, y = crps_mean, color = group)) + 
#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~newsite) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = model_name, y = crps_mean), width=.2, height = 0, alpha = .3, size=4) + 
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)




# Compare means by plot
scored_hindcasts_plot <- scored_hindcasts %>% group_by(pretty_group, group, plotID, model_name, newsite) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))
scored_hindcasts_site <- scored_hindcasts %>% group_by(pretty_group, group, siteID, model_name, newsite) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))
scored_hindcasts_mean <- scored_hindcasts %>% group_by(pretty_group, group, uncert, model_name, newsite) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))

ggplot(scored_hindcasts_plot) + 
		geom_violin(aes(x = model_name, y = crps_mean, color = group), 
								draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~newsite) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = model_name, y = crps_mean, color = group), #width=.2, #height = 0, 
							alpha = .3, size=4, 
							position=position_jitterdodge(dodge.width = 1, jitter.width = .1)) + 
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)


# Do T-tests

# F vs. B difference between observed site forecasts
# plot means
oldsites_allcov <- scored_hindcasts_plot %>% filter(newsite=="Observed site" & model_name == "all_covariates")
oldsites_bac <- oldsites_allcov[oldsites_allcov$group=="16S",]
oldsites_fun <- oldsites_allcov[oldsites_allcov$group=="ITS",]
t.test(oldsites_bac$crps_mean, oldsites_fun$crps_mean)
# all values
oldsites_allcov <- scored_hindcasts %>% filter(newsite=="Observed site" & model_name == "all_covariates")
oldsites_bac <- oldsites_allcov[oldsites_allcov$group=="16S",]
oldsites_fun <- oldsites_allcov[oldsites_allcov$group=="ITS",]
t.test(oldsites_bac$crps, oldsites_fun$crps)



# F vs. B difference between new site forecasts
# plot means
newsites_allcov_plot <- scored_hindcasts_plot %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov_plot[newsites_allcov_plot$group=="16S",]
newsites_fun <- newsites_allcov_plot[newsites_allcov_plot$group=="ITS",]
t.test(newsites_bac$crps_mean, newsites_fun$crps_mean)
# site means
newsites_allcov_site <- scored_hindcasts_site %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov_site[newsites_allcov_site$group=="16S",]
newsites_fun <- newsites_allcov_site[newsites_allcov_site$group=="ITS",]
t.test(newsites_bac$crps_mean, newsites_fun$crps_mean)
# all values
newsites_allcov <- scored_hindcasts %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov %>% filter(group=="16S")
newsites_fun <- newsites_allcov %>% filter(group=="ITS")
t.test(newsites_bac$crps, newsites_fun$crps)



ggplot(skill_score, aes(x = pretty_group, y = skill_score, color = pretty_group)) + 
	ggtitle("Skill score: forecasting to new sites") +
	facet_wrap(~model_name) +  
	geom_point(size=4, show.legend = F) + 
	#ylim(c( -.25, 0.)) +
	ylab("Skill score (decrease in mean CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)






its_bart_fcast <- hindcast_data_oldsites[hindcast_data_oldsites$plotID %in% c("CPER_001","CPER_002","CPER_004"),] %>% 
	filter(scenario == "full_uncertainty_ITS")

plot_summaries <- summaries$plot_est %>% filter(time_period == "calibration" & scenario == "full_uncertainty_ITS")
its_bart_est <- plot_summaries %>% filter(plotID %in% c("CPER_001","CPER_002","CPER_004"))  %>% arrange(siteID,plotID,dateID) #,"BART_002","BART_004"),]



ggplot() +
	facet_grid(rows=vars(plotID), 
						 cols=vars(model_name), 
						 drop=T, scales="free", space="free") +
	geom_line(data = its_bart_fcast, aes(x = dates, y = `50%`), show.legend = F, na.rm = T) +
	geom_ribbon(data = its_bart_fcast, aes(x = dates, ymin = `2.5%`, ymax = `97.5%`), 
							alpha=0.7,  na.rm = T, show.legend = F)  + 
	geom_point(data = its_bart_fcast, aes(x = dates, y = as.numeric(truth)))  +
	scale_x_date(date_breaks = "6 month", 
							 labels=date_format("%b-%Y"),
							 limits = as.Date(c('2014-01-01','2021-01-01'))) + #ylim(c(-1.1,1.1)) + 
	geom_line(data = its_bart_est, aes(x = dates, y = `50%`), show.legend = F, color = 2, na.rm = T) +
	geom_ribbon(data = its_bart_est, aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill = 2,
							alpha=0.7,  na.rm = T, show.legend = F)  


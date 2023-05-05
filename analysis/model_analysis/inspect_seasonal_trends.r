# read in seasonal values
library(lubridate)
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")


seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
converged = scores_list$converged_list
converged_strict = scores_list$converged_strict_list
seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% converged)
seas_vals_wide <- seas_in[[6]] %>% filter(model_id %in% converged)

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>% 
	select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)
cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")


sum.in <- readRDS(here("data", "summary/logit_beta_regression_summaries.rds"))
plot_estimates = sum.in$plot_est %>% filter(model_name != "all_covariates")

plot_estimates$month = lubridate::month(plot_estimates$dates)
plot_estimates$year = lubridate::year(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cov_est = plot_estimates %>% filter(grepl("env_cov",model_name))


pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
pheno_categories_long = pheno_categories_in[[2]] %>% filter(ID %in% unique(cycl_only_est$siteID))


max_dates = seas_vals_short %>% mutate(dates = max_y_date) %>% filter(time_period == "2015-11_2018-01")
#max_dates$sampling_season =  apply(max_dates, 1, assign_pheno_category)
#max_dates_cycl_only$sampling_season = factor(max_dates_cycl_only$sampling_season, ordered = T, levels = c("dormancy","dormancy_greenup", "greenup","greenup_peak", "peak", "greendown_peak", "greendown","greendown_dormancy"))

max_cycl <- max_dates %>% 
	mutate(dates = ymd(paste0("2014-", max_dates$max_month, "-15")))# %>% 
	#mutate(site_cat = assign_pheno_date(dates))

# Subset model estimates to years 2016-2018, which have data for most sites
new_max_abun = plot_estimates %>% 
	filter(year %in% c("2016","2017")) %>% 
	filter(time_period == "2015-11_2018-01") %>% 
	group_by(model_id, siteID, year, month, dates) %>% 
	# Use median plot estimates; take site/date mean of these values
	summarize(mean_modeled_abun = mean(`50%`, na.rm=T)) %>%
	ungroup %>% 
	group_by(model_id, siteID, year) %>% 
	# Only retain highest value per site-year
	filter(mean_modeled_abun == max(mean_modeled_abun, na.rm=T))

# Approx 60k rows; takes a few minute to assign phenophases to each site/date
max_abun <- new_max_abun %>% 
	mutate(dates = ymd(paste0(year, "-", month, "-15"))) %>% 
	filter(siteID %in% pheno_categories_long$ID) %>% 
	mutate(site_cat = assign_pheno_site_date(siteID, dates))
# About 6 sites are missing all MODIS phenophase data; these "NA" values will be excluded
# About 7000 additional dates have some missing data 
max_abun$sampling_season = factor(max_abun$site_cat, ordered = T, levels = c("dormancy",#"dormancy_greenup", 
																																						 "greenup",#"greenup_peak", 
																																						 "peak", #"greendown_peak", 
																																						 "greendown"#,"dormancy_greendown"
))

# temp <- max_cycl
# temp$rowname = 1:nrow(temp)
# cat_list2 = temp %>%
# 	split( .$rowname) %>%
# 	map(assign_pheno_date, dates)
# site_cats2 = as.matrix(cat_list2)
# max_abun <- cbind.data.frame(new_max_abun, site_cats2)

site_cats2 = lapply(max_cycl$dates, assign_pheno_date) %>% as.matrix() %>% unlist()
max_cycl <- cbind.data.frame(max_cycl, site_cats2)
max_cycl$sampling_season = factor(max_cycl$site_cats2, ordered = T, levels = c("dormancy","dormancy_greenup", 
																																						 "greenup","greenup_peak", 
																																						 "peak", "greendown_peak", 
																																						 "greendown","dormancy_greendown"
))
max_cycl <- merge(max_cycl, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")], all.x=T)


# Get most frequently-assigned phenophase per group	
seasonality_mode = max_abun %>% 
	filter(!is.na(sampling_season)) %>% 
	group_by(model_id, sampling_season) %>% 
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>% 
	select(-n) #%>% 
#pivot_wider(names_from="sampling_season", values_from = "freq", values_fill = 0)
seasonality_mode2 = seasonality_mode %>% group_by(model_id) %>% filter(freq == max(freq))
seasonality_mode2 <- merge(seasonality_mode2, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")], all.x=T)

seasonality_mode3 = seasonality_mode %>% group_by(model_id) %>% filter(freq > .5)
seasonality_mode3 <- merge(seasonality_mode3, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")], all.x=T)

max_abun_to_plot <- merge(max_abun, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")])

seasonality_mode_to_plot <- merge(seasonality_mode, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")])

saveRDS(list(seasonality_mode2, max_abun), here("data/clean/group_peak_phenophases.rds"))


phenophase_in = readRDS(here("data/clean/group_peak_phenophases.rds"))
seasonality_mode2 = phenophase_in[[1]]
max_abun = phenophase_in[[2]]

library(ggrepel)
ggplot(max_abun %>% filter(fcast_type=="Functional"),
			 aes(color=pretty_group, x = sampling_season, y = amplitude)) +
	geom_point( alpha=.2, 
							position=position_jitter(height=0, width=.2)) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(model_name~pretty_group, scales="free") +
	geom_label_repel(aes(label=taxon), size=3)

ggplot(max_abun_to_plot %>% filter(fcast_type=="Functional" & model_name=="cycl_only"),
			 aes(fill=pretty_group, color = pretty_group, 
			 		x = sampling_season, group=model_id)) +
	geom_histogram(stat="count") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~pretty_group, scales="free") 

ggplot(max_abun_to_plot  %>% filter(fcast_type=="Functional" & 
																			model_name=="cycl_only"), aes(sampling_season, colour = pretty_group, group=model_id)) +
	geom_freqpoly(stat="count", show.legend = F, alpha=.5) +
	theme_bw(base_size = 20)+
	facet_grid(~pretty_group, scales="free") + 
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase")


out = ggplot(seasonality_mode_to_plot %>% filter(#fcast_type != "Functional" & 
																					 	model_name=="cycl_only"),
			 aes(fill=pretty_group, color = pretty_group, y = freq,
			 		x = sampling_season, group=model_id)) +
	geom_line(alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~pretty_group, scales="free") + 
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase")
tag_facet(out, tag_pool = LETTERS)

ggplot(seasonality_mode_to_plot %>% filter(fcast_type != "Functional" & 
																					 	model_name=="cycl_only"),
			 aes(fill=pretty_group, color = pretty_group, y = freq,
			 		x = sampling_season)) +
	geom_smooth(alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~pretty_group, scales="free") + 
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase")


# CHeck against observed trends for speciifc taxa
# ggplot(plot_estimates %>% filter(taxon=="glucose_simple"),# %>% filter(siteID %in% c("HARV","BART","DSNY")),
ggplot(plot_estimates %>% filter(taxon=="alphaproteobacteria"),# %>% filter(siteID %in% c("HARV","BART","DSNY")),
			 aes(fill=species, x = as.numeric(month), y = `50%`)) +
	geom_jitter(aes(y = `50%`), show.legend = F, color="red", alpha=.1) +
	geom_smooth() +
	#geom_smooth(aes(y = `truth`)) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_jitter(aes(y = as.numeric(truth)), alpha=.2, height = 0, width=.2) + xlab(NULL) + labs(fill='') +
	facet_grid(~model_name, scales="free")




ggplot(cycl_only_vals,
			 aes(x=max_y_date, y=pretty_group)) +
	geom_jitter(aes(colour = pretty_group), width = 0.1, height=.1, alpha=.5, show.legend = F) +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Peak month of seasonal trend")) +
	xlab(NULL) + ylab(NULL)
#facet_grid(rows=vars(taxon_name)) +

seasonal_cycl = ggplot(seas_vals_long %>% filter(grepl("cycl_only",model_name)) %>%
											 	filter(grepl("ectomycorrhizal", taxon)),
											 aes(x=dates, y=y_cycl)) +
	geom_line(colour = "red") +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Seasonal trend in abundances")) +
	ylab("Cyclic component") +
	xlab(NULL) +
	facet_grid(rows=vars(model_id)) +
	#facet_wrap(~model_id) +
	scale_x_date(date_labels = "%b")
seasonal_cycl

pathogen_est = cycl_only_est  %>%
	filter(grepl("plant_pathogen", taxon))

ecto_est = cycl_only_est  %>%
	filter(grepl("ecto", taxon))
ecto_est_env_cov = env_cov_est  %>%
	filter(grepl("ecto", taxon))
ecto_est_env_cycl = env_cycl_est  %>%
	filter(grepl("ecto", taxon))




ggplot(max_dates_cycl_only,
			 aes(color=fcast_type, x = sampling_season)) +
	geom_point(aes(y = amplitude), alpha=.2, 
						 position=position_jitter(height=0, width=.2)) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(pretty_group~fcast_type, scales="free")




new_max_abun = plot_estimates %>% 
	filter(time_period == "2015-11_2018-01") %>% 
	group_by(model_id, siteID, month) %>% 
	summarize(mean_modeled_abun = mean(Mean, na.rm=T)) %>%
	filter(mean_modeled_abun == max(mean_modeled_abun, na.rm=T)) %>% 
	mutate(dates =  ymd(paste0("2014-",month, "-01")),)



# %>% group_by(fcast_type, pretty_group, model_id) %>% 
# 	filter(mean_modeled_abun == max(mean_modeled_abun, na.rm=T)) %>% 
# 	mutate(dates =  ymd(paste0("2014-",month, "-01")),)

cat_list=list()
for (i in 1:nrow(new_max_abun)){
	cat_list[[i]] <- assign_pheno_site_date(site = new_max_abun[i,]$siteID, date =  new_max_abun[i,]$dates)
}
site_cats = as.matrix(cat_list)
max_abun <- cbind.data.frame(new_max_abun, site_cats)


ggplot(seasonality_mode2 %>% filter(fcast_type=="Functional" & model_name != "env_cov"),
			 aes(color=pretty_group, x = sampling_season, y = amplitude)) +
	geom_point( alpha=.2, 
							position=position_jitter(height=0, width=.1)) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	# theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
	# 			legend.position = "bottom",legend.title = element_text(NULL),
	# 			plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
	ylab("Seasonal amplitude") + 
	facet_grid(model_name~pretty_group, scales="free") +
	geom_label_repel(aes(label=taxon))

ggplot(new_max_abun  %>% filter(rank_only=="phylum" & model_name != "env_cov"),
			 aes(color=pretty_group, x = sampling_season, y = amplitude)) +
	geom_point( alpha=.2, 
							position=position_jitter(height=0, width=.1)) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	# theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
	# 			legend.position = "bottom",legend.title = element_text(NULL),
	# 			plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
	ylab("Seasonal amplitude") + 
	facet_grid(model_name~pretty_group, scales="free") 

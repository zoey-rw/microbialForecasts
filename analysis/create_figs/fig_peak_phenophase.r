library(lubridate)
library(ggallin)
library(ggrepel)

source("source.R")

phenophase_in = readRDS(here("data/clean/group_peak_phenophases.rds"))
seasonality_mode_max = phenophase_in[[1]]
seasonality_mode_half = phenophase_in[[3]]
max_abun = phenophase_in[[4]] %>% 
	filter(significant_sin == 1 | significant_cos)
seasonality_mode_to_plot = phenophase_in[[5]] %>% 
	filter(significant_sin == 1 | significant_cos)



phenophase_in_legacy = readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
seasonality_mode_max = phenophase_in_legacy[[1]]
seasonality_mode_half = phenophase_in_legacy[[3]]
seasonality_mode_to_plot = phenophase_in_legacy[[5]]

max_abun = phenophase_in_legacy[[4]] %>% filter(siteID %in% c("CPER","DSNY","HARV","OSBS","STER"))
# Get most frequently-assigned phenophase per group	
seasonality_mode_filtered = max_abun %>% 
	filter(!is.na(sampling_season)) %>% 
	group_by(model_id, sampling_season) %>% 
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>% 
	select(-n) 
seasonality_mode_filtered <- merge(seasonality_mode_filtered, unique(max_abun[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")]), all.x=T)

seasonality_mode_to_plot = seasonality_mode_filtered



freq_phenophase = seasonality_mode_to_plot %>% 
	pivot_wider(values_from=freq, names_from=sampling_season, values_fill = 0)
freq_phenophase$greenup_peak = freq_phenophase$greenup + freq_phenophase$peak
freq_phenophase$greenup_peak_greendown = freq_phenophase$greenup + freq_phenophase$peak + freq_phenophase$greendown
freq_phenophase$max_phase_green = ifelse(freq_phenophase$greenup_peak > .5, 1, 0)
freq_phenophase$max_phase_growing = ifelse(freq_phenophase$greenup_peak_greendown > .5, 1, 0)
freq_phenophase_long = freq_phenophase %>% 
	pivot_longer(cols=c(greenup_peak, greendown, dormancy), 
							 names_to = "sampling_season", values_to="freq") 
freq_phenophase_long$sampling_season = ordered(freq_phenophase_long$sampling_season, levels = c("greenup_peak","greendown","dormancy"))


phenophase_proportion = freq_phenophase %>% 
	group_by(model_name, pretty_group, rank_only) %>% 
	add_tally(name = "total") %>% 
	ungroup() %>% 
	group_by(model_name, pretty_group, rank_only, total, max_phase_growing) %>% 
	tally(name = "count") %>% 
	filter(max_phase_growing == 1 ) %>% 
	mutate(proportion_growing_season = count/total) %>% select(-c(max_phase_growing, count))


ggplot(freq_phenophase_long %>% filter(#fcast_type == "Functional" & 
																			 	model_name=="cycl_only"),
			 aes(fill=pretty_group, y = freq,
			 		x = sampling_season, group=model_id)) +
	geom_line(aes(color = pretty_group), alpha=.5, show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~pretty_group, scales="free") + 
	ylab("Frequency of max abundance in phenophase") +
	xlab("Phenophase") +
	geom_label_repel(aes(label=taxon), size=4)


ggplot(max_abun,# %>% filter(fcast_type=="Functional"),
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

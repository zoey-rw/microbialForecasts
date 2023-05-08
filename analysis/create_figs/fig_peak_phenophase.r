library(lubridate)
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")


phenophase_in = readRDS(here("data/clean/group_peak_phenophases.rds"))
seasonality_mode2 = phenophase_in[[1]]
max_abun = phenophase_in[[2]]
freq_phenophase = phenophase_in[[5]] %>% pivot_wider(values_from=freq, names_from=sampling_season)

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

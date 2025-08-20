source("source.R")


sum.in <- readRDS(here("data", "summary/pheno_summaries.rds"))

plot_estimates = sum.in$plot_est %>% filter(model_id %in% sum.in$keep_list)
plot_estimates$month = lubridate::month(plot_estimates$dates)
plot_estimates$month_date = as.Date(paste0(plot_estimates$month, "-01-2016"), format = "%m-%d-%Y")
plot_estimates$fg_category <- microbialForecast:::assign_fg_categories(plot_estimates$taxon)

cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))


pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
bart_harv_pheno = pheno_categories_in[[1]] %>% filter(ID %in% c("BART","HARV") & year == "2016") %>% ungroup


northern_sites <- c("HARV",
										"BART")

plant_associated = env_cycl_est  %>%
	filter(taxon %in% c("plant_pathogen","oligotroph","heat_stress","copiotroph","saprotroph","lignolytic"))

plant_associated2 = cycl_only_est  %>%
	filter(taxon %in% c("plant_pathogen","oligotroph","heat_stress","copiotroph","saprotroph","lignolytic"))

#cycl_only_est$pretty_fg_names <- recode(cycl_only_est$taxon, !!!pretty_names)



ggplot(plant_associated %>% filter(siteID %in% northern_sites) %>% filter(taxon != "heat_stress"),
			 aes(x = month_date, color=siteID)) +
	
	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_classic()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(~taxon, ncol=1, scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site")



data_rich_sites = c("CPER","DSNY","HARV","OSBS","STER")
ggplot(plant_associated %>% 
			 	filter(siteID %in% data_rich_sites),
			 aes(x = as.numeric(month))) +
	geom_point(aes(y = `50%`, color = siteID), show.legend = F, #color="red", 
						 alpha=.1, position=position_jitter(height=0)) +
	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .5, 
						 position=position_jitter(height=0)) + 
	geom_smooth(aes(y = `50%`), se = F, color=1) +
	theme_bw( base_size = 22) +
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(fg_category~taxon, scales="free") + ggtitle("Full models, estimates plotted by month ")




ggplot(plant_associated %>% filter(siteID %in% data_rich_sites) %>% filter(taxon != "heat_stress"),
			 aes(x = month_date, color=siteID)) +
	
	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_classic()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(~taxon, ncol=1, scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site") + 	
	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .3,
						 position=position_jitter(height=0, width=.5))

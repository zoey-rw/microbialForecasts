# read in seasonal values
library(lubridate)
source("source.R")

scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list
#converged_strict = scores_list$converged_strict_list

sum.in <- readRDS(here("data", paste0("summary/logit_beta_regression_summaries.rds")))
plot_estimates = sum.in$plot_est %>% filter(model_id %in% converged)
plot_estimates$month = lubridate::month(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
cycl_only_est$fg_category <- microbialForecast:::assign_fg_categories(cycl_only_est$taxon)


env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cycl_est$fg_category <- microbialForecast:::assign_fg_categories(env_cycl_est$taxon)

seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))


seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% converged)
seas_vals_wide <- seas_in[[6]] %>% filter(model_id %in% converged)

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>% select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)

cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")

winter_groups = cycl_only_vals %>% filter(
	time_period=="2015-11_2018-01") %>%
	filter(max_month %in% c(1,2,12))
unique(winter_groups$taxon)


fg_seasonal = cycl_only_vals %>% filter(fcast_type=="Functional")
fg_seasonal$fg_category <- microbialForecast:::assign_fg_categories(fg_seasonal$taxon)

ggplot(fg_seasonal,
			 aes(x=max_y_date, y=fg_category)) +
	geom_jitter(aes(colour = pretty_group), width = 0.1, height=.1, alpha=.5, show.legend = F) +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Peak month of seasonal trend")) +
	xlab(NULL) + ylab(NULL)
#facet_grid(rows=vars(taxon_name)) +



cycl_only_est$month_date = as.Date(paste0(cycl_only_est$month, "-01-2016"), format = "%m-%d-%Y")

# "pretty_names" is from globalVariables.r but isn't exporting...
cycl_only_est$pretty_fg_names <- recode(cycl_only_est$taxon, !!!microbialForecast:::pretty_names)

simple_sugars = cycl_only_est %>% filter(fg_category=="Simple substrates" & fcast_type=="Functional")
complex_sugars = cycl_only_est %>% filter(fg_category=="Complex substrates" & fcast_type=="Functional")

simple_complex = cycl_only_est %>% filter(fg_category %in% c("Simple substrates","Complex substrates") & fcast_type=="Functional") %>%
	filter(!taxon %in% c("chitinolytic","cellulolytic","lignolytic"))



copio_oligo_est = cycl_only_est  %>%
	filter(grepl("otroph", taxon) &  time_period=="2015-11_2018-01")



pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
bart_harv_pheno = pheno_categories_in[[1]] %>% filter(ID %in% c("BART","HARV") & year == "2016") %>% ungroup

northern_sites <- c("ABBY", "BART", "BONA",  "DEJU", "HARV",
										"HEAL", "NIWO", "ONAQ", "RMNP", "STEI", "TOOL",
										"TREE", "UNDE", "WREF", "YELL")

northern_sites <- c("BART", "BONA", "HARV",
										"HEAL", "STEI", "TOOL",
										"TREE", "UNDE", "WREF", "YELL")

northern_sites <- c("HARV",
										"BART")

plant_associated = cycl_only_est  %>%
	filter(time_period=="2015-11_2018-01")  %>%
	filter(taxon %in% c("plant_pathogen","oligotroph","heat_stress","copiotroph","saprotroph"))


ggplot(plant_associated %>% filter(siteID %in% northern_sites) %>% filter(taxon != "heat_stress"),
			 aes(x = month_date, color=siteID)) +

	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_classic()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(~pretty_fg_names, ncol=1, scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site")



ggplot(simple_complex %>% filter(siteID %in% northern_sites) %>% filter(taxon != "heat_stress"),
			 aes(x = month_date, color=siteID)) +

	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_classic()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(fg_category~pretty_fg_names,  scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site")



ggplot(simple_sugars %>% filter(siteID %in% northern_sites),
			 aes(fill=species, x = as.numeric(month))) +
	#geom_point(aes(y = `Mean`), show.legend = F, color="red",
	#					 alpha=.1, position=position_jitter(height=0)) +
	#geom_point(aes(y = as.numeric(truth)), alpha = .3, position=position_jitter(height=0)) +
	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(~taxon, scales="free")

ggplot(complex_sugars %>% # filter(taxon %in% c("sucrose_complex",
												#											 "cellulose_complex")) %>%
			 	filter(siteID %in% northern_sites),
			 aes(fill=species, x = as.numeric(month))) +
	geom_point(aes(y = `Mean`), show.legend = F, color="red",
						 alpha=.1, position=position_jitter(height=0)) +
	geom_point(aes(y = as.numeric(truth)), alpha = .5, position=position_jitter(height=0)) +
	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	xlab(NULL) + labs(fill='') +
	facet_wrap(~taxon,nrow=2, scales="free")


copio_oligo_env_cycl_est = env_cycl_est  %>%
	filter(grepl("otroph", taxon) | taxon %in% c("glucose_simple",
																							 "cellulose_complex","acetate_simple") &  time_period=="2015-11_2018-01")

ggplot(copio_oligo_env_cycl_est %>% filter(siteID %in% northern_sites),
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
	facet_wrap(fg_category~taxon, scales="free")#, ncol=1)


all_decomposers <- cycl_only_est %>% filter(fcast_type=="Functional") %>%
	filter(fg_category %in% c("Complex substrates", "Simple substrates") | taxon %in% c("saprotroph","copiotroph","oligotroph"))

ggplot(all_decomposers %>%
			 	filter(siteID %in% northern_sites),
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
	facet_wrap(fg_category~taxon, scales="free")#, ncol=1)


all_decomposers_env <- env_cycl_est %>% filter(fcast_type=="Functional") %>%
	filter(fg_category %in% c("Complex substrates", "Simple substrates") | taxon %in% c("saprotroph","copiotroph","oligotroph"))

ggplot(all_decomposers_env %>%
			 	filter(siteID %in% northern_sites),
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


ggplot(all_decomposers_env %>%
			 	filter(siteID %in% c("HARV","BART")),
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

remotes::install_github("giocomai/ganttrify")

remotes::install_github("giocomai/ganttrify")
p_load(ganttrify)

data_in <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))


model_date_summaries <- data_in$plot_est %>% 
	group_by(time_period) %>% 
	filter(!is.na(truth)) %>% 
	summarise(start_date = min(dates, na.rm = T),
						end_date = max(dates, na.rm = T)) %>% 
	mutate(wp = recode(as.character(time_period), !!!microbialForecast:::calibration_label),
				 activity = recode(as.character(time_period), !!!microbialForecast:::calibration_label))

ganttrify(project = model_date_summaries,
					month_number_label = F,
					by_date = T,
					month_breaks = 6,
					font_family = "Roboto Condensed")

heat_stress
plant_pathogen
russula

cal_data <- data_in$plot_est %>%  filter(model_id %in% keep_list) %>% 
mutate(calibration_label= 
			 recode(as.character(time_period), `2015-11_2018-01` = "Partial dataset for calibration", 
			 																	 `2015-11_2020-01` = "Full dataset", 
			 			 `20130601_20151101` = "Provisional data"))
cal_data$calibration_label = factor(cal_data$calibration_label, ordered=T, levels = c("Provisional data",
																																											"Partial dataset for calibration","Full dataset"))

to_plot = cal_data %>% filter(#plotID=="HARV_001" & 
	taxon=="plant_pathogen" & model_name == "cycl_only" & plotID %in% c("CPER_004"))
# View actual models & hindcasts for a single site
ggplot(to_plot, aes(group = calibration_label)) + 
	# facet_grid(cols=vars(plotID), 
	# 					 rows = vars(time_period), drop=T, scales="free") +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`, fill = calibration_label), alpha=0.3) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`, fill = calibration_label), alpha=0.4) +
	theme_bw()+
	scale_fill_brewer(palette = "Set1") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab("Modeled abundance") +
	scale_y_sqrt() +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') + 
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) #+ 
#	ggtitle("Provisional dataset is quantitatively ")


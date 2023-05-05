# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce, ggrepel)
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(sp)
library(rgdal)
library(cowplot)
library(ggrepel)

# Read in hindcast data
hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
scores_list$converged_list

hindcast=hindcast_in %>%
	mutate(timepoint = date_num)

hindcast$month <- lubridate::month(hindcast$dates)
hindcast$month_label <- lubridate::month(hindcast$dates, label = T)

# Add duplicated data row so that plotting ribbons are continuous.
max_cal_date <- hindcast %>%
	group_by(taxon_name, plotID) %>%
	filter(fcast_period=="calibration") %>%
	filter(timepoint == max(timepoint, na.rm=T)) %>%
	mutate(fcast_period="hindcast")
hindcast <- rbind.data.frame(hindcast, max_cal_date)

all_cov_hindcast = hindcast %>% filter(is.na(site_prediction) |
																			 	site_prediction != "New time x site (random effect)")

select_plots <- c("BART_001","CPER_001","BLAN_001","CLBJ_001")
select_plots <- c("HARV_033","CPER_004")
taxon_names =  c("ectomycorrhizal")
taxon_names =  c("basidiomycota","mortierellomycota")
taxon_names =  c("denitrification","actinobacteria")
model_name = "env_cycl"


select_hindcasts <- all_cov_hindcast %>% filter(model_name == !!model_name & taxon %in% taxon_names)
select_hindcasts_select_plots <- select_hindcasts[select_hindcasts$plotID %in% select_plots,]


# COlor by species
ggplot(select_hindcasts_select_plots, aes(fill=species, x = dates, y =med, group=plotID)) +
	facet_grid(cols=vars(siteID),
						 rows=vars(species),
						 #cols=vars(factor(siteID, levels=c('KONZ','HARV'))),
						 drop=T, scales="free"	) +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo, ymax = hi), alpha=0.2) +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo_25, ymax = hi_75), alpha=.5) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo_25, ymax = hi_75), alpha=0.3) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo, ymax = hi), alpha=0.1) +
	geom_line(data = ~filter(.x, fcast_period=="calibration"),
						alpha=0.8) +
	geom_line(data = ~filter(.x, fcast_period=="hindcast"),
						aes(x = dates, y = med), alpha=0.3) +
	geom_point(aes(y = as.numeric(truth)), position = position_jitter()) +
	xlab(NULL) + labs(fill='') +
	geom_label_repel(data = ~filter(.x, timepoint==40),
									 aes(x=dates, y=med, label = plotID),
									 box.padding = 3, point.padding = 0) +
	scale_fill_brewer(palette = "Set2") +
	scale_color_brewer(palette = "Set2") +
	theme(panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	ggtitle("Hindcasts at 4 plots") +
	theme_minimal(base_size = 20) +
	scale_y_sqrt() +
	theme(legend.position = "none")


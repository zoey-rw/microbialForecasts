# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce, ggrepel, terra)
source("source.R")
library(sp)
library(cowplot)

#devtools::install_github("eliocamp/tagger")
library(tagger)


# Read in hindcast data
hindcast_in <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))
scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list

# Remove first timepoint from calibration
# hindcast = hindcast_in %>%
# 	filter(!(fcast_period=="calibration" & is_any_start_date)) %>%
# 	mutate(timepoint = date_num)

hindcast=hindcast_in %>%
	mutate(timepoint = date_num,
				 month = lubridate::month(dates),
				 month_label = lubridate::month(dates, label = T))

# Add duplicated data row so that plotting ribbons are continuous.
max_cal_date <- hindcast %>%
	group_by(taxon_name, plotID) %>%
	filter(fcast_period=="calibration") %>%
	filter(timepoint == max(timepoint, na.rm=T)) %>%
	mutate(fcast_period="hindcast")
hindcast <- rbind.data.frame(hindcast, max_cal_date)

all_cov_hindcast = hindcast %>% filter(is.na(site_prediction) | grepl("observed", site_prediction))
																			 # is.na(site_prediction) |
																			 # 	site_prediction != "New time x site (random effect)")
# select_hindcasts <- all_cov_hindcast %>% filter(model_name == "all_covariates" & taxon %in% c("cellulolytic","nitrification","heat_stress","plant_pathogen","copiotroph","acidobacteriota"))


#hindcast_in$fcast_period <- ifelse(hindcast_in$dates > "2018-01-01", "hindcast", "calibration")



select_plots <- c("HARV_033","OSBS_026","WOOD_044","KONZ_001")
select_plots <- c("HARV_033","WOOD_044")
select_plots <- c("HARV_033","KONZ_001")
select_plots <- c("HARV_033","HARV_004","HARV_001","HARV_034")
select_plots <- c("HARV_033","HARV_004","KONZ_001","KONZ_002")
select_plots <- c("HARV_033","HARV_004","KONZ_001","KONZ_002")

#select_plots <- c("BART_001","CPER_001","BLAN_001","CLBJ_001")

select_plots<- c("CPER_006", #good for plant path
"KONZ_043") #good for plant path

select_hindcasts <- hindcast %>% filter(model_name == "env_cycl" & taxon %in% c("dissim_nitrate_reduction","plant_pathogen"))
# select_hindcasts <- hindcast %>% filter(model_name == "all_covariates" & taxon %in% c("actinobacteriota","acidobacteriota"))

select_hindcasts_select_plots <- select_hindcasts[select_hindcasts$plotID %in% select_plots,]


# New facet label names for supp variable
tax.labs <- c("Bacteria (chitin-enriched)", "Fungi (Order: Mortierellales)", "Fungi (ectomycorrhizae)", "Fungi (plant pathogen)","Bacteria\n(dissimilatory-nitrate-reducer)")
names(tax.labs) <- c("chitin_complex", "mortierellales", "ectomycorrhizal","plant_pathogen","dissim_nitrate_reduction")
# New facet label names for supp variable
site.labs <- c("Harvard Forest, MA", "Konza Prairie, KS","Central Plains, CO")
names(site.labs) <- c("HARV", "KONZ", "CPER")

select_hindcasts_select_plots$siteID = ordered(select_hindcasts_select_plots$siteID, levels=c("CPER",'KONZ','HARV'))

# COlor by species
examples = ggplot(select_hindcasts_select_plots, aes(fill=species, x = dates, y =med, group=plotID)) +
	facet_grid(#cols=vars(siteID),
						 rows=vars(species),
						 cols=vars(siteID),
						 drop=T,
						 #space="free",
						 scales="free",
						 labeller = labeller(species = tax.labs, siteID = site.labs)
	) +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo, ymax = hi), alpha=0.35) +
	geom_ribbon(data = ~filter(.x, fcast_period=="calibration"),
							aes(x = dates,ymin = lo_25, ymax = hi_75), alpha=.65) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo_25, ymax = hi_75), alpha=0.2) +
	geom_ribbon(data = ~filter(.x, fcast_period=="hindcast"),
							aes(x = dates, ymin = lo, ymax = hi), alpha=0.08) +
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
	theme(panel.spacing = unit(.1, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.1, .2, 2, 0), "cm")) + ylab("Abundance") +
	#ggtitle("Hindcasts at 4 plots") +
	theme_classic(base_size = 18) +
	#scale_y_sqrt() +
	theme(legend.position = "none") +
	tag_facets(tag_pool = c("B","C","D","E")) +
#	tag_facets(position = list("x"=0,"y"=.8,"hjust"=0,"vjust" = 0)) +
	theme(tagger.panel.tag.text = element_text(size = 16),
				tagger.panel.tag.background = element_rect(fill = "white", color = "white"))

examples


options(stringsAsFactors=F)

neonDomains <- terra::vect("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/NEON_Domains.shp", layer="NEON_Domains")

neonSites = read_csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
# neonMercator <- spTransform(neonDomains,
# 														CRS("+proj=merc")) # requires newer sqlite version
neonSites <- neonSites %>% filter(grepl("Terrestrial", field_site_type)) %>%
	filter(field_latitude < 60)
neonDomains <- sf::st_read("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/NEON_Domains.shp", layer="NEON_Domains")

neonMap = ggplot(neonDomains) +
	geom_sf( fill = "white", color =1) +
	#geom_polygon(aes(x=long, y = lat, group=group), 	fill = "white", color =1) +
																		#fill = neonDomains@data$DomainID)) +
	theme_no_axes() +
	geom_point(data=neonSites, aes(x=field_longitude, y=field_latitude), color ="black", size=4, alpha=.6)  +
	geom_point(data=neonSites %>% filter(field_site_id %in% c("CPER", "KONZ")),
																			 aes(x=field_longitude, y=field_latitude), color ="red", size=5) +
	geom_label_repel(data = neonSites %>% filter(field_site_id %in% c("CPER", "KONZ")),
									 aes(x=field_longitude, y=field_latitude, label = word(field_site_name,1, 2)), size=5,
								 box.padding = 4, point.padding = .1) +
	ylim(c(25,50)) + xlim(c(-125,-65)) +
	theme(
				panel.border = element_blank(),
		panel.background = element_rect(fill='transparent'), #transparent panel bg
		plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
		panel.grid.major = element_blank(), #remove major gridlines
		panel.grid.minor = element_blank() #remove minor gridlines
	)  + 	theme(plot.margin = unit(c(0,0,0,0), "cm"))

# Map next to time-series
fig_to_save <- ggarrange(neonMap, examples, nrow=2,
												 #widths = c(1, 2)),
													heights = c(1, 2), labels=c("A)", NULL))

# Map within time-series
# fig_to_save = ggdraw() +
# 	draw_plot(examples) +
# 	draw_plot(neonMap, x = 0.25, y = .39, width = .5, height = .35)
# fig_to_save

png(here("figures","example_hindcasts.png"), width = 1000, height=1400)
print(fig_to_save)

dev.off()




ggplot(hindcast_in %>% filter(model_name == "all_covariates" & siteID=="SJER" & taxon=="nitrification")) +
	facet_grid(rows=vars(plotID), drop=T, scales="free") +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = med), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +ggtitle("'Good' forecast",  "lowest CRPS score, for site x taxon combination")

ggplot(hindcast_in %>% filter(model_name == "all_covariates" & siteID=="OSBS" & taxon=="agaricomycetes")) +
	facet_grid(rows=vars(plotID), drop=T, scales="free") +
	#rows = vars(fcast_type), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +ggtitle("'Bad' forecast",  "highest CRPS score, for site x taxon combination")

# one plot for testing
ggplot(good_hindcasts_top_plots %>% filter(plotID=="CPER_004")) +
	facet_grid(rows=vars(taxon), drop=T, scales="free") +
		#rows = vars(fcast_type), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')

good_hindcasts_top_plots[good_hindcasts_top_plots$taxon_name == "armatimonadota",]

pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_good_example.pdf", height = 6)
for(i in 1:length(top_plots)){
	print(ggplot(good_hindcasts_top_plots) +
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
					facet_grid_paginate(
						taxon~plotID,
						drop=T, scales="free",
						ncol = 1, nrow = 3, page = i))
}
dev.off()


pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_bad_example.pdf", height = 6)
for(i in 1:length(top_plots)){
	print(ggplot(bad_hindcasts_top_plots) +
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
					facet_grid_paginate(
						taxon~plotID,
						drop=T, scales="free",
						ncol = 1, nrow = 3, page = i))
}
dev.off()




seasonality= readRDS(here("data/summary/seasonal_amplitude.rds"))
all_cov_vals = seasonality[[1]]
cycl_vals = seasonality[[2]]
all_cov_vals_refit = seasonality[[3]]
cycl_vals_refit = seasonality[[4]]

taxon_name = "mortierellales"
taxon_name = "chitin_complex"
taxon_names = c("mortierellales", "chitin_complex")

taxon_names = unique(cycl_vals$taxon)


input_dateID = c("201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201412")
input_months = rep(seq(1:12),2)
dates = fixDate(input_dateID)
dates_month = lubridate::month(dates)
input_date_df = data.frame(x = dates_month,
													 dates = dates,
													 dates_month_label = lubridate::month(dates_month, label = T))

plot_list = list()
out_seas_vals = list()
for (taxon_name in taxon_names) {

alpha = cycl_vals[cycl_vals$taxon==taxon_name,]$sin
beta = cycl_vals[cycl_vals$taxon==taxon_name,]$cos
alpha2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$sin
beta2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$cos
alpha_refit = cycl_vals_refit[cycl_vals_refit$taxon==taxon_name,]$sin
beta_refit = cycl_vals_refit[cycl_vals_refit$taxon==taxon_name,]$cos
alpha2_refit = all_cov_vals_refit[all_cov_vals_refit$taxon==taxon_name,]$sin
beta2_refit = all_cov_vals_refit[all_cov_vals_refit$taxon==taxon_name,]$cos

# df = data.frame(x = input_date_df$month,
# 								dates = input_date_df$dates)
df = input_date_df
df$y_cycl = alpha * sin(2*pi*df$x/12) + beta * cos(2*pi*df$x/12)
# df$y_allcov = alpha2 * sin(2*pi*df$x/12) + beta2 * cos(2*pi*df$x/12)
# df$y_allcov_refit = alpha2_refit * sin(2*pi*df$x/12) + beta2_refit * cos(2*pi*df$x/12)
# df$y_cycl_refit = alpha_refit * sin(2*pi*df$x/12) + beta_refit * cos(2*pi*df$x/12)

out_seas_vals[[taxon_name]] = df %>% mutate(taxon_name = !!taxon_name)

seasonal_cycl = ggplot(df, aes(x=dates, y=y_cycl)) +
	geom_line(colour = "red") +
	theme_bw(base_size = 20) +
	#ggtitle(paste0("Seasonal trend in abundance: ", taxon_name)) +
	ylab("Cyclic component") +
	#ylim(c(-.2,.2))
ylim(c(-.1,.1))
#
# seasonal_cycl_refit = ggplot(df, aes(x=dates, y=y_cycl_refit)) +
# 	geom_line(colour = "red") +
# 	theme_bw(base_size = 20) +
# 	ggtitle(paste0("Seasonal trend in abundance: ", taxon_name)) + ylab("Cyclic component") + ylim(c(-.2,.2))

plot_list[[taxon_name]] <- seasonal_cycl
# seasonal_allcov = ggplot(df, aes(x=dates, y=y_allcov)) +
# 	geom_line(colour = "red") +
# 	theme_bw(base_size = 16) +
# 	ggtitle("Residual seasonal trend in abundance") + ylab("Abundance")

}

seas_vals_to_plot = rbindlist(out_seas_vals, fill=T)

seasonal_cycl = ggplot(seas_vals_to_plot, aes(x=dates, y=y_cycl)) +
	geom_line(colour = "red") +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Seasonal trend in abundances")) +
	ylab("Cyclic component") +
	xlab(NULL) +
	#facet_grid(rows=vars(taxon_name)) +
	facet_wrap(~taxon_name) +
	scale_x_date(date_labels = "%b")
seasonal_cycl

ggarrange(plot_list[[1]],plot_list[[2]],nrow=1, common.legend = T)
ggarrange(plot_list[[1]],
					plot_list[[2]],
					plot_list[[3]],
					plot_list[[4]],
					plot_list[[5]],
					plot_list[[6]],
					plot_list[[7]],
					ncol=1, common.legend = T)
ggarrange(plot_list[1:7],ncol=1, common.legend = T)



hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))
hindcast_in$month <- lubridate::month(hindcast_in$dates)
hindcast_in$month_label <- lubridate::month(hindcast_in$dates, label = T)
hindcast_in <- hindcast_in %>% mutate(time_period = recode(time_period, !!!microbialForecast:::date_recode))


seasonality= readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_vals = seasonality[[1]]
#cycl_vals = seasonality[[2]]
seas_vals = seasonality[[6]]



ggplot(cycl_vals %>% filter(grepl("cycl",model_name)),
			 aes(x=max_y_date, y=pretty_group)) +
	geom_point(aes(colour = pretty_group)) +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Seasonal trend in abundances")) +
	xlab(NULL)

# Merge seasonality with hindcast data
# All cov model
hindcast_seas <- merge(hindcast_in %>% filter(is.na(site_prediction) | site_prediction != "New time x site (random effect)"),
											 seas_vals, all=T)

cycl_only <-  hindcast_seas %>%
	filter(model_name == "cycl_only") %>%
	distinct()


seasonal_cycl = ggplot(seas_vals_long %>% filter(grepl("cycl_only",model_name)),
											 aes(x=dates, y=y_cycl)) +
	geom_line(aes(color = taxon)) +
	theme_bw(base_size = 20) +
	ggtitle(paste0("Seasonal trend in abundances")) +
	ylab("Cyclic component") +
	xlab(NULL) +
	#facet_grid(rows=vars(taxon_name)) +
	facet_wrap(~pretty_group, nrow=2) +
	scale_x_date(date_labels = "%b")
seasonal_cycl


phylum_cycl_only <-  hindcast_seas %>%
	filter(pretty_name == "phylum" & model_name == "cycl_only") %>%
	distinct()
phylum_env_cov <-  hindcast_seas %>%
	filter(pretty_name == "phylum" & model_name == "env_cov") %>%
	distinct()
phylum_env_cycl <-  hindcast_seas %>%
	filter(pretty_name == "phylum" & model_name == "env_cycl") %>%
	distinct()

phylum_cycl_only[phylum_cycl_only$amplitude > .05,]
table(seas_vals[seas_vals$amplitude > .05,]$taxon)

taxon %in% c("chitin_complex","mortierellales")

select_hindcasts <- cycl_only %>% filter(taxon %in% c("ascomycota","saprotroph","lignolytic"))
select_hindcasts <- cycl_only %>% filter(taxon %in% c("russula","saprotroph","ectomycorrhizal"))
select_hindcasts <- cycl_only %>% filter(taxon %in% c("basidiomycota","copiotroph","lignolytic"))
select_hindcasts <- hindcast_seas %>% filter(taxon %in% c("ectomycorrhizal","chitinolytic"))
select_plots <- c("HARV_033","OSBS_026","WOOD_044","KONZ_001")
select_plots <- c("BART_001","CLBJ_001")

ggplot(select_hindcasts  %>% filter(plotID %in% select_plots), 
			 aes(fill=taxon, x = dates, y =med, group=plotID)) +
	facet_grid(siteID~ species+model_name,
						 #cols=vars(factor(siteID, levels=c('KONZ','HARV'))),
						 drop=T, scales="free"    ) +
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
	scale_fill_brewer(palette = "Set2") +
	scale_color_brewer(palette = "Set2") +
	theme(panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	ggtitle("Hindcasts at 4 plots") +
	theme_minimal(base_size = 20) +
	scale_y_sqrt() +
	theme(legend.position = "none")


ggplot(select_hindcasts  %>% filter(plotID %in% select_plots), 
			 aes(fill=taxon, x = dates, y =med, group=plotID)) +
	facet_grid(cols=vars(siteID),
						 rows=vars(species),
						 #cols=vars(factor(siteID, levels=c('KONZ','HARV'))),
						 drop=T, scales="free"    ) +
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
	scale_fill_brewer(palette = "Set2") +
	scale_color_brewer(palette = "Set2") +
	theme(panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	ggtitle("Hindcasts at 4 plots") +
	theme_minimal(base_size = 20) +
	scale_y_sqrt() +
	theme(legend.position = "none")


b <- ggplot(phylum_cycl_only %>% filter(taxon %in% c("chloroflexi","chytridiomycota") &
																	 	siteID %in% c("HARV","OSBS","CPER"))) +
	geom_smooth(aes(x = month, y = `50%`,
									color = siteID),
							size=1, alpha = .2, na.rm = T, se = F) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends") +
	theme_bw(base_size = 18) +
	#scale_y_log10()
	scale_y_sqrt()


a <- ggplot(phylum_all_cov %>% filter(taxon %in% c("chloroflexi","chytridiomycota") &
																	 	siteID %in% c("HARV","OSBS","CPER"))) +
	geom_jitter(aes(x = month, y = `50%`,
								 color = siteID),
						 size=1, alpha = .4, na.rm = T, se = F) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Modeled abundances") +
	theme_bw(base_size = 18) +
	ylim(c(0,.1)) +
	#scale_y_log10()
scale_y_sqrt()


c <- ggplot(phylum_all_cov %>% filter(taxon %in% c("chloroflexi","chytridiomycota") &
																				siteID %in% c("HARV","OSBS","CPER"))) +
	geom_smooth(aes(x = month, y = `50%`,
									color = siteID),
							size=1, alpha = .4, na.rm = T, se = F) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Modeled abundances") +
	theme_bw(base_size = 18) +
	ylim(c(0,.1)) +
	#scale_y_log10()
	scale_y_sqrt()

ggarrange(a,b, nrow=2, common.legend = T)
ggarrange(a,b, c,nrow=3, common.legend = T)




# Couldn't figure out how to vectorize.
vals$max <- NA
vals$amplitude <- NA
for (i in 1:nrow(vals)){
	vals$max[[i]] <- getMaxMin(vals$sin[[i]],
														 vals$cos[[i]])
	vals$amplitude[[i]] <- sqrt(vals$sin[[i]]^2 + vals$cos[[i]]^2)
}
vals <- as.data.frame(vals)
vals$max <- unlist(vals$max)
vals$amplitude <- unlist(vals$amplitude)

cycl_vals$plot_sin <- sin((2*pi*cycl_vals$max)/12)
cycl_vals$plot_cos <- cos((2*pi*cycl_vals$max)/12)



all_cov_vals$plot_sin <- sin((2*pi*all_cov_vals$max)/12)
all_cov_vals$plot_cos <- cos((2*pi*all_cov_vals$max)/12)


taxon_name = "chytridiomycota"
taxon_name = "mortierellales"
x<-seq(from = 1,to=12,by = 0.1)

alpha = cycl_vals[cycl_vals$taxon==taxon_name,]$sin
beta = cycl_vals[cycl_vals$taxon==taxon_name,]$cos
seasonal_cycl = curve(alpha * sin(x) + beta * cos(x), lwd=3, lty=3, col="Red")


alpha2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$sin
beta2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$cos
seasonal_allcov = curve(alpha2 * sin(2*pi*x)/12 + beta2 * cos(2*pi*x)/12, lwd=3, lty=3, col="Red")
seasonal_allcov = curve(exp(alpha2 * sin(x) + beta2 * cos(x)), lwd=3, lty=3, col="Red")


abun = ggplot(hindcast_all_cov %>% filter( taxon %in% c("chitin_complex","mortierellales") &
																	 	siteID %in% c("HARV","KONZ"))) +
	geom_smooth(aes(x = month, y = mean,
									color = siteID),
							size=1, alpha = .2, na.rm = T, se = F) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Chytridiomycota: modeled seasonal abundances") +
	theme_bw(base_size = 16) +
	#scale_y_log10()
	scale_y_sqrt()


df = data.frame(x = seq(from = 1,to=12,by = 0.1))
seasonal_cycl = ggplot(df, aes(x)) +
	geom_function(fun = function(x) alpha * sin(2*pi*x/12) + beta * cos(2*pi*x/12),, colour = "red") +
	theme_bw(base_size = 16) +
	ggtitle("Seasonality")
seasonal_allcov = ggplot(df, aes(x)) +
	geom_function(fun = function(x) alpha2 * sin(2*pi*x/12) + beta2 * cos(2*pi*x/12), colour = "red") +
	theme_bw(base_size = 16) +
	ggtitle("Seasonality unexplained by environmental predictors")

ggarrange(abun,seasonal_cycl,seasonal_allcov,nrow=3, common.legend = T)


seasonal_allcov = curve(alpha2 * sin(2*pi*x/12) + beta2 * cos(2*pi*x/12), lwd=3, lty=3, col="Red")

month_fraction = 1:12/12
seasonal_allcov = curve(alpha2 *  sin(2*pi*x/12) + beta2 * cos(2*pi*x/12), lwd=3, lty=3, col="Red")


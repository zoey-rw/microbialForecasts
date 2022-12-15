
hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds"))
hindcast_data$month <- lubridate::month(hindcast_data$dates)
hindcast_data$month_label <- lubridate::month(hindcast_data$dates, label = T)


seasonality= readRDS(here("data/summary/seasonalAmplitude.rds"))
all_cov_vals = seasonality[[1]]
cycl_vals = seasonality[[2]]

# Merge seasonality with hindcast data
# All cov model
hindcast_all_cov <- merge(hindcast_data %>% filter(model_name == "all_covariates") %>% select(-rank_only),
													all_cov_vals, by = c("taxon", "model_name", "fcast_type",
																							 "pretty_name", "pretty_group",
																							 "only_rank"))
phylum_all_cov <-  hindcast_all_cov %>%
	filter(pretty_name == "Phylum") %>%
	distinct()

# Cycl model
hindcast_cycl_only <- merge(hindcast_data %>% filter(model_name == "cycl_only") %>% select(-rank_only),
														cycl_vals, by = c("taxon", "model_name", "fcast_type",
																							"pretty_name", "pretty_group",
																							"only_rank"))
phylum_cycl_only <-  hindcast_cycl_only %>%
	filter(pretty_name == "Phylum") %>%
	distinct()

cycl_vals[cycl_vals$amplitude > .1 & cycl_vals$pretty_name=="Phylum",]


ggplot(phylum_all_cov %>% filter(taxon %in% c("chloroflexi","chytridiomycota"))) +
	geom_smooth(aes(x = month, y = `50%`,
									group = siteID),
							size=1, alpha = .2, na.rm = T, se = F) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~taxon, scales="free") +
	ggtitle("Seasonal trends (all-covariate model)") +
	theme_bw(base_size = 18)


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
x<-seq(from = 1,to=12,by = 0.1)

alpha = cycl_vals[cycl_vals$taxon==taxon_name,]$sin
beta = cycl_vals[cycl_vals$taxon==taxon_name,]$cos
seasonal_cycl = curve(alpha * sin(x) + beta * cos(x), lwd=3, lty=3, col="Red")


alpha2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$sin
beta2 = all_cov_vals[all_cov_vals$taxon==taxon_name,]$cos
seasonal_allcov = curve(alpha2 * sin(x) + beta2 * cos(x), lwd=3, lty=3, col="Red")


abun = ggplot(phylum_all_cov %>% filter(taxon %in% c("chytridiomycota") &
																	 	siteID %in% c("HARV","OSBS","CPER","STER"))) +
	geom_smooth(aes(x = month, y = `50%`,
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
	geom_function(fun = function(x) alpha * sin(x) + beta * cos(x), colour = "red") +
	theme_bw(base_size = 16) +
	ggtitle("Seasonality")
seasonal_allcov = ggplot(df, aes(x)) +
	geom_function(fun = function(x) alpha2 * sin(x) + beta2 * cos(x), colour = "red") +
	theme_bw(base_size = 16) +
	ggtitle("Seasonality unexplained by environmental predictors")

ggarrange(abun,seasonal_cycl,seasonal_allcov,nrow=3, common.legend = T)

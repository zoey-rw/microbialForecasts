# Combine effect size estimates (beta covariates) from all workflows
# This script must be run before other analysis scripts

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups
sum.in <- readRDS(here("data", paste0("summary/logit_beta_regression_summaries.rds")))
sum.all <- sum.in$summary_df  %>% mutate(tax_rank = rank)


df <- data.table::rbindlist(list(sum.all), fill=TRUE) %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))


# Add prettier data values
df$pretty_name <- recode(df$rank, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$tax_rank, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

df$tax_rank <- ordered(df$tax_rank, levels = c("genus_bac","genus_fun",
																														 "family_bac","family_fun",
																														 "order_bac", "order_fun",
																														 "class_bac", "class_fun",
																														 "phylum_bac","phylum_fun",
																														 "functional_group", "diversity_16S", "diversity_ITS"))


# For saving: filter by effect type

# Linear model beta (covariate) effects
beta_effects <- df %>% filter(grepl("beta", rowname))
beta_effects$beta <- ordered(beta_effects$beta, levels = c("sin", "cos",
																			 "Ectomycorrhizal trees",
																			 "LAI",
																			 "pC",
																			 "pH",
																			 "Temperature",
																			 "Moisture","rho"))
levels(beta_effects$beta)[levels(beta_effects$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
saveRDS(beta_effects, here("data", "summary/predictor_effects.rds"))

# Site effects
site_effects <- df %>% filter(grepl("site", rowname))
saveRDS(site_effects, here("data", "summary/site_effects.rds"))


# Linear model beta (covariate) effects
rho_effects <- df %>% filter(grepl("rho", rowname))
saveRDS(rho_effects, here("data", "summary/rho_effects.rds"))

# Seasonality (cos/sin) effects
seas_params <- df %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
seas_vals <- seas_params %>% pivot_wider(id_cols = c("taxon","model_name","fcast_type","time_period",
																												"pretty_name","pretty_group","rank","only_rank"),
																						names_from = beta,
																						values_from = "Mean")

# Convert to amplitude and max
# Didn't vectorize this function, oops
out <- list()
for (i in 1:nrow(seas_vals)) {
	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
}
out <- rbindlist(out)
seas_vals <- cbind.data.frame(seas_vals, out)


cycl_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
																		model_name == "cycl_only")
all_cov_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
																		 	model_name == "all_covariates")

cycl_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
																					model_name == "cycl_only")
all_cov_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
																					 	model_name == "all_covariates")
saveRDS(list(all_cov_vals, cycl_vals, all_cov_vals_refit, cycl_vals_refit), here("data/summary/seasonalAmplitude.rds"))





ggplot(rho_effects %>% filter(model_name=="all_covariates"),
			 aes(x = pretty_group,y = abs(Mean), color=pretty_group)) +
	geom_boxplot() +
	geom_jitter( size = 5, height = 0, width=.4, alpha = .3,
							 shape = 16, show.legend = F) +
	ylab(NULL) + stat_compare_means()




input_coefficients = list(sin = -0.0101773497, cos = -0.019152772)

input_dates = c("201301","201302","201303","201304","201305","201306","201307","201308","201309","201310","201311","201312","201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201412")
input_months = rep(seq(1:12),2)

	plot_seasonal_trend = function(input_coefficients, input_dates) {


		asDate = fixDate(input_dates)
		dates_month = lubridate::month(asDate)
		dates_month_label <- lubridate::month(dates_month, label = T)
		x = dates_month
		#x<-seq(from = 1,to=12,by = 0.1)

		sin_coeff = input_coefficients$sin
		cos_coeff = input_coefficients$cos

		sin_vals = sin(2*pi*x/12) %>% scale
		cos_vals = cos(2*pi*x/12) %>% scale

		seas_fun = function(x) sin_coeff * sin_vals + cos_coeff * cos_vals

		y = seas_fun(x)
		df = data.frame(input_dates, asDate, dates_month_label, x, y)

		seas_plot = ggplot(df, aes(x=asDate, y = y)) +
			geom_line(colour = "red") +
			theme_bw(base_size = 16) +
			ggtitle("Seasonality")
return(seas_plot)
	}

	taxon_name="basidiomycota"
	taxon_name="basidiomycota"
	taxon_name="ascomycota"

	taxon_calibration= sum.in$plot_est %>% filter(taxon==taxon_name & siteID %in% c("HARV"))
	taxon_calibration= sum.in$plot_est %>% filter(taxon==taxon_name & plotID %in% select_plots) #c("HARV_033","HARV_004"))
	input_dates = taxon_calibration$dateID

	taxon_row = cycl_vals %>% filter(taxon==!!taxon_name)
	input_coefficients = data.frame(sin = taxon_row$sin,
												 cos = taxon_row$cos)
	plot_seasonal_trend(input_coefficients, input_dates)



	examples = ggplot(taxon_calibration, aes(fill=species, x = dates,group=plotID)) +
		#rows = vars(fcast_type), drop=T, scales="free") +
		geom_ribbon(
								aes(x = dates,ymin = `2.5%`, ymax = `97.5%`), alpha=0.2) +
		geom_ribbon(
								aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=.5) +
		geom_line(aes(x = dates, y = `50%`), alpha=0.3) +
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
	examples

	seasonal_cycl = ggplot(df, aes(x=dates, y=y_cycl)) +
		geom_line(colour = "red") +
		theme_bw(base_size = 20) +
		ggtitle(paste0("Seasonal trend in abundance: ", taxon_name)) + ylab("Cyclic component") +
		#ylim(c(-.2,.2))
		ylim(c(-.1,.1))



	# Check the predicted~observed time series
	ggplot(sum.in$plot_est) +
		geom_point(aes( x = Mean, y = as.numeric(truth),
										color = siteID), show.legend = F, alpha=.5)  +
		ylab("Observed") + xlab("Predicted") +
		theme_minimal(base_size = 18) +
		ggtitle("Overestimating plot means") +
		geom_abline(slope=1, intercept = 0) + facet_wrap(~taxon)

	# Check the predicted~observed time series
	ggplot(sum.in$plot_est %>% filter(taxon %in% c("basidiomycota"))) +
		geom_point(aes( x = Mean, y = as.numeric(truth),
										color = siteID), show.legend = F, alpha=.5)  +
		ylab("Observed") + xlab("Predicted") +
		theme_minimal(base_size = 18) +
		ggtitle("Overestimating plot means") +
		geom_abline(slope=1, intercept = 0)




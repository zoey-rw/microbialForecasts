# Combine & separately save model parameter and effect size estimates (beta covariates) from all models

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups
sum.in <- readRDS(here("data", paste0("summary/clr_regression_summaries.rds")))
sum.all <- sum.in$summary_df  %>% filter(model_name != "all_covariates") %>% 
	mutate(tax_rank = rank,
				 time_period = recode(time_period, !!!microbialForecast:::date_recode))
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# df$tax_rank <- ordered(df$tax_rank, levels = c("genus_bac","genus_fun",
# 																														 "family_bac","family_fun",
# 																														 "order_bac", "order_fun",
# 																														 "class_bac", "class_fun",
# 																														 "phylum_bac","phylum_fun",
# 																														 "functional_group", "diversity_16S", "diversity_ITS"))


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
saveRDS(beta_effects, here("data", "summary/clr_predictor_effects.rds"))

# # Site effects
# site_effects <- df %>% filter(grepl("site", rowname))
# saveRDS(site_effects, here("data", "summary/site_effects.rds"))
# 
# # Linear model beta (covariate) effects
# rho_effects <- df %>% filter(grepl("rho", rowname) | grepl("core_sd", rowname))
# saveRDS(rho_effects, here("data", "summary/rho_core_sd_effects.rds"))
# 
# # Seasonality (cos/sin) effects
# seas_params <- df %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
# seas_vals <- seas_params %>% pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","time_period",
# 																										 #"pretty_name","
# 																										 "pretty_group","rank","rank_only"),
# 																				 names_from = beta,
# 																				 values_from = c("Mean","significant")) %>% rename(sin="Mean_sin", cos = "Mean_cos")
# 
# # Convert to amplitude and max
# # Didn't vectorize this function, oops
# out <- list()
# for (i in 1:nrow(seas_vals)) {
# 	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
# }
# out <- rbindlist(out)
# seas_vals <- cbind.data.frame(seas_vals, out)
# 
# 
# cycl_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
# 																		model_name == "cycl_only")
# all_cov_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
# 																		 	model_name == "all_covariates|env_cov")
# 
# cycl_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
# 																					model_name == "cycl_only")
# all_cov_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
# 																					 	model_name == "all_covariates|env_cov")
# 
# 
# 
# 
# input_dateID = c("201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201412")
# dates = fixDate(input_dateID)
# input_date_df = data.frame(x = lubridate::month(dates),
# 													 dates = dates)
# out_seas_vals =list()
# 
# seas_vals_only = seas_vals %>% filter(grepl("cycl",model_name))
# for (row in 1:nrow(seas_vals_only)){
# 	
# 	alpha = seas_vals_only[row,]$sin
# 	beta = seas_vals_only[row,]$cos
# 	df = input_date_df
# 	y_cycl = alpha * sin(2*pi*input_date_df$x/12) + beta * cos(2*pi*input_date_df$x/12)
# 	names(y_cycl) = input_date_df$dates
# 	out_seas_vals[[row]] <- data.frame(y_cycl) %>% t %>% as.data.frame()
# }
# seas_vals_to_plot = rbindlist(out_seas_vals, fill=T)
# seas_vals_to_plot <- cbind.data.frame(seas_vals_only, seas_vals_to_plot)
# seas_vals_long = seas_vals_to_plot %>% 
# 	pivot_longer(cols = c(16:26), names_to = "dates", values_to = "y_cycl") %>%
# 	mutate(dates = as.Date(dates))
# 
# max_vals =	seas_vals_long %>% group_by(model_name, taxon, time_period, model_id) %>%
# 	filter(y_cycl == max(y_cycl, na.rm=T)) %>%
# 	mutate(max_y_date = dates) %>% select(-c(dates, y_cycl))
# seas_vals_long <- merge(seas_vals_long, max_vals, all=T)
# 
# saveRDS(list(seas_vals_long, all_cov_vals, cycl_vals, all_cov_vals_refit, cycl_vals_refit, seas_vals), here("data/summary/seasonal_amplitude.rds"))

#
# ggplot(rho_effects %>% filter(model_name=="all_covariates"),
# 			 aes(x = pretty_group,y = abs(Mean), color=pretty_group)) +
# 	geom_boxplot() +
# 	geom_jitter( size = 5, height = 0, width=.4, alpha = .3,
# 							 shape = 16, show.legend = F) +
# 	ylab(NULL) + stat_compare_means()
#
#
#
# 	select_plots <- c("HARV_033","OSBS_026","WOOD_044","KONZ_001")
#
# 	examples = ggplot(sum.in$plot_est %>% filter(taxon %in% c("basidiomycota") & plotID %in% select_plots),
# 																							 aes(fill=species, x = dates, group=plotID)) +
# 		#rows = vars(fcast_type), drop=T, scales="free") +
# 		geom_ribbon(
# 								aes(x = dates,ymin = `2.5%`, ymax = `97.5%`), alpha=0.2) +
# 		geom_ribbon(
# 								aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=.5) +
# 		geom_line(aes(x = dates, y = `50%`), alpha=0.3) +
# 		geom_point(aes(y = as.numeric(truth)), position = position_jitter()) +
# 		xlab(NULL) + labs(fill='') +
# 		scale_fill_brewer(palette = "Set2") +
# 		scale_color_brewer(palette = "Set2") +
# 		theme(panel.spacing = unit(.2, "cm"),
# 					legend.position = "bottom",legend.title = element_text(NULL),
# 					plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 		ggtitle("Hindcasts at 4 plots") +
# 		theme_minimal(base_size = 20) +
# 		scale_y_sqrt() +
# 		theme(legend.position = "none") + facet_grid(plotID~model_name, scales="free")
# 	examples
#
# 	seasonal_cycl = ggplot(df, aes(x=dates, y=y_cycl)) +
# 		geom_line(colour = "red") +
# 		theme_bw(base_size = 20) +
# 		ggtitle(paste0("Seasonal trend in abundance: ", taxon_name)) + ylab("Cyclic component") +
# 		#ylim(c(-.2,.2))
# 		ylim(c(-.1,.1))
#
#
# 	# Check the predicted~observed time series
# 	ggplot(sum.in$plot_est) +
# 		geom_point(aes( x = Mean, y = as.numeric(truth),
# 										color = siteID), show.legend = F, alpha=.5)  +
# 		ylab("Observed") + xlab("Predicted") +
# 		theme_minimal(base_size = 18) +
# 		ggtitle("Overestimating plot means") +
# 		geom_abline(slope=1, intercept = 0) + facet_wrap(~taxon)




seas_vals =readRDS(here("data/summary/seasonal_amplitude.rds"))
all_seas_vals = seas_vals[[5]]


taxon_name="ectomycorrhizal_ectomycorrhizal"
input_dates = unique(select_hindcasts$dateID)
input_dates = c("201601", "201602", "201603", "201605",
								"201606", "201607", "201608", "201609", "201610", "201611", "201612")

taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "2015-11_2018-01")
input_coefficients = data.frame(sin = taxon_row$sin,
																cos = taxon_row$cos)
trend_201511_201801 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2015-11_2018-01")

taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "2015-11_2020-01")
input_coefficients = data.frame(sin = taxon_row$sin,
																cos = taxon_row$cos)
trend_201511_202001 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2015-11_2020-01")

taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "20130601_20151101")
input_coefficients = data.frame(sin = taxon_row$sin,
																cos = taxon_row$cos)
trend_20130601_20151101 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2013-06_2015-11")
ggarrange(trend_20130601_20151101, trend_201511_201801, trend_201511_202001, nrow = 3)


summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
plot_est=summaries$plot_est %>% filter(model_name=="env_cycl")
plot_est_cycl=summaries$plot_est %>% filter(model_name=="cycl_only")

plot_est <- plot_est[!(plot_est$time_period %in% c("2015-11_2018-01","2015-11_2020-01") & plot_est$date_num == 5),]
plot_est_cycl <- plot_est_cycl[!(plot_est_cycl$time_period %in% c("2015-11_2018-01","2015-11_2020-01") & plot_est_cycl$date_num == 5),]



ggplot(plot_est_cycl %>% filter(plotID %in% c("HARV_004","HARV_013"),
													 taxon %in% c("ectomycorrhizal") &
													 	model_name == "cycl_only" &
													 	time_period %in% c("20130601_20151101","2015-11_2018-01","2015-11_2020-01")),
			 aes(fill=species, x = dates, y = `50%`, group=plotID)) +
	geom_line(data = ~filter(.x, time_period=="20130601_20151101"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="20130601_20151101"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="red") +
	geom_ribbon(data = ~filter(.x, time_period=="20130601_20151101"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="red") +
	geom_line(data = ~filter(.x, time_period=="2015-11_2018-01"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2018-01"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="blue") +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2018-01"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="blue") +
	geom_line(data = ~filter(.x, time_period=="2015-11_2020-01"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2020-01"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="green") +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2020-01"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="green") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
	facet_grid(~plotID)



ggplot(plot_est %>% filter(plotID %in% c("HARV_004","HARV_013"),
																taxon %in% c("ectomycorrhizal") &
																	model_name == "all_covariates" &
																	time_period %in% c("20130601_20151101","2015-11_2018-01","2015-11_2020-01")),
			 aes(fill=species, x = dates, y = `50%`, group=plotID)) +
	geom_line(data = ~filter(.x, time_period=="20130601_20151101"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="20130601_20151101"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="red") +
	geom_ribbon(data = ~filter(.x, time_period=="20130601_20151101"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="red") +
	geom_line(data = ~filter(.x, time_period=="2015-11_2018-01"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2018-01"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="blue") +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2018-01"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="blue") +
	geom_line(data = ~filter(.x, time_period=="2015-11_2020-01"),
						aes(x = dates, y =  `50%`), show.legend = F) +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2020-01"),
							aes(x = dates,ymin=`2.5%`, ymax = `97.5%`), alpha=0.1, fill="green") +
	geom_ribbon(data = ~filter(.x, time_period=="2015-11_2020-01"),
							aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=0.2, fill="green") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
	facet_grid(~plotID)



# Get best plots for examples (most calibration AND validation points)
legacy_not_na <- plot_est %>% filter(time_period == "20130601_20151101" & !is.na(truth))
not_na <- plot_est %>% filter(!is.na(truth) & plotID %in% legacy_not_na$plotID)

plot_hindcast_counts <- sort(table(not_na$plotID))



cycl_calibration = seas_vals[[2]]
all_cov_calibration = seas_vals[[1]]
ggplot(cycl_calibration) + geom_point(aes(y = rank, x = max)) + xlab("Month in which peak seasonal trend is observed")
ggplot(all_cov_calibration) + geom_point(aes(y = rank, x = max)) + xlab("Month in which peak residual seasonal trend is observed")

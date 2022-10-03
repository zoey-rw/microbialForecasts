
# Create forecasts for Shannon diversity
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("./functions/prepDiversityData.r")
source("./functions/forecastDiversity.r")

pacman::p_load(Rfast, data.table, doParallel)

summaries <- readRDS("./data/summary/div_summaries.rds")

# for testing
model_name = "cycl_only"
model_name = "all_covariates"
group = "ITS"
group = "16S"

#Run for multiple groups at once, in parallel (via PSOCK)
cl <- makeCluster(2, type="FORK", outfile="")
registerDoParallel(cl)

# Loop through each group
output.list = foreach(group = c("ITS","16S"), 
											.errorhandling = 'pass') %dopar% {															
												message("Beginning forecast loop for: ", group)
												time_period = "2015-11_2018-01"
												#for (time_period in c("2015-11_2018-01", "2016-01_2020-01")){
												message("Calibration time period: ", time_period)
												# Read in observations and covariate data
												div_in = switch(group,
																				"ITS" = readRDS("./data/clean/alpha_div_ITS.rds"),
																				"16S" = readRDS("./data/clean/alpha_div_16S.rds"))
												
												if (time_period %in% c("2015-11_2018-01", "2016-01_2020-01")) {
													rank.df = div_in$recent # since values are center-scaled 
													if (time_period == "2015-11_2018-01") {
														min.date = "20151101"
														max.date = "20180101"		
													} else if (time_period == "2016-01_2020-01") {
														min.date = "20160101"
														max.date = "20200101" }
												} else {
													rank.df = div_in$full
													min.date = "20130601" }
												# Get forecasting time series
												model.inputs <- prepDivData(rank.df = rank.df, min.prev = 3, 
																										min.date = min.date, max.date = "20200101", 
																										full_timeseries = T)
												
												# Loop through linear model versions
												model_output_list <- list()
												for (model_name in c("all_covariates", "cycl_only")){
													message("Forecasting with model: ", model_name)
													
													# Filter model estimates for each plot abundance
													plot_summary <- summaries$plot_est %>% filter(time_period == !!time_period & 
																																					model_name == !!model_name &
																																					group == !!group)
													
													# Get actual model MCMC samples & inputs
													f <- paste0("./data/model_outputs/diversity/", model_name, "/samples_div_", group,"_", min.date, "_", max.date, ".rds")
													read_in <- readRDS(f)
													param_samples <- as.data.frame(as.matrix(read_in$samples))
													truth.plot.long <- model.dat <- read_in$metadata$model_data
													plot_site_key <- model.dat %>% select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
													site_list <- unique(plot_site_key$siteID)
													
													# Use new model inputs for full date, site, and plot keys
													new_plot_site_key <- model.inputs$truth.plot.long %>% 
														select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% 
														distinct() %>% filter(!siteID %in% plot_site_key$siteID)
													new_site_list <- unique(new_plot_site_key$siteID)
													
													# Forecast at both observed and unobserved sites
													full_site_list <- c(site_list, new_site_list)
													#full_site_list <- site_list[1:2] #testing
													site_output_list <- list()
													for (siteID in full_site_list){
														message("SiteID: ", siteID)
														newsite <- siteID %in% new_plot_site_key$siteID
														plot_key <- if (newsite) new_plot_site_key else plot_site_key
														plot_key <- plot_key %>% filter(siteID == !!siteID)
														plot_list <- unique(plot_key$plotID)
														
														plot_output_list <- list()
														for (plotID in plot_list){
															message("PlotID: ", plotID)
															#go for it!!!
															hindcast.plot <- diversity_fcast(plotID, model.inputs, 
																															 param_samples, group,
																															 truth.plot.long, 
																															 Nmc = 5000,
																															 plot_summary) %>% 
																mutate(model_name = !!model_name,
																			 time_period = !!time_period)
															plot_output_list[[plotID]] <- hindcast.plot
														}
														site_output_list[[siteID]] <- rbindlist(plot_output_list)	
													}
													model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)	
												}
												group_output <- rbindlist(model_output_list)	
												return(group_output)
											}				
out <- rbindlist(output.list)	
out$fcast_period <- ifelse(out$dates < "2018-01-01", "calibration", "hindcast")
out$fcast_type <- "Diversity"
saveRDS(out, "./data/summary/hindcast_div.rds")

out <- readRDS("./data/summary/hindcast_div.rds")



#### Fr
# Merging in some truth values that weren't added correctly. Not ideal but it works.
group <- "16S"
div_in <- readRDS("./data/clean/alpha_div_16S.rds")
rank.df = div_in$full %>% filter(!siteID %in% c("ABBY","LAJA"))
model.inputs_16S <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101", full_timeseries = T)
truth_16S <- model.inputs_16S$truth.plot.long  %>% select(-c(plot_num, site_num, timepoint))
truth_16S$scenario = "full_uncertainty_16S"

group <- "ITS"
div_in <- readRDS("./data/clean/alpha_div_ITS.rds")
rank.df = div_in$full %>% filter(!siteID %in% c("ABBY","LAJA"))
model.inputs_ITS <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101", full_timeseries = T)
truth_ITS <- model.inputs_ITS$truth.plot.long %>% select(-c(plot_num, site_num, timepoint))
truth_ITS$scenario = "full_uncertainty_ITS"

out <- out %>% select(-c(Shannon_orig, Shannon_scale_site, truth))
truth <- rbind(truth_16S, truth_ITS)
out_fixed <- left_join(out, truth)
saveRDS(out_fixed, "./data/summary/hindcast_div.rds")

# View example output
ggplot(out %>% filter(plotID=="BART_002")) + 
	facet_grid(#rows=vars(taxon), 
		cols = vars(model_name), drop=T, scales="free") +
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
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ylim(c(-2,2))

old_hindcast <- readRDS("./data/summary/hindcast_div_legacy_incl.rds")
# View example output
ggplot(old_hindcast %>% filter(plotID=="BART_002" & group=="ITS" & model_name == "all_covariates")) + 
	facet_grid(#rows=vars(taxon), 
		cols = vars(model_name), drop=T, scales="free") +
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
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ylim(c(-2,2))


ggplot(hindcast.plot) + 
	facet_grid(#rows=vars(taxon), 
		cols = vars(model_name), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') + 
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01")))# + ylim(c(-2,2))

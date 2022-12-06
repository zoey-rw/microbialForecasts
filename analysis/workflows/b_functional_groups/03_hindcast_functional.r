
# Create forecasts for functional groups, using structure from SOBOL code
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

pacman::p_load(readxl, rjags, Rfast, moments, scales, data.table, doParallel)

# Read in microbial abundances
cal <- c(readRDS("./data/clean/cal_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("./data/clean/val_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/val_groupAbundances_ITS_2021.rds"))

# Read in model outputs
summaries <- readRDS(here("data/summary/bychain_beta_fg_summaries_20151101_20180101.rds"))

# Read in predicted site effects
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))

# R-ead in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS("./data/clean/all_predictor_data.rds")


# Loop through each group
fcast_ranks = microbialForecast:::keep_fg_names


# for testing
#testing = T
#if (testing) fcast_ranks = fcast_ranks[7:8]
#if (testing)
#fcast_ranks = tail(fcast_ranks, 28)
fcast_ranks = head(fcast_ranks, 12)
#i = 1
model_name = "all_covariates"
time_period = "2015-11_2018-01"

#Run for multiple groups at once, in parallel (via PSOCK)
cl <- makeCluster(28, type="FORK", outfile="")
registerDoParallel(cl)

output.list = foreach(i=1:length(fcast_ranks),
											#.errorhandling = 'pass',
											.verbose = T
											#outfile=""
											) %dopar% {

												pacman::p_load(microbialForecast)

												source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

	testing = F
	print(testing)
	time_period = "2015-11_2018-01"
	min.date = "20151101"
	max.date = "20180101"
	rank.name = fcast_ranks[[i]]

	#rank.name = "ectomycorrhizal"
	#rank.name = "oligotroph"
	message("Beginning forecast loop for: ", rank.name)
	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	# Prep validation data
	# model.inputs <- prepFunctionalData(rank.df = rank.df, min.prev = 3,
	# 																	 min.date = "20151101", max.date = "20200101",
	# 																	 full_timeseries = T)

	# Testing WITHOUT the prevalence filter
	model.inputs <- prepFunctionalData(rank.df = rank.df, min.prev = 0,
																		 min.date = "20151101", max.date = "20200101",
																		 full_timeseries = T)


	model_output_list <- list()
	#for (model_name in c("all_covariates", "cycl_only")){
		for (model_name in c("all_covariates")){
			message("Forecasting with model: ", model_name)

		# Filter model estimates for each plot abundance
			plot_summary <- summaries$plot_est %>%
			filter(time_period == !!time_period &
						 	model_name == !!model_name &
						 	species == !!rank.name)

			model.inputs$truth.plot.long <- model.inputs$truth.plot.long %>%
				mutate(species = ifelse(species=="other", paste0(species, "_", !!rank.name), species))


		# Get model outputs
		f <- paste0("./data/model_outputs/functional_groups/", model_name, "/beta_samples_", rank.name,"_", min.date, "_", max.date, ".rds")

		if (!file.exists(f)) next

		read_in <- readRDS(f)
		param_samples <- as.data.frame(as.matrix(read_in$samples))
		truth.plot.long <- model.dat <- read_in$metadata$model_data
		plot_site_key <- model.dat %>%
			select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
			distinct()
		site_list <- unique(plot_site_key$siteID)

		# Use new model inputs for full date, site, and plot keys
		new_plot_site_key <- model.inputs$truth.plot.long %>%
			select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
			distinct() %>%
			filter(!siteID %in% plot_site_key$siteID)
		new_site_list <- unique(new_plot_site_key$siteID)

		# Forecast at both observed and unobserved sites
		full_site_list <- c(site_list, new_site_list)
		site_output_list <- list()

		#if (testing) full_site_list = full_site_list[1:2]
		siteID = "HEAL"

		for (siteID in full_site_list){
			message("SiteID: ", siteID)

			if (siteID %in% all_predictors$site_skip) next()

			newsite <- siteID %in% new_plot_site_key$siteID
			plot_key <- if (newsite) new_plot_site_key else plot_site_key
			plot_key <- plot_key %>% filter(siteID == !!siteID)
			plot_list <- unique(plot_key$plotID)

			plot_output_list <- list()

			if (testing) plot_list = head(plot_list, 3)

			plotID = plot_list[[1]]

			for (plotID in plot_list){
				message("PlotID: ", plotID)
				#go for it!!!

				# Forecast with random site effect
				hindcast.plot <- #microbialForecast::
					fg_fcast_beta(plotID,
				model.inputs,
				param_samples,
				truth.plot.long,
				plot_summary,
				Nmc = 1000,
				predict_site_effects = NULL)  %>% mutate(model_name = !!model_name,
																time_period = !!time_period,
																# species = !!rank.name,
																# taxon = !!rank.name,
																fcast_type = "Functional group", predicted_site_effect=F,
																newsite = !!newsite)

					# Forecast with estimated site effect, if available
				if (newsite){
					if (!siteID %in% unobs_sites$siteID) next()
				hindcast.plot_pred_site_eff <- #microbialForecast::
					fg_fcast_beta(plotID,
									model.inputs,
									param_samples,
									truth.plot.long,
									plot_summary,
									Nmc = 1000,
									predict_site_effects = unobs_sites)  %>%
					mutate(model_name = !!model_name,
									time_period = !!time_period,
									# species = !!rank.name,
									# taxon = !!rank.name,
									fcast_type = "Functional group", predicted_site_effect=T,
								 newsite = !!newsite)
				hindcast.plot <- rbind(hindcast.plot,
															 hindcast.plot_pred_site_eff)
				}

				plot_output_list[[plotID]] <- hindcast.plot
			}
			site_output_list[[siteID]] <- rbindlist(plot_output_list, fill = T)
		}
		model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)
	}
	rank_output <- rbindlist(model_output_list)
	#rank_output_list[[rank.name]] = rank_output
	saveRDS(rank_output, paste0("./data/summary/beta_hindcast_fg_", rank.name, "_", time_period, ".rds"))

	return(rank_output)
											}

stopCluster(cl)

# table(output.list$taxon_name)
# output.list

file.list <- list.files(path = "./data/summary/",
												recursive = T,
												pattern = "beta_hindcast_fg_[a-z]",
												full.names = T)

#file_path = file.list[2:41]

hindcast_list <- purrr::map(file.list, readRDS)
out <- rbindlist(hindcast_list,fill=TRUE)	%>%
	mutate(fcast_period = ifelse(dates < "2018-01-01", "calibration", "hindcast"),
	category = assign_fg_categories(taxon),
	group = assign_fg_kingdoms(category))
saveRDS(out, paste0("./data/summary/beta_hindcast_fg_2015-11_2018-01.rds"))
#
#
#
# out <- readRDS(here("data", "summary/beta_hindcast_fg_2015-11_2018-01.rds"))
#
#
# to_plot <- hindcast.plot %>% filter(!grepl("other", taxon))
# # View example output
# ggplot(to_plot %>% filter(plotID=="YELL_002")) +
# 	facet_grid(rows=vars(plotID),
# 						 cols = vars(predicted_site_effect), drop=T, scales="free") +
# 	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
# 	# geom_line(aes(x = dates, y = `50%`), show.legend = F) +
# 	# #geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
# 	# geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`), fill=2, alpha=0.2, na.rm=T) +
# 	# geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`), fill=2, alpha=0.5, na.rm=T) +
# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), fill=4, alpha=0.2, na.rm=T) +
# 	geom_ribbon(aes(x = dates, ymin = lo_25, ymax = lo_75), fill=4, alpha=0.5, na.rm=T) +
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) +
# 	labs(fill='') #+ ylim(c(0,.2))


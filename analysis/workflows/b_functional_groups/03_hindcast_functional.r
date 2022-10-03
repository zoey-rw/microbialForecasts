
# Create forecasts for functional groups, using structure from SOBOL code
#source("./source.R")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

pacman::p_load(readxl, rjags, Rfast, moments, scales, data.table, doParallel)

# Set some global parameters
# Nmc_AB <- 5000 # Number of samples for subset
# Nmc <- 5000
# N.beta = 8
# k = 1

# Read in microbial abundances
cal <- c(readRDS("./data/clean/cal_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("./data/clean/val_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/val_groupAbundances_ITS_2021.rds"))

#summaries <- readRDS("./data/summary/fg_summaries.rds")
summaries <- readRDS(here("data", paste0("summary/beta_fg_summaries_20151101_20180101.rds")))


# # for testing
# model_name = "cycl_only"
# model_name = "all_covariates"
# scenario = "full_uncertainty"
# rank.name = microbialForecast:::keep_fg_names[37]
# rank.name = microbialForecast:::keep_fg_names[12]

#Run for multiple groups at once, in parallel (via PSOCK)
cl <- makeCluster(28, type="FORK", outfile="")
registerDoParallel(cl)


testing = F
# Loop through each group
fcast_ranks = microbialForecast:::keep_fg_names

if (testing) fcast_ranks = fcast_ranks[7:8]

output.list = foreach(i=1:length(fcast_ranks),
											.errorhandling = 'pass') %dopar% {

												pacman::p_load(microbialForecast)


	testing = T
	time_period = "2015-11_2018-01"
	min.date = "20151101"
	max.date = "20180101"
	rank.name = fcast_ranks[[i]]

	message("Beginning forecast loop for: ", rank.name)
	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	# Prep validation data
	model.inputs <- prepFunctionalData(rank.df = rank.df, min.prev = 3,
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

		for (siteID in full_site_list){
			message("SiteID: ", siteID)
			newsite <- siteID %in% new_plot_site_key$siteID
			plot_key <- if (newsite) new_plot_site_key else plot_site_key
			plot_key <- plot_key %>% filter(siteID == !!siteID)
			plot_list <- unique(plot_key$plotID)

			plot_output_list <- list()

			#if (testing) plot_list = plot_list[1:2]

			for (plotID in plot_list){
				message("PlotID: ", plotID)
				#go for it!!!
				hindcast.plot <- microbialForecast::fg_fcast_beta(plotID,
				model.inputs,
				param_samples,
				truth.plot.long,
				plot_summary,
				Nmc = 1000)  %>% mutate(model_name = !!model_name,
																time_period = !!time_period,
																# species = !!rank.name,
																# taxon = !!rank.name,
																fcast_type = "Functional group")

				plot_output_list[[plotID]] <- hindcast.plot
			}
			site_output_list[[siteID]] <- rbindlist(plot_output_list)
		}
		model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)
	}
	rank_output <- rbindlist(model_output_list)
	#rank_output_list[[rank.name]] = rank_output
	return(rank_output)
											}
output.list

out <- rbindlist(output.list)	%>%
	mutate(fcast_period = ifelse(dates < "2018-01-01", "calibration", "hindcast"),
	category = assign_fg_categories(taxon),
	group = assign_fg_kingdoms(category))
saveRDS(out, paste0("./data/summary/beta_hindcast_fg_", time_period, ".rds"))

out <- readRDS("./data/summary/beta_hindcast_fg.rds")

# View example output
ggplot(out %>% filter(plotID=="BART_002" & taxon == "oligotroph")) +
	facet_grid(rows=vars(plotID),
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
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')


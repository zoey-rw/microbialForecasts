#### Create forecasts for individual taxonomic groups

#### Commands for running as a batch job through OGE cluster ####

# Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else k = 1

#### Reading in files ####

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Read in microbial abundances
cal <- c(readRDS("./data/clean/cal_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("./data/clean/val_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/val_groupAbundances_ITS_2021.rds"))

# Read in model outputs to grab parameter estimates
single_tax_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_201801.rds")

# Read in predicted site effects
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))

# Read in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS("./data/clean/all_predictor_data.rds")

# Read in list of taxa that passed convergence with simpler model.
keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")


#### Set some parameters for testing ####
# #k = 1
# #model_name = "cycl_only"
model_name = "all_covariates"
time_period = "2015-11_2018-01"
min.date = "20151101"
max.date = "20180101"
testing = F
#taxon = keep_names[1]


#rank_output_list = foreach(k=1:10, .errorhandling = 'pass') %do% {
#run_single_taxon_fcast <- function(k = 1) {

	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

	rank.name <- microbialForecast:::tax_names[[k]]
	message("Beginning forecast loop for: ", rank.name)

	# Filter model estimates for each plot abundance
	rank_plot_summary <- single_tax_summaries$plot_est %>%
		filter(time_period == !!time_period &
					 	rank_name == !!rank.name)

	#rm(single_tax_summaries); gc()

	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	# Grab names of taxa to keep
	keep_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]

	#Run for multiple ranks, in parallel
	cl <- makeCluster(27, type="FORK", outfile="")
	registerDoParallel(cl)

	tax_output_list = foreach(taxon=keep_names, .errorhandling = 'remove') %dopar% {

		message("Forecasting for taxon: ", taxon)

		keep_vec <- c(taxon, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")

		# Prep validation data
		model.inputs <- prepTaxonomicData(rank.df = rank.df,
																			min.prev = 0,
																			min.date = min.date,
																			max.date = "20200101",
																			full_timeseries = T,
																			keep_vec = keep_vec)

		# library(profvis)
		# profvis({
		model_output_list <- list()
		for (model_name in c("all_covariates", "cycl_only")){


			message("Forecasting with model: ", model_name)
			#rank_taxon = paste0(rank_name, "_", )
			# Filter model estimates for each plot abundance
			plot_summary <- rank_plot_summary %>%
				#filter(taxon == !!taxon)
				filter(rank_name == !!rank.name & grepl(!!taxon, rank) &
							 model_name == !!model_name)

			# Get model outputs
			f <- file.path("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/single_taxon", model_name,  paste0("samples_", rank.name, "_", taxon, "_", min.date, "_", max.date, ".rds"))

			if(!file.exists(f)) next()

			read_in <- readRDS(f)

			# trim.jags <- mcmc.list(read_in$samples[[1]], read_in$samples[[2]])
			# param_samples <- as.data.frame(as.matrix(trim.jags))
			#
			param_samples <- as.data.frame(as.matrix(read_in$samples))
			#truth.plot.long <- model.dat <- read_in$metadata$model_data


			rank.df_spec <- rank.df %>%
				select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!taxon)
			rank.df_spec$other <- 1-rank.df_spec[[taxon]]
			new_model_inputs <- prepTaxonomicData(rank.df = rank.df_spec,
																				min.prev = 3,
																				min.date = min.date,
																				max.date = "20180101")

			truth.plot.long <- model.dat <- new_model_inputs$truth.plot.long

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

			if (testing == T) {
				full_site_list <- c(head(site_list,2), head(new_site_list,2))
				siteID = new_site_list[[2]]
			} else {
				full_site_list <- c(site_list, new_site_list)
			}
			#full_site_list <- c(site_list[[1]], new_site_list[[1]])

			# Loop through both observed and unobserved sites
			site_output_list <- list()
			for (siteID in full_site_list) {
				message("SiteID: ", siteID)

				if (siteID %in% all_predictors$site_skip) next()

				#Sample covariate data
				# covar <- create_covariate_samples(model.inputs, siteID,
				# 																	Nmc_large, Nmc)

				newsite <- siteID %in% new_plot_site_key$siteID
				plot_key <- if (newsite) new_plot_site_key else plot_site_key
				plot_key <- plot_key %>% filter(siteID == !!siteID)
				plot_list <- unique(plot_key$plotID)

				plot_output_list <- list()

				if (testing == T) {
				plotID <- plot_list[[1]] #testing
				plot_list <- plot_list[[1]]
				}

				for (plotID in plot_list){
					message("PlotID: ", plotID)

					# Forecast with known or random site effect
					#go for it!!!
					hindcast.plot <- #microbialForecast::
						fg_fcast_beta(plotID,
													model.inputs,
													param_samples,
													truth.plot.long,
													plot_summary,
													Nmc = 1000,
													predict_site_effects = NULL,
													rank.name = rank.name)  %>% mutate(model_name = !!model_name,
																																	 time_period = !!time_period,
																																	 species = !!taxon,
																																	 rank_name = !!rank.name,
																																	 fcast_type = "Taxonomic group",
																																	 predicted_site_effect=F,
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
														predict_site_effects = unobs_sites,
														rank.name = rank.name)  %>%
							mutate(model_name = !!model_name,
										 time_period = !!time_period,
										 # species = !!rank.name,
										 # taxon = !!rank.name,
										 species = !!taxon,
										 rank_name = !!rank.name,
										 fcast_type = "Taxonomic group",
										 predicted_site_effect=T,
										 newsite = !!newsite)
						hindcast.plot <- rbind(hindcast.plot,
																	 hindcast.plot_pred_site_eff)
					}

	print(tail(hindcast.plot, 1))
					plot_output_list[[plotID]] <- hindcast.plot
				}
				site_output_list[[siteID]] <- rbindlist(plot_output_list)
			}
			model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)
		}
			#tax_output_list[[taxon]] <- rbindlist(model_output_list)
			tax_output <- rbindlist(model_output_list, fill = T)
			return(tax_output)
		#})
	}
	rank_output <- rbindlist(tax_output_list, fill = T)
# 	return(rank_output)
# }

#rank_out <- run_single_taxon_fcast(k)
#rank_out <- rbindlist(rank_output_list, fill = T)
	rank_out <- rank_output
rank.name <- microbialForecast:::tax_names[[k]]

print("output:")
print(dim(rank_out))
out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_single_tax_", rank.name, ".rds")
saveRDS(rank_out, out.path)
#
# all_out <- rbindlist(rank_output_list, fill = T)
# all_out$group <- ifelse(grepl("_bac", all_out$rank, fixed = T), "16S", "ITS")
# saveRDS(all_out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_single_tax_test.rds")


stopCluster(cl)

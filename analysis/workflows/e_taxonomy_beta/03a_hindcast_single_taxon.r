# hindcast single taxa

#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else k = 1

# Create forecasts for individual taxonomic groups
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

pacman::p_load(Rfast, moments, scales)

# Read in microbial abundances
cal <- c(readRDS("./data/clean/cal_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("./data/clean/val_groupAbundances_16S_2021.rds"),
				 readRDS("./data/clean/val_groupAbundances_ITS_2021.rds"))

single_tax_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_201801.rds")

# Loop through each model

# Set some parameters for testing
# #k = 1
# #model_name = "cycl_only"
model_name = "all_covariates"
time_period = "2015-11_2018-01"
min.date = "20151101"
max.date = "20180101"

testing = F
keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")


#rank_output_list = foreach(k=1:10, .errorhandling = 'pass') %do% {
#run_single_taxon_fcast <- function(k = 1) {


	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

	rank.name <- microbialForecast:::tax_names[[k]]
	message("Beginning forecast loop for: ", rank.name)

	# Filter model estimates for each plot abundance
	rank_plot_summary <- single_tax_summaries$plot_est %>%
		filter(time_period == !!time_period &
					 	model_name == !!model_name &
					 	rank_name == !!rank.name)

	#rm(single_tax_summaries)

	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	keep_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]

	# Grab names of taxa to keep
	# keep_names <- keep_list[[rank.name]]$taxon.name
	#
	# print(keep_names)
	taxon = keep_names[1]
	# taxon = keep_names[2]
	# taxon = keep_names[4]

	#Run for multiple ranks, in parallel
	time_period = "2015-11_2018-01"

	cl <- makeCluster(28, type="FORK", outfile="")
	registerDoParallel(cl)

	tax_output_list = foreach(taxon=keep_names, .errorhandling = 'remove') %dopar% {
	#for (taxon in keep_names) {

		message("Forecasting for taxon: ", taxon)

		keep_vec <- c(taxon, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")

		# Prep validation data
		model.inputs <- prepTaxonomicData(rank.df = rank.df, min.prev = 3,
																			min.date = min.date,
																			max.date = "20200101",
																			full_timeseries = T,
																			keep_vec = keep_vec)

		# library(profvis)
		# profvis({
		model_output_list <- list()
		for (model_name in c("all_covariates", "cycl_only")){


			message("Forecasting with model: ", model_name)

			# Filter model estimates for each plot abundance
			plot_summary <- rank_plot_summary %>%
				filter(taxon == !!taxon)

			# Get model outputs
			f <- file.path("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/single_taxon", model_name,  paste0("samples_", rank.name, "_", taxon, "_", min.date, "_", max.date, ".rds"))

			if(!file.exists(f)) next()

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
			siteID = full_site_list[[1]]
			#full_site_list <- c(site_list[[1]], new_site_list[[1]])
			for (siteID in full_site_list) {
				message("SiteID: ", siteID)

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

					#go for it!!!
					hindcast.plot <- taxa_fcast(plotID,
																			model.inputs,
																			param_samples,
																			truth.plot.long,
																			plot_summary,
																			Nmc = 1000)
						hindcast.plot <- hindcast.plot %>% mutate(model_name = !!model_name,
																														 time_period = !!time_period,
																														 species = !!taxon,
																														 rank_name = !!rank.name,
																														 fcast_type = "Taxonomic")
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
	#stopCluster(cl)
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

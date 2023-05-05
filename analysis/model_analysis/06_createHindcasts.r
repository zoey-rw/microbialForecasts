#### Create forecasts for individual taxonomic groups

testing=F
#### Reading in files ####

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/run_hindcast.r")

# Read in microbial abundances
all.ranks <- c(readRDS("./data/clean/groupAbundances_16S_2023.rds"),
							 readRDS("./data/clean/groupAbundances_ITS_2023.rds"))

# Read in model outputs to grab parameter estimates
calibration_model_summaries <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))

# Read in predicted site effects
unobs_sites_env_cov <- readRDS(here("data/summary/site_effects_unobserved_env_cov.rds"))
unobs_sites_cycl_only <- readRDS(here("data/summary/site_effects_unobserved_cycl_only.rds")) %>% 
	mutate(model_name="cycl_only")
unobs_sites_env_cycl <- readRDS(here("data/summary/site_effects_unobserved_env_cycl.rds"))
pred_effects <- list(unobs_sites_env_cov, unobs_sites_cycl_only, unobs_sites_env_cycl) %>% rbindlist()

# Read in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS("./data/clean/all_predictor_data.rds")

# Read in list of taxa that passed convergence
keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))

model_id_list = unique(pred_effects$model_id)

#model_id_list = model_id_list[model_id_list %in% keep_list] %>% tail(100)

# model_id_list = c( "cycl_only_chitinolytic_20151101_20180101",  "cycl_only_lignolytic_20151101_20180101")
# model_id_list = c( "env_cov_actinobacteriota_20151101_20180101",  "env_cycl_acidobacteriota_20151101_20180101")
# model_id_list = keep_list
# k = 1
min.date = "20151101"

#Run for multiple ranks, in parallel
cl <- makeCluster(28, type="FORK", outfile="")
registerDoParallel(cl)

# for testing
#tax_output_list = foreach(k = 1:100, .errorhandling = 'pass') %dopar% {
tax_output_list = foreach(k = 1:length(model_id_list), .errorhandling = 'remove') %dopar% {
# tax_output_list = list()
# 	for(k in 1:length(model_id_list)){
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/run_hindcast.r")



	model_id=model_id_list[k]
		message("Beginning forecast loop for: ", model_id)
		parsed_id = parse_model_id(model_id)
		rank.name = parsed_id[[1]]
		time_period = parsed_id[[2]]
		taxon = species = parsed_id[[4]]
		fcast_type = parsed_id[[7]]
		model_name = parsed_id[[6]]

			# Filter model estimates for each plot abundance
		plot_summary <- rank_plot_summary <- calibration_model_summaries$plot_est %>%
				filter(model_id == !!model_id)

			keep_vec <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date",taxon)
			rank.df = all.ranks[[as.character(rank.name)]]
			rank.df_spec <- rank.df %>%
				select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!taxon)
			rank.df_spec$other <- 1-rank.df_spec[[taxon]]
			# Prep validation data using full time series
			full.ts.model.inputs <- prepBetaRegData(rank.df = rank.df,
																							min.prev = 0,
																							min.date = min.date,
																							max.date = "20200101",
																							full_timeseries = T,
																							keep_vec = keep_vec)

				message("Forecasting with model: ", model_name)

				# Get model outputs
				f <- here(file.path("data/model_outputs/logit_beta_regression/", model_name,  paste0("samples_", model_id, ".rds")))
				if(!file.exists(f)) next()
				read_in <- readRDS(f)
				param_samples <- as.data.frame(as.matrix(read_in$samples))
				truth.plot.long <- model.dat <- read_in$metadata$model_data

				plot_site_key <- model.dat %>%
					select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
					distinct()
				site_list <- unique(plot_site_key$siteID)

				# Use new model inputs for full date, site, and plot keys
				new_plot_site_key <- full.ts.model.inputs$truth.plot.long %>%
					select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
					distinct() %>%
					filter(!siteID %in% plot_site_key$siteID)
				new_site_list <- unique(new_plot_site_key$siteID)

				if (testing == T) {
					full_site_list <- c(head(site_list,5), head(new_site_list,5))
					#full_site_list <- new_site_list
					#siteID = site_list[[1]]
					#
					 siteID = new_site_list[[2]]
				} else {
					full_site_list <- c(site_list, new_site_list)
				}

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
						pred_rank = ifelse(fcast_type=="Functional", rank, taxon)

						# Forecast with known or random site effect
						hindcast.plot <-
							fcast_logit_beta(plotID,
															 full.ts.model.inputs,
															 param_samples,
															 truth.plot.long,
															 plot_summary,
															 Nmc = 1000,
															 predict_site_effects = NULL,
															 rank.name = pred_rank,
															 model_id=model_id)  %>% mutate(model_name = !!model_name,
															 																	 time_period = !!time_period,
															 																	 species = !!taxon,
															 																	 rank_name = !!rank.name,
															 																	 predicted_site_effect=F,
															 																	 newsite = !!newsite,
															 															 model_id=!!model_id)


						# Forecast with estimated site effect, if available
						if (siteID %in% pred_effects$siteID) {

							hindcast.plot_pred_site_eff <-
								fcast_logit_beta(plotID,
																 full.ts.model.inputs,
																 param_samples,
																 truth.plot.long,
																 plot_summary,
																 Nmc = 1000,
																 predict_site_effects = pred_effects,
																 rank.name = pred_rank,
																 model_id=model_id)  %>%
								mutate(model_name = !!model_name,
											 time_period = !!time_period,
											 species = !!taxon,
											 rank_name = !!rank.name,
											 predicted_site_effect=T,
											 newsite = !!newsite,
											 model_id=!!model_id)
							hindcast.plot <- rbindlist(list(hindcast.plot,
																		 hindcast.plot_pred_site_eff), fill=T)
						}

						print(tail(hindcast.plot, 1))
						plot_output_list[[plotID]] <- hindcast.plot
					}
					site_output_list[[siteID]] <- rbindlist(plot_output_list, fill=T)
				}
			tax_output <- rbindlist(site_output_list, fill=T)
			#tax_output_list[[taxon]] <- tax_output
			return(tax_output)
		}

		#rank_out <- run_single_taxon_fcast(k)
out.fcast <- rbindlist(tax_output_list, fill = T)

		message("output nrow: ", nrow(out.fcast))
		out.path <- here(paste0("data/summary/all_hindcasts_raw.rds"))
		saveRDS(out.fcast, out.path)
		stopCluster(cl)

#
# 	rank_output = out.fcast
#
# 	rank_output <- rank_output %>% arrange(model_name, plotID, dateID)
# 	rank_output$dates <- fixDate(rank_output$dateID)
# 	ggplot(rank_output %>% filter(plotID %in% c("BONA_001","BART_001"))
# 	) +
# 		#facet_grid(rows=vars(species), drop=T, scales="free") +
# 		geom_line(aes(x = dates, y = med), show.legend = F, na.rm = T) +
# 		geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue", na.rm = T) +
# 		geom_ribbon(aes(x = dates, ymin = lo_25, ymax = hi_75),fill="red", alpha=0.6, na.rm = T) +
# 		theme_bw()+
# 		scale_fill_brewer(palette = "Paired") +
# 		theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 					legend.position = "bottom",legend.title = element_text(NULL),
# 					plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 		facet_grid(model_id~plotID + predicted_site_effect) +
# 		geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
#
#
	# ggplot(tax_output %>% 
	# 			 	filter(#plotID %in% c("BART_001") &
	# 			 				 	model_id %in% c("env_cycl_acidobacteriales_20151101_20180101"))
	# ) +
	# 	#facet_grid(rows=vars(species), drop=T, scales="free") +
	# 	geom_line(aes(x = dates, y = med), show.legend = F, na.rm = T) +
	# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue", na.rm = T) +
	# 	geom_ribbon(aes(x = dates, ymin = lo_25, ymax = hi_75),fill="red", alpha=0.6, na.rm = T) +
	# 	theme_bw()+
	# 	scale_fill_brewer(palette = "Paired") +
	# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
	# 				legend.position = "bottom",legend.title = element_text(NULL),
	# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	# 	facet_grid(model_id+plotID~predicted_site_effect) +
	# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
	# 
	# 	
	# 	ggplot(out.fcast %>% filter(!is.na(crps) & newsite==T)) +
	# 		theme_bw()+
	# 		geom_boxplot(aes(x = predicted_site_effect, y = crps)) +
	# 		geom_jitter(aes(x = predicted_site_effect, y = crps)) 

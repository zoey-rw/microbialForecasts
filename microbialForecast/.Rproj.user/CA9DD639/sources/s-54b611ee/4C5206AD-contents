#' @title run_MCMC_dirichlet
#' @description run_MCMC_dirichlet
#' @export

# iter <- 1000
# burnin <- 500
# thin <- 1
# test = F
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"
# nchains=3
# model_name = "cycl_only"
# model_name = "all_covariates"
# thin = 10
# iter_per_chunk = 1000
# thin = 5
# init_iter = 5000
# chain_no = 1
# max.date = "20180101"
# min.date = "20151101"


#' @title run_MCMC_taxa
#' @description run_MCMC_taxa for entire rank (more than one taxon)
#'
#' @export
run_MCMC_dirichlet <- function(k = 1, iter = 1000, init_iter = 100000, iter_per_chunk = 10000, burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 s = "acidobacteriota",
										 model_name = "all_covariates",
										 #time_period = "calibration",
										 chain_no = 1,
										 max.date = "20180101",
										 min.date = "20151101",
										 ...) {
	#set.seed(chain_no)
	#pacman::p_load(reshape2, parallel)

	# Subset to one rank.
	rank.name <- microbialForecast:::tax_names[k]

	# Read in microbial abundances
	cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
	val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))

	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	# Grab names of taxa to keep
	keep_list <- readRDS("./data/summary/converged_taxa_list.rds")
	keep_names <- keep_list[[rank.name]]$taxon.name
	# keep_names <- c("acidobacteriota","actinobacteriota","bacteroidota","proteobacteria","verrucomicrobiota","firmicutes")
	# keep_names <- c("acidobacteriota","actinobacteriota","proteobacteria")
	keep_vec <- c(keep_names, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")


	# Prep model inputs/outputs.
	cat(paste0("\nPreparing model data for ", model_name, ", ", scenario, ", ", rank.name, "\n"))


	out.path <-
		paste0("./data/model_outputs/taxa/",
					 model_name, "/samples_", rank.name,"_",min.date,"_",max.date, ".rds")


		# Reduce size for testing
		if (test == T) {
			rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
			rank.df <- rank.df %>% filter(siteID %in% c("HARV","CPER"))
#			rank.df = rank.df[1:500,]
		}
		model.dat <- prepTaxonomicData(rank.df = rank.df,
																	 min.prev = 3,
																	 min.date = min.date,
																	 max.date = max.date,
																	 keep_vec = keep_vec)

		constants <- model.dat[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index",
														 "plot_num", "plot_site_num",
														 "N.plot", "N.spp", "N.core", "N.site", "N.date",
														 "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd",
														 "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]
	cat(paste0("\nCompleted model data for ", model_name, ", ", scenario, ", ", rank.name, "\n"))

	if (model_name == "cycl_only"){
		constants$N.beta = 2
		Nimble_model = nimbleModTaxa_cycl_only
		constants <- constants[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index",
														 "plot_num", "plot_site_num",
														 "N.plot", "N.spp", "N.core", "N.site", "N.date", "sin_mo", "cos_mo","N.beta")]
	} else if(model_name == "all_covariates"){
		constants$N.beta = 8
		Nimble_model = nimbleModTaxa
	} else message("Missing specification of Nimble model.")

	constants$omega <- 0.0001 * diag(constants$N.spp)
	constants$zeros = rep(0, constants$N.spp)
	if (constants$N.spp < 8) {
		constants$omega <- 0.0001 * diag(8)
		constants$zeros = rep(0, 8)
	}
	if (constants$N.spp < 8) {
		constants$omega <- 0.0001 * diag(8)
		constants$zeros = rep(0, 8)
	}


	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long)



	inits <- initsFun(constants, type = "tax")
	cat("\nInits created\n")

	#inits
	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants, data = list(y=model.dat$y),
												inits = inits)
	cat("\nModel created")

	# Compile model
	cModel <- compileNimble(Rmodel)
	cat("\nModel compiled")

#	nimbleOptions(multivariateNodesAsScalars = TRUE)
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
																								 "sig","intercept", #"omega",
																								 "rho"),
														monitors2 = c("plot_rel"), thin2 = 25,
														useConjugacy = T)

	cat("\nMCMC configured")

	#mcmcConf$removeSamplers(c('sigma', 'intercept', 'rho'))
	#mcmcConf$addSampler(target = c('sigma', 'intercept', 'rho'), type = 'RW_block')

	myMCMC <- buildMCMC(mcmcConf)
	cat("\nModel built")

	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	cat("\nMCMC compiled")

	out.path2 <- gsub(".rds", paste0("_chain", chain_no, ".rds"), out.path)
	#iter_per_chunk <- 50000

	compiled$run(niter=init_iter, thin=thin, thin2 = 25, nburnin = burnin)
	cat(paste0("\nInitial run finished for chain", chain_no))

	out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
					out.path2)

	continue <- check_continue(out.run, min_eff_size = 10)
	while (continue){
		cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
		compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
		out.run<-as.mcmc(as.matrix(compiled$mvSamples))
		out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
		saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
						out.path2)
		continue <- check_continue(out.run, min_eff_size = 10)
	}

	# 	continue <- check_continue(out.run)
	# 	if (continue){
	# 		cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
	# 		compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
	# 		out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	# 		out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	# 		saveRDS(list(out.run, out.run2), out.path2)
	#
	# 		continue <- check_continue(out.run)
	# 		if (continue){
	# 			cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
	# 			compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
	# 			out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	# 			out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	# 			saveRDS(list(out.run, out.run2), out.path2)
	# 		}
	# 	}
	# }


	message("Finished sampling for fit for run: ", model_name, ", ",
					scenario, ", ", rank.name, ", chain: ", chain_no, "\n"); gc()
	return("Success!")
	}

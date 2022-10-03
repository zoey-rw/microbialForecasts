#
# iter <- 1000
# burnin <- 500
# thin <- 1
# test = F
# k = 38
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"
# nchains=3
# model_name = "cycl_only"
# model_name = "all_covariates"
# time_period = "calibration"
# max.date = "20180101"
# min.date = "20151101"

#' @title run_MCMC_functional
#' @description run_MCMC_functional
#' @export

run_MCMC_functional <- function(k = 17,
																iter = 1000,
																burnin = 500,
																thin = 1,
																test = F,
																temporalDriverUncertainty = TRUE,
																spatialDriverUncertainty = TRUE,
																scenario = NULL,
																model_name = "cycl_only",
																time_period = "calibration",
																max.date = "20200101",
																min.date = "20160101",
																...) {
	#pacman::p_load(reshape2, parallel)
	source(here::here("source.R"))
	source(here::here("functions", "prepFunctionalData.r"))

	ranks.keep <- keep_fg_names
	rank.name <- ranks.keep[k]

	# Read in microbial abundances
	cal <- c(readRDS(here("data","/clean/cal_groupAbundances_16S_2021.rds")),
					 readRDS(here("data","/clean/cal_groupAbundances_ITS_2021.rds")))
	val <- c(readRDS(here("data","/clean/val_groupAbundances_16S_2021.rds")),
					 readRDS(here("data","/clean/val_groupAbundances_ITS_2021.rds")))

	# Only keep functional groups; subset to one rank.
	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	if (test == T) rank.df = rank.df[1:500,]

	# Prep model inputs/outputs.
	model.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3,
																	min.date = min.date,max.date = max.date)
	constants <- model.dat[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index",
													 "plot_num", "plot_site_num",
													 "N.plot", "N.spp", "N.core", "N.site", "N.date",
													 "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd",
													 "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]
	# for output.
	metadata <- list(rank.name = rank.name,
									 niter = iter,
									 nburnin = burnin,
									 thin = thin,
									 model_data = model.dat$truth.plot.long,
									 model_name = model_name)
	if (!is.null(scenario)) metadata$scenario = scenario

	if (model_name == "cycl_only"){
		constants$N.beta = 2
		Nimble_model = nimbleModFunctional_cycl_only
	} else if(model_name == "all_covariates"){
		constants$N.beta = 8
		Nimble_model = nimbleModFunctional_trunc
	} else message("Missing specification of Nimble model.")


	# Reduce size for testing
	if (test == T) {
		out.path <-
			paste0("./data/model_outputs/functional_groups/",
						 model_name, "/test_samples_", rank.name,"_",min.date,"_",max.date, ".rds")
	} else {

		out.path <-
			paste0("./data/model_outputs/functional_groups/",
						 model_name, "/samples_", rank.name,"_",min.date,"_",max.date, ".rds")
	}


	initsList <- initsFun(constants, type = "fg")

	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants,
												data = list(y=model.dat$y),
												inits = initsList) # Set initial values using function
	# Rmodel$plot_mu <- initsList$plot_mu
	# Rmodel$Ex <- initsList$Ex

	#	Rmodel$checkConjugacy()
	nimbleOptions(MCMCjointlySamplePredictiveBranches = FALSE)
	cModel <- compileNimble(Rmodel)
	# Configure & compile MCMC
	myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","rho","sig","core_sd","intercept",
																									"site_effect"),
											monitors2 = c("plot_mu"), thin2 = 25)
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)

	# Sample from MCMC and save, just in case
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin,
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	cat(paste0("Finished sampling for ", rank.name))

	out <- list(samples = samples, metadata = metadata)
	saveRDS(out, out.path, compress = FALSE)
	cat(paste0("Saved raw samples for ", rank.name))

	# Calculate summary and save output.
	param_summary <- fast.summary.mcmc(samples$samples)
	plot_summary <- fast.summary.mcmc(samples$samples2)
	cat(paste0("Summarized output for ", rank.name))

	out <- list(samples = samples$samples,
							param_summary = param_summary,
							metadata = metadata,
							plot_summary = plot_summary)
	saveRDS(out, out.path, compress = FALSE)

	message("\n Functional group samples saved for fit for run: ", rank.name,
					"\n Model: ", model_name,
					"\n Scenario: ", scenario)
	return("ok")
}

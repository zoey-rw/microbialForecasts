source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
iter <- 1000
burnin <- 500
thin <- 1

# iter <- 50000
# burnin <- 10000
# thin <- 2
test = T
test = F
n.chains = 3
group = "ITS"
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
model_name = "all_covariates"
model_name = "cycl_only"
min.date = "20151101"
max.date = "20180101"

run_MCMC <- function(group = "ITS", 	
										 iter = 1000,
										 burnin = 500,
										 thin = 1,
										 test = F,
										 n.chains = 3,
										 temporalDriverUncertainty = TRUE,
										 spatialDriverUncertainty = TRUE,
										 scenario = NULL,
										 min.date = "20160101",
										 max.date = "20180101",
										 model_name = "cycl_only",
										 time_period = "calibration") {
	message("\nBeginning model fit for run: ", group, 
					"\nModel: ", model_name)
	
	pacman::p_load(reshape2, parallel, lubridate, nimble, coda, tidyverse) 
	
	source("./source.R")
	source("./functions/prepDiversityData.r")
	
	div_in = switch(group,
									"ITS" = readRDS("./data/clean/alpha_div_ITS.rds"),
									"16S" = readRDS("./data/clean/alpha_div_16S.rds"))
	if (min.date == "20151101") {
		rank.df = div_in$recent # since values are center-scaled 
	} else rank.df = div_in$full
	

		out.path <- 
			paste0("./data/model_outputs/diversity/", 
						 model_name, "/samples_div_", group,"_",min.date,"_",max.date, ".rds")
		
		# Reduce size if testing
		if (test == T){ 
			rank.df = rank.df %>% filter(siteID %in% c("BART","HARV","WREF"))
			out.path <- 
				paste0("./data/model_outputs/diversity/", 
							 model_name, "/test_div_", group,"_",min.date,"_",max.date, ".rds")
		}
		
		model.dat <- prepDivData(rank.df = rank.df, min.prev = 3, min.date = "20151101",max.date = "20180101")
		constants <- model.dat[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index", 
																 "plot_num", "plot_site_num", 
																 "N.plot", "N.spp", "N.core", "N.site", "N.date",
																 "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd", 
																 "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]
		
		# for output.
		metadata <- list(niter = iter,
										 nburnin = burnin,
										 thin = thin,
										 model_data = model.dat$truth.plot.long)
		if (!is.null(scenario)) metadata$scenario = scenario
		
	if (model_name == "cycl_only"){
		constants$N.beta = 2
		Nimble_model = nimbleMod_shannon_cycl_only
	} else if(model_name == "all_covariates"){
		constants$N.beta = 8
		Nimble_model = nimbleMod_shannon
	} else message("Missing specification of Nimble model.")
	
	# Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants, data = list(y=model.dat$y),
												inits = initsFun(constants))
	
	# Remove large objects due to data usage
	#rm(model.dat, div_in, rank.df, constants); gc()
	
	cModel <- compileNimble(Rmodel)
	
	# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
	# nimbleOptions(MCMCsaveHistory = TRUE)
	nimbleOptions(multivariateNodesAsScalars = TRUE)
	
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect","intercept", 
																								 "rho","sig","core_sd"),
														monitors2 = c("plot_mu"))
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	samples.out <- runMCMC(compiled, niter = iter,
												 nchains = n.chains, nburnin = burnin,
												 samplesAsCodaMCMC = T, thin = thin, thin2 = 20)
	
	message("\nFinished sampling for run: ", group, 
					"\nModel: ", model_name)
	
	# Save just in case something goes wrong during the model summary
	out <- list(samples = samples.out, metadata = metadata)
	saveRDS(out, out.path)
	message("\nDiversity raw samples saved for fit for run: ", group, 
						 "\nModel: ", model_name)
	
	# Calculate summary and save this output instead
	param_summary <- fast.summary.mcmc(samples.out$samples)
	plot_summary <- fast.summary.mcmc(samples.out$samples2)
	out <- list(samples = samples.out$samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	saveRDS(out, out.path)
	
	message("Diversity output saved for fit for run: ", group, 
						 "\nModel: ", model_name)
	return("Success!")
}

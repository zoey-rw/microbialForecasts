
iter <- 1000
burnin <- 500
thin <- 1

iter <- 50000
burnin <- 10000
thin <- 2
test = T
test = F
n.chains = 3
group = "ITS"
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario = "full_uncertainty_ITS"
model_name = "cycl_only"
model_name = "all_covariates"
time_period = "refit"

run_MCMC <- function(group = "ITS", 	
										 iter = 1000,
										 burnin = 500,
										 thin = 1,
										 test = F,
										 n.chains = 3,
										 temporalDriverUncertainty = TRUE,
										 spatialDriverUncertainty = TRUE,
										 scenario = NULL,
										 model_name = "cycl_only",
										 time_period = "calibration") {
	
	pacman::p_load(reshape2, parallel, lubridate, nimble, coda, tidyverse) 
	
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/cyclical_model_source.r")
	
	
	div_in = switch(group,
									"ITS" = 
										readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds"),
									"16S" = 
										readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds"))
	
	
	
	if (time_period == "refit"){
		
		out.path <- 
			paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/refit", 
						 model_name, "/samples_div_", scenario, ".rds")
		
		rank.df = div_in$full
		rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),] #Remove sites missing key covariates

		if (test == T) rank.df = rank.df[1:500,]
		
		
		# Custom function for organizing model data.
		# model.cal.dat <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20170101")
		# # Only keep plots that were in calibration data
		# model.dat <- prepDivData(rank.df = rank.df[rank.df$plotID %in% model.cal.dat$plotID,], min.prev = 3, max.date = "20200101")
		
		
		model.dat <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101")
		
	} else if (time_period == "calibration") {
		
		rank.df = div_in$cal
		
		if (test == T) rank.df = rank.df[1:500,]
		
		
		# Custom function for organizing model data.
		model.dat <- prepDivData(rank.df = rank.df, min.prev = 3)
	} else cat("Missing specification of time period.")
	
	# Reduce size for testing
	if (test == T) {
		out.path <- 
			paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/", 
						 model_name,"/test_",time_period, "_samples_div_", scenario, ".rds")
	} else {
		out.path <- 
			paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/", 
						 model_name,"/", time_period, "_samples_div_", scenario, ".rds")
		
	}
	
	
	if (model_name == "cycl_only"){
		n.beta = 2
		Nimble_model = nimbleMod_shannon_cycl_only
	} else if(model_name == "all_covariates"){
		n.beta = 8
		Nimble_model = nimbleMod_shannon
	} else cat("Missing specification of Nimble model.")
	

	
	constants <- list(N.plot =  length(unique(model.dat$plotID)), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(unique(model.dat$siteID)),
										timepoint = model.dat$timepoint,
										mois = model.dat[["mois"]],
										temp = model.dat[["temp"]],
										mois_sd = model.dat[["mois_sd"]],
										temp_sd = model.dat[["temp_sd"]],
										pH = model.dat[["pH"]],
										pC = model.dat[["pC"]],
										pH_sd = model.dat[["pH_sd"]],
										pC_sd = model.dat[["pH_sd"]],
										relEM = model.dat[["relEM"]],
										nspp = model.dat[["nspp"]],
										LAI = model.dat[["LAI"]],
										rc_grass = model.dat[["rc_grass"]],
										sin_mo = model.dat[["y_sin"]],
										cos_mo = model.dat[["y_cos"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = n.beta)
	truth <- model.dat$truth.plot.long # for output.
	
	# Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants, data = list(y=model.dat$y),
												inits = initsFun(constants))
	
	# Remove large objects due to data usage
	rm(model.dat, div_in, rank.df, constants); gc()
	
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
	cat(paste0("\nFinished sampling for run: ", group))
	
	metadata <- list(niter = iter,
									 nburnin = burnin,
									 thin = thin,
									 model_data = truth)
	
	out <- list(samples = samples.out, metadata = metadata)
	saveRDS(out, out.path)
	cat(paste0("\nDiversity raw samples saved for fit for run: ", group, 
						 "\n Model: ", model_name, 
						 "\n Scenario: ", scenario)) 
	
	# Calculate summary and save output.
	param_summary <- fast.summary.mcmc(samples.out$samples)
	plot_summary <- fast.summary.mcmc(samples.out$samples2)
	
	# Create outputs and clear compiled code
	
	out <- list(samples = samples.out$samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	saveRDS(out, out.path)
	
	cat(paste0("Diversity output saved for fit for run: ", group, 
						 "\n Model: ", model_name, 
						 "\n Scenario: ", scenario)) 
	return(out)
}


iter <- 1000
burnin <- 500
thin <- 1
test = F
k = 38
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario <- "full_uncertainty"
nchains=3
model_name = "cycl_only"
model_name = "all_covariates"
time_period = "calibration"

run_MCMC <- function(k = 17, 	
										 iter = 1000,
										 burnin = 500,
										 thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE,
										 spatialDriverUncertainty = TRUE,
										 scenario = NULL,
										 model_name = "cycl_only",
										 time_period = "calibration",
										 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepFunctionalData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/cyclical_model_source.r")
	
	
	ranks.keep <- keep_fg_names
	rank.name <- ranks.keep[k]
	
	if (time_period == "refit"){
		
		# Read in microbial abundances
		cal <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
						 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
		val <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds"), 
						 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds"))
		
		# Only keep functional groups; subset to one rank.
		cal.rank.df <- cal[[rank.name]] 
		val.rank.df <- val[[rank.name]] 
		rank.df <- rbind(cal.rank.df, val.rank.df)	
		rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),]
		
		if (test == T) rank.df = rank.df[1:500,]
		
		model.cal.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3, max.date = "20170101")
		# Only keep plots that were in calibration data
		model.dat <- prepFunctionalData(rank.df = rank.df[rank.df$plotID %in% model.cal.dat$plotID,], min.prev = 3, max.date = "20200101")
		
	} else if (time_period == "calibration") {
		
		# Read in microbial abundances
		d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
					 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
		
		# Only keep functional groups; subset to one rank.
		rank.df <- d[[rank.name]] 
		if (test == T) rank.df = rank.df[1:500,]
		
		# Prep model inputs/outputs.
		model.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3)
	
		} else cat("Missing specification of time period.")
	
	
	
	# Reduce size for testing
	if (test == T) {
		out.path <- 
			paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups/", 
						 model_name,"/test_",time_period, "_samples_", rank.name,"_", scenario, ".rds")
	} else {
		out.path <- 
			paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups/", 
						 model_name,"/", time_period,  "_samples_", rank.name,"_",scenario, ".rds")
	}
	
	mo <- month(as.Date(paste0(colnames(model.dat$mois), "01"), format="%Y%m%d"))
	y_sin = sin((2*pi*mo)/12)
	y_cos = cos((2*pi*mo)/12)
	
	if (model_name == "cycl_only"){
		n.beta = 2
		Nimble_model = nimbleModFunctional_cycl_only
	} else {
		n.beta = 8
		Nimble_model = nimbleModFunctional_trunc
	}
	
	# create data object
	constants <- list(N.plot =  length(model.dat[["plot_start"]]), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(model.dat$site_start),
										timepoint = model.dat$timepoint,
										sin_mo = y_sin,
										cos_mo = y_cos,
										mois    = model.dat[["mois"]],
										temp    = model.dat[["temp"]],
										mois_sd = model.dat[["mois_sd"]],
										temp_sd = model.dat[["temp_sd"]],
										pH      = model.dat[["pH"]],
										pC      = model.dat[["pC"]],
										pH_sd   = model.dat[["pH_sd"]],
										pC_sd   = model.dat[["pC_sd"]],
										nspp    = model.dat[["nspp"]],
										#										rc_grass = model.dat[["rc_grass"]],
										LAI = model.dat[["LAI"]],
										relEM = model.dat[["relEM"]],
										#										plotID = model.dat$plotID,
										#										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = n.beta)
	# for output.
	metadata <- list(c("rank.name" = rank.name, "niter" = iter,"nburnin" = burnin, "thin" = thin),
									 "model_data" = model.dat$truth.plot.long, "model_name" = model_name)
	
	initsList <- initsFun(constants, type = "fg")
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model, 
												#Rmodel <- nimbleModel(code = nimbleModFunctional, 
												constants = constants, 
												data = list(y=model.dat$y), 
												inits = initsList) # Set initial values using function
	# Rmodel$plot_mu <- initsList$plot_mu
	# Rmodel$Ex <- initsList$Ex
	
	# Remove large objects due to data usage
	rm(model.dat, d, rank.df, initsList, constants); gc()
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

	cat(paste0("\n Functional group samples saved for fit for run: ", rank.name, 
						 "\n Model: ", model_name, 
						 "\n Scenario: ", scenario,
						 "\n Time period: ", time_period)) 
	return("ok")
}

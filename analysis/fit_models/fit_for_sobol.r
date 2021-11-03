# Running shannon diversity model with different NA values
# 

# iter <- 1000
# burnin <- 500
# thin <- 1
# test = T
# n.chains = 3
# group = "ITS"
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE

run_MCMC <- function(group = "ITS", 	
										 iter = 1000,
										 burnin = 500,
										 thin = 1,
										 test = F,
										 n.chains = 3,
										 temporalDriverUncertainty = TRUE,
										 spatialDriverUncertainty = TRUE) {
	pacman::p_load(reshape2, parallel, lubridate, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
	
	div_in = switch(group,
									"ITS" = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds"),
									"16S" = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds"))
	
	rank.df = rbind(div_in$cal, div_in$val)
	
	# Leaving out two sites because they're missing data, plus Harvard Forest because it's our hindcast test site.
	rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA","HARV")),]
	rank.df$Shannon <- scale(rank.df$Shannon, scale = F)
	if (test == T) {
		rank.df = rank.df[1:500,]
	}
	# Custom function for organizing model data.
	model.dat <- prepDivData(rank.df = rank.df, min.prev = 3,max.date = "20200101")
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
										rc_grass = model.dat[["rc_grass"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 6)
	truth <- model.dat$truth.plot.long # for output.
	
	# Configure & compile model
	Rmodel <- nimbleModel(code = nimbleMod_shannon,
												constants = constants, data = list(y=model.dat$y),
												inits = initsFun(constants))
	
	# Remove large objects due to data usage
	rm(model.dat, div_in, rank.df, constants); gc()
	
	cModel <- compileNimble(Rmodel)
	
	# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
	# nimbleOptions(MCMCsaveHistory = TRUE)
	nimbleOptions(multivariateNodesAsScalars = TRUE)
	
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",#"intercept", 
																								 "rho","sig","core_sd"),
														monitors2 = c("plot_mu"))
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	samples.out <- runMCMC(compiled, niter = iter,
												 nchains = n.chains, nburnin = burnin,
												 samplesAsCodaMCMC = T, thin = thin)
	samples2 <- rm.NA.mcmc(samples.out$samples2)
	samples <- rm.NA.mcmc(samples.out$samples)
	# plot(samples)
	param_summary <- summary(samples)
	plot_summary <- summary(samples2)
	
	# Create outputs and clear compiled code
	metadata <- list(niter = iter,
									 nburnin = burnin,
									 thin = thin,
									 model_data = truth)
	
	cat(paste0("Diversity model fit for run: ", group)) 
	out <- list(samples = samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	return(out)
}

 
pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 

out <- run_MCMC(group = "16S", iter = 50000, burnin = 20000, thin = 3, 
									test=F, 
									temporalDriverUncertainty = TRUE, 
									spatialDriverUncertainty = TRUE)

saveRDS(out, paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/16S_sobol_calibration.rds"))

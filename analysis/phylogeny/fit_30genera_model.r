# Running dirichlet model on 30 soil genera (for phylogenetic analysis of environmental sensitivities)

# iter <- 1000
# burnin <- 500
# thin <- 1
# test = T
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE


run_MCMC <- function(iter = 1000,
										 burnin = 500,
										 thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE,
										 spatialDriverUncertainty = TRUE,
										 min_prev) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Set model run parameters
	# iter <- 5000
	# burnin <- 2000
	# thin <- 2

	
	# Read in dataset
	rank.df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_phylo_30tax.rds")
	
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:1000,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for bacteria (phylogeny)"))
	model.dat <- prepModelData(rank.df = rank.df, min.prev = min_prev)
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_phylo_min", min_prev, ".rds")
	print(paste0("Completed model data for phylo"))
	
	
	# Set up constants for NIMBLE model
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
										nspp = model.dat[["nspp"]],
										rc_grass = model.dat[["rc_grass"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,									
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 6,
										alpha0 =  rep(0, ncol(model.dat$y)))
	truth <- model.dat$truth.plot.long # for output.
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModLong,
												constants = constants, data = list(y=model.dat$y),
												inits = initsFun(constants))
	# Compile model
	cModel <- compileNimble(Rmodel)
	
	# Configure & compile MCMC
	myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","site_effect","intercept",
																									"rho"),
											monitors2 = c("plot_rel","SIGMA"), useConjugacy = F)
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)

	# Remove large objects due to data usage
	rm(model.dat, d, prepModelData, rank.df)
	gc()
	
	# Sample from MCMC
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	samples2 <- rm.NA.mcmc(samples$samples2)
	samples <- rm.NA.mcmc(samples$samples)
	
	# Calculate summary and save output.
	param_summary <- summary(samples)
	plot_summary <- summary(samples2)
	metadata <- c("niter" = iter,
								"nburnin" = burnin,
								"thin" = thin,
								"model_data" = truth)
	out <- list(samples = samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	saveRDS(out, out.path)
	return(out)
}


run_MCMC(min_prev = 5, iter = 100000, burnin = 80000, thin = 5, 
				 test=T,
				 temporalDriverUncertainty = T,
				 spatialDriverUncertainty = T)


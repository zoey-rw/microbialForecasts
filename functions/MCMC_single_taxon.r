
run_MCMC <- function(k = 1, iter = 1000,  burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 s = "acidobacteriota",
										 min.date = "20151101",
										 max.date = "20180101",
										 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source(here("functions","prepTaxonomicData.r"))
	
	
	# Subset to one rank.
	rank.name <- tax_names[k]
	
	# Read in microbial abundances
	cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
	val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")), 
					 readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))
	
	cal.rank.df <- cal[[rank.name]] 
	val.rank.df <- val[[rank.name]] 
	rank.df <- rbind(cal.rank.df, val.rank.df)	
	
	
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))
	
	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")] 
	rank.df_spec <- rank.df %>% 
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s) 
	rank.df_spec$other <- 1-rank.df_spec[[s]]
	
	model.dat <- prepTaxonomicData(rank.df = rank.df_spec, 
																 min.prev = 3, 
																 min.date = min.date,
																 max.date = max.date)
	
	out.path <- here("data",paste0("model_outputs/convergence_testing/samples_", rank.name, "_", s, "_",min.date,"_",max.date, ".rds"))
	
	
	print(paste("Completed model data for", rank.name, ", group:", s))
	
	constants <- model.dat[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index", 
													 "plot_num", "plot_site_num", 
													 "N.plot", "N.spp", "N.core", "N.site", "N.date", "sin_mo", "cos_mo")]
	constants$N.beta = 2
	
	constants$omega <- 0.0001 * diag(constants$N.spp)
	constants$zeros = rep(0, constants$N.spp)
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
	
	inits <- inits[c("y", "plot_mu", "intercept", "sig", "beta", "rho", "sigma", 
									 "plot_rel", "site_effect")]
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModTaxa_cycl_only,
												constants = constants, data = list(y=model.dat$y),
												inits = inits)
	#Rmodel$Ex <- inits$Ex
	#Rmodel$site_effect <- rep(0, constants$N.site)
	#Rmodel$y <- inits$y
	# Compile model
	cModel <- compileNimble(Rmodel)
	nimbleOptions(multivariateNodesAsScalars = TRUE)
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
																								 "sig","intercept", 
																								 "rho"),
														monitors2 = c("plot_rel"), thin2 = 25,
														useConjugacy = T)
	
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	
	# Sample from MCMCv
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	
	cat(paste("Finished sampling for ", rank.name, ", group: ", s)); gc()
	
	# Calculate summary and save output.
	param_summary <- fast.summary.mcmc(samples$samples)
	plot_summary <- fast.summary.mcmc(samples$samples2)
	
	out <- list(samples = samples$samples,
							param_summary = param_summary,
							metadata = metadata,
							plot_summary = plot_summary,
							gelman = gelman.diag(samples$samples))
	saveRDS(out, out.path, compress = F)
	cat(paste("Saved summary output for ", rank.name, ", group: ", s, "to: \n", out.path))
	return("ok")
}

# Refit functional group models using full dataset (incl 2018)

iter <- 1000
burnin <- 500
thin <- 1
test = T
k = 7
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
nchains=3
scenario = "full_uncertainty"

run_MCMC <- function(k = 17, 	
														 iter = 1000,
														 burnin = 500,
														 thin = 1,
														 test = F,
														 temporalDriverUncertainty = TRUE,
														 spatialDriverUncertainty = TRUE,
														scenario = NULL,
														 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepFunctionalData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Read in microbial abundances
	cal <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
					 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	val <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds"), 
					 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds"))
	
	# Only keep functional groups; subset to one rank.
	rank.name <- keep_fg_names[k]
	cal.rank.df <- cal[[rank.name]] 
	val.rank.df <- val[[rank.name]] 
	rank.df <- rbind(cal.rank.df, val.rank.df)	
	rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),]
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_samples_", rank.name, ".rds")
	print(paste0("Completed model data for ", rank.name))
	
	model.cal.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3, max.date = "20170101")
	model.dat <- prepFunctionalData(rank.df = rank.df[rank.df$plotID %in% model.cal.dat$plotID,], min.prev = 3, max.date = "20200101")
	
	# create data object
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
										pC_sd = model.dat[["pC_sd"]],
										nspp = model.dat[["nspp"]],
										rc_grass = model.dat[["rc_grass"]],
										relEM = model.dat[["relEM"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 6)
	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long)
	
	initsList <- initsFun(constants, type = "fg")
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModFunctional_trunc, 
												constants = constants, 
												data = list(y=model.dat$y), 
												inits = initsList) # Set initial values using function
	
#	Rmodel$checkConjugacy()
	cModel <- compileNimble(Rmodel)
	# Configure & compile MCMC
	myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","rho","sig","intercept","core_sd",
																									"site_effect"), monitors2 = c("plot_mu"), thin2 = 20)
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)
	# Remove large objects due to data usage
	# rm(model.dat, d, prepModelData, rank.df); gc()
	# Sample from MCMC
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	out <- list(samples = samples, 
							metadata = metadata)
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
	
	saveRDS(out, out.path, compress = F)
	cat(paste0("Saved summarized output for ", rank.name))
	return()
}

#### Run on all groups ----

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Running functional group model using fill dataset
pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 


# Create parameters to pass	
n.groups <- length(keep_fg_names)
params = data.frame(group = 1:n.groups,
										scenario = rep("full_uncertainty", n.groups),
										temporalDriverUncertainty = rep(TRUE, n.groups),
										spatialDriverUncertainty =  rep(TRUE, n.groups), 
										ranks = keep_fg_names)

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	out <- run_MCMC(k = params$group[[j]], iter = 300000, burnin = 150000, thin = 10, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
									scenario = params$scenario[[j]])
	return(out)
}


library(doParallel)
cl <- makeCluster(10, type="PSOCK", outfile="")
registerDoParallel(cl)


# j <- 7
# out <- run_MCMC(k = params$group[[j]], iter = 5000, burnin = 1000, thin = 2, 
# 								test=F, 
# 								temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
# 								spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
# 								scenario = params$scenario[[j]])


# # # get missing scenarios for re-running
# missing_scenarios <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/analysis/missing_scenarios_refit_fg.csv")
# 
# to_rerun <- as.numeric(rownames(params[params$ranks %in% missing_scenarios$scenario,]))
# 
# output.list = foreach(j=to_rerun,
# 											.errorhandling = 'pass') %dopar% {
# 	run_scenarios(j)
# 	return()
# }

output.list = foreach(j=1:40,
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j)
return()
										}

stopCluster(cl)

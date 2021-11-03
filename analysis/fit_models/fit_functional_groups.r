#Get arguments from the command line
# argv <- commandArgs(TRUE)
# # Check if the command line is not empty and convert values to numerical values
# if (length(argv) > 0){
# 	k <- as.numeric( argv[1] )
# } else {
# 	k=1
# }


iter <- 1000
burnin <- 500
thin <- 1
test = T
k = 38
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario <- "full_uncertainty"
nchains=3

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
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Only keep functional groups; subset to one rank.
	ranks.keep <- keep_fg_names
	rank.name <- ranks.keep[k]
	rank.df <- d[[rank.name]] 
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:700,]
		out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/test_samples_", rank.name, "_", scenario,".rds")
	} else {
		out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario,".rds")
	}
	
	# Prep model inputs/outputs.
	model.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3)

	# create data object
	constants <- list(N.plot =  length(model.dat[["plot_start"]]), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(model.dat$site_start),
										timepoint = model.dat$timepoint,
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
										relEM = model.dat[["relEM"]],
#										plotID = model.dat$plotID,
#										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 6)
	# for output.
	metadata <- list(c("rank.name" = rank.name, "niter" = iter,"nburnin" = burnin, "thin" = thin),
									 "model_data" = model.dat$truth.plot.long)
	
	initsList <- initsFun(constants, type = "fg")
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModFunctional_trunc, 
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
	cat(paste0("Saved output for ", rank.name))
	return("ok")
}

#### Run on all groups ----


# Running dirichlet model 
pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Create parameters to pass	
n.groups <- length(keep_fg_names)
params = data.frame(group = rep(1:n.groups),
										scenario = c(rep("no_uncertainty", n.groups), 
																 rep("spatial_uncertainty", n.groups), 
																 rep("temporal_uncertainty", n.groups), 
																 rep("full_uncertainty", n.groups)),
										temporalDriverUncertainty = c(rep(FALSE, n.groups*2),
																									rep(TRUE, n.groups*2)),
										spatialDriverUncertainty = c(rep(FALSE, n.groups),
																								 rep(TRUE, n.groups),
																								 rep(FALSE, n.groups),
																								 rep(TRUE, n.groups)),
										ranks = rep(keep_fg_names, 4))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	#run_MCMC(k = params$group[[j]], iter = 60000, burnin = 20000, thin = 3, 
	run_MCMC(k = params$group[[j]], iter = 600000, burnin = 300000, thin = 20, 

																		test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
									scenario = params$scenario[[j]])
	return()
}



library(doParallel)
#cl <- makeCluster(16, type="PSOCK", outfile="")
cl <- makeCluster(8, type="PSOCK", outfile="")
registerDoParallel(cl)

output.list = foreach(j=120:160,
											#.export=c("run_scenarios","params","run_MCMC"),
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j)
	return()
											}
stopCluster(cl)






# # get missing scenarios for re-running
# missing_scenarios <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/analysis/missing_scenarios.csv")
# params$specific <- paste(params$ranks, params$scenario)
# missing_scenarios$specific <- paste(missing_scenarios$scenario, missing_scenarios$rank.name)
# to_rerun <- as.numeric(rownames(params[params$specific %in% missing_scenarios$specific,]))
# 
# fixed_fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction",
# 										"assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction",
# 										"nitrification", "denitrification", "chitinolytic", "lignolytic",
# 										"methanotroph", "copiotroph", "oligotroph", "glucose_simple", "glycine_simple",
# 										"streptomycin_antibiotic","fsfeso4_anaerobic", "ironcitrate_anaerobic",
# 										"light_stress", "endophyte", "plant_pathogen",
# 										"animal_pathogen", "ectomycorrhizal", "lichenized", "wood_saprotroph",
# 										"soil_saprotroph", "litter_saprotroph", "saprotroph")
# to_rerun <- as.numeric(rownames(params[params$scenario=="full_uncertainty" & !params$ranks %in% fixed_fg_names,]))
# output.list = foreach(j=to_rerun,
# 											.errorhandling = 'pass') %dopar% {
# 	run_scenarios(j)
# 	return()
# }

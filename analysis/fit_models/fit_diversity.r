# Running shannon diversity model with different NA values


iter <- 1000
burnin <- 500
thin <- 1
test = F
group = "16S"
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE

run_MCMC <- function(group = "ITS", 	
														 iter = 1000,
														 burnin = 500,
														 thin = 1,
														 test = F,
														 out.path = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min5.rds",
														 temporalDriverUncertainty = TRUE,
														 spatialDriverUncertainty = TRUE,
										 scenario = NULL,
										 ...
														 ) {

	
	pacman::p_load(reshape2, parallel, lubridate, nimble, coda, tidyverse) 
	
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
	
	if (group == "ITS"){
		div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
	} else if (group == "16S") {
			div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
	}
	
	rank.df = div_in$cal
	
	if (test == T) {
	rank.df = rank.df[1:500,]
	}
	
	# Custom function for organizing model data.
	model.dat <- prepDivData(rank.df = rank.df, min.prev = 3)
	
	# mo <- month(as.Date(paste0(colnames(model.dat$mois), "01"), format="%Y%m%d"))
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
										
										N.beta = 6)
										# y_sin = scale(sin((2*pi*mo)/12), center = F), # specific to diversity models
										# y_cos = scale(cos((2*pi*mo)/12), center = F), # specific to diversity models
										# N.beta = 8)
	truth <- model.dat$truth.plot.long # for output.
	
	# Configure & compile model
	Rmodel <- nimbleModel(code = nimbleMod_shannon,
												constants = constants, data = list(y=model.dat$y),
												inits = initsFun(constants))
	cModel <- compileNimble(Rmodel)
	
	# Configure & compile MCMC
	myMCMC <- buildMCMC(cModel,monitors = c("beta","sigma","site_effect","intercept","rho"),
											 monitors2 = c("plot_mu"), useConjugacy = F, 
											 control = list(#scale = 1, adaptInterval=50, 
											 							 adaptive=F))
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)
	
	
	# Remove large objects due to data usage
	# rm(model.dat, div_in, rank.df)
	# gc()
	
	# Run MCMC
	samples.out <- runMCMC(compiled, niter = iter,
												 nchains = 3, nburnin = burnin,
												 samplesAsCodaMCMC = T, thin = thin)
	
	# Process and summarize outputs
	samples2 <- rm.NA.mcmc(samples.out$samples2)
	samples <- rm.NA.mcmc(samples.out$samples)
	
	plot(samples)
	param_summary <- summary(samples)
	plot_summary <- summary(samples2)

	# Create outputs and clear compiled code
	metadata <- list(niter = iter,
									 nburnin = burnin,
									 thin = thin,
									 model_data = truth,
									 forecast_issue_time = Sys.Date())
	
	cat(paste0("Diversity model fit for run: ", group)) 
	out <- list(samples = samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	return(out)
}


pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
# Read in the microbial abundance data and covariates
# Create parameters to pass	
params = data.frame(index = 1:8,
										scenario = c("no_uncertainty_ITS", "spatial_uncertainty_ITS",
																 "temporal_uncertainty_ITS", "full_uncertainty_ITS",
																 "no_uncertainty_16S", "spatial_uncertainty_16S",
																 "temporal_uncertainty_16S", "full_uncertainty_16S"),
										group = c(rep("ITS", 4),rep("16S", 4)),
										temporalDriverUncertainty = c(F, F, T, T, F, F, T, T),
										spatialDriverUncertainty = c(F, T, F, T, F, T, F, T))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	out <- run_MCMC(group = params$group[[j]], iter = 150000, burnin = 100000, thin = 10, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]])
	return(out)
}

# Create cluster and pass it everything in the workspace
# cl <- makeCluster(8, outfile="")
# clusterExport(cl, ls())
# 
# # Run for all functional groups, in parallel (via PSOCK)
# output.list <- parLapply(cl, 1:8, run_scenarios)
out.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds"

library(doParallel)
cl <- makeCluster(4, type="PSOCK", outfile="")
registerDoParallel(cl)

# Only running the full driver uncertainty and zero driver uncertainty
output.list = foreach(j=c(1,4,5,8), #.export=c("run_scenarios","params","run_MCMC"), 
											.errorhandling = 'pass') %dopar% {
	print(params[j,])
	run_scenarios(j)
}
saveRDS(output.list, out.path)

sample.list <- lapply(output.list, "[[", 1)
param.summary.list <- lapply(output.list, "[[", 2)
metadata.list <- lapply(output.list, "[[", 3)
plot.summary.list <- lapply(output.list, "[[", 4)

names(sample.list) <- params$scenario[c(1,4,5,8)]
names(param.summary.list) <- params$scenario[c(1,4,5,8)]
names(metadata.list) <- params$scenario[c(1,4,5,8)]
names(plot.summary.list) <- params$scenario[c(1,4,5,8)]

saveRDS(list(sample.list = sample.list,
						 param.summary.list = param.summary.list, 
						 metadata.list = metadata.list,
						 plot.summary.list = plot.summary.list),
						 				out.path)

stopCluster(cl)




# no_uncertainty_ITS <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=F, 
# 													 temporalDriverUncertainty = F, 
# 													 spatialDriverUncertainty = F)
# spatial_uncertainty_ITS <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=F, 
# 															 temporalDriverUncertainty = F, 
# 															 spatialDriverUncertainty = T)
# full_uncertainty_ITS <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=F, 
# 																		temporalDriverUncertainty = T, 
# 																		spatialDriverUncertainty = T)
# saveRDS(list(no_uncertainty_ITS = no_uncertainty_ITS,
# 						 spatial_uncertainty_ITS = spatial_uncertainty_ITS,
# 						 full_uncertainty_ITS = full_uncertainty_ITS),
# 				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_ITS_samples_min3.rds")
# 
# no_uncertainty_16S <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=T, 
# 															 temporalDriverUncertainty = F, 
# 															 spatialDriverUncertainty = F)
# 
# spatial_uncertainty_16S <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=T, 
# 																		temporalDriverUncertainty = F, 
# 																		spatialDriverUncertainty = T)
# full_uncertainty_16S <- run_MCMC(group = "ITS", iter = 1000, burnin = 500, test=T, 
# 																 temporalDriverUncertainty = T, 
# 																 spatialDriverUncertainty = T)





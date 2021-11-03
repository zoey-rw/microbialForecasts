# Running dirichlet regression on 5 taxonomic ranks for soil fungi and bacteria


#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	j <- as.numeric( argv[1] )
} #else {
	#j=1
#}

# # # For testing
# iter <- 500
# burnin <- 200
# thin <- 1
# test = F=T
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"

run_MCMC <- function(k = 17, iter = 1000,  burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Subset to one rank.
	rank.name <- tax_names[k]
	rank.df <- d[[rank.name]] 
	# rank.df <- rank.df[which(rank.df$siteID %in% c("NIWO")),]
	# rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),]
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))
	model.dat <- prepTaxonomicData(rank.df = rank.df, min.prev = 3)
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario, ".rds")
	print(paste0("Completed model data for ", rank.name))

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
									relEM = model.dat[["relEM"]],
									plotID = model.dat$plotID,
									plot_site = model.dat$plot_site,
									plot_num = model.dat$plot_num,
									plot_site_num = model.dat$plot_site_num,									
									plot_start = model.dat[["plot_start"]],
									plot_index = model.dat[["plot_index"]],
									site_start = model.dat[["site_start"]],
									N.beta = 6
									)

# for output.
metadata <- list("rank.name" = rank.name,
								 "niter" = iter,
								 "nburnin" = burnin,
								 "thin" = thin,
								 "model_data" = model.dat$truth.plot.long)

inits <- initsFun(constants, type = "tax")
## Configure & compile model
Rmodel <- nimbleModel(code = nimbleModTaxa,
											constants = constants, data = list(y=model.dat$y),
											inits = inits)
# Compile model
cModel <- compileNimble(Rmodel)
nimbleOptions(multivariateNodesAsScalars = TRUE)
# Configure & compile MCMC
mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect","sig","intercept", 
																							 "rho"),
													monitors2 = c("plot_rel"), thin2 = 25,
													useConjugacy = T#, control = c(scale = .1)
													)
# mcmcConf$removeSamplers(c('site_effect'))
# mcmcConf$addSampler(target = c('site_effect'), type = 'AF_slice')
myMCMC <- buildMCMC(mcmcConf)
#mcmcConf
# myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","site_effect","intercept", "rho"),
# 													monitors2 = c("plot_rel"), useConjugacy = F)
compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
# Remove large objects due to data usage
 rm(model.dat, d, rank.df, constants); gc()

# Sample from MCMC
samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
									 nchains = 3, samplesAsCodaMCMC = T, thin = thin)

cat(paste0("Finished sampling for ", rank.name)); gc()

out <- list(samples = samples$samples, samples2 = samples$samples2,
						metadata = metadata)
saveRDS(out, out.path, compress = FALSE)

cat(paste0("Saved raw samples for ", rank.name))

# Calculate summary and save output.
param_summary <- fast.summary.mcmc(samples$samples)
plot_summary <- fast.summary.mcmc(samples$samples2)

out <- list(samples = samples$samples,
						param_summary = param_summary,
						metadata = metadata,
						plot_summary = plot_summary)
saveRDS(out, out.path, compress = F)
cat(paste0("Saved summary output for ", rank.name," to ", out.path))
return("ok")
}


# Running dirichlet model 
pacman::p_load(doParallel, reshape2, nimble, coda, tidyverse) 


# Create parameters to pass	
n.groups <- 10 # Number of taxonomic ranks
params = data.frame(group = rep(1:n.groups, 2),
										scenario = c(rep("no_uncertainty", n.groups), 
																 rep("full_uncertainty", n.groups)),
										temporalDriverUncertainty = c(rep(FALSE, n.groups),
																									rep(TRUE, n.groups)),
										spatialDriverUncertainty = c(rep(FALSE, n.groups),
																								 rep(TRUE, n.groups)))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	out <- run_MCMC(k = params$group[[j]], iter = 1000000, burnin = 500000, thin = 15, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]],
									scenario = params$scenario[[j]])
	return(out)
}

# cl <- makeCluster(4, type="PSOCK", outfile="fit_taxa.log")
# registerDoParallel(cl)
# #Run for all 10 taxonomic ranks, in parallel (via PSOCK)
# output.list = foreach(j=c(14:20),
# 											#.export=c("run_scenarios","params","run_MCMC"),
# 											.errorhandling = 'pass') %dopar% {
# 												run_scenarios(j)
# 											}
# 
# stopCluster(cl)
# 
# # NOT PARALLEL
# for(j in c(11:20)){
# 	print(params[j,])
# 	run_scenarios(j)
#}

# Run for a single rank using command-line arguments 
 run_scenarios(j)

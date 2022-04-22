# iter <- 1000
# burnin <- 500
# thin <- 1
# test = F
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"
# nchains=3
# model_name = "cycl_only"
# #model_name = "all_covariates"
# time_period = "calibration"
# #time_period = "refit"
# thin = 10
# iter_per_chunk = 1000
# thin = 5
# init_iter = 5000
# chain_no = 1
run_MCMC <- function(k = 1, iter = 1000, init_iter = 100000, iter_per_chunk = 10000, burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 s = "acidobacteriota",
										 model_name = "all_covariates",
										 time_period = "calibration", chain_no = 1, 
										 ...) {
	#set.seed(chain_no)
	

	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/cyclical_model_source.r")
	
	# Read in microbial abundances
	cal <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
					 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	val <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds"), 
					 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds"))
	
	# Subset to one rank.
	rank.name <- tax_names[k]
	
	# Grab names of taxa to keep
	keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
	keep_names <- keep_list[[rank.name]]$taxon.name
	keep_vec <- c(keep_names, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
	
	
	# Prep model inputs/outputs.
	cat(paste0("\nPreparing model data for ", model_name, ", ", time_period, ", ", rank.name, "\n"))
	
	out.path <- 
		paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", 
					 model_name,"/", time_period,  "_samples_", rank.name,"_",scenario, ".rds")
	
	if (time_period == "refit"){
		
		
		cal.rank.df <- cal[[rank.name]] 
		val.rank.df <- val[[rank.name]] 
		rank.df <- rbind(cal.rank.df, val.rank.df)	
		rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),] #Remove sites missing key covariates
		
		if (test == T) rank.df = rank.df[1:500,]
		
		colnames(rank.df) <- lapply(strsplit(colnames(rank.df), "\\."), "[[", 1)
		rank.df_spec <- rank.df[,colnames(rank.df) %in% keep_vec] 
		tots <- rowSums(rank.df_spec[,7:ncol(rank.df_spec)])
		rank.df_spec$other <- 1-tots
		# Remove samples where more than 99% of reads are "Other"
		rank.df_spec <- rank.df_spec %>% filter(tots > .01)
		model.dat <- prepTaxonomicData(rank.df = rank.df_spec, min.prev = 3, max.date = "20200101")
		
	} else if (time_period == "calibration") {
		
		rank.df <- cal[[rank.name]] 
		rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),]
		
		# Reduce size for testing
		if (test == T) {
			rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
			rank.df = rank.df[1:500,]
		}
		
		colnames(rank.df) <- lapply(strsplit(colnames(rank.df), "\\."), "[[", 1)
		rank.df_spec <- rank.df[,colnames(rank.df) %in% keep_vec] 
		tots <- rowSums(rank.df_spec[,7:ncol(rank.df_spec)])
		rank.df_spec$other <- 1-tots
		# Remove samples where more than 99% of reads are "Other"
		rank.df_spec <- rank.df_spec %>% filter(tots > .01)
		model.dat <- prepTaxonomicData(rank.df = rank.df_spec, min.prev = 3)
		
	} else cat("Missing specification of time period.")
	
	
	cat(paste0("\nCompleted model data for ", model_name, ", ", time_period, ", ", rank.name, "\n"))
	
	mo <- month(as.Date(paste0(colnames(model.dat$mois), "01"), format="%Y%m%d"))
	y_sin = sin((2*pi*mo)/12)
	y_cos = cos((2*pi*mo)/12)
	
	
	if (model_name == "cycl_only"){
		n.beta = 2
		Nimble_model = nimbleModTaxa_cycl_only
	} else {
		n.beta = 8
		Nimble_model = nimbleModTaxa
	}
	
	constants <- list(N.plot =  length(unique(model.dat$plotID)), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(unique(model.dat$siteID)),
										timepoint = model.dat$timepoint,
										sin_mo = y_sin,
										cos_mo = y_cos,
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
										LAI = model.dat[["LAI"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,									
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = n.beta
	)
	
	constants$omega <- 0.0001 * diag(constants$N.spp)
	constants$zeros = rep(0, constants$N.spp)
	if (constants$N.spp < 8) {
		constants$omega <- 0.0001 * diag(8)
		constants$zeros = rep(0, 8)
	}
	
	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = init_iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long)
	
	inits <- initsFun(constants, type = "tax")
	cat("\nInits created\n")
	
	#inits
	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants, data = list(y=model.dat$y),
												inits = inits)
	
	#return(list(rnorm(1,0,1), Rmodel$plot_mu[1,1,30]))
	
	cat("\nModel compiled")
	
	# Compile model
	cModel <- compileNimble(Rmodel)
	
	
	cat("\nModel compiled")
	
#	nimbleOptions(multivariateNodesAsScalars = TRUE)
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
																								 "sig","intercept", 
																								 "rho"),
														monitors2 = c("plot_rel"), thin2 = 25,
														useConjugacy = T, )
	
	cat("\nMCMC configured")
	
	#mcmcConf$removeSamplers(c('sigma', 'intercept', 'rho'))
	#mcmcConf$addSampler(target = c('sigma', 'intercept', 'rho'), type = 'RW_block')
	
	myMCMC <- buildMCMC(mcmcConf)
	cat("\nModel built")
	
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	cat("\nMCMC compiled")
	
	#out.path2 <- gsub(".rds", paste0("_", format(Sys.time(), "%F_%H-%M-%S"), ".rds"), out.path)
	out.path2 <- gsub(".rds", paste0("_chain", chain_no, ".rds"), out.path)
	#iter_per_chunk <- 50000
	
	compiled$run(niter=init_iter, thin=thin, nburnin = burnin)
	cat(paste0("\nInitial run finished for chain", chain_no))
	
	out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata), 
					out.path2)

	continue <- check_continue(out.run, min_eff_size = 10) 
	while (continue){
		cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
		compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
		out.run<-as.mcmc(as.matrix(compiled$mvSamples))
		out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
		saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata), 
						out.path2)
		continue <- check_continue(out.run, min_eff_size = 10)
	}
		
	# 	continue <- check_continue(out.run)
	# 	if (continue){
	# 		cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
	# 		compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
	# 		out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	# 		out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	# 		saveRDS(list(out.run, out.run2), out.path2)
	# 
	# 		continue <- check_continue(out.run)
	# 		if (continue){
	# 			cat(paste0("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n"))
	# 			compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
	# 			out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	# 			out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	# 			saveRDS(list(out.run, out.run2), out.path2)
	# 		}
	# 	}
	# }
	
	# # Sample from MCMC
	# samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
	# 									 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	
	cat(paste0("\nFinished sampling for ", model_name, ", ", time_period, ", ", rank.name, ", chain: ", chain_no, "\n")); gc()
	
	# 
	# # Calculate summary and save output.
	# param_summary <- fast.summary.mcmc(samples$samples)
	# plot_summary <- fast.summary.mcmc(samples$samples2)
	# 
	# out <- list(samples = samples$samples,
	# 						param_summary = param_summary,
	# 						metadata = metadata,
	# 						plot_summary = plot_summary,
	# 						gelman = gelman.diag(samples$samples))
	# saveRDS(out, out.path, compress = F)
	# cat(paste0("\nSaved summary output for ", model_name, ", ", time_period, ", ", 
	# 					 rank.name, " to: \n", out.path, "\n"))
	
	return("ok")
}

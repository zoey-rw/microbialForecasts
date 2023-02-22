#' @title run_MCMC_bychain
#' @description run_MCMC_bychain
#' @export
#
# # For testing
# iter <- 500
# burnin <- 200
# thin <- 1
# test =F
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"
# s = "acidobacteriota"
# k = 6
# s = "ascomycota"
# model_name = "all_covariates"
# chain_no = 1
# min.date = "20151101"
# max.date = "20200101"
# max.date = "20180101"
# iter_per_chunk = 500
# init_iter = 1000
# model_name = "cycl_only"
# min_eff_size_perchain = 10
# max_loops=100
# thin2 = 20

run_MCMC_bychain <- function(k = 1,
														 iter = 1000,  burnin = 500, thin = 1,

														 thin2 = 20,
														 test = F, chain_no = 1,
														 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
														 scenario=NULL,
														 s = "acidobacteriota",
														 min.date = "20151101",
														 max.date = "20180101",
														 model_name = "cycl_only",

														 iter_per_chunk = 5000,
														 init_iter = 10000,
														 min_eff_size_perchain = 5,
														 max_loops=100,
														 ...) {


	## ################################################
	## A DISTRIBUTION GIVING THE LOGIT OF A BETA DISTRIBUTION ##
	## ################################################
	dLogitBeta <- nimbleFunction (
		## Returns density of x where
		##                    y ~ Beta(a1,a2)
		##                    x = logit(y)
		run = function(x = double(0),
									 shape1=double(0, default=1.0),
									 shape2=double(0, default=1.0),
									 log = integer(0, default=0)) {
			returnType(double(0))
			y = ilogit(x)
			logProbX = log(y) + log(1 - y) + dbeta(y, shape1=shape1, shape2=shape2, log=TRUE) ## Via change of variables
			if (log)
				return(logProbX)
			return(exp(logProbX))
		}
	)

	rLogitBeta <- nimbleFunction (
		## Generates y ~ Beta(a1,a2)
		## Returns   x = logit(y)
		run = function(n = integer(0, default=1),
									 shape1 = double(0, default=1.0),
									 shape2 = double(0, default=1.0)) {
			returnType(double(0))
			if(n != 1)
				nimPrint("Warning: rLogitBeta only allows n = 1; Using n = 1.\n")
			y <- rbeta(1, shape1=shape1, shape2=shape2)
			x <- logit(y)
			return(x)
		}
	)

	registerDistributions(list(dLogitBeta = list(
		BUGSdist = "dLogitBeta(shape1, shape2)",
		discrete = FALSE,
		types    = c("value=double(0)"), ## , "para=double(0)"
		pqAvail  = FALSE)))

	assign('dLogitBeta', dLogitBeta, envir = .GlobalEnv)
	assign('rLogitBeta', rLogitBeta, envir = .GlobalEnv)


#####
	# to consider: https://groups.google.com/g/nimble-users/c/t1ArfNDPcdg/m/1bT8qXecEAAJ
	# nimbleModBeta_cycl_only <- nimble::nimbleCode({
	#
	# 	# Loop through core observations ----
	# 	for (i in 1:N.core) {
	# 		y[i, 1] ~ dbeta(mean = plot_mu[plot_num[i], timepoint[i]],
	# 										sd = core_sd)
	# 	}
	#
	# 	# Plot-level process model ----
	# 	for (p in 1:N.plot) {
	#
	# 		for (t in plot_start[p]) {
	#
	# 			# Plot means for first date
	# 			Ex[p, t] ~ dbeta(mean = X_init, sd = .2)
	#
	# 			# Add process error (sigma)
	# 			plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
	# 			#plot_mu[p, t] ~ T(dnorm(mean = Ex[p, t], sd = sigma), 0, 1)
	# 		}
	#
	# 		# Second date onwards
	# 		for (t in plot_index[p]:N.date) {
	# 			# Previous value * rho
	# 			logit(Ex[p, t]) <-
	# 				rho * logit(plot_mu[p, t - 1]) +
	# 				beta[1] * sin_mo[t] +
	# 				beta[2] * cos_mo[t] +
	# 				site_effect[plot_site_num[p]] +
	# 				intercept
	#
	# 			# Add process error (sigma)
	# 			plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
	# 		}
	# 	}
	#
	# 	# Priors for site effects----
	# 	for (k in 1:N.site) {
	# 		site_effect[k] ~ dnorm(0,  sd = sig)
	# 	}
	#
	# 	core_sd ~ dinvgamma(5,.5)
	# 	sigma ~ dinvgamma(5,.5)
	# 	sig ~ dinvgamma(3,1)
	#
	# 	intercept ~ dnorm(0, sd = 1)
	# 	rho ~ dnorm(0, sd = 1) # Autocorrelation
	#
	# 	X_init <- .1 # Initial plot mean
	#
	#
	# 	beta[1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
	#
	# }) #end NIMBLE model.
#####

	nimbleModBeta_cycl_only <- nimble::nimbleCode({


		# Loop through core observations ----
		for (i in 1:N.core) {
			# y[i, 1] ~ dbeta(mean = Ex[plot_num[i], timepoint[i]],
			# 								sd = core_sd)
			y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]], sd = core_sd), 0, 1)
		}


		# Plot-level process model ----
		for (p in 1:N.plot) {

			logit(Ex[p, 1]) ~ dnorm(mean = X_init, sd = .2)
			plot_mu[p, 1] ~ dnorm(mean = Ex[p, 1] , sd = .01)

			for (t in 2:N.date) {

				# Dynamic linear model
				logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
					beta[1] * sin_mo[t] +
					beta[2] * cos_mo[t] +
					site_effect[plot_site_num[p]] +
					intercept

				# Define shape parameters for process error (sigma) model
				shape1[p, t] <- Ex[p, t] * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)
				shape2[p, t] <- (1 - Ex[p, t]) * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)

				# Add process error (sigma)
				logit(plot_mu[p, t]) ~ dLogitBeta(shape1[p, t], shape2[p, t])
			}
		}



		# Priors for site effects----
		for (k in 1:N.site) {
			site_effect[k] ~ dnorm(0,  sd = sig)
		}
		# Priors for everything else ----
		core_sd ~ dinvgamma(5,.5)
		sigma ~ dinvgamma(5,.5)
		sig ~ dinvgamma(3,1)

		intercept ~ dnorm(0, sd = 1)
		rho ~ dnorm(0, sd = 1) # Autocorrelation

		X_init <- .5 # Initial plot mean
		beta[1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])

	}) #end NIMBLE model.

	nimbleModBeta_all_covariates <- nimble::nimbleCode({


		# Loop through core observations ----
		for (i in 1:N.core) {
			# y[i, 1] ~ dbeta(mean = plot_mu[plot_num[i], timepoint[i]],
			# 								sd = core_sd)
			y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]], sd = core_sd), 0, 1)
		}

		# Plot-level process model ----
		for (p in 1:N.plot) {

			for (t in plot_start[p]) {

				# Plot means for first date
				Ex[p, t] ~ dbeta(mean = X_init, sd = .2)

				# Add process error (sigma)
				plot_mu[p, t] ~ dbeta(mean = Ex[p, t], sd = sigma)
				#plot_mu[p, t] ~ T(dnorm(mean = Ex[p, t], sd = sigma), 0, 1)
			}

			# Second date onwards
			for (t in plot_index[p]:N.date) {
				# Dynamic linear model
				logit(Ex[p, t]) <- rho * logit(plot_mu[p, t - 1]) +
					beta[1] * temp_est[plot_site_num[p], t] +
					beta[2] * mois_est[plot_site_num[p], t] +
					beta[3] * pH_est[p, plot_start[p]] +
					beta[4] * pC_est[p, plot_start[p]] +
					#beta[5]*nspp[p,t] +
					beta[5] * relEM[p, t] +
					beta[6] * LAI[plot_site_num[p], t] +
					beta[7] * sin_mo[t] +
					beta[8] * cos_mo[t] +
					site_effect[plot_site_num[p]] +
					intercept

				# Define shape parameters for process error (sigma) model
				shape1[p, t] <- Ex[p, t] * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)
				shape2[p, t] <- (1 - Ex[p, t]) * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)

				# Add process error (sigma)
				logit(plot_mu[p, t]) ~ dLogitBeta(shape1[p, t], shape2[p, t])
			}
		}

		# Priors for site effects----
		for (k in 1:N.site) {
			site_effect[k] ~ dnorm(0,  sd = sig)
		}

		core_sd ~ dinvgamma(5,.5)
		sigma ~ dinvgamma(5,.5)
		sig ~ dinvgamma(3,1)

		intercept ~ dnorm(0, sd = 1)
		rho ~ dnorm(0, sd = 1) # Autocorrelation

		X_init <- .1 # Initial plot mean

		beta[1:8] ~ dmnorm(zeros[1:8], omega[1:8, 1:8])



		# Add driver uncertainty if desired ----
		if (temporalDriverUncertainty) {
			for (k in 1:N.site) {
				for (t in site_start[k]:N.date) {
					mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
					temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
				}
			}
		} else {
			for (k in 1:N.site) {
				for (t in site_start[k]:N.date) {
					mois_est[k, t] <- mois[k, t]
					temp_est[k, t] <- temp[k, t]
				}
			}
		}

		# Using 40th time point (values are constant over time)
		if (spatialDriverUncertainty) {
			for (p in 1:N.plot) {
				pH_est[p, plot_start[p]] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
				pC_est[p, plot_start[p]] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
			}
		} else {
			for (p in 1:N.plot) {
				pH_est[p, plot_start[p]] <- pH[p, plot_start[p]]
				pC_est[p, plot_start[p]] <- pC[p, plot_start[p]]
			}
		}


	}) #end NIMBLE model.


	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse)

	# Subset to one rank.
	rank.name <- microbialForecast:::tax_names[k]

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

	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID",
																														"dateID", "sampleID",
																														"dates", "plot_date","other")]
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s)
	rank.df_spec$other <- 1-rank.df_spec[[s]]

	model.dat <- prepBetaRegData(rank.df = rank.df_spec,
															 min.prev = 3,
															 min.date = min.date,
															 max.date = max.date)

	out.path <- here("data",paste0("model_outputs/logit_beta_regression/",model_name,"/samples_", rank.name, "_", s, "_",min.date,"_",max.date, ".rds"))
	out.path2 <- gsub(".rds", paste0("_chain", chain_no, ".rds"), out.path)


	print(paste("Completed model data for", rank.name, ", group:", s))
	constants <- model.dat[c("plotID",  "timepoint","plot_site",
													 "site_start", "plot_start", "plot_index",
													 "plot_num", "plot_site_num",
													 "N.plot", "N.spp", "N.core", "N.site", "N.date",
													 "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd",
													 "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]


	if (model_name == "cycl_only"){
		constants$N.beta = 2
		Nimble_model = nimbleModBeta_cycl_only
		constants <- constants[c("plotID",  "timepoint","plot_site",
														 "site_start", "plot_start", "plot_index",
														 "plot_num", "plot_site_num",
														 "N.plot", "N.spp", "N.core", "N.site",
														 "N.date", "N.beta","sin_mo", "cos_mo")]
	} else if(model_name == "all_covariates"){
		constants$N.beta = 8
		Nimble_model = nimbleModBeta_all_covariates
	} else message("Missing specification of Nimble model.")

	# Some model hyperparameters
	constants$omega <- 0.0001 * diag(constants$N.beta)
	constants$zeros = rep(0, constants$N.beta)

	# This date scaling should be done within the prepBetaRegData function instead
	constants$sin_mo = scale(constants$sin_mo, center = F) %>% as.numeric()
	constants$cos_mo = scale(constants$cos_mo, center = F) %>% as.numeric()


	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long,
									 "sample_values" = model.dat$sample_values)

	inits <- initsFun(constants, type = "fg")
	#
	# inits <- inits[c("y", "plot_mu", "intercept", "sig", "beta", "rho", "sigma",
	# 								 "plot_rel", "site_effect")]
	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
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
																								 "rho","core_sd"),
														monitors2 = c("plot_mu"), thin2 = thin2,
														useConjugacy = T)

	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)


	compiled$run(niter=init_iter, thin=thin, thin2 = thin2, nburnin = burnin)
	cat(paste0("\nInitial run finished for chain", chain_no))

	# Sample from MCMC

	out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	out.run2 <- na.omit(out.run2)
	# Remove timepoints that had no sampling at all (only zeros)
	out.run2 <- mcmc(out.run2[, which(colSums(out.run2) != 0)])

	saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
					out.path2)

	continue <- check_continue(out.run, min_eff_size = min_eff_size_perchain)
	loop_counter = 0
	while (continue){
		message("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n")
		if (loop_counter < max_loops) {

			compiled$run(niter=iter_per_chunk, thin=thin, reset=F)

			# Shorten if more than 10k samples have accumulated
			out.run = window_chain(compiled$mvSamples, max_size = 20000)
			out.run2 <- window_chain(compiled$mvSamples2, max_size = 20000)
			out.run2 <- na.omit(out.run2)
			# Remove timepoints that had no sampling at all (only zeros)
			out.run2 <- mcmc(out.run2[, which(colSums(out.run2) != 0)])

			metadata$niteration = iter_per_chunk * loop_counter

			saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
							out.path2)
			continue <- check_continue(out.run, min_eff_size = min_eff_size_perchain)
			loop_counter = loop_counter + 1
		} else {
			message("Exceeded the max number of loops (max_loops). No more sampling!")
			continue <- FALSE
		}
	}

	return("ok")
}

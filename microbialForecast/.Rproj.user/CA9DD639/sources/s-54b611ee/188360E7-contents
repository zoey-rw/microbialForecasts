
#' @title nimbleMod_shannon
#' @description nimbleMod_shannon
#' @export
nimbleMod_shannon <- nimbleCode({

	# Observation model (cores ~ plot means)
	for(i in 1:N.core){
		y[i,1] ~ dnorm(plot_mu[plot_num[i],timepoint[i]],  core_sd)
	}

	# Process model
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dgamma(2, 1) # Plot means for first date
		}

		for (t in plot_index[p]:N.date) { # Starts from second date
			# Previous value * rho + covariates
			Ex[p,t] <- rho * plot_mu[p,t-1] +
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,plot_start[p]] +
				beta[4]*pC_est[p,plot_start[p]] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
			# Add process error, sigma
			plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
		}
	}

	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		}
	}

	# Add spatial uncertainty (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] ~ dnorm(pH[p,plot_start[p]], sd = pH_sd[p,plot_start[p]])
			pC_est[p,plot_start[p]] ~ dnorm(pC[p,plot_start[p]], sd = pC_sd[p,plot_start[p]])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] <- pH[p,plot_start[p]]
			pC_est[p,plot_start[p]] <- pC[p,plot_start[p]]
		}
	}

	rho ~ dnorm(0, sd = 1)
	#core_sd ~ dinvgamma(3, .5)
	core_sd ~ dgamma(.1, 1)
	sigma ~ dgamma(.5, .1)
	#	intercept ~ dgamma(.1, .1)
	intercept ~ dnorm(0, sd = 1)

	# Priors for covariates:
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}
	# Priors for site effect variance ----
	sig ~ dgamma(.5,1)

}) #end NIMBLE model.



#' @title nimbleModTaxa
#' @description nimbleModTaxa
#' @export
nimbleModTaxa <- nimbleCode({

	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
	}

	# Plot-level process model ----
	for(s in 1:N.spp){
		for(p in 1:N.plot){
			for (t in plot_start[p]) {
				plot_mu[p,s,t] ~ dgamma(0.5, 1) # Plot means for first date
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}

			for (t in plot_index[p]:N.date) {
				# Previous value * rho
				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) +
					beta[s,1]*temp_est[plot_site_num[p],t] +
					beta[s,2]*mois_est[plot_site_num[p],t] +
					beta[s,3]*pH_est[p,plot_start[p]] +
					beta[s,4]*pC_est[p,plot_start[p]] +
					beta[s,5]*relEM[p,t] +
					beta[s,6]*LAI[plot_site_num[p],t] +
					beta[s,7]*sin_mo[t] + beta[s,8]*cos_mo[t] +
					site_effect[plot_site_num[p],s] +
					intercept[s]
				# Add process error (sigma)
				#plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sigma[s])
				plot_mu[p, s, t] ~ T(dnorm(mean = Ex[p, s, t],
																	 sigma[s]), 0, Inf)
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}
		}
	}

	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		}
	}

	# Using 40th time point (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] ~ dnorm(pH[p,plot_start[p]], sd = pH_sd[p,plot_start[p]])
			pC_est[p,plot_start[p]] ~ dnorm(pC[p,plot_start[p]], sd = pC_sd[p,plot_start[p]])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] <- pH[p,plot_start[p]]
			pC_est[p,plot_start[p]] <- pC[p,plot_start[p]]
		}
	}

	# Priors for site effect covariance matrix ----
	sig ~ dgamma(3,1)
	# SIGMA[1:N.spp,1:N.spp] <- diag(rep(sig^2, N.spp))
	#
	# # Priors for site random effects:
	# for(k in 1:N.site){
	# 	site_effect[k,1:N.spp] ~ dmnorm(alpha0[1:N.spp], # vector of zeros
	# 																	cov = SIGMA[1:N.spp,1:N.spp])
	# }

	# Priors for site random effects:
	for(s in 1:N.spp){
		for(k in 1:N.site){
			site_effect[k,s] ~ dnorm(0, sig)
		}
	}



	# Priors for everything else ----
	for (s in 1:N.spp){
		sigma[s] ~ dgamma(.1, .1)
		beta[s,1:N.beta] ~ dmnorm(zeros[1:N.beta], omega[1:N.beta, 1:N.beta])
	}
	rho[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	intercept[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])


}) #end NIMBLE model.




#' @title nimbleModFunctional
#' @description nimbleModFunctional
#' @export
nimbleModFunctional <- nimbleCode({

	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1] ~ dbeta(mean = plot_mu[plot_num[i],timepoint[i]],
									 sd = core_sd)
	}

	# Plot-level process model ----
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dbeta(mean=.3, sd=.3) # Plot means for first date
		}

		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) +
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,plot_start[p]] +
				beta[4]*pC_est[p,plot_start[p]] +
				#beta[5]*nspp[p,t] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ dnorm(mean = Ex[p,t], sigma)
			#plot_mu[p,t] <- Ex[p,t]
		}
	}


	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		}
	}

	# Using 40th time point (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] ~ dnorm(pH[p,plot_start[p]], sd = pH_sd[p,plot_start[p]])
			pC_est[p,plot_start[p]] ~ dnorm(pC[p,plot_start[p]], sd = pC_sd[p,plot_start[p]])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] <- pH[p,plot_start[p]]
			pC_est[p,plot_start[p]] <- pC[p,plot_start[p]]
		}
	}


	# Priors for site effect covariance matrix ----
	sig ~ dgamma(.5,1)
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}

	# Priors for everything else ----
	core_sd ~ dgamma(.5,1)
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	#intercept ~ dgamma(.1, .1)
	intercept ~ dnorm(0, sd = 1)

	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.




#' @title nimbleModFunctional_trunc
#' @description Model version with truncated Normal distribution to prevent negative values
#' @export
nimbleModFunctional_trunc <- nimbleCode({

	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1] ~ dbeta(mean = plot_mu[plot_num[i],timepoint[i]],
									 sd = core_sd)
	}

	# Plot-level process model ----
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dbeta(mean=.3, sd=.3) # Plot means for first date
		}

		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) +
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,plot_start[p]] +
				beta[4]*pC_est[p,plot_start[p]] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ T(dnorm(mean = Ex[p,t], sigma), 0, Inf)
			#plot_mu[p,t] <- Ex[p,t]
		}
	}


	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		}
	}

	# Using 40th time point (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] ~ dnorm(pH[p,plot_start[p]], sd = pH_sd[p,plot_start[p]])
			pC_est[p,plot_start[p]] ~ dnorm(pC[p,plot_start[p]], sd = pC_sd[p,plot_start[p]])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,plot_start[p]] <- pH[p,plot_start[p]]
			pC_est[p,plot_start[p]] <- pC[p,plot_start[p]]
		}
	}


	# Priors for site effect covariance matrix ----
	sig ~ dgamma(.5,1)
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}

	# Priors for everything else ----
	core_sd ~ dgamma(.5,1)
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	intercept ~ dnorm(0, sd = 1)

	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.






#' @title nimbleModFunctional_cycl_only
#' @description nimbleModFunctional_cycl_only
#' @export
nimbleModFunctional_cycl_only <- nimbleCode({

	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1] ~ dbeta(mean = plot_mu[plot_num[i],timepoint[i]],
									 sd = core_sd)
	}

	# Plot-level process model ----
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dbeta(mean=.3, sd=.3) # Plot means for first date
		}

		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) +
				beta[1]*sin_mo[t] + beta[2]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept

			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ T(dnorm(mean = Ex[p,t], sigma), 0, Inf)
			#plot_mu[p,t] <- Ex[p,t]
		}
	}

	# Priors for site effect covariance matrix ----
	sig ~ dgamma(.5,1)
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}

	# Priors for everything else ----
	core_sd ~ dgamma(.5,1)
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	intercept ~ dnorm(0, sd = 1)

	for (n in 1:2){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.





#' @title nimbleMod_shannon_cycl_only
#' @description nimbleMod_shannon_cycl_only
#' @export
nimbleMod_shannon_cycl_only <- nimbleCode({

	# Observation model (cores ~ plot means)
	for(i in 1:N.core){
		y[i,1] ~ dnorm(plot_mu[plot_num[i],timepoint[i]],  core_sd)
	}

	# Process model
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dgamma(2, 1) # Plot means for first date
		}

		for (t in plot_index[p]:N.date) { # Starts from second date
			# Previous value * rho + covariates
			Ex[p,t] <- rho * plot_mu[p,t-1] +
				beta[1]*sin_mo[t] + beta[2]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
			# Add process error, sigma
			plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
		}
	}

	rho ~ dnorm(0, sd = 1)
	#core_sd ~ dinvgamma(3, .5)
	core_sd ~ dgamma(.1, 1)
	sigma ~ dgamma(.5, .1)
	#	intercept ~ dgamma(.1, .1)
	intercept ~ dnorm(0, sd = 1)

	# Priors for covariates:
	for (n in 1:2){
		beta[n] ~ dnorm(0, sd = 1)
	}
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}
	# Priors for site effect variance ----
	sig ~ dgamma(.5,1)

}) #end NIMBLE model.





#' @title nimbleModTaxa_cycl_only
#' @description nimbleModTaxa_cycl_only
#' @export
nimbleModTaxa_cycl_only <- nimbleCode({

	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
	}

	# Plot-level process model ----
	for(s in 1:N.spp){
		for(p in 1:N.plot){
			for (t in plot_start[p]) {
				plot_mu[p,s,t] ~ dgamma(0.5, 1) # Plot means for first date
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}

			for (t in plot_index[p]:N.date) {
				# Previous value * rho
				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) +
					beta[s,1]*sin_mo[t] + beta[s,2]*cos_mo[t] +
					site_effect[plot_site_num[p],s] +
					intercept[s]
				# Add process error (sigma)
				plot_mu[p,s,t] ~ T(dnorm(mean = Ex[p,s,t], sigma[s]), 0, Inf)
				#plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sigma[s])
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}
		}
	}

	# Priors for site effect covariance matrix ----
	sig ~ dgamma(3,1)

	# Priors for site random effects:
	for(s in 1:N.spp){
		for(k in 1:N.site){
			site_effect[k,s] ~ dnorm(0, sig)
		}
	}


	# Priors for everything else ----
	for (s in 1:N.spp){
		sigma[s] ~ dgamma(.1, .1)
		beta[s,1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
	}
	rho[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	intercept[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])

}) #end NIMBLE model.

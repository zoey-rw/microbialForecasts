

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
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
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
	
	#for (n in 1:2){
		beta[1:2] ~ dnorm(0, sd = 1)
	#}
}) #end NIMBLE model.




#####
# nimbleModTaxa_cycl_only <- nimbleCode({ 
# 	
# 	# Loop through core observations ----
# 	for(i in 1:N.core){
# 		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
# 	}
# 	
# 	# Plot-level process model ----
# 	for(s in 1:N.spp){
# 		for(p in 1:N.plot){
# 			for (t in plot_start[p]) {
# 				plot_mu[p,s,t] ~ dgamma(0.5, 1) # Plot means for first date
# 				# Convert back to relative abundance
# 				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
# 			}
# 			
# 			for (t in plot_index[p]:N.date) {
# 				# Previous value * rho
# 				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) + 
# 					beta[s,1]*sin_mo[t] + beta[s,2]*cos_mo[t] +
# 					site_effect[plot_site_num[p],s] +
# 					intercept[s]
# 				# Add process error (sigma)
# 				plot_mu[p,s,t] ~ T(dnorm(mean = Ex[p,s,t], sigma[s]), 0, Inf)
# 				#plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sigma[s])
# 				# Convert back to relative abundance
# 				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
# 			}
# 		}
# 	}
# 	
# 	# Add driver uncertainty if desired ----
# 	if(temporalDriverUncertainty) {
# 		for(k in 1:N.site){
# 			for (t in site_start[k]:N.date) {
# 				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
# 				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
# 			}
# 		}
# 	} else {
# 		for(k in 1:N.site){
# 			for (t in site_start[k]:N.date) {
# 				mois_est[k,t] <- mois[k,t]
# 				temp_est[k,t] <- temp[k,t]
# 			}
# 		} 
# 	}
# 	
# 	# Using 40th time point (values are constant over time)
# 	if(spatialDriverUncertainty) {
# 		for(p in 1:N.plot){
# 			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
# 			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
# 		}
# 	} else {
# 		for(p in 1:N.plot){
# 			pH_est[p,1] <- pH[p,1]
# 			pC_est[p,1] <- pC[p,1]
# 		} 
# 	}
# 	
# 	# Priors for site effect covariance matrix ----
# 	sig ~ dgamma(3,1)
# 	# SIGMA[1:N.spp,1:N.spp] <- diag(rep(sig^2, N.spp))		
# 	# 
# 	# # Priors for site random effects:
# 	# for(k in 1:N.site){
# 	# 	site_effect[k,1:N.spp] ~ dmnorm(alpha0[1:N.spp], # vector of zeros
# 	# 																	cov = SIGMA[1:N.spp,1:N.spp])
# 	# }
# 	
# 	# Priors for site random effects:
# 	for(s in 1:N.spp){
# 		for(k in 1:N.site){
# 			site_effect[k,s] ~ dnorm(0, sig)
# 		}
# 	}
# 	
# 	
# 	# Priors for everything else ----
# 	for (s in 1:N.spp){
# 		rho[s] ~ dnorm(0, sd = 1)
# 		sigma[s] ~ dgamma(.1, .1)
# 		intercept[s] ~ dnorm(0, sd = 1)
# 		for (n in 1:2){
# 			beta[s,n] ~ dnorm(0, sd = 1)
# 		}
# 	}
# 	
# 	
# }) #end NIMBLE model.
#####



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
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
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
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
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
		beta[s,1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
	}
	rho[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	intercept[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	
}) #end NIMBLE model.


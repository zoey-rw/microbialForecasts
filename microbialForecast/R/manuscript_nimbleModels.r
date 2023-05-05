# source(here("source.R"))

# ggplot() + geom_density(aes(x = rinvgamma(100000,5,.5))) + theme_minimal(base_size = 14) + ggtitle("Inverse gamma with shape parameters .5, .5")
#
# ggplot() + geom_density(aes(x = rinvgamma(100000,3,1))) + theme_minimal(base_size = 14) + ggtitle("Inverse gamma with shape parameters 3, 1")

nimbleModBeta_all_covariates <- nimble::nimbleCode({


	# Loop through core observations (column 1 of input df "y") ----
	for (i in 1:N.core) {
		# y[i, 1] ~ dbeta(mean = plot_mu[plot_num[i], timepoint[i]],
		# 								sd = core_sd)
		y[i, 1] ~ T(dnorm(mean = plot_mu[plot_num[i], timepoint[i]],
											sd = core_sd), 0, 1)
	}

	# Plot-level process model ----
	for (p in 1:N.plot) {

		for (t in plot_start[p]) {

			# Plot means for first date
			logit(Ex[p, t]) ~ dnorm(mean = X_init, sd = .2)
			plot_mu[p, t] ~ T(dnorm(mean = Ex[p, t], sd = sd_init), 0, 1)

			# # Plot means for first date - for some reason this approach leads to huge CI,
			# with multimodal densities by 1
			# logit(Ex[p, t]) ~ dnorm(mean = X_init, sd = .2)
			# shape1[p, t] <- Ex[p, t] * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)
			# shape2[p, t] <- (1 - Ex[p, t]) * ((Ex[p, t] * (1-Ex[p, t]))/sigma^2 - 1)
			# logit(plot_mu[p, t]) ~ dLogitBeta(shape1[p, t], shape2[p, t])
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
	X_init <- .5 # Initial plot mean - will be constrained by data
	sd_init <- .01 # Initial process error
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

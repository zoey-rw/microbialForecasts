# Forecasting/uncertainty partitioning function
# model_outputs = model.outputs
# model.dat = model.dat
# plot_est = model_plot_est
# model_samples = model_samples
# plot_num = 107
# include = c("E","D","P")
# group = "ITS"
# N.beta = 6

forecast_fn <- function(model_outputs = NULL,
												plot_est = NULL,
												model.dat = NULL,
												model_samples = NULL,
												species_num = NULL, 
												plot_num = 1,
												group = "ITS",
												include = c("E"),
												N.beta = 6
){
	require(dplyr)
	if (is.null(species_num)) {
		species_num <- 1
	}
	
	
	single_plot_est <- plot_est %>% dplyr::filter(plot_num == !!plot_num)
	# Determine site number
	site_num <- single_plot_est %>% filter(date_num==max(date_num)) %>% select(site_num) %>% unlist()
	
	# Get last model estimate
	IC <- plot_est %>% dplyr::filter(plot_num == !!plot_num) %>% 
			filter(date_num==max(date_num, na.rm=T)) %>% select(`50%`) %>% unlist()
	
	### Initial condition uncertainty???
	# No: Initialize w/ last model value 
	if (!("I" %in% include)) {
		x <- IC
	} else {
	# Yes: Full distribution of estimates of last model value 
		# TO DO: write this code
	}
	
	# Parameters to decide forecast length
	last.obs <- model.dat$N.date
	NT <- ncol(model.dat$nspp) 
	NT.fcast <- NT - last.obs
	
	### Driver uncertainty???
	# No: Drivers are just means (with rows repeated as iterations)
	newdata <- array(NA, dim = c(Nmc, N.beta, NT))
	if (!("D" %in% include)) {
		for (time in last.obs:NT) {
			for (i in 1:Nmc) {
				newdata[i, ,time] <- c(model.dat$temp[site_num, time],
															 model.dat$mois[site_num, time],
															 model.dat$pH[plot_num],
															 model.dat$pC[plot_num],
															 model.dat$nspp[plot_num, time],
															 model.dat$rc_grass[plot_num, time])
			}
		}
	} else {
		# Yes: Drivers are sampled from means & SD's
		for (time in last.obs:NT) {
			newdata[, ,time] <- c(rnorm(Nmc, model.dat$temp[site_num, time], 
																	model.dat$temp_sd[site_num, time]),
														rnorm(Nmc, model.dat$mois[site_num, time], 
																	model.dat$mois_sd[site_num, time]),
														rnorm(Nmc, model.dat$pH[plot_num], 
																	model.dat$pH_sd[plot_num]),
														rnorm(Nmc, model.dat$pC[plot_num], 
																	model.dat$pC_sd[plot_num]),
														rep(model.dat$nspp[plot_num, time], Nmc),
														rep(model.dat$rc_grass[plot_num, time], Nmc))
		}
		
	}
	# Create random index so that all samples come from the same iteration
	random.index <- sample(nrow(newdata), replace=T)
	
	
	### Parameter uncertainty???
	# No: Parameters are just posterior medians (with rows repeated as iterations)
	if (!("P" %in% include)) {
		rho <- model_outputs %>% dplyr::filter(beta == "rho") %>% select(`50%`)
		beta <- model_outputs %>% dplyr::filter(!is.na(beta_num)) %>% select(`50%`)
		site_effect <- model_outputs %>% 
			dplyr::filter(!is.na(site_num)) %>% 
			dplyr::filter(site_num == !!site_num) %>% select(`50%`)
		# intercept <- model_outputs %>% dplyr::filter(grepl("intercept", rowname)) %>% 
		# 	#dplyr::filter(species_num == !!site_num) %>% 
		# 	select(`50%`)
		
		# Add dimensions
		# intercept <- t(rbind(rep(intercept[,1], Nmc)))
		site_effect <- t(rbind(rep(site_effect[,1], Nmc)))
		beta <- t(data.frame(rep(beta, Nmc)) )
		rho <- t(rbind(rep(rho[,1], Nmc)))
	} else {
		# Yes: full distribution of samples
		model_samples <- model_samples
		# intercept <- do.call(rbind, model_samples[,"intercept",drop=F])
		rho <- do.call(rbind, model_samples[,"rho",drop=F])
		beta <- do.call(rbind, model_samples[,grepl("beta", colnames(model_samples[[1]]))])
		site_effect <- do.call(rbind, model_samples[,grepl(paste0("effect[",site_num,"]"), colnames(model_samples[[1]]), fixed=T), drop=F])
	}

	### Process error???
	# No: add 0 process error
	if (!"E" %in% include) {
		sigma <- t(rbind(rep(0, Nmc)))
	} else {
	# Yes: full distribution of samples
		sigma <- do.call(rbind, model_samples[,grepl("sigma", colnames(model_samples[[1]])),drop=F])
		sigma <- unlist(lapply(sigma, function(y) lapply(y, function(x) 1/sqrt(x))))
		#sigma <- sigma^2
	}
	
	#### MAKE PREDICTIONS!!! ####
	
	## set up storage
	predict <- matrix(NA, Nmc, NT)
	## simulate
	for (time in (last.obs+1):NT) {
		Z  <- newdata[random.index, ,time]
		beta_sample  <- beta[random.index, ]
		
		# before everything was sampled
		# mu <- rho * log(x) + sum(Z[1,] * beta) + intercept + site_effect
		
		mu <- rho[random.index, ] * x + apply(Z * beta[random.index, ], 1, sum) + 
			#intercept[random.index, ] + 
			site_effect[random.index, ]
		x  <- rnorm(Nmc, mu, sigma)
		predict[, time] <- x
	}
	
	ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
	ci$mean <- apply(predict, 2, mean, na.rm=T)
	ci$sd <- apply(predict, 2, sd, na.rm=T)
	ci$date_num <- as.numeric(1:NT)
	colnames(ci)[1:5] <- paste0(colnames(ci)[1:5], "_", paste0(include, collapse="")) 
	return(ci)
}

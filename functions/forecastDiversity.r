# Function to forecast functional groups at all NEON sites, using parameters estimated from models (from summary script) 
# 
# taxon_name <- "copiotroph"
# 

# library(profvis)
# group = "ITS"
# Nmc <- 15000
# IC = .01
# NT = 80
# model_outputs = model_outputs
# scenario = "full_uncertainty_ITS"
# model_name = "all_covariates"
# test = TRUE
#time_period = "calibration"

fcast_all_plots <- function(model_val,
														model_outputs,
														scenario = "full_uncertainty_ITS",
														group = "ITS",
														model_name = "all_covariates",
														NT = 80,
														time_period = "calibration",
														taxon_name = NULL,
														Nmc = 15000,
														IC = .3){
	

	require(tidyverse)
	require(nimble)
	require(lubridate)
	
	# Create key of all available plots for given taxon
	#truth.plot.long1 <- model_cal$truth.plot.long %>% mutate(dates = fixDate(dateID))
	truth.plot.long <- model_val$truth.plot.long %>% mutate(dates = fixDate(dateID))
	#truth.plot.long <- cbind()
	val_key <- truth.plot.long %>% 
		select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
	val_plot_key <- val_key %>% select(-c(dateID, date_num)) %>% distinct() #%>% 
		#filter(plotID %in% model_cal$plotID)
	
	# val_key <- model_cal$truth.plot.long %>% mutate(dates = fixDate(dateID)) %>% 
	# 	select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct() 
	# val_plot_key <- val_key %>% select(-c(dateID, date_num)) %>% distinct()
		
	group_name <- truth.plot.long %>% select(name) %>% distinct() %>% unlist()
	cat(paste0("\nGroup: ", group_name, "\nScenario: ", scenario, "\nModel: ", model_name))
	
	
	if (model_name == "cycl_only"){
		N.beta = 2
	} else if(model_name == "all_covariates"){
		N.beta = 8
	} else cat("Missing specification of Nimble model.")
	
	
	model_summary <- model_outputs$summary_df 
	model_plot_est <- model_outputs$plot_est 
	#model_samples <- model_outputs$samples
	
	# Filter to scenario of interest
	model_out <- model_summary %>% 
		filter(scenario == !!scenario & time_period == !!time_period & 
					 	model_name == !!model_name) 

	
	ci_list <- list()
	# Loop through all plots
	for (i in 1:nrow(val_plot_key)) {
#	for (plot_num in val_plot_key$plot_num) {
	#	for (plot_num in c(1,20, 30, 50, 85, 87, 100, 150)){
		plotID <- val_plot_key[i,]$plotID 
		plot_num <- val_plot_key[i,]$plot_num 
		site_num <- val_plot_key[i,]$site_num 
		siteID <- val_plot_key[i,]$siteID 
		
		# 
		# 
		# site_num <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
		# 	select(site_num) %>% unique() %>% unlist()
		# siteID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
		# 	select(siteID) %>% unique() %>% unlist()
		# plotID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
		# 	select(plotID) %>% unique() %>% unlist()
		cat(paste0("\nForecasting for plot: ", plotID))#, "\nGroup: ", group, "\nScenario: ", scenario))
		
		
		
		# Parameters to decide forecast length
		start_date <- model_val$site_start[siteID]
		
		# Get covariates for site/plot/date
		covar <- array(NA, dim = c(Nmc, N.beta, NT))
		
		# Desperately tried to vectorize this loop but couldn't figure it out.
		# time_vec <- start_date:NT
		# covar[,,time_vec] <- c(rnorm(Nmc, model_val$temp[site_num, time_vec], 
		# 												 model_val$temp_sd[site_num, time_vec]),
		# 									 rnorm(Nmc, model_val$mois[site_num, time_vec], 
		# 									 			model_val$mois_sd[site_num, time_vec]),
		# 									 rnorm(Nmc, model_val$pH[plot_num], 
		# 									 			model_val$pH_sd[plot_num]),
		# 									 rnorm(Nmc, model_val$pC[plot_num], 
		# 									 			model_val$pC_sd[plot_num]),
		# 									 rep(model_val$nspp[plot_num, 80], Nmc),
		# 									 rep(model_val$relEM[plot_num, 80], Nmc))
		
		set.seed(1)
		
		
		mo <- month(as.Date(paste0(colnames(model_val$mois), "01"), format="%Y%m%d"))
		y_sin = sin((2*pi*mo)/12)
		y_cos = cos((2*pi*mo)/12)
		
		for (time in start_date:NT) {
			
			if (model_name == "all_covariates") {
			covar[,,time] <- c(rnorm(Nmc, model_val$temp[site_num, time],
															 model_val$temp_sd[site_num, time]),
												 rnorm(Nmc, model_val$mois[site_num, time],
												 			model_val$mois_sd[site_num, time]),
												 rnorm(Nmc, model_val$pH[plot_num],
												 			model_val$pH_sd[plot_num]),
												 rnorm(Nmc, model_val$pC[plot_num],
												 			model_val$pC_sd[plot_num]),
												 rep(model_val$relEM[plot_num, 80], Nmc),
												 rep(model_val$LAI[site_num, time], Nmc),
												 rep(y_sin[time], Nmc),
												 rep(y_cos[time], Nmc))
			} else if (model_name == "cycl_only") {
				covar[,,time] <- c(rep(y_sin[time], Nmc),
													 rep(y_cos[time], Nmc))
			}
		}
		
		# Check whether there's already an estimated site effect. If not, we'll sample!
		is_new_site <- ifelse(siteID %in% model_out$siteID, FALSE, TRUE)
		if (!is_new_site) {
			site_effect <- model_out %>% filter(model_out$siteID == !!siteID & 
																						grepl("site_effect", rowname) &
																										time_period == !!time_period & 
																						model_name == !!model_name) #%>% 
			site_effect_samp <- rnorm(Nmc, site_effect$Mean, site_effect$SD)
			#select(`50%`) %>% unlist()
			plot_est <- model_plot_est %>% 
				filter(group==!!group & 
							 	scenario == !!scenario & 
							 	plotID == !!plotID &
							 	time_period == !!time_period & 
							 	model_name == !!model_name) %>% 
				select(-c(date_num,plot_num, site_num, uncert))
			
		} else {
			# Sample from site effect variance
			site_effect_tau <- model_out %>% filter(grepl("sig$", rowname))
			# Convert precision to SD
			site_effect_tau <- unlist(lapply(site_effect_tau$Mean, 
																			 function(y) lapply(y, function(x) 1/sqrt(x))))
			site_tau <- mean(site_effect_tau)
			new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
			site_effect_samp <- unlist(new_site_effect)
			plot_est <- structure(list(siteID = siteID, 
																 plotID = plotID, dateID = NA, name = NA,
																  truth = NA, fcast_type = NA, fcast_period = NA,
																 `2.5%` = NA, `25%` = NA, `50%` = NA, `75%` = NA, `97.5%` = NA, 
																  scenario = scenario, group = group), row.names = 1L, class = "data.frame")
		}
		
		### Get other parameter estimates
		
		rho <- model_out[model_out$rowname=="rho",]
		rho_samp <-  rnorm(Nmc, rho$Mean, rho$SD)
		
		beta <- model_out[grepl("beta", model_out$rowname),]
		beta_samp <- apply(beta, 1, function(x) {
			rnorm(Nmc, as.numeric(x[["Mean"]]), as.numeric(x[["SD"]]))
		})
		
		int <- model_out[model_out$rowname=="intercept",]
		if (nrow(int) == 0) { # temporary till intercept is monitored
			int_samp <- rep(0, Nmc) 
			} else {
				int_samp <-  rnorm(Nmc, int$Mean, int$SD) 
			}
		
		sigma <- model_out[model_out$rowname=="sigma",]
		sigma_samp <-  rnorm(Nmc, sigma$Mean, sigma$SD)
		sigma_samp <- suppressWarnings(1/sqrt(sigma_samp))
#		sigma_samp <- suppressWarnings(unlist(lapply(sigma_samp, function(y) lapply(y, function(x) 1/sqrt(x)))))
		# Replace any NAs
		to_replace <- length(sigma_samp[is.na(sigma_samp)])
		sigma_samp[is.na(sigma_samp)] <- sample(na.omit(sigma_samp), to_replace)
		
		
		#cat(paste0("\nIC: ", IC, " sigma: ", sigma[1], " beta[1]: ", beta[1], " rho: ", rho, " site_effect: ", site_effect[1]))
		#### MAKE PREDICTIONS!!! ####
		### Initial condition uncertainty??? # Yes: input as argument
		x <- IC 
		## set up storage
		predict <- matrix(NA, Nmc, NT)
		## simulate
		for (time in (start_date):NT) {
		#for (time in (start_date+1):NT) {
			Z  <- covar[, ,time]
			#mu <- unlist(rho)  * logit(x) + apply(Z * beta[, ], 1, sum) + site_effect
			
			
			if (model_name == "cycl_only"){
				
				mu <- rho_samp * x + 
					beta_samp[,1]*Z[,1] + 
					beta_samp[,2]*Z[,2] + 
					site_effect_samp + int_samp
				
			} else if(model_name == "all_covariates"){

			mu <- rho_samp * x +
				beta_samp[,1]*Z[,1] +
				beta_samp[,2]*Z[,2] +
				beta_samp[,3]*Z[,3] +
				beta_samp[,4]*Z[,4] +
				beta_samp[,5]*Z[,5] +
				beta_samp[,6]*Z[,6] +
				beta_samp[,7]*Z[,7] +
				beta_samp[,8]*Z[,8] +
				site_effect_samp + int_samp

			} else cat("Missing specification of Nimble model.")
			# } else if(model_name == "all_covariates"){
			# 
			# mu <- rho_samp * x +
			# 	beta_samp[1,1]*Z[,1] +
			# 	beta_samp[1,2]*Z[,2] +
			# 	beta_samp[1,3]*Z[,3] +
			# 	beta_samp[1,4]*Z[,4] +
			# 	beta_samp[1,5]*Z[,5] +
			# 	beta_samp[1,6]*Z[,6] +
			# 	beta_samp[1,7]*Z[,7] +
			# 	beta_samp[1,8]*Z[,8] +
			# 	site_effect_samp[1] + int_samp
			# 
			# } else cat("Missing specification of Nimble model.")
				
				
			#at("\nmu: "); print(head(mu))
			# Add process error 
			#	x <- rep(expit(unlist(mu)), Nmc)
			#x <- expit(unlist(mu))
			#x  <- rnorm(Nmc, expit(mu), sigma_samp)
			x  <- rnorm(Nmc, mu, sigma_samp)
			predict[, time] <- x
			#cat("\nx: "); cat(x[1:2])
		}
		ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
		ci$mean <- apply(predict, 2, mean, na.rm=T)
		ci$sd <- apply(predict, 2, sd, na.rm=T)
		ci$date_num <- as.numeric(1:NT)
		#colnames(ci)[1:5] <- paste0(colnames(ci)[1:5], "_", paste0(include, collapse="")) 
		colnames(ci)[1:3] <- c("lo","med","hi")
		ci$plotID <- plotID
		ci$siteID <- siteID
		ci$group <- group
		ci$scenario <- scenario
		ci$model_name <- model_name
		ci$new_site <- ifelse(is_new_site, T, F) 
		
		ci <- left_join(ci, truth.plot.long, by = c("date_num", "plotID", "siteID"))
		plot_est$truth <- NULL
		ci$timepoint<- as.numeric(ci$timepoint)
		ci <- suppressWarnings(left_join(ci, plot_est))
		ci$time_period <- time_period
		# Check concurrence between model and formula estimates
		# plot(ci$med, ci$`50%`); abline(0,1)
		
#		print(tail(ci, 2))
		ci_list[[plotID]] <- ci
	}
	ci_allplots <- data.table::rbindlist(ci_list, fill = T)
	#plot(ci_allplots$med, ci_allplots$`50%`); abline(0,1)
	return(ci_allplots)
}



# 
# 
# ggplot(ci) +
# 	geom_line(aes(x = dates, y = mean), show.legend = F, linetype=2) +
# 	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')

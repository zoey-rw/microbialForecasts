# Function to forecast functional groups at all NEON sites, using parameters estimated from models (from summary script) 
# 
# taxon_name <- "copiotroph"
# # 
# N.beta = 6
# Nmc <- 15000
# IC = .01
# NT = 80
# model.dat <- model_val
# model_outputs = data_in
# scenario = "full_uncertainty"
# test = TRUE

fcast_all_plots <- function(model.dat,
														model_outputs,
														scenario = "full_uncertainty",
														NT = 80,
														taxon_name = NULL,
														Nmc = 15000,test = FALSE,
														N.beta = 6, IC = .3){
	require(tidyverse)
	#require(tidyr)
	require(nimble)
	
	# Create key of all available plots for given taxon
	truth.plot.long <- model.dat$truth.plot.long %>% mutate(dates = fixDate(dateID))
	val_key <- truth.plot.long %>% 
		select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
	val_plot_key <- val_key %>% select(-c(dateID, date_num)) %>% distinct()
	taxon_name <- truth.plot.long %>% select(name) %>% distinct() %>% unlist()
	cat(paste0("\nGroup: ", taxon_name, "\nScenario: ", scenario))
	
	
	
	model_summary <- model_outputs[[scenario]]$summary_df 
	model_plot_est <- model_outputs[[scenario]]$plot_est 
	model_samples <- model_outputs[[scenario]]$samples
	
	# Filter to taxon of interest
	model_out <- model_summary %>% 
		filter(taxon==!!taxon_name & scenario == !!scenario) 
	# Reshape
	params <- model_out %>% pivot_wider(id_cols = taxon, names_from = "rowname", values_from = c(Mean, SD))
	
	ci_list <- list()
	# Loop through all plots
	to_loop <- 1:nrow(val_plot_key)
	
	if(test) to_loop <- c(1:5,10,20,30,40,50,60,70,80,90)
	print(to_loop)
	for (plot_num in to_loop) {
		print(plot_num)
		site_num <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(site_num) %>% unique() %>% unlist()
		siteID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(siteID) %>% unique() %>% unlist()
		plotID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(plotID) %>% unique() %>% unlist()

		cat(paste0("\nForecasting for plot: ", plotID, "\nGroup: ", taxon_name, "\nScenario: ", scenario))

		
		
		# Parameters to decide forecast length
		# last.obs <- model.cal$N.date
		# NT <- ncol(model.cal$nspp) 
		# NT.fcast <- NT - last.obs
		#NT <- model.dat$N.date
		start_date <- model.dat$site_start[siteID]
		
		set.seed(1)
	
		# Get covariates for site/plot/date
		covar <- array(NA, dim = c(Nmc, N.beta, NT))
		for (time in start_date:NT) {
			covar[,,time] <- c(rnorm(Nmc, model_val$temp[site_num, time],
															 model_val$temp_sd[site_num, time]),
												 rnorm(Nmc, model_val$mois[site_num, time],
												 			model_val$mois_sd[site_num, time]),
												 rnorm(Nmc, model_val$pH[plot_num],
												 			model_val$pH_sd[plot_num]),
												 rnorm(Nmc, model_val$pC[plot_num],
												 			model_val$pC_sd[plot_num]),
												 rep(model_val$nspp[plot_num, 80], Nmc),
												 rep(model_val$relEM[plot_num, 80], Nmc))
		}
		
		
		# Check whether there's already an estimated site effect. If not, we'll sample!
		is_new_site <- ifelse(siteID %in% model_out$siteID, FALSE, TRUE)
		if (!is_new_site) {
			site_effect <- model_out %>% filter(model_out$siteID == !!siteID & 
																						grepl("site_effect", rowname)) #%>% 
			site_effect_samp <- rnorm(Nmc, site_effect$Mean, site_effect$SD)
				#select(`50%`) %>% unlist()
			plot_est <- model_plot_est %>% 
				filter(taxon==!!taxon_name & scenario == !!scenario & plotID == !!plotID) %>% select(-c(date_num,plot_num, site_num,taxon)) %>% rename(taxon = name)
			# 
			#site_effect_samp <- site_effect %>% select(Mean) %>% unlist()

		} else {
			# Sample from site effect variance
			site_effect_tau <- model_out %>% filter(grepl("sig$", rowname))
			# Convert precision to SD
			site_effect_tau <- unlist(lapply(site_effect_tau$Mean, 
																			 function(y) lapply(y, function(x) 1/sqrt(x))))
			site_tau <- mean(site_effect_tau)
			new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
			site_effect_samp <- unlist(new_site_effect)
			# plot_est <- structure(list(plot_num = plot_num, timepoint = NA, siteID = siteID, 
			# 													 plotID = plotID, dateID = NA, other = NA, 
			# 													 date_num = NA, name = NA, truth = NA, 
			# 													 site_num = site_num, dates = NA, 
			# 													 `2.5%` = NA, `25%` = NA, `50%` = NA, `75%` = NA, `97.5%` = NA, 
			# 													 rank = NA, scenario = scenario, 
			# 													 taxon = taxon_name), row.names = 1L, class = "data.frame")
		}
		
		### Get other parameter estimates
		
		rho <- model_out[model_out$rowname=="rho",]
		rho_samp <-  rnorm(Nmc, rho$Mean, rho$SD)

		beta <- model_out[grepl("beta", model_out$rowname),]
		beta_samp <- apply(beta, 1, function(x) {
			rnorm(Nmc, as.numeric(x[["Mean"]]), as.numeric(x[["SD"]]))
			})
		
		sigma <- model_out[model_out$rowname=="sigma",]
		sigma_samp <-  rnorm(Nmc, sigma$Mean, sigma$SD)
#		sigma_samp <- truncnorm::rtruncnorm(Nmc, mean = sigma$Mean, sd = sigma$SD, a = 0)
		
		
		int <- model_out[model_out$rowname=="intercept",]
		int_samp <-  rnorm(Nmc, int$Mean, int$SD)
		sigma_samp <- suppressWarnings(1/sqrt(sigma_samp))
		
		# Remove parameter uncertainty to better match model CI
		sigma_samp <- rep(1/sqrt(sigma$Mean), Nmc)
		beta_samp <- matrix(rep(beta$Mean, Nmc), ncol = 6, byrow = T)
		rho_samp <- rep(rho$Mean, Nmc)
		if(length(int$Mean)==0) { int_samp <- rep(0, Nmc) 
		} else int_samp <- rep(int$Mean, Nmc)
		
		
		#cat(paste0("\nIC: ", IC, " sigma: ", sigma[1], " beta[1]: ", beta[1], " rho: ", rho, " site_effect: ", site_effect[1]))
		#### MAKE PREDICTIONS!!! ####
		### Initial condition uncertainty??? # Yes: input as argument
		x <- mean(as.numeric(truth.plot.long$truth), na.rm=T)
		x <- as.numeric(head(truth.plot.long$truth[!is.na(truth.plot.long$truth)],1))
		if(is.na(x)) x <- IC
		## set up storage
		predict <- matrix(NA, Nmc, NT)
		## simulate
		#for (time in (start_date):NT) {
			for (time in (start_date+1):NT) {
			Z  <- covar[, ,time]
			#mu <- unlist(rho)  * logit(x) + apply(Z * beta[, ], 1, sum) + site_effect
			
			 
			mu <- rho_samp * logit(x) + beta_samp[,1]*Z[,1] +
				beta_samp[,2]*Z[,2] +
	beta_samp[,3]*Z[,3] +
	beta_samp[,4]*Z[,4] +
	beta_samp[,5]*Z[,5] +
	beta_samp[,6]*Z[,6] + site_effect_samp + int_samp
			#at("\nmu: "); print(head(mu))
			# Add process error 
			#	x <- rep(expit(unlist(mu)), Nmc)
			#x <- expit(unlist(mu))
			# x  <- rnorm(Nmc, expit(mu), sigma_samp)
			# x[x < 0] <- 0.00001
			x <- truncnorm::rtruncnorm(Nmc, mean = expit(mu), sd = sigma_samp, a = 0)
			#x  <- expit(rnorm(Nmc, mu, sigma_samp))
			predict[, time] <- x
			#cat("\nx: "); cat(x[1:2])
		}
		ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
		ci <- ci %>% mutate(mean = apply(predict, 2, mean, na.rm=T),
									sd = apply(predict, 2, sd, na.rm=T),
									date_num = as.numeric(1:NT),
									plotID = plotID,
									siteID = siteID,
									taxon_name = taxon_name,
									name = taxon_name,
									taxon = taxon_name,
									scenario = scenario,
									new_site = ifelse(is_new_site, T, F))
		colnames(ci)[1:3] <- c("lo","med","hi")
		
		# Merge with model estimates and observed values
		
		if (!is_new_site) {
			plot_est <- model_plot_est %>% 
				filter(taxon==!!taxon_name & plotID == !!plotID) %>% 
				select(-c(plot_num, site_num, truth, dateID, dates, other)) 
			ci <- left_join(ci, plot_est, by = intersect(colnames(ci), colnames(plot_est))) %>% 
				mutate(timepoint = NULL,
							 truth = NULL)
		}
		truth.plot.long$taxon_name <- truth.plot.long$name
		
		ci$timepoint <- NULL
		ci$truth <- NULL
		ci <- left_join(ci, truth.plot.long, by = c("date_num", "plotID", "siteID", "taxon_name"))
		
		# Check concurrence between model and formula estimates
	# plot(ci$med, ci$`50%`, col = ci$timepoint); abline(0,1)
	# plot(ci$med, ci$truth, col = ci$timepoint); abline(0,1)
	# plot(ci$truth, ci$`50%`, col = ci$timepoint); abline(0,1)
										
		
		#print(tail(ci, 2))
		ci_list[[plotID]] <- ci
	}
	ci_allplots <- data.table::rbindlist(ci_list, fill = T)
	return(ci_allplots)
}


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

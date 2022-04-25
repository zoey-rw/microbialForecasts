# Function to forecast functional groups at all NEON sites, using parameters estimated from model samples
# 
##### Use Nmc samples to make predictions, returns a dataframe with Nmc rows 
#

fg_fcast <- function(
	...
	# plotID = plotID,
	# covar = covar,
	# param_samples = param_samples,
	# ic = ic,
	# truth.plot.long = truth.plot.long,
	# model.inputs = model.inputs,
	# Nmc = 5000, 
	# plot_summary = plot_summary,
	# plot_start_date = plot_start_date
) {
	
	siteID <- substr(plotID, 1, 4)
	taxon_name <- truth.plot.long %>% select(name) %>% distinct() %>% unlist()
	
	print(paste0("Forecasting for group: ", taxon_name))
	
	# Check whether there's already an estimated site effect. If not, we'll sample!
	is_new_site <- ifelse(siteID %in% truth.plot.long$siteID, FALSE, TRUE)
	if (!is_new_site) {
		
		# (model.inputs has the full time-series)
		plot_obs <- model.inputs$truth.plot.long %>% filter(plotID==!!plotID) %>% select(-c(plot_num,site_num)) %>% rename(species = name) 
		site_num <- unique(truth.plot.long[truth.plot.long$siteID==siteID,]$site_num)
		site_param <- paste0("site_effect[", site_num, "]")
		site_effect <- 	param_samples[row_samples,] %>% select(!!site_param) %>% unlist()
		
		plot_est <- plot_summary %>% 
			filter(plotID == !!plotID) %>% 
			select(-c(plot_num, site_num, dateID)) %>% rename(species = name) 
		
		# Take initial condition & start forecast from last observed value if possible
		last_obs <- plot_est %>% filter(timepoint==max(timepoint)) 
		
		#last_obs <- truth.plot.long %>% filter(!is.na(truth)) %>% tail(1)
		plot_start_date <- last_obs$timepoint
		ic <- last_obs$`50%`
		
	} else {
		
		plot_obs <- model.inputs$truth.plot.long %>% filter(plotID==!!plotID) %>% select(-c(plot_num,site_num)) %>% rename(species = name) 
		# Sample from site effect variance
		site_effect_tau <- param_samples[row_samples,] %>% select(grep("sig$", colnames(.))) %>% unlist()
		# Convert precision to SD
		site_effect_tau <- unlist(lapply(site_effect_tau, 
																		 function(y) lapply(y, function(x) 1/sqrt(x))))
		site_tau <- mean(site_effect_tau)
		new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
		site_effect <- unlist(new_site_effect)
		
		# Take initial condition & start forecast from mean observed value if possible
		plot_start_date <- model.inputs$plot_index[plotID]
		ic <- mean(as.numeric(plot_obs$truth), na.rm = T)

	}
	
	### Get other parameter estimates
	### Rho
	rho <- param_samples[row_samples,] %>% select(grep("rho", colnames(.))) %>% unlist()
	### Betas
	betas <- param_samples[row_samples,] %>% select(grep("beta", colnames(.)))
	### Intercept
	intercept <- param_samples[row_samples,] %>% select(grep("intercept", colnames(.))) %>% unlist()
	### Process error 
	sigma_samp <- param_samples[row_samples,] %>% select(grep("sigma", colnames(.))) %>% unlist()
	sigma <- lapply(sigma_samp, function(y) lapply(y, function(x) 1/sqrt(x))) %>% unlist()
	
	sig_mean <- mean(sigma)
	
	# If the model only had sin/cosine, remove the other covariate data
	if (ncol(betas)==2) {
		if(ncol(covar)==8) {
			covar <- covar[,c(7,8),]
		}
	}
	
	predict <- matrix(NA, Nmc, NT)
	## simulate

	# In case initial condition wasn't set
	if(is.na(ic)) ic <- .0001
	
	x <- ic
	for (time in (plot_start_date+1):NT) {
		Z  <- covar[, ,time]
		mu <- rho * logit(x) + apply(Z * betas, 1, sum) + site_effect + intercept
		
		
		# Truncated to prevent negative values
		#x <- truncnorm::rtruncnorm(Nmc, mean = expit(mu), sd = sigma, a = 0)
		x <- truncnorm::rtruncnorm(Nmc, mean = expit(mu), sd = sigma, a = 0, b = 1)
		# Save to array
		predict[, time] <- x
	}
	
	ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
	ci <- ci %>% mutate(mean = apply(predict, 2, mean, na.rm=T),
											sd = apply(predict, 2, sd, na.rm=T),
											date_num = as.numeric(1:NT),
											plotID = plotID,
											siteID = siteID,
											taxon_name = taxon_name,
											rank = "functional_group",
											species = taxon_name,
											taxon = taxon_name,
											taxon_name = taxon_name,
											new_site = ifelse(is_new_site, T, F))
	colnames(ci)[1:3] <- c("lo","med","hi")
	
	if (!is_new_site) {
		plot_est_join <- plot_est %>% 
			select(-c(truth, timepoint, other)) 
		ci <- left_join(ci, plot_est_join, by = intersect(colnames(ci), colnames(plot_est_join)))
	}
	ci <- left_join(ci, date_key, by=c("date_num"))
	ci$dates <- fixDate(ci$dateID)

	ci <- left_join(ci, plot_obs, by = c("date_num", "plotID", "siteID", "species", "dateID"))

	return(ci)
}



# 
# ggplot(ci) +
# 	facet_grid(rows=vars(species), drop=T, scales="free") +
# 	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
# 	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')

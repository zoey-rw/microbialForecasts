# Function to forecast taxonomic groups at all NEON sites, using parameters estimated from models (from summary script) 
# 

# N.beta = 6
# Nmc <- 15000
# IC = .01
# model_outputs = data_in
# rank.name = "phylum_bac"
# NT = 80

fcast_all_plots <- function(model_val,
														model_outputs,
														rank.name = NULL,
														taxon_name = NULL,
														Nmc = 15000,
														NT = 80,
														test = F,
														N.beta = 6, IC = .3){
	require(dplyr)
	require(tidyr)
	require(nimble)
	
	# Create key of all available plots for given taxon
	truth.plot.long <- model_val$truth.plot.long
	val_key <- truth.plot.long %>% 
		select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
	val_plot_key <- val_key %>% select(-c(dateID, date_num)) %>% distinct()
	taxon_names <- truth.plot.long %>% select(species) %>% distinct() %>% unlist()
	cat(paste0("\nGroup: ", rank.name))
	
	# Filter to rank of interest, exclude site effects
	model_plot_est <- model_outputs$plot_est %>% filter(rank.name==!!rank.name) 
	model_out <- model_outputs$summary_df %>% filter(rank.name==!!rank.name) 
	
	ci_list <- list()
	# Loop through all plots
		if (test){
			to_loop <- c(1,20, 30, 50, 84, 85, 87, 100, 150) 
		} else {
			to_loop <-  1:nrow(val_plot_key)
		}
 for (plot_num in to_loop){
	#for (plot_num in c(1:2)){
			site_num <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(site_num) %>% unique() %>% unlist()
		siteID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(siteID) %>% unique() %>% unlist()
		plotID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% 
			select(plotID) %>% unique() %>% unlist()
		cat(paste0("\nForecasting for plot: ", plotID))
		
		
		# Parameters to decide forecast length
		start_date <- model_val$site_start[siteID]
		
		# Get covariates for site/plot/date
		covar <- array(NA, dim = c(Nmc, N.beta, NT))
		set.seed(1)
		for (time in start_date:NT) {
			covar[,,time] <- c(rnorm(Nmc, model_val$temp[site_num, time], 
															 model_val$temp_sd[site_num, time]),
												 rnorm(Nmc, model_val$mois[site_num, time], 
												 			model_val$mois_sd[site_num, time]),
												 rnorm(Nmc, model_val$pH[plot_num], 
												 			model_val$pH_sd[plot_num]),
												 rnorm(Nmc, model_val$pC[plot_num], 
												 			model_val$pC_sd[plot_num]),
												 rep(model_val$nspp[plot_num, time], Nmc),
												 rep(model_val$relEM[plot_num, time], Nmc))
		}
		# TODO: add this uncertainty back
		# for (time in start_date:NT) {
		# 	covar[,,time] <- c(rep(model_val$temp[site_num, time], Nmc),
		# 										 rep(model_val$mois[site_num, time], Nmc),
		# 										 rep(model_val$pH[plot_num], Nmc),
		# 										 rep(model_val$pC[plot_num], Nmc),
		# 										 rep(model_val$nspp[plot_num, time], Nmc),
		# 										 rep(model_val$relEM[plot_num, time], Nmc))
		# }
		# 
		# Check whether there's already an estimated site effect. If not, we'll sample!
		is_new_site <- ifelse(siteID %in% model_out$siteID, FALSE, TRUE)
		
		n.taxa <- length(taxon_names)
		all_tax_abs <- array(dim = c(Nmc, n.taxa, NT))
		
		for (i in 1:n.taxa){
			taxon_name <- taxon_names[i]
			print(paste0("Forecasting for taxon: ", taxon_name))
			taxon_summary <- model_out %>% filter(taxon == !!taxon_name)
			
		if (!is_new_site) {
			site_effect <- taxon_summary %>% filter(taxon_summary$siteID == !!siteID & 
																						grepl("site_effect", rowname)) #%>% 
			#site_effect_samp <- rnorm(Nmc, site_effect$Mean, site_effect$SD)
			
			# TODO: change back
			site_effect_samp <- site_effect %>% select(Mean) %>% unlist()
			
			plot_est <- model_plot_est %>% 
				filter(taxon==!!taxon_name & plotID == !!plotID & rank == rank.name) %>% 
				select(-c(plot_num, site_num, dateID)) #%>% rename(taxon = name)
			# set initial condition to first observed value
			IC <- as.numeric(na.omit(plot_est$truth)[1])
		} else {
			# Sample from site effect variance
			site_effect_tau <- model_out %>% filter(grepl("sig$", rowname))
			# Convert precision to SD
			site_effect_tau <- unlist(lapply(site_effect_tau$Mean, 
																			 function(y) lapply(y, function(x) 1/sqrt(x))))
			site_tau <- mean(site_effect_tau)
			new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
			site_effect <- unlist(new_site_effect)
		}
		
		### Get other parameter estimates
		rho <- taxon_summary[grep("rho", taxon_summary$rowname),]
		rho_samp <-  rnorm(Nmc, rho$Mean, rho$SD)
		beta <- taxon_summary[grepl("beta", taxon_summary$rowname),]
		int <- taxon_summary[grepl("int", taxon_summary$rowname),]
		# TODO: change back
		# beta_samp <- apply(beta, 1, function(x) {
		# 	rnorm(Nmc, as.numeric(x[["Mean"]]), as.numeric(x[["SD"]]))
		# })
		sigma <- taxon_summary[grep("sigma", taxon_summary$rowname),]
		
		sigma_samp <- rep(1/sqrt(sigma$Mean), Nmc)
		beta_samp <- matrix(rep(beta$Mean, Nmc), ncol = 6, byrow = T)
		rho_samp <- rep(rho$Mean, Nmc)
		int_samp <- rep(int$Mean, Nmc)
		
		# sigma_samp <-  rnorm(Nmc, sigma$Mean, sigma$SD)
		# # Convert tau to SD
		# sigma_samp <- suppressWarnings(1/sqrt(sigma_samp))
		# # Replace any NAs
		# to_replace <- length(sigma_samp[is.na(sigma_samp)])
		# sigma_samp[is.na(sigma_samp)] <- sample(na.omit(sigma_samp), to_replace)
		# just seeing if this helps..
		#sigma_samp <- 0
		
		#### MAKE PREDICTIONS!!! ####
		## set up storage
		predict <- matrix(NA, Nmc, NT)
		## simulate
		
		#for (time in (start_date):NT) {
#		time <- 13
		### Initial condition uncertainty??? # Yes: input as argument
		x <- IC
		
		for (time in (start_date+1):NT) {
			Z  <- covar[, ,time]
			
			mu <- rho_samp * log(x) + beta_samp[,1]*Z[,1] +
				beta_samp[,2]*Z[,2] +
				beta_samp[,3]*Z[,3] +
				beta_samp[,4]*Z[,4] +
				beta_samp[,5]*Z[,5] +
				beta_samp[,6]*Z[,6] + site_effect_samp + int_samp
			
			# Truncated to prevent negative values
			x <- truncnorm::rtruncnorm(Nmc, mean = exp(mu), sd = sigma_samp, a = 0, b = Inf)
		
			# Save to array
			all_tax_abs[,i,time] <- x
			}
		}
		
		# Make abundances relative within each MCMC sample
		all_tax_rel <- array(dim = c(Nmc, n.taxa, NT))
		for (time in (start_date+1):NT) {
			all_tax_rel[,,time] <- all_tax_abs[,,time]/rowSums(all_tax_abs[,,time])
		}
		colnames(all_tax_rel) <- taxon_names
		
		all_tax_single_plot <- list() 
		for (i in 1:n.taxa){
			taxon_name <- taxon_names[i]
			
			predict <- all_tax_rel[,i,]
		ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
		ci <- ci %>% mutate(mean = apply(predict, 2, mean, na.rm=T),
												sd = apply(predict, 2, sd, na.rm=T),
												date_num = as.numeric(1:NT),
												plotID = plotID,
												siteID = siteID,
												taxon_name = taxon_name,
												rank = rank.name,
												species = taxon_name,
												taxon = taxon_name,
												taxon_name = taxon_name,
												new_site = ifelse(is_new_site, T, F))
		colnames(ci)[1:3] <- c("lo","med","hi")
		
		
		if (!is_new_site) {
		plot_est <- model_plot_est %>% 
			filter(taxon==!!taxon_name & plotID == !!plotID & rank == rank.name) %>% 
			select(-c(plot_num, site_num, truth, dateID)) #%>% rename(taxon = name)
		ci <- left_join(ci, plot_est, by = intersect(colnames(ci), colnames(plot_est)))
		}
		ci$timepoint <- NULL
		ci$truth <- NULL
		ci <- left_join(ci, truth.plot.long, by = c("date_num", "plotID", "siteID", "species"))
		# Check concurrence between model and formula estimates
		# plot(ci$med, ci$`50%`); abline(0,1)
		#print(tail(ci))
		all_tax_single_plot[[i]] <- ci
		}
		ci_single_plot <-  do.call(rbind, all_tax_single_plot)
		ci_single_plot$dates <- fixDate(ci_single_plot$dateID)
		
		ci_list[[plotID]] <- ci_single_plot
		# plot(ci_single_plot$med, ci_single_plot$truth, col = as.factor(ci_single_plot$species)); abline(0,1)
		# plot(ci_single_plot$mean, ci_single_plot$truth, col = as.factor(ci_single_plot$species)); abline(0,1)
		
	}
	ci_allplots <- do.call(plyr::rbind.fill, ci_list)
	return(ci_allplots)
}


# 
# ggplot(ci_single_plot) +
# 	facet_grid(rows=vars(species), drop=T, scales="free") +
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

# Forecasting/uncertainty partitioning function for NEW SITES
# 
# plot_num = 11
# species_num = NULL
# model_outputs = model.outputs
# model.dat = model_val
# plot_est = model_plot_est
# model_samples = model_samples
# plot_num = 107
# group = "16S"
# N.beta = 6
# Nmc = 1000
# include = c("E","D","P")

newsite_forecast_fn <- function(model_outputs = NULL,
												model.dat = NULL,
												model.cal = NULL,
												model_samples = NULL,
												species_num = NULL, 
												plot_num = 1,
												group = "ITS",
												include = c("E"),
												val_key,
												N.beta = 6){
	require(dplyr)
	if (is.null(species_num)) species_num <- 1
	# Arbitrary starting point
	IC <- 0
	site_num <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% select(site_num) %>% unique() %>% unlist()
	siteID <- val_key %>% dplyr::filter(plot_num == !!plot_num) %>% select(siteID) %>% unique() %>% unlist()
	
	### Initial condition uncertainty???
	# No: Initialize w/ last model value 
	if (!("I" %in% include)) {
		x <- IC
	} else {
		# Yes: Full distribution of estimates of last model value 
		# TO DO: write this code
	}
	
	# Parameters to decide forecast length
	last.obs <- model.cal$N.date
	NT <- ncol(model.cal$nspp) 
	NT.fcast <- NT - last.obs
	
	start_date <- model.dat$site_start[siteID]
	
	### Driver uncertainty???
	# No: Drivers are just means (with rows repeated as iterations)
	newdata <- array(NA, dim = c(Nmc, N.beta, NT))
	if (!("D" %in% include)) {
		for (time in start_date:NT) {
			for (i in 1:Nmc) {
				newdata[i, ,time] <- c(model.dat$temp[site_num, time],
															 model.dat$mois[site_num, time],
															 model.dat$pH[plot_num],
															 model.dat$pC[plot_num],
															 model.dat$nspp[plot_num, time],
															 model.dat$rc_grass[plot_num, time])
			}}
	} else {
		# Yes: Drivers are sampled from means & SD's
		for (time in start_date:NT) {
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
		}}
	# Create random index so that all samples come from the same iteration
	random.index <- sample(nrow(newdata), replace=T)
	
	# Sample from site effect variance
	site_effect_tau <- do.call(rbind, model_samples[,grepl("sig$", colnames(model_samples[[1]])), drop=F])
	# Convert precision to SD
	site_effect_tau <- unlist(lapply(site_effect_tau, function(y) lapply(y, function(x) 1/sqrt(x))))
	site_tau <- mean(site_effect_tau)
	new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
	
	### Parameter uncertainty???
		# Yes: full distribution of samples
		model_samples <- model_samples
		rho <- do.call(rbind, model_samples[,"rho",drop=F])
		beta <- do.call(rbind, model_samples[,grepl("beta", colnames(model_samples[[1]]))])
		site_effect <- new_site_effect 

	### Process error???
	# No: add 0 process error
	if (!"E" %in% include) {
		sigma <- t(rbind(rep(0, Nmc)))
	} else {
		# Yes: full distribution of samples
		sigma <- do.call(rbind, model_samples[,grepl("sigma", colnames(model_samples[[1]])),drop=F])
		sigma <- unlist(lapply(sigma, function(y) lapply(y, function(x) 1/sqrt(x))))
	}
	
	#### MAKE PREDICTIONS!!! ####
	## set up storage
	predict <- matrix(NA, Nmc, NT)
	## simulate
	for (time in (start_date+1):NT) {
		Z  <- newdata[random.index, ,time]
		beta_sample  <- beta[random.index, ]
		mu <- rho[random.index, ] * x + apply(Z * beta[random.index, ], 1, sum) + site_effect[random.index, ]
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













library(coda)
library(forestplot)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)
library(dplyr)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/forecast_fn.r")

group = "ITS"
group = "16S"

data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")

if (group == "ITS"){
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
	model_samples <- data_in$samples$full_uncertainty_ITS
	
} else if (group == "16S") {
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
	model_samples <- data_in$samples$full_uncertainty_16S
}
# Prep model outputs
model_summary <- data_in$summary_df %>% dplyr::filter(grepl("full",  data_in$summary_df$scenario) & group == !!group)
model_plot_est <- data_in$plot_est %>% dplyr::filter(grepl("full", data_in$plot_est$scenario) & group == !!group) 


# Model truth/inputs
cal.rank.df <- div_in$cal
cal.rank.df$Shannon <- scale(cal.rank.df$Shannon, scale = F)
val.rank.df <- div_in$val
val.rank.df$Shannon <- scale(val.rank.df$Shannon, scale = F)
val.rank.df <- val.rank.df[which(!val.rank.df$siteID %in% c("ABBY","LAJA")),]
# Custom function for organizing model data.
model_dat <- prepDivData(rank.df = cal.rank.df, min.prev = 3, max.date = "20170101")
# Prep validation data (2017-onward)
model_val <- prepDivData(rank.df = rbind(cal.rank.df,val.rank.df), min.prev = 3, max.date = "20200801")



val_plot_key <- cbind.data.frame(plot_num = model_val$plot_num, 
																 plotID = model_val$plotID, 
																 siteID = model_val$siteID)
new_sites <- setdiff(model_val$siteID, model_dat$siteID)
val_plot_key <- val_plot_key %>% filter(siteID %in% new_sites)

truth.plot.long <- model_val$truth.plot.long
val_key <- truth.plot.long %>% select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
val_key <- val_key %>% filter(siteID %in% new_sites)
new_plots <- unique(val_key$plot_num)

# plot_num = 9
# species_num = NULL
model_outputs = model.outputs
# group = "16S"
# N.beta = 6
# Nmc = 1000

# plot list to work with. have to skip 1 and 2 because no chemistry data for ABBY.
# include = c("E","D","P")

out.data.list <- list()
for (p in new_plots) {
	print(p)
	### Function inputs (for testing)
	Nmc = 2000
	plot_num <- p
	# NT <- ncol(model_dat$nspp) 
	model.outputs = model_summary
	# model.dat = model_dat
	# plot_est = model_plot_est
	# model.samples = model_samples
	# # group = "ITS"
	# # include = c("E","D","P")
	# # include = c("E")
	# # N.beta = 6
	plotID <- unique(val_key[which(val_key$plot_num == plot_num),]$plotID)
	print(plotID)
	
	# Prep validation data (2017-onward) - specify plot
	full_obs <- model_val$truth.plot.long 
	full_obs$group <- group
	plot_obs <- full_obs %>% filter(plotID==!!plotID)


	ci_DEP <- newsite_forecast_fn(model_outputs = model.outputs,
																model.dat = model_val,
																model.cal = model_dat,
												model_samples = model_samples,
												plot_num = p,
												val_key = val_key,
												include = c("E","D","P"),
												group = group)
	ci_EP <- newsite_forecast_fn(model_outputs = model.outputs,
															 model.dat = model_val,
															 model.cal = model_dat,
											 model_samples = model_samples,
											 val_key = val_key,
											 group = group,
											 plot_num = p,
											 include = c("E","P"))
	# Add CI estimates
	plotting_data <- merge(plot_obs, ci_DEP, all=T)
	plotting_data <- merge(plotting_data, ci_EP, all=T) 
	plotting_data$plot_num <- plot_num
	
	
	tail(plotting_data, 10); dim(plotting_data)
	out.data.list[[plotID]] <- plotting_data
}


#### MUST RUN TWICE AND SAVE OUTPUT FOR EACH ####
out.data.bacteria <- do.call(rbind, out.data.list)
out.data.fungi <- do.call(rbind, out.data.list)

out.data <- rbind(out.data.fungi, out.data.bacteria)
out.data <- out.data[!is.na(out.data$date_num),]
out.data$fcast_period <- "Hindcast"
out.data$pretty_group <- ifelse(out.data$group=="16S","Bacteria","Fungi")
out.data$lo <- ifelse(is.na(out.data$`2.5%`), out.data$`2.5%_EDP`, out.data$`2.5%`)
out.data$hi <- ifelse(is.na(out.data$`97.5%`), out.data$`97.5%_EDP`, out.data$`97.5%`)
out.data$med <- ifelse(is.na(out.data$`50%`), out.data$`50%_EDP`, out.data$`50%`)

saveRDS(out.data, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data_newsites.rds")

out.data <- out.data[out.data$plotID %in% c("UKFS_031","NOGP_003","KONZ_024"),]

# Simple CI
output.plot <-  ggplot(out.data) +
	#	ggplot(out.data[out.data$group=="16S",]) +
	
	facet_grid(rows=vars(plotID), cols=vars(pretty_group), drop=T, scales="free_x", space="free") +
	geom_line(aes(x = date_num, y = `50%_EDP`), show.legend = F) + 
	geom_ribbon(aes(x = date_num, ymin = `2.5%_EDP`, ymax = `97.5%_EDP`, fill = pretty_group), 
							#fill = "darkblue", 
							alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Shannon diversity hindcasts at 3 new plots")) + #scale_x_date() +	
	theme_minimal(base_size=20) +
	#scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Shannon diversity anomaly") + 
	geom_point(aes(x = date_num, y = as.numeric(truth))) + xlim(c(25,65)) + xlab(NULL)

output.plot

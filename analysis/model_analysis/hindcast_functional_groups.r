# Create forecasts for functional groups
library(tidyverse)
model_outputs <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepFunctionalData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/forecastFunctional.r")

# Read in microbial abundances
cal <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds"))

# Alternative to parallel
#output.list <- list()
#for (rank in ranks.keep){
	
library(doParallel)
cl <- makeCluster(28,  outfile="")
registerDoParallel(cl)
# Run in parallel 
# k <- 12
output.list = foreach(k=1:length(keep_fg_names),
										 # .export=c("fcast_all_plots","ranks.keep","data_in"), .packages = c("tidyr","tidyverse","nimble"),
											.errorhandling = 'pass') %dopar% {
												
												library(tidyverse)
												
												# Only functional groups; subset to one rank.
												rank.name <- keep_fg_names[k]
												cal.rank.df <- cal[[rank.name]] 
												val.rank.df <- val[[rank.name]] 
												rank.df <- rbind(cal.rank.df, val.rank.df)	
												rank.df <- rank.df[!rank.df$siteID %in% c("ABBY","LAJA"),]
												
												# Model truth/inputs
												# Prep validation data
												model.cal.dat <- prepFunctionalData(rank.df = rank.df, 
																														min.prev = 3, max.date = "20170101")
												model_val <- prepFunctionalData(rank.df = 
																													rank.df,#[rank.df$plotID %in% model.cal.dat$plotID,], 
																												min.prev = 3, max.date = "20200101", full_timeseries = T)
												
												
												scenario_list <- list()
												for (scenario in c("full_uncertainty")){
													
													
													model_list <- list()
										
													for (model_name in c("cycl_only", "all_covariates")){
													
												#for (scenario_name in names(model_outputs)){
												out <- fcast_all_plots(taxon_name = rank.name, 
																																		 model_val = model_val, 
																																		 scenario = scenario,
																				model_outputs = model_outputs,
																				test=F,
																				model_name = model_name,
																				time_period = "calibration",
																				Nmc = 10000, IC = .01)
												
												model_list[[model_name]] <- out
												
													}
													out_scenario <- do.call(rbind, model_list)
													scenario_list[[scenario]] <- out_scenario
												}
												
												scenarios <- do.call(rbind, scenario_list)
												return(scenarios)
											}

out <- do.call(rbind, output.list)
out$fcast_period <- ifelse(out$dates < "2017-01-01", "calibration", "hindcast")
saveRDS(out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_fg.rds")


stopCluster(cl)



out.data_all <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_fg.rds")
#out.data <- scenarios


out.data <- out.data_all %>% filter(taxon_name == "copiotroph")
out.data <- out.data_all %>% filter(taxon_name == "ectomycorrhizal")

#out.data <- out.data[out.data$time_period != "refit",]
out.data <- out.data[out.data$plotID %in% c("HARV_001","BART_001","DSNY_001"),]
out.data <- out.data[out.data$plotID %in% c("UKFS_031","NOGP_003","KONZ_024"),]
out.data <- out.data[out.data$plotID %in% c("BART_001","CPER_002","HARV_004"),]

# Simple CI
output.plot <- ggplot(out.data) +
					facet_grid(rows=vars(plotID), drop=T, scales="free") +
	geom_line(data = out.data[out.data$fcast_period =="calibration",], 
						aes(x = dates, y = `50%`), show.legend = F) + 
	geom_line(data = out.data[out.data$fcast_period =="hindcast",], 
						aes(x = dates, y = med), show.legend = F) + 
	geom_ribbon(data = out.data[out.data$fcast_period =="calibration",],
							aes(x = dates, ymin = `2.5%`, ymax = `97.5%`, fill = fcast_period), alpha=0.6,  na.rm = F) +
	geom_ribbon(data = out.data[out.data$fcast_period =="hindcast",],
							aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6,  na.rm = F) +
	labs(title = "Example hindcasts (calibration: 2013-2016)") + 
				theme_bw()+
				scale_fill_brewer(palette = "Paired") + 
				theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"), 
							legend.position = "bottom",legend.title = element_text(NULL),
							plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) + 
				geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') 
				# facet_grid_paginate(#cols=vars(taxon), rows=vars(plotID), 
				# 	plotID~taxon,
				# 	ncol = 1, nrow = 8, page = i)

output.plot



ggplot(out.data) +
	facet_grid(rows=vars(plotID), cols=vars(scenario), drop=T, scales="free", space="free") +
	geom_line(aes(x = dates, y = `med`), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = `lo`, ymax = `hi`, fill = taxon_name), 
							fill = "darkblue",
							alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	#	ggtitle(paste0("Functional group hindcasts at 3 new plots")) + #scale_x_date() +	
	theme_minimal(base_size=20) +
	#scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Abundance") + 
	geom_point(aes(x = dates, y = as.numeric(truth))) #+ xlim(c(25,65)) + xlab(NULL)




out.data$calibration <- ifelse(out.data$dates < "2016-12-31", T, F)


fungi <- out.data %>% filter(taxon_name %in% c("wood_saprotroph", "plant_pathogen", "animal_pathogen","litter_saprotroph", "soil_saprotroph", "litter_saprotroph", "ectomycorrhizal") & scenario == "full_uncertainty")
ggplot(fungi[fungi$plotID=="CPER_001",]) +
	
	facet_grid(rows=vars(taxon_name), drop=T, scales="free"#, space="free"
	) +
	geom_line(aes(x = dates, y = `med`), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = `lo`, ymax = `hi`, fill = taxon_name), 
							#fill = "darkblue",
							alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Functional group hindcasts at plot CPER_001")) + #scale_x_date() +	
	#theme_minimal(base_size=20) +
	#scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Abundance") +  
	geom_point(aes(x = dates, y = as.numeric(truth))) + 
	scale_shape_manual(values = c(16, 21), name = NULL, labels = c("Calibration","Validation")) 

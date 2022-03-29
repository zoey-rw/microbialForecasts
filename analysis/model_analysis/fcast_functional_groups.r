# Create forecasts for functional groups
library(tidyverse)

k <- 38

model_outputs <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/refit_fg_summaries.rds")
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
output.list = foreach(k=1:lengths(keep_fg_names),
										 # .export=c("fcast_all_plots","ranks.keep","data_in"), .packages = c("tidyr","tidyverse","nimble"),
											.errorhandling = 'pass') %dopar% {
												
												rank.name <- keep_fg_names[k]
												cal.rank.df <- cal[[rank.name]] 
												val.rank.df <- val[[rank.name]] 
												rank.df <- rbind(cal.rank.df, val.rank.df)	
												rank.df <- rank.df[!rank.df$siteID %in% c("ABBY","LAJA"),]
												
												# Model truth/inputs
												# Prep validation data
												model_val <- prepFunctionalData(rank.df = rank.df, min.prev = 3, max.date = "20200801", full_timeseries = T)
												#scenario_name  <- "full_uncertainty"
												scenario_list <- list()
												for (scenario_name in names(data_in)){
												scenario_list[[scenario_name]] <- fcast_all_plots(taxon_name = rank.name, 
																																		 model.dat = model_val, 
																																		 scenario = scenario_name,
																				model_outputs = model_outputs,
																				N.beta = 6, Nmc = 10000, IC = .01)
												}
												scenarios <- do.call(rbind, scenario_list)
												return(scenarios)
											}

out <- do.call(rbind, output.list)
saveRDS(out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fcast_fg.rds")


stopCluster(cl)




fcast_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fcast_fg.rds")


out.data <- fcast_data %>% filter(taxon_name == "copiotroph")
out.data <- out.data[out.data$plotID %in% c("HARV_001","BART_001","KONZ_024"),]
out.data <- out.data[out.data$plotID %in% c("UKFS_031","NOGP_003","KONZ_024"),]
out.data <- out.data[out.data$plotID %in% c("BART_001","CPER_001","DSNY_001"),]

# All points
output.plot <-  ggplot(out.data) +
	facet_grid(rows=vars(plotID), cols=vars(scenario), drop=T, scales="free_x", space="free") +
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

output.plot



ggplot(out.data[out.data$siteID=="HARV",]) +
	facet_grid(rows=vars(plotID), drop=T, scales="free", space="free") +
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








out.data <- fcast_data %>% filter(taxon_name %in% c("nitrification", "cellulose_complex", "cellobiose_complex","copiotroph", "heat_stress", "litter_saprotroph", "ectomycorrhizal"))

out.data <- fcast_data %>% filter(taxon_name %in% c("wood_saprotroph", "plant_pathogen", "animal_pathogen","litter_saprotroph", "soil_saprotroph", "litter_saprotroph", "ectomycorrhizal"))

#ggplot(out.data[out.data$siteID=="HARV",]) +
ggplot(out.data[out.data$plotID=="CPER_001",]) +
	
	facet_grid(rows=vars(taxon_name), drop=T, scales="free"#, space="free"
						 ) +
	geom_line(aes(x = dates, y = `med`), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = `lo`, ymax = `hi`, fill = taxon_name), 
							#fill = "darkblue",
							alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Functional group models at CPER_001")) + #scale_x_date() +	
	#theme_minimal(base_size=20) +
	#scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Abundance") + 
	geom_point(aes(x = dates, y = as.numeric(truth)), shape = 21)



#Just looking at the real values for nitrification
ggplot(plot_est) + geom_jitter(aes(x = siteID, y = as.numeric(truth))) + ggtitle("nitrification abundances")
not_na <- out.data[!is.na(out.data$truth),]
sort(table(not_na$plotID))

harv <- out.data[out.data$siteID=="HARV",]

harv_nitr <- harv %>% filter(taxon_name %in% c("nitrification"))




plot_est <- data_in$full_uncertainty$plot_est %>% filter(taxon %in% c("nitrification"))
sum_nitr <- data_in$full_uncertainty$summary_df %>% filter(taxon %in% c("nitrification")) 

nitr_site <- sum_nitr %>% dplyr::filter(grepl("site", sum_nitr$rowname))

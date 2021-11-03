# Run on node with 10+ cores

library(tidyverse)
library(doParallel)

model_outputs <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/allrank_summaries.rds")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/forecastTaxonomic.r")
#source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/forecastFunctional.r")


# Read in microbial abundances
cal <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds"))

# Subset to one rank.
# k <- 3
# rank.name = "order_bac"

rank.list <- list()
cl <- makeCluster(10,  outfile="")
registerDoParallel(cl)
# Run in parallel 
output.list = foreach(k=1:10, .errorhandling = 'pass') %dopar% {

rank.name <- tax_names[[k]]
cal.rank.df <- cal[[rank.name]] 
val.rank.df <- val[[rank.name]] 
rank.df <- rbind(cal.rank.df, val.rank.df)	
rank.df <- rank.df[!rank.df$siteID %in% c("ABBY","LAJA"),]
# Model truth/inputs
# Prep validation data
model_val <- prepTaxonomicData(rank.df = rank.df, min.prev = 3, max.date = "20200801", 	full_timeseries = T)

out <- fcast_all_plots(model_val = model_val, 
																					model_outputs = model_outputs,
																					rank.name = rank.name, test = F,
																					#scenario = scenario, 
																					Nmc = 5000, IC = .01)
return(out)
}
all.ranks <- data.table::rbindlist(output.list, fill = T)

all.ranks$fcast_period <- ifelse(all.ranks$dates < "2017-01-01", "calibration", "hindcast")
all.ranks$fcast_period <- ifelse(all.ranks$timepoint < 44, "calibration", "hindcast")

saveRDS(all.ranks, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_tax.rds")


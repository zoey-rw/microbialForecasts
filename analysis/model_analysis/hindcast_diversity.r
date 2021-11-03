# Create hindcasts for diversity scenarios
library(tidyverse)
library(hrbrthemes)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/forecastDiversity.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

group = "ITS"
group = "16S"

data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")

# 
scenario_list <- list()

for (scenario in div_scenarios[1:8]){

if (grepl("ITS", scenario)){
	group <- "ITS"
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
} else {
	group <- "16S"
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
}

# Model truth/inputs
cal.rank.df <- div_in$cal
cal.rank.df$Shannon <- scale(cal.rank.df$Shannon, scale = F)
val.rank.df <- div_in$val
val.rank.df$Shannon <- scale(val.rank.df$Shannon, scale = F)
val.rank.df <- val.rank.df[which(!val.rank.df$siteID %in% c("ABBY","LAJA")),]

# Custom function for organizing model data.
model_val <- prepDivData(rank.df = rbind(cal.rank.df,val.rank.df), min.prev = 3, max.date = "20200801", full_timeseries = T)

scenario_list[[scenario]] <- fcast_all_plots(model_val = model_val, 
																						 model_outputs = data_in,
																										scenario = scenario, 
																						 group = group,
																						 Nmc = 5000, IC = .01)
}
scenarios <- do.call(rbind, scenario_list)




scenarios$fcast_period <- ifelse(scenarios$dates < "2017-01-01", "calibration", "hindcast")
scenarios$pretty_group <- ifelse(scenarios$group=="16S","Bacteria","Fungi")

saveRDS(scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_div.rds")


scenarios <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_div.rds")
out.data <- scenarios[scenarios$plotID %in% c("UKFS_031","NOGP_003","KONZ_024"),]
out.data <- scenarios[scenarios$plotID %in% c("BART_001","CPER_002","HARV_004"),]

# Simple CI
#output.plot <-  ggplot(out.data) +
		ggplot(out.data[out.data$group=="16S",]) +
	
	facet_grid(rows=vars(plotID), cols=vars(scenario), drop=T, scales="free", space="free") +
	geom_line(aes(x = date_num, y = med), show.legend = F) + 
	geom_ribbon(aes(x = date_num, ymin = lo, ymax = hi, fill = pretty_group), 
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

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

group_out <- list()
for (group in c("16S","ITS")){

data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")

if (group == "ITS"){
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
	model_samples <- data_in$samples$full_uncertainty_ITS
} else if (group == "16S") {
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
	model_samples <- data_in$samples$full_uncertainty_16S
	
}
# Prep model outputs
model_summary <- data_in$summary_df %>% 
 dplyr::filter(grepl("full",  data_in$summary_df$scenario) & group == !!group)


	model_plot_est <- data_in$plot_est %>% 
	 dplyr::filter(grepl("full",  data_in$plot_est$scenario) & group == !!group)
	
cal.rank.df <- div_in$cal
cal.rank.df$Shannon <- scale(cal.rank.df$Shannon, scale = F)
val.rank.df <- div_in$val
val.rank.df$Shannon <- scale(val.rank.df$Shannon, scale = F)
val.rank.df <- val.rank.df[which(!val.rank.df$siteID %in% c("ABBY","LAJA")),]

# Custom function for organizing model data.
model_dat <- prepDivData(rank.df = cal.rank.df, min.prev = 3, max.date = "20170101")
# Prep validation data (2017-onward)
model_val <- prepDivData(rank.df = rbind(cal.rank.df,val.rank.df), min.prev = 3, max.date = "20200801")
# Create date/plot keys
plot_key <- model_plot_est[,c("plot_num","plotID")] %>% distinct()
date_key <- model_plot_est %>% pad(end_val = as.Date("2020-09-01")) %>% select(date_num, dates) %>% distinct()
date_key$date_num <- 1:nrow(date_key)
date_key <- date_key[1:80,]
date_key$dateID <- gsub("-","",substr(date_key$dates, 1, 7))

out.data.list <- list()
for (p in unique(plot_key$plot_num)) {
		print(p)
Nmc = 2000
plot_num <- p
model.outputs = model_summary
model.dat = model_dat
plot_est = model_plot_est
model.samples = model_samples
plotID <- plot_key[which(plot_key$plot_num == plot_num),]$plotID

# Prep validation data (2017-onward) - specify plot
full_obs <- model_val$truth.plot.long 
full_obs$group <- group
plot_obs <- full_obs %>% filter(plotID==!!plotID)
plot_obs$plot_num <- plot_key[match(plot_obs$plotID, plot_key$plotID),]$plot_num
plot_obs$date_num <- date_key[match(plot_obs$dateID, date_key$dateID),]$date_num

# Combine model predictions and truths, for later plotting
single_plot_est <- model_plot_est %>% dplyr::filter(plot_num == !!plot_num)
single_plot_est$date_num <- as.integer(single_plot_est$timepoint)
single_plot_est <- merge(single_plot_est, date_key, all.y=T)
single_plot_est$plot_num <- plot_num
plot_est_truth <- merge(single_plot_est, plot_obs[,c("date_num","plot_num","truth")], all=T, by=c("date_num","plot_num")) 
plot_est_truth$group <- group


ci_DEP <- forecast_fn(model_outputs = model.outputs,
											model.dat = model.dat,
											plot_est = model_plot_est,
											model_samples = model_samples,
											plot_num = p,
											include = c("E","D","P"),
											group = group)
ci_EP <- forecast_fn(model_outputs = model.outputs,
										 model.dat = model.dat,
										 plot_est = model_plot_est,
										 model_samples = model_samples,
										 group = group,
										 plot_num = p,
										 include = c("E","P"))
# Add CI estimates
plotting_data <- merge(plot_est_truth, ci_DEP, all=T)
plotting_data <- merge(plotting_data, ci_EP, all=T) 
plotting_data$plot_num <- plot_num

plotting_data$truth <- ifelse(is.na(plotting_data$truth.x), 
															plotting_data$truth.y, plotting_data$truth.x)
plotting_data <- plotting_data %>% filter(plot_num==!!plot_num)
tail(plotting_data, 40); dim(plotting_data)
plotting_data$plotID <- plotID
out.data.list[[plotID]] <- plotting_data
}
group_out[[group]] <- do.call(rbind, out.data.list)
}

#### MUST RUN TWICE AND SAVE OUTPUT FOR EACH ####
# out.data.bacteria <- do.call(rbind, out.data.list)
# out.data.fungi <- do.call(rbind, out.data.list)
out.data <- rbind(group_out[["ITS"]], group_out[["16S"]])
# out.data <- rbind(out.data.fungi, out.data.bacteria)
out.data <- out.data[!is.na(out.data$date_num),]
out.data[which(out.data$`50%`==0),]$`50%` <- NA
out.data[which(out.data$`2.5%`==0),]$`2.5%` <- NA
out.data[which(out.data$`97.5%`==0),]$`97.5%` <- NA
out.data$fcast_period <- NA
out.data[which(out.data$dates < "2017-01-01"),]$fcast_period <- "Calibration"
out.data[which(out.data$dates >= "2017-01-01"),]$fcast_period <- "Hindcast"
out.data$pretty_group <- ifelse(out.data$group=="16S","Bacteria","Fungi")
out.data$lo <- ifelse(is.na(out.data$`2.5%`), out.data$`2.5%_EDP`, out.data$`2.5%`)
out.data$hi <- ifelse(is.na(out.data$`97.5%`), out.data$`97.5%_EDP`, out.data$`97.5%`)
out.data$med <- ifelse(is.na(out.data$`50%`), out.data$`50%_EDP`, out.data$`50%`)

saveRDS(out.data, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data.rds")

out.data <- out.data[out.data$plotID %in% c("DSNY_001",
																						"HARV_001","ORNL_001"#,
																					#	"TALL_003","WOOD_003"
																						),]

# Simple CI
output.plot <-  ggplot(out.data) +
	#	ggplot(out.data[out.data$group=="16S",]) +
	
	facet_grid(rows=vars(plotID), cols=vars(pretty_group), drop=T, scales="free", space="free") +
	geom_line(aes(x = dates, y = med), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, 
									fill = pretty_group, alpha = fcast_period
									), 
							#fill = "darkblue", 
							#alpha=0.4, 
							show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Shannon diversity hindcasts at 3 plots")) + #scale_x_date() +	
	scale_alpha_discrete(range = c(.7, .3)) +
	theme_minimal(base_size=20) +
#	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Shannon diversity anomaly") + 
	geom_point(aes(x = dates, y = as.numeric(truth))) +
	geom_text(aes(x = as.Date("2017-01-01"), y = 2, label = "Hindcast"), 
						size = 5, hjust = "left") +
	geom_text(aes(x = as.Date("2015-01-01"), y = 2, label = "Calibration"), 
						size = 5, hjust = "left") + xlab(NULL)

output.plot


output.plot <-  ggplot(out.data) +
#	ggplot(out.data[out.data$group=="16S",]) +
	
	facet_grid(rows=vars(plotID), cols=vars(pretty_group), drop=T, scales="free", space="free") +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`), 
							fill = "darkblue", alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Shannon diversity hindcasts for ")) + #scale_x_date() +	
		theme_minimal(base_size=18) +
	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
		ylab("Shannon diversity") + 
	geom_line(aes(x = dates, y = `50%_EDP`)) +
	geom_ribbon(aes(x = dates, ymin = `2.5%_EDP`, ymax = `97.5%_EDP`), fill = "darkblue", alpha=0.4, show.legend = F) +
	#geom_ribbon(aes(x = dates, ymin = `2.5%_EP`, ymax = `97.5%_EP`), fill = "darkblue", alpha=0.4, show.legend = F) +
	geom_point(aes(x = dates, y = as.numeric(truth))) +
		geom_text(aes(x = as.Date("2017-01-01"), y = 3, label = "Hindcast"), 
							size = 5, hjust = "left") +
		
		geom_text(aes(x = as.Date("2015-01-01"), y = 3, label = "Calibration"), 
							size = 5, hjust = "left")

output.plot


# run twice
output.plots_fungi <- output.plot
output.plots_bacteria <- output.plot

library(grid)

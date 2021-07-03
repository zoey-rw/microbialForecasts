# Summarize outputs from diversity models
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(coda)
library(dplyr)
library(tibble)

# # Read in samples for visualization
read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")

plot_est_df_all <- list()
summary_df_all <- list()
for (scenario in names(read_in$sample.list)){
	
	cat(paste0("\nSummarizing ", scenario, "..."))
	samples <- read_in$sample.list[[scenario]]
	param_summary <- read_in$param.summary.list[[scenario]]
	plot_summary <- read_in$plot.summary.list[[scenario]]
	truth.plot.long <- read_in$metadata.list[[scenario]]$model_data 

	
	# Calculate plot means and actual values per rank
	plot_est <- plot_summary[[2]]
	pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)),
					 rank = "bacterial_richness")
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d"),
					 timepoint = as.integer(timepoint))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	
	means <- param_summary[[2]]
	beta_out <- means[grep("beta", rownames(means)),]
	rho_out <- means[grep("rho", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% mutate(beta = "rho")
	intercept_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	
	
	# Get beta sizes per rank
	beta_out <-  beta_out %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(beta_num,
												 "1" = "Temperature",
												 "2" = "Moisture",
												 "3" = "pH",
												 "4" = "pC",
												 "5" = "Plant species richness",
												 "6" = "% grasses",
												 "7" = "sin(month)",
												 "8" = "cosin(month)"))
	beta_out <- plyr::rbind.fill(beta_out, rho_out)
	beta_out$significant <- ifelse(beta_out$`2.5%` < 0 & beta_out$`97.5%` < 0 |
																 	beta_out$`2.5%` > -0 & beta_out$`97.5%` > -0,
																 1, 0)
	beta_out$effSize <- abs(beta_out$`50%`)
	
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID)
	
	summary_df <- plyr::rbind.fill(site_eff_out, beta_out, intercept_out, sigma_out)
	summary_df$scenario <- scenario
	summary_df_all[[scenario]] <- summary_df
	
	allplots$scenario <- scenario
	plot_est_df_all[[scenario]] <- allplots
}

plot_est_df_all <- do.call(rbind, plot_est_df_all)
summary_df_all <- do.call(rbind, summary_df_all)
summary_df_all <- summary_df_all %>% tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"))
summary_df_all$scenario <- paste(summary_df_all$uncert1, summary_df_all$uncert2)

saveRDS(list(plot_est = plot_est_df_all,
						 summary_df = summary_df_all,
						 samples = read_in$sample.list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")




data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")
# ## GGPLOT
sum.all <- data_in$summary_df

beta_out <- sum.all[which(!is.na(sum.all$beta)),]
# By rank with every taxon - cluttered
ggplot(data=beta_out,
			 aes(x = reorder(beta, effSize),y = effSize)) +
	facet_grid(rows = vars(group), cols = vars(scenario), drop = T) +
	geom_point(aes(shape = as.factor(significant), color = beta), size = 4) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Not significant","Significant")) 

# View fcast w/ confidence intervals for no and full uncertainties
allplots <- plot_est_df_all %>% filter(scenario %in% c("no_uncertainty_ITS", "full_uncertainty_ITS","spatial_uncertainty_ITS","temporal_uncertainty_ITS"))
plot_data <- allplots %>% filter(plotID == "HARV_001")
plot_data$observed <- ifelse(is.na(plot_data$truth), "Estimated", "Observed")
# Fcast
ggplot(pred.plot1, aes(x = dates)) +
	#facet_wrap(~scenario) + 
	geom_line(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	#geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, show.legend = F) +
	geom_ribbon(data = pred.plot1[pred.plot1$scenario == "no_uncertainty_ITS",], 
							aes(ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.4, show.legend = F) +
	geom_ribbon(data = pred.plot1[pred.plot1$scenario == "full_uncertainty_ITS",], ,
							aes(ymin = `25%`, ymax = `75%`),fill = "darkblue", alpha=0.4, show.legend = F) +
	geom_point(aes(y = as.numeric(truth))) + theme_ipsum(base_size = 14, strip_text_size = 22) + ggtitle(paste0("Shannon diversity at ", unique(pred.plot1$plotID))) + scale_x_date() +	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) 


full_uncert <- plot_data %>% filter(scenario == "full_uncertainty_ITS")
no_uncert <- plot_data %>% filter(scenario == "no_uncertainty_ITS")
spat_uncert <- plot_data %>% filter(scenario == "spatial_uncertainty_ITS")
temp_uncert <- plot_data %>% filter(scenario == "temporal_uncertainty_ITS")

ggplot() +
	#facet_wrap(~scenario) + 
	geom_line(data = pred.plot1[pred.plot1$scenario == "no_uncertainty_ITS",],
						aes(x = dates, y = `50%`), show.legend = F) + 
	geom_line(data = pred.plot1[pred.plot1$scenario == "full_uncertainty_ITS",],
						aes(x = dates, y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	#geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, show.legend = F) +
	geom_ribbon(data = pred.plot1[pred.plot1$scenario == "no_uncertainty_ITS",], 
							aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.4, show.legend = F) +
	geom_ribbon(data = pred.plot1[pred.plot1$scenario == "full_uncertainty_ITS",],
							aes(x = dates, ymin = `25%`, ymax = `75%`),fill = "darkblue", alpha=0.4, show.legend = F) +
	geom_point(data = pred.plot1[pred.plot1$scenario == "full_uncertainty_ITS",],
						 aes(x = dates, y = as.numeric(truth))) + 
	theme_ipsum(base_size = 14, strip_text_size = 22) + ggtitle(paste0("Shannon diversity at ", unique(pred.plot1$plotID))) + scale_x_date() +	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) 



# Side by side...
ggplot(pred.plot1, aes(x = dates)) + 
	facet_wrap(~scenario) + 
	geom_line(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, show.legend = F) +
	
	geom_ribbon(
		aes(ymin = `25%`, ymax = `75%`),fill = "darkblue", alpha=0.4, show.legend = F) +
	geom_point(aes(y = as.numeric(truth))) + theme_ipsum(base_size = 14, strip_text_size = 22) + ggtitle(paste0("Shannon diversity at ", unique(pred.plot1$plotID))) + scale_x_date() +	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm"))  + ylim(c(4,6))










output.plot <-  ggplot() +
	geom_line(data = no_uncert,
						aes(x = dates, y = `50%`), show.legend = F) + 
	geom_ribbon(data = no_uncert, 
							aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.4) +
	geom_ribbon(data = spat_uncert, 
							aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "green", alpha=0.4) +
	geom_point(data = no_uncert,
						 aes(x = dates, y = as.numeric(no_uncert$truth))) + 
	theme_ipsum(base_size = 14, strip_text_size = 22) + 
	ggtitle(paste0("Shannon diversity at CPER_004")) + scale_x_date() +	
	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) + ylim(c(4,6)) +
	geom_ribbon(data = full_uncert,
							aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "darkblue", alpha=0.4) +
	scale_fill_brewer(palette = "Paired")
output.plot

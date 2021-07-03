library(gganimate)
library(reshape2)
library(scales)
library(viridis)
library(hrbrthemes)
library(coda)
library(tidyverse)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data_construction/microbe/prep_diversity_fn.r")


anim_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/animation_samples.rds")
div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S_full.rds")
model.dat <- prepDivData(rank.df = div_in, min.prev = 5)


full_summary.list <- anim_data[[2]]
animation <- list()
for (i in 1:6){
	
	means <- full_summary.list[[i]]$quantiles
	plot_mu_out <- means[grep("plot_mu", rownames(means)),]
	beta_out <- means[grep("beta", rownames(means)),]
	
	pred.plot <- plot_mu_out %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","date_num")) %>%
		mutate(plot_num = as.numeric(gsub("plot_mu\\[", "", plot_num)),
					 date_num = as.numeric(gsub("\\]", "", date_num)))
	# 
	truth.plot.long <- model.dat$plot.truth %>% as.data.frame() %>%
		mutate(truth = as.numeric(as.character(truth)))
	
	date_limit <- switch(i,
											 "1" = 34,
											 "2" = 39,
											 "3" = 44,
											 "4" = 49,
											 "5" = 54,
											 "6" = 68)
	
	truth.plot.long$dates <- as.Date(paste0(truth.plot.long$dateID, "01"), "%Y%m%d")
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","date_num"), all=T)
	allplots[allplots$date_num > date_limit,]$truth <- NA

	real_date_limit <- unique(allplots[allplots$date_num == date_limit,]$dates)
	print(real_date_limit)
	allplots$date_limit <- real_date_limit
	# 
	# 
	allplots$run <- i
	
	animation[[i]] <- allplots  %>% arrange(plotID, dateID)
	
}

# Remove zero estimates (dumb method)
animation.all <- do.call(rbind, animation)
zero_est <- which(animation.all$`97.5%`==0 & animation.all$`2.5%`==0)
animation.all <- animation.all[!rownames(animation.all) %in% zero_est,] 

# Choose plots for visualizing (most sampling times)
not_na <- animation.all[!is.na(animation.all$truth),] %>% distinct(plotID, dateID)
top_plots <- sort(table(not_na$plotID), decreasing = T)[1:9]
top_plots_animate <- animation.all[animation.all$plotID %in% names(top_plots),]

pred.plot1 <- animation.all %>% filter(plotID == "CPER_001")
p <- ggplot(top_plots_animate, aes(x = dates)) +
	geom_line(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = siteID), alpha=0.2, show.legend = F) +
	geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = siteID), alpha=0.4, show.legend = F) +
	geom_point(aes(y = truth)) + 
	#theme_ipsum() + 
	coord_cartesian(ylim = c(2, 10), # This focuses the x-axis on the range of interest
									clip = 'off') +  
	theme_ipsum(base_size = 14, strip_text_size = 16) + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) + 
	#ggtitle(unique(pred.plot1$plotID)) + 
	geom_vline(aes(xintercept = date_limit), linetype = 2) + 
scale_x_date() + facet_wrap(~plotID) +
	geom_text(aes(x = date_limit + 120, y = 6.5, label = "Forecast"), 
						size = 5, hjust = "left", angle = 90)
	 

p2 <- p + transition_time(run) + 
	ggtitle(paste0("Shannon diversity at plot", unique(pred.plot1$plotID)),
					subtitle = 'Frame {frame} of {nframes}')


p3 <- p +  transition_states(run,
														 transition_length = .01,
														 state_length = .005) + 
	ggtitle("Shannon diversity at NEON plots")  + enter_fade()
p3


animate(p3, height = 4, width = 5, units = "in", res = 150)
anim_save(animation = p3, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/anim_diversity1.gif", width = 1000, height = 800)

anim_save("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/anim_diversity.gif")
animate(p2, fps=25)

p2 +
	transition_reveal(dates)


# speed up or slow down frames per second
animate(p2, fps=20)

p <- ggplot(pred.plot1, aes(x = dates)) +
	geom_point(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	# geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = siteID), alpha=0.2, show.legend = F) +
	# geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = siteID), alpha=0.4, show.legend = F) +
	geom_point(aes(y = truth)) + theme_ipsum(base_size = 14, strip_text_size = 22) + 
	ggtitle(unique(pred.plot1$plotID)) + scale_x_date()



p + transition_time(run) 



# Just one plot forecast
ggplot(pred.plot1, aes(x = dates)) +
	geom_line(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = siteID), alpha=0.2, show.legend = F) +
	geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = siteID), alpha=0.4, show.legend = F) +
	geom_point(aes(y = truth)) + 
	#theme_ipsum() + 
	coord_cartesian(ylim = c(4, 6), # This focuses the x-axis on the range of interest
									clip = 'off') +  
	theme_ipsum(base_size = 14, strip_text_size = 16) + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) + 
	ggtitle(paste0("Shannon diversity at: ", unique(pred.plot1$plotID))) + 
	scale_x_date() 

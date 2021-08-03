# Predictive loss for diversity
# Testing CRPS in scoringRules

group <- "16S"

pacman::p_load(scoringRules, tidyr, dplyr, reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/")

group = "16S"


newsites_est <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data_newsites.rds")
oldsites_est <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_forecast_data.rds")

newsites_est <- newsites_est %>% filter(!is.na(truth)) %>% mutate(truth = as.numeric(truth))
oldsites_est <- oldsites_est %>% filter(!is.na(truth)) %>% mutate(truth = as.numeric(truth))

scored_newsites <- newsites_est %>% group_by(group) %>% 
	mutate(crps = crps_norm(truth, mean_EDP, sd_EDP))# %>% summarize(mean_score = mean(crps, na.rm=T))
scored_newsites$newsites <- "New sites"
scored_oldsites <- oldsites_est %>% group_by(group) %>% filter(fcast_period=="Hindcast") %>% 
	mutate(crps = crps_norm(truth, mean_EDP, sd_EDP)) 
# %>% summarize(mean_score = mean(crps, 
scored_oldsites$newsites <- "Observed sites"
scores_out <- plyr::rbind.fill(scored_newsites, scored_oldsites)
## Compare overall scores for out of sample bacteria vs fungi
p <- ggplot(scores_out, aes(x = pretty_group, y = crps)) + 
	geom_violin(aes(fill = pretty_group), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~newsites) +   coord_trans(y = "log10") +
	geom_jitter(aes(x = pretty_group, y = crps), width=.2, height = 0, alpha = .3) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)
p

# Do T-tests
newsites <- scores_out[scores_out$newsites=="New sites",]
newsites_bac <- newsites[newsites$group=="16S",]
newsites_fun <- newsites[newsites$group=="ITS",]
oldsites <- scores_out[scores_out$newsites=="Observed sites",]
oldsites_bac <- oldsites[oldsites$group=="16S",]
oldsites_fun <- oldsites[oldsites$group=="ITS",]
t.test(oldsites_bac$crps, oldsites_fun$crps)
t.test(newsites_bac$crps, newsites_fun$crps)


library(ggpubr)
p + stat_compare_means(method = "t.test", size = 5)


data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")

scores.list <- data_in$scores
to_score <- scores.list %>% #filter(!is.na(scores.list$truth)) %>% 
	select(plot_num, date_num, site_num, siteID, plotID, dateID, truth, Mean, SD, scenario, group) %>% mutate(truth = as.numeric(truth))

scored <- to_score %>% group_by(group, scenario) %>% mutate(crps = crps_norm(truth, Mean, SD)) %>% summarize(mean_score = mean(crps, na.rm=T))
scored <- to_score %>% group_by(group, scenario) %>% mutate(crps = crps_norm(truth, Mean, SD))

to_score <- to_score[to_score$scenario=="no uncertainty" & to_score$group == group,]
#scored <- scoringRules::crps_norm(to_score$truth, mean = to_score$Mean, sd = to_score$SD)
to_score$crps <- scoringRules::crps_norm(to_score$truth, mean = to_score$Mean, sd = to_score$SD)

to_score[to_score$Mean==0,]$Mean <- NA
to_score[to_score$SD==0,]$SD <- NA
ggplot(scored_newsites[scored_newsites$plotID=="CLBJ_003",]) +
	#geom_line(aes(x = date_num, y = Mean), show.legend = F) + 
	geom_ribbon(aes(x = date_num, ymin = mean_EDP-sd_EDP, ymax = mean_EDP+sd_EDP), 
							fill = "darkblue", alpha=0.4, show.legend = F, na.rm=T) +
	geom_point(aes(x = date_num, y = as.numeric(truth))) +
	geom_point(aes(x = date_num, y = 1, size = as.numeric(crps)), color= "red")



if (group == "ITS"){
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
} else if (group == "16S") {
	div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
}

rank.df = div_in$cal

# Custom function for organizing model data.
model.cal <- prepDivData(rank.df = div_in$cal, min.prev = 3)
model.val <- prepDivData(rank.df = rbind( div_in$cal,
																					div_in$val), min.prev = 3)



data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")

plot_est <- data_in$plot_est


allrank_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/allrank_summaries.rds")
rank.name <- "phylum_bac"

cal.pl.out <- list()
pl.out <- list()

plot_est <- allrank_summaries$plot_est
rank.name <- "bac_diversity"
for (rank.name in c("full_uncertainty_16S","full_uncertainty_ITS")){
	
	
	
	
	## CALIBRATION
	plot_est_rank <- plot_est[plot_est$rank==rank.name,]
	obs.df_rank <- plot_est_rank[!is.na(plot_est_rank$truth),]
	
		obs <- obs.df$truth
		
		# Calculate predictive loss
		npred.real <- nrow(obs.df)
		## Residual variance
		ybar <- obs.df$`50%`
		G <- sum((obs-ybar)^2, na.rm=T)/npred.real
		## Predictive variance 
		## TODO: change from range -> actual variance (need full samples)
		P <- sum(range(obs.df[,c("2.5%","97.5%")], na.rm=T))/npred.real
		#P <- sum(apply(ypred,2,var))/npred.real
		Dpl <- G + P
		PL <- c(G,P,Dpl)
		names(PL) <- c("Residual variance","Predictive variance","Total variance")
		
		cal.pl.out[[rank.name]] <- as.data.frame(t(PL))
		
}
cal.pl.df <- rbind.named.dfs(cal.pl.out)


cal.pl.df <- cal.pl.df %>% tidyr::separate(name, sep = "_", into = c("rank", "group"))
cal.pl.df$rank <- ordered(cal.pl.df$rank, levels = c("genus",
																										 "family",
																										 "order", 
																										 "class",
																										 "phylum"))
ggplot(cal.pl.df) + geom_point(aes(x = rank, y = `Total variance`)) + facet_grid(~group)



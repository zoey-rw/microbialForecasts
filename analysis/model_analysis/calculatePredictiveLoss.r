# Predictive loss for diversity

group <- "16S"

pacman::p_load(tidyr, dplyr, reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDiversityData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/")


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
rank.name <- "phylum_bac"
for (rank.name in c("phylum_bac","class_bac","family_bac","order_bac","genus_bac",
										"phylum_fun",
										"class_fun","family_fun","order_fun","genus_fun"
)){
	
	## CALIBRATION
	plot_est_rank <- plot_est[plot_est$rank==rank.name,]
	obs.df_rank <- plot_est_rank[!is.na(plot_est_rank$truth),]
	
	for(s in 1:11){
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
}
cal.pl.df <- rbind.named.dfs(cal.pl.out)


cal.pl.df <- cal.pl.df %>% tidyr::separate(name, sep = "_", into = c("rank", "group"))
cal.pl.df$rank <- ordered(cal.pl.df$rank, levels = c("genus",
																								 "family",
																								 "order", 
																								 "class",
																								 "phylum"))
ggplot(cal.pl.df) + geom_point(aes(x = rank, y = `Total variance`)) + facet_grid(~group)










	div_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
	div_16S <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
	
	div_ITS <- rbind(div_ITS$cal, div_ITS$val)
	div_16S <- rbind(div_16S$cal, div_16S$val)
	
	div_ITS <- div_ITS %>% rename(Shannon_ITS = Shannon)
	div_16S <- div_16S %>% rename(Shannon_16S = Shannon)

	div_both <- merge(div_16S, div_ITS, by = c("sampleID","siteID","plotID"))
	
	summary(lm(Shannon_ITS ~ Shannon_16S, div_both))
	
	ggplot(div_both) + geom_point(aes(x = Shannon_16S, y = Shannon_ITS)) + facet_grid(~siteID) 
	

# Calculate and visualize CRPS for functional group forecasts
pacman::p_load(scoringRules, tidyr, dplyr, reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

hindcast_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_tax.rds")
hindcast_data$dates <- fixDate(hindcast_data$dateID)
hindcast_data <- hindcast_data %>% filter(dates > "2016-12-31"  & !is.na(hindcast_data$truth))


hindcast_data$pretty_name <- recode(hindcast_data$rank,
												 "genus_bac" = "Genus",
												 "family_bac" = "Family",
												 "order_bac" = "Order", 
												 "class_bac" = "Class", 
												 "phylum_bac" = "Phylum",
												 "genus_fun" = "Genus",
												 "family_fun" = "Family",
												 "order_fun" = "Order", 
												 "class_fun" = "Class", 
												 "phylum_fun" = "Phylum",
												 "functional_group" = "Functional group")

hindcast_data$group <- sapply(str_split(hindcast_data$rank, "_",  n = 2), `[`, 2)
hindcast_data$group <- ifelse(hindcast_data$group == "bac", "16S", "ITS")
hindcast_data$pretty_group <- ifelse(hindcast_data$group == "16S", "Bacteria", "Fungi")
hindcast_data$only_rank <- sapply(str_split(hindcast_data$rank, "_",  n = 2), `[`, 1)
hindcast_data$only_rank <- ordered(hindcast_data$only_rank, levels = c("genus",
																								 "family",
																								 "order", 
																								 "class",
																								 "phylum", "functional"))
hindcast_data$rank <- ordered(hindcast_data$rank, levels = c("genus_bac","genus_fun",
																			 "family_bac","family_fun",
																			 "order_bac", "order_fun",
																			 "class_bac", "class_fun",
																			 "phylum_bac","phylum_fun",
																			 "functional_group"))


scored_hindcasts <- hindcast_data %>% mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd),
																																									 fcast_type = "Taxonomic")
saveRDS(scored_hindcasts, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_tax.rds")




# COMPARE CRPS BY ERRORS-IN-VARIABLES UNCERTAINTY
# All values
ggplot(scored_hindcasts, aes(x = pretty_name, y = crps)) + 
	geom_violin(aes(fill = pretty_name), draw_quantiles = c(0.5), show.legend=F) + 
	facet_grid(rows=vars(group), drop=T, scales = "free") +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = pretty_name, y = crps), width=.2, height = 0, alpha = .3) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

# Mean values
scored_hindcasts_summary <- scored_hindcasts %>% group_by(group, pretty_name) %>% summarize(mean_score = mean(crps, na.rm=T))
ggplot(scored_hindcasts_summary, aes(x = pretty_name, y = mean_score)) + 
	geom_violin(aes(fill = pretty_name), draw_quantiles = c(0.5), show.legend=F) + 
	facet_grid(rows=vars(group), drop=T, scales = "free") +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = pretty_name, y = mean_score), width=.2, height = 0, alpha = .8, show.legend = F) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)







# Do T-tests
newsites <- scores_out[scores_out$newsites=="New sites",]
newsites_bac <- newsites[newsites$group=="16S",]
newsites_fun <- newsites[newsites$group=="ITS",]
oldsites <- scores_out[scores_out$newsites=="Observed sites",]
oldsites_bac <- oldsites[oldsites$group=="16S",]
oldsites_fun <- oldsites[oldsites$group=="ITS",]
t.test(oldsites_bac$crps, oldsites_fun$crps)
t.test(newsites_bac$crps, newsites_fun$crps)


tukey <- function(df, x, y, extra_info = NULL){
	new.df <- cbind.data.frame(x = df[,x], y = df[,y])
	abs_max <- max(new.df[,"y"])
	maxs <- new.df %>%
		group_by(x) %>%
		summarise(tot=max(y)+ 0.3 * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>% 
		dplyr::select(-y) %>%
		left_join(maxs, by="x") 
	if (!is.null(extra_info)){
		Tukey_test <- cbind.data.frame(Tukey_test, extra_info)
	}
	return(Tukey_test)
}
tukey(df = scored_hindcasts, x = "rank", y = "crps")


x = "rank"
y = "crps"
new.df <- cbind.data.frame(x = df$rank, y = df$crps)
abs_max <- max(new.df[,"y"], na.rm = T)
maxs <- new.df %>%
	group_by(x) %>%
	summarise(tot=max(y, na.rm=T)+ 0.6 * abs_max)
Tukey_test <- aov(y ~ x, data=new.df) %>%
	agricolae::HSD.test("x", group=TRUE) %>%
	.$groups %>%
	as_tibble(rownames="x") %>%
	rename("Letters_Tukey"="groups") %>% 
	dplyr::select(-y) %>%
	left_join(maxs, by="x") %>% 
	rename("rank"="x") 

Tukey_test$group <- sapply(str_split(Tukey_test$rank, "_",  n = 2), `[`, 2)
Tukey_test$group <- ifelse(Tukey_test$group == "bac", "Bacteria", "Fungi")
Tukey_test$pretty_name <- recode(Tukey_test$rank,
																		"genus_bac" = "Genus",
																		"family_bac" = "Family",
																		"order_bac" = "Order", 
																		"class_bac" = "Class", 
																		"phylum_bac" = "Phylum",
																		"genus_fun" = "Genus",
																		"family_fun" = "Family",
																		"order_fun" = "Order", 
																		"class_fun" = "Class", 
																		"phylum_fun" = "Phylum",
																		"functional_group" = "Functional group")
scored_hindcasts <- left_join(scored_hindcasts, Tukey_test)
ggplot(scored_hindcasts, aes(x = pretty_name, y = crps)) + 
	geom_violin(aes(fill = pretty_name), color = "black", draw_quantiles = c(.25, 0.5, .75), show.legend=F) + 
	facet_grid(rows=vars(group), drop=T, scales = "free") +  
	geom_text(data = Tukey_test, aes(x = pretty_name, y = tot, label = Letters_Tukey)) + 
	coord_trans(y = "log2") +
	geom_jitter(aes(x = pretty_name, y = crps#, color = pretty_name
									), width=.2, height = 0, alpha = .2) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)


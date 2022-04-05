pacman::p_load(scoringRules, tidyr, dplyr, reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 


crps_div <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_div.rds")
crps_div$taxon_name <- paste(crps_div$fcast_type, "_", crps_div$pretty_group)
crps_fg <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_fg.rds")
crps_fg$only_rank <- "fg"
crps_fg$pretty_name <- "Functional group"
crps_fg$pretty_group <- "Bacteria"
crps_fg[grep("sapr|path|arbusc|ecto|endo|lichen", crps_fg$taxon_name),]$pretty_group <- "Fungi"
crps_fg$rank <- paste0(crps_fg$only_rank,"_", crps_fg$pretty_group)

crps_fg <- crps_fg %>% filter(!taxon_name %in% c("litter_saprotroph","soil_saprotroph","wood_saprotroph"))

crps_tax <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_tax.rds")
crps_tax <- crps_tax %>% filter(species != "other")
crps_all <- data.table::rbindlist(list(crps_div, crps_fg, crps_tax), fill=T)
crps_all$pretty_name <- ordered(crps_all$pretty_name, levels = c("Genus",
																																			 "Family",
																																			 "Order", 
																																			 "Class",
																																			 "Phylum", "Functional group"))


# Compare taxonomic ranks and functional groups at observed sites, F and B
crps_fg_tax <- crps_all[crps_all$fcast_type != "Diversity",]
crps_fg_tax <- crps_fg_tax[!crps_fg_tax$new_site,]

df <- crps_fg_tax
x = "rank"
y = "crps"
new.df <- cbind.data.frame(x = df$rank, y = df$crps)
abs_max <- max(new.df[,"y"], na.rm = T)
maxs <- new.df %>%
	group_by(x) %>%
	summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
Tukey_test <- aov(y ~ x, data=new.df) %>%
	agricolae::HSD.test("x", group=TRUE) %>%
	.$groups %>%
	as_tibble(rownames="x") %>%
	rename("Letters_Tukey"="groups") %>% 
	dplyr::select(-y) %>%
	left_join(maxs, by="x") %>% 
	rename("rank"="x") 

Tukey_test$group <- sapply(str_split(Tukey_test$rank, "_",  n = 2), `[`, 2)
Tukey_test$pretty_group <- ifelse(Tukey_test$group %in% c("bac","Bacteria"), "Bacteria", "Fungi")
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
																 "fg_Fungi" = "Functional group",
																 "fg_Bacteria" = "Functional group")
Tukey_test$pretty_name <- ordered(Tukey_test$pretty_name, levels = c("Genus",
																																				 "Family",
																																				 "Order", 
																																				 "Class",
																																				 "Phylum", "Functional group"))
brk_y <- c(0.0001,0.001, 0.01, 0.05, 0.1, 0.5, 1, 1.5)
lab_y <- c( "< 0.0001", "0.001","0.01", "0.05", "0.1", "0.5", "1", "1.5")
#### FIGURE: within-site CRPS of F vs B taxonomic ranks and functional groups 
ggplot(crps_fg_tax, aes(x = pretty_name, y = crps, fill = pretty_group)) + 
	geom_violin(color = "black", draw_quantiles = c(.25, 0.5, .75), show.legend=F, scale = "width") + 
	geom_text(data = Tukey_test, aes(x = pretty_name, y = tot, label = Letters_Tukey), show.legend = F) + 
	facet_grid(rows=vars(pretty_group), drop=T, scales = "free") +  
	coord_trans(y = "log10") +  scale_y_continuous(breaks = brk_y, labels = lab_y) +
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + theme(axis.text.x = element_text(hjust=0,vjust = 1,angle=-40))

crps_fg <- crps_fg_tax[crps_fg_tax$fcast_type=="Functional group",]




library(tidyverse)

scored_hindcasts_summary <- crps_all %>% #filter(fcast_type != "Diversity") %>% 
	group_by(pretty_group, fcast_type, pretty_name, taxon_name, new_site) %>% 
	summarize(mean_score = mean(crps, na.rm=T),
						median_score = median(crps, na.rm=T)) %>% ungroup #%>% 
#	group_by(pretty_group, fcast_type) %>% summarize(decrease_in_accuracy = )
summary_mean <- scored_hindcasts_summary %>% select(-c(median_score)) %>% pivot_wider(names_from = new_site, values_from = mean_score, names_prefix = "new_site_")
#summary_median <- scored_hindcasts_summary %>% select(-c(mean_score)) %>% pivot_wider(names_from = new_site, values_from = median_score, names_prefix = "new_site_")
summary_mean <- summary_mean %>% mutate(
	diff_in_accuracy = new_site_TRUE - new_site_FALSE,
	pct_decrease = diff_in_accuracy/new_site_FALSE,
	decrease_in_accuracy = 1 - new_site_TRUE/new_site_FALSE
)

#summary_median$decrease_in_accuracy <- summary_median$new_site_TRUE/summary_median$new_site_FALSE

# Compare decrease in accuracy for fungi vs bacteria
ggplot(summary_mean[summary_mean$fcast_type != "Diversity",], 
			 aes(x = pretty_name, y = decrease_in_accuracy, color = pretty_group)) + 
	# facet_grid(rows=vars(fcast_type), #cols = vars(pretty_group),
	# 					 drop=T, scales = "free_x") +
	#coord_trans(y = "log10") +
	geom_jitter(aes(fill = pretty_group), shape=21, color="black", width=.1, height = 0, size=4) + 
	ylab("Relative accuracy at new sites") + xlab(NULL) + 
	theme_bw(base_size=18)+ scale_fill_manual(values = c("grey30","grey90")) + geom_hline(yintercept = 0, linetype=2) 


ggplot(summary_mean, aes(x = pretty_group, y = decrease_in_accuracy)) + 
	facet_grid(rows=vars(fcast_type), #cols = vars(pretty_group),
						 drop=T, scales = "free_x") +
	#coord_trans(y = "log10") +
	geom_jitter(aes(fill = pretty_group), shape=21, color="black", width=.1, height = 0, size=4) + 
		ylab("Relative accuracy at new sites") + xlab(NULL) + 
	theme_bw(base_size=18) + scale_fill_manual(values = c("grey30","grey90")) + geom_hline(yintercept = 0, linetype=2) 














df <- summary_mean# %>% filter(pretty_name != "Functional group")
new.df <- cbind.data.frame(x = df$pretty_group, y = df$decrease_in_accuracy)
abs_max <- max(new.df[,"y"], na.rm = T)
maxs <- new.df %>%
	group_by(x) %>%
	summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
Tukey_test <- aov(y ~ x, data=new.df) %>%
	agricolae::HSD.test("x", group=TRUE) %>%
	.$groups %>%
	as_tibble(rownames="x") %>%
	rename("Letters_Tukey"="groups") %>% 
	dplyr::select(-y) %>%
	left_join(maxs, by="x") %>% 
	rename("pretty_group"="x") 


# FIGURE: Compare decrease in accuracy across taxonomic ranks
ggplot(summary_mean[summary_mean$fcast_type != "Diversity",], aes(x = pretty_name, y = decrease_in_accuracy, color = pretty_group)) + 
	# geom_violin(color = "black", draw_quantiles = c(.25, 0.5, .75), show.legend=F, scale = "width") + 
	facet_grid(rows=vars(pretty_group), #cols = vars(pretty_group),
						 drop=T, scales = "free") +
	#coord_trans(y = "log10") +
	geom_jitter( width=.1, height = 0, alpha = .7, size = 4, show.legend = F) + 
	ylab("Relative accuracy at new sites") + xlab(NULL) + 
	theme_bw(base_size=20) + theme(axis.text.x = element_text(hjust=0,vjust = 1,angle=-40))


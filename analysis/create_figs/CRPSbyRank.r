pacman::p_load(scoringRules,reshape2, parallel, lubridate) 
source("source.R")

crps_in <- readRDS(here("data/summary/CRPS_hindcasts.rds"))

#crps_all <- crps_all %>% filter(!taxon_name %in% c("litter_saprotroph","soil_saprotroph","wood_saprotroph"))

crps_all <- crps_in$scored_hindcasts # CRPS for each observation
# Compare taxonomic ranks and functional groups at observed sites, F and B
crps_fg_tax <- crps_all[crps_all$fcast_type != "Diversity",]
crps_fg <- crps_fg_tax[crps_fg_tax$fcast_type=="Functional group",]


crps_fg_tax <- crps_fg_tax[!crps_fg_tax$new_site,]


scored_hindcasts_plot <- crps_in$scored_hindcasts_plot

skill_score <- crps_in$skill_score_rank
skill_score_taxon <- crps_in$skill_score_taxon

##### Compare CRPS for all in-site predictions by rank (fg/tax) -----
Tukey_test <- tukey(crps_fg_tax$rank, crps_fg_tax$crps, y.offset = .1) 
colnames(Tukey_test)[[1]] <- c("rank")
Tukey_test$group <- sapply(str_split(Tukey_test$rank, "_",  n = 2), `[`, 2)
Tukey_test$pretty_group <- ifelse(Tukey_test$group %in% c("bac","Bacteria"), "Bacteria", "Fungi")
Tukey_test$pretty_name <- recode(Tukey_test$rank, !!!pretty_rank_names)
Tukey_test$pretty_name <- ordered(Tukey_test$pretty_name, levels = c("Genus",
																																		 "Family",
																																		 "Order", 
																																		 "Class",
																																		 "Phylum", "Functional group"))
brk_y <- c(0.0001,0.001, 0.01, 0.05, 0.1, 0.5, 1)
lab_y <- c( "< 0.0001", "0.001","0.01", "0.05", "0.1", "0.5", "1")
#### FIGURE: within-site CRPS of taxonomic ranks and functional groups, sep by F/B
in_site_rank <- ggplot(crps_in$scored_hindcasts_plot %>% filter(fcast_type != "Diversity")) + 
	geom_point(aes(x = pretty_name, y = crps_mean, color = pretty_group), 
						 alpha = .1, size=4,
						 position=position_jitterdodge(dodge.width = 1)) + 
	geom_violin(aes(x = pretty_name, y = crps_mean, color = pretty_group), 
							draw_quantiles = c(0.5), show.legend=F) + 
	coord_trans(y = "log10") +
	scale_y_continuous(breaks = brk_y, labels = lab_y) +
	facet_grid(rows=vars(pretty_group), drop=T, scales = "free") +  
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=20) + ggtitle("Predictability at observed sites") + 
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1)) + 
	guides(color=guide_legend(title=NULL))
in_site_rank <- in_site_rank + 	
	geom_text(data = Tukey_test, aes(x = pretty_name, y = tot, label = Letters_Tukey), col = 2, size = 6, 
						show.legend = F) 

#####
in_site_rank


##### Skill score by each taxon -----
pseudoLog <- scales::pseudo_log_trans(base = 10)

skill_taxon <- skill_score_taxon %>% filter(fcast_type != "Diversity")
skill_score_taxon_plot <- ggplot(skill_taxon) + 
	geom_point(aes(x = pretty_name, y = skill_score, color = pretty_group), 
						 alpha = .3, size=4,
						 position=position_jitterdodge(dodge.width = 1)) + 
	geom_violin(aes(x = pretty_name, y = skill_score, color = pretty_group), 
							draw_quantiles = c(0.5), show.legend=F) + 
	facet_grid(rows=vars(pretty_group), drop=T, scales = "free") +  
	ylab("Skill score (% change in CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") + 
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1)) + 
	guides(color=guide_legend(title=NULL)) + 
	coord_trans(y=pseudoLog)
Tukey_test <- tukey(skill_taxon$pretty_name, skill_taxon$skill_score, y.offset = .1) 
colnames(Tukey_test)[[1]] <- c("pretty_name")
skill_score_taxon_plot <- skill_score_taxon_plot + 	
	geom_text(data = Tukey_test, aes(x = pretty_name, y = tot, label = Letters_Tukey), col = 2, size = 6, 
						show.legend = F) 
skill_score_taxon_plot
##### 

##### Compare skill score for ranks ------

skill_score_rank_allcov <- skill_score %>% filter(model_name == "all_covariates")
skill_score_rank <- ggplot(skill_score_rank_allcov, 
			 aes(x = pretty_name, y = skill_score, 
			 		color = pretty_group, shape = fcast_type)) + 
	#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_grid(~model_name, drop = T) +  
	geom_jitter(aes(x = pretty_name, y = skill_score), width=.1, 
							height = 0, alpha = .8, size=4) + 
	#coord_trans(y = "sqrt") +  
	ylab("Skill score (% change in CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") + 
	theme(text = element_text(size = 20),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=22)) + 
	guides(shape=guide_legend(title=NULL),
				 color=guide_legend(title=NULL)) +
	coord_trans(y=pseudoLog)
Tukey_test <- tukey(skill_taxon$pretty_name, skill_taxon$skill_score, y.offset = .1) 
colnames(Tukey_test)[[1]] <- c("pretty_name")
skill_score_rank <- skill_score_rank + 	
	geom_text(data = Tukey_test, aes(x = pretty_name, y = tot, label = Letters_Tukey), col = 2, size = 6, 
						show.legend = F) 
skill_score_rank
#####



##### Compare skill score for F vs B ------
skill_score_f_vs_b <- ggplot(skill_score %>% filter(model_name == "all_covariates"), 
														 aes(x = pretty_group, y = skill_score, 
														 		color = pretty_name, shape = fcast_type)) + 
	#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_grid(~model_name, drop = T) +  
	geom_jitter(aes(x = pretty_group, y = skill_score), width=.1, 
							height = 0, alpha = .8, size=4) + 
	ylab("Skill score (% change in CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") + 
	theme(text = element_text(size = 20),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=22)) + 
	guides(shape=guide_legend(title=NULL),
				 color=guide_legend(title=NULL)) +
	coord_trans(y=pseudoLog)
skill_score_f_vs_b
#####

crps_in$scored_hindcasts_plot
in_site_f_vs_b <- ggplot(crps_in$scored_hindcasts_taxon %>% filter(model_name == "all_covariates",
																																	 newsite == "Observed site"), 
														 aes(x = pretty_group, y = crps_mean, 
														 		color = pretty_name, shape = fcast_type)) + 
	#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_grid(~model_name, drop = T) +  
	geom_jitter(aes(x = pretty_group, y = crps_mean), width=.1, 
							height = 0, alpha = .8, size=4) + 
	ylab("CRPS") + xlab(NULL) + 
	theme_minimal(base_size=18) + ggtitle("Predictability at observed sites") + 
	theme(text = element_text(size = 20),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=22)) + 
	guides(shape=guide_legend(title=NULL),
				 color=guide_legend(title=NULL))
in_site_f_vs_b



scored_hindcasts_summary <- crps_fg_tax %>% #filter(fcast_type != "Diversity") %>% 
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


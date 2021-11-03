# Calculate and visualize CRPS for functional group forecasts
pacman::p_load(scoringRules, tidyr, dplyr, reshape2, parallel, lubridate, nimble, coda, tidyverse, runjags) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

hindcast_data <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_fg.rds") %>% 
	mutate(fg_cat = assign_fg_categories(taxon_name)) %>% filter(taxon_name %in% keep_fg_names)
hindcast_data$scenario <- factor(hindcast_data$scenario, levels = c("no_uncertainty","temporal_uncertainty","spatial_uncertainty","full_uncertainty"))
hindcast_data_filt <- hindcast_data %>% filter(dates > "2016-12-31" & !is.na(hindcast_data$truth))

scored_hindcasts <- hindcast_data_filt %>% mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd),
																																									 fcast_type = "Functional group")

saveRDS(scored_hindcasts, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_fg.rds")





scored_hindcasts <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/CRPS_fg.rds")

# Compare the functional groups
scored_hindcasts_summary <- scored_hindcasts %>% 
	mutate(fg_cat = assign_fg_categories(taxon_name)) %>% 
	group_by(fg_cat, taxon_name, new_site) %>% 
	summarize(mean_score = mean(crps, na.rm=T))

# Plot every group
ggplot(scored_hindcasts_summary,
			 aes(x = new_site, y = mean_score)) + 
	#facet_grid(~fg_cat, space="free", scales = "free_x") +  
	coord_trans(y = "log10") +
	geom_jitter(aes(color = fg_cat), alpha = .7, size = 3, width = .03) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=16) + theme(axis.text.x = element_text()) + 
	scale_x_discrete(labels=c("FALSE" = "Within-site", "TRUE" = "Out-of-site")) +
	scale_color_discrete(name = "Relative accuracy at new sites") +
	geom_text(aes(label=ifelse(mean_score>.06,as.character(taxon_name),'')),vjust = .1, hjust=-0.2)






summary_mean <- scored_hindcasts_summary %>% pivot_wider(names_from = new_site, values_from = mean_score, names_prefix = "new_site_")
#summary_median <- scored_hindcasts_summary %>% select(-c(mean_score)) %>% pivot_wider(names_from = new_site, values_from = median_score, names_prefix = "new_site_")
summary_mean <- summary_mean %>% mutate(
	diff_in_accuracy = new_site_TRUE - new_site_FALSE,
	pct_decrease = diff_in_accuracy/new_site_FALSE,
	decrease_in_accuracy = 1 - new_site_TRUE/new_site_FALSE
)

# Plot every group
ggplot(summary_mean,
			 aes(x = fg_cat, y = decrease_in_accuracy)) + 
	#facet_grid(~fg_cat, space="free", scales = "free_x") +  
	#coord_trans(y = "log10") +
	geom_jitter(aes(color = fg_cat), alpha = .7, size = 3, width = .03, show.legend = F) + 
	ylab("Relative accuracy at new sites") + xlab(NULL) + 
	theme_minimal(base_size=16) + theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) + 
	ggrepel::geom_text_repel(aes(label=ifelse(decrease_in_accuracy>-.4,as.character(taxon_name),'')),vjust = .1, hjust=-0.2) + geom_hline(yintercept = 0, linetype=2) 






# COMPARE CRPS BY ERRORS-IN-VARIABLES UNCERTAINTY
# All values
ggplot(scored_hindcasts, aes(x = scenario, y = crps)) + 
	geom_violin(aes(fill = scenario), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_wrap(~scenario) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = scenario, y = crps), width=.2, height = 0, alpha = .3) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

# Mean values
scored_hindcasts_summary <- scored_hindcasts %>% group_by(name, scenario) %>% summarize(mean_score = mean(crps, na.rm=T))
ggplot(scored_hindcasts_summary, aes(x = scenario, y = mean_score)) + 
	geom_violin(aes(fill = scenario), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_wrap(~scenario) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = scenario, y = mean_score), width=.2, height = 0, alpha = .8, show.legend = F) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

# COMPARE FORECASTS AT NEW VS OBSERVED SITES

scored_hindcasts_summary <- scored_hindcasts %>% group_by(name, scenario, new_site) %>% summarize(mean_score = mean(crps, na.rm=T))
ggplot(scored_hindcasts_summary, aes(x = scenario, y = mean_score)) + 
	geom_violin(aes(fill = scenario), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~new_site) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = scenario, y = mean_score), width=.2, height = 0, alpha = .8, show.legend = F) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)



ggplot(scored_hindcasts_summary, aes(x = name, y = mean_score)) + 
	geom_violin(aes(fill = name), draw_quantiles = c(0.5), show.legend=F) + 
	facet_grid(rows = vars(scenario), cols = vars(new_site)) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = name, y = mean_score), width=.2, height = 0, alpha = .8, show.legend = F) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)



# Not tested (from old code
)
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






scored_hindcasts$fg_cat <- NA
scored_hindcasts[grep("simple", scored_hindcasts$name),]$fg_cat <- "Simple substrates"
scored_hindcasts[grep("complex", scored_hindcasts$name),]$fg_cat <- "Complex substrates"
scored_hindcasts[grep("stress", scored_hindcasts$name),]$fg_cat <- "Stresses"
scored_hindcasts[grep("antibiotic", scored_hindcasts$name),]$fg_cat <- "Antibiotic resistance"
scored_hindcasts[grep("anaerobic", scored_hindcasts$name),]$fg_cat <- "Anaerobic"
scored_hindcasts[grep("nitr|fixa", scored_hindcasts$name),]$fg_cat <- "N-cycling"
scored_hindcasts[grep("sapr|path|arbusc|ecto|endo|lichen", scored_hindcasts$name),]$fg_cat <- "Trophic guild"
scored_hindcasts[grep("copio|oligo", scored_hindcasts$name),]$fg_cat <- "Life-history"
scored_hindcasts[is.na(scored_hindcasts$fg_cat),]$fg_cat <- "Other"

scored_hindcasts$group <- ifelse(grepl("Troph", scored_hindcasts$fg_cat), "ITS", "16S")

full_uncert <- scored_hindcasts %>% filter(scenario == "full_uncertainty")
scored_hindcasts_summary <- full_uncert %>% group_by(group, name, new_site) %>% summarize(mean_score = mean(crps, na.rm=T))
ggplot(scored_hindcasts_summary, aes(x = group, y = mean_score)) + 
	geom_violin(aes(fill = group), draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~new_site) +  
	coord_trans(y = "log10") +
	geom_jitter(aes(x = group, y = mean_score), width=.2, height = 0, alpha = .8, show.legend = F) + ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

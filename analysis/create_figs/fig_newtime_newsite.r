
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
#p_load(forestmangr)
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
#scores_list = readRDS(here("data/summary/scoring_metrics_cv_jan18.rds"))

library(ggallin)

#keep = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
unconverged = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
converged = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
#keep = keep[keep$median_gbr < 3,]
unconverged = unconverged[unconverged$median_gbr > 7,]
#unconverged = unconverged[unconverged$median_gbr > 3,]
#converged = converged[converged$median_gbr < 4,]
all_cov = converged[converged$model_name=="all_covariates",]



# Get mean values for plotting
site_time_scores_to_plot = scores_list$scoring_metrics %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	#	filter(!taxon_model_rank %in% unconverged$taxon_model_rank)%>%
	#filter(taxon_model_rank %in% converged$taxon_model_rank)%>%
	filter(model_name=="all_covariates" &
				 	fcast_type !="Diversity" &
				 	!grepl("random", site_prediction))


rank_vals = site_time_scores_to_plot %>% filter(grepl("observed", site_prediction))
rank_means = rank_vals %>% ungroup() %>%
	group_by(fcast_type, pretty_group, model_name, pretty_name, rank) %>%
	summarize(mean_RSQ=mean(RSQ, na.rm=T),
						sd_RSQ = sd(RSQ, na.rm=T)) %>%
	mutate(ymax = mean_RSQ + 1.96*sd_RSQ,
				 ymin = mean_RSQ - 1.96*sd_RSQ) %>% arrange(pretty_group, pretty_name)



# Get mean values for plotting
skill_score_vals = scores_list$skill_score_taxon %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity")
skill_score_to_plot = skill_score_vals %>% ungroup() %>%
	group_by(fcast_type, pretty_group, model_name, pretty_name, rank) %>%
	summarize(mean_skill_score=mean(skill_score, na.rm=T),
						sd_skill_score = sd(skill_score, na.rm=T)) %>%
	mutate(ymax = mean_skill_score + 1.96*sd_skill_score,
				 ymin = mean_skill_score - 1.96*sd_skill_score) %>% arrange(pretty_group, pretty_name)



tukey2 = function (x, y, extra_info = NULL, y.offset = 0.3)
{
	new.df <- cbind.data.frame(x = x, y = y)
	abs_max <- max(new.df[, 2], na.rm=T)
	maxs <- new.df %>% group_by(x) %>% summarise(tot = max(y, na.rm=T) +
																							 	y.offset * abs_max)
	Tukey_test <- aov(y ~ x, data = new.df) %>% agricolae::HSD.test("x",
																																	group = TRUE) %>% .$groups %>% as_tibble(rownames = "x") %>%
		rename(Letters_Tukey = "groups") %>% dplyr::select(-y) %>%
		left_join(maxs, by = "x")
	if (!is.null(extra_info)) {
		Tukey_test <- cbind.data.frame(Tukey_test)
	}
	return(Tukey_test)
}

# Calculate Tukey groups
tukey_list <- list()
for (pretty_group in c("Bacteria","Fungi")){
	group_df <- rank_vals %>% filter(pretty_group == !!pretty_group)
	out <- tukey2(group_df$pretty_name, group_df$RSQ	, y.offset = .1)
	tukey_list[[pretty_group]] <- out %>% mutate(pretty_group = !!pretty_group)
}
tukey_list_newtime <- data.table::rbindlist(tukey_list)

tukey_list <- list()
for (pretty_group in c("Bacteria","Fungi")){
	group_df <- skill_score_vals %>% filter(pretty_group == !!pretty_group)
	out <- tukey(group_df$pretty_name, group_df$skill_score	, y.offset = .1)
	tukey_list[[pretty_group]] <- out %>% mutate(pretty_group = !!pretty_group)
}
tukey_list_newsite <- data.table::rbindlist(tukey_list)


# Compare new-time only across taxonomic ranks and functional groups
newtime_by_rank = ggplot(rank_means,
												 aes(x = pretty_name, y = mean_RSQ, color = pretty_group)) +
	facet_grid(cols=vars(pretty_group), scales="free") +
	theme_minimal(base_size = 18) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ at new times") + xlab("Rank") +
	#scale_color_discrete(name=NULL, labels=c("Within-site","New site")) +
	#scale_y_continuous(trans=scales::pseudo_log_trans(base = 2)) +
	geom_pointrange(aes(ymin = ymin, ymax = ymax, color = pretty_group), show.legend = F) + ggtitle(NULL)  +
	stat_compare_means(data = rank_vals,
										 aes(x = pretty_name, y = RSQ),
										 method = "anova", inherit.aes = F, size=5) + # Add global p-value
	theme(plot.margin = margin(1,2,.5,1, "cm")) +
	geom_rug(data = rank_vals %>% filter(pretty_group == "Bacteria"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="r", show.legend = F, color = "#F8766D") +
	geom_rug(data = rank_vals %>% filter(pretty_group == "Fungi"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="l", show.legend = F, color = "#00BFC4") #+
#geom_rug(data = rank_vals, aes(x = pretty_name, y = RSQ, color = pretty_group), alpha=.5, sides="l", show.legend = F)
newtime_by_rank +
	geom_point(data = rank_vals,
						 aes(x = pretty_name, y = RSQ, color = pretty_group), size=3, alpha=.1, show.legend = F) +
	geom_text(data = tukey_list_newtime,
																																																																	 aes(x = x, y = tot + .1,
																																																																	 		label = Letters_Tukey), show.legend = F,
																																																																	 color = 1, size =6)






# Compare new-time only across taxonomic ranks and functional groups
newsite_by_rank = ggplot(skill_score_to_plot,
												 aes(x = pretty_name, y = mean_skill_score, color = pretty_group)) +
	facet_grid(cols=vars(pretty_group), scales="free") +
	theme_minimal(base_size = 18) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("Relative skill at new sites") + xlab("Rank") +
	#scale_color_discrete(name=NULL, labels=c("Within-site","New site")) +
	#scale_y_continuous(trans=scales::pseudo_log_trans(base = 2)) +
	geom_pointrange(aes(ymin = ymin, ymax = ymax, color = pretty_group), show.legend = F) + ggtitle(NULL)  +
	stat_compare_means(data = skill_score_vals,
										 aes(x = pretty_name, y = skill_score),
										 method = "anova", inherit.aes = F, size=5) + # Add global p-value
	theme(plot.margin = margin(1,2,.5,1, "cm")) +
	geom_rug(data = skill_score_vals %>% filter(pretty_group == "Bacteria"),
					 aes(x = pretty_name, y = skill_score), alpha=.5, sides="r", show.legend = F, color = "#F8766D") +
	geom_rug(data = skill_score_vals %>% filter(pretty_group == "Fungi"),
					 aes(x = pretty_name, y = skill_score), alpha=.5, sides="l", show.legend = F, color = "#00BFC4") #+
#geom_rug(data = rank_vals, aes(x = pretty_name, y = RSQ, color = pretty_group), alpha=.5, sides="l", show.legend = F)
newsite_by_rank +
	# geom_point(data = skill_score_vals,
	# 														aes(x = pretty_name, y = skill_score, color = pretty_group), size=3, alpha=.1, show.legend = F) +
	scale_y_continuous(trans = pseudolog10_trans)  + geom_text(data = tukey_list_newsite,
																														 aes(x = x, y = tot + 1.5,
																														 		label = Letters_Tukey), show.legend = F,
																														 color = 1, size =6)




cal_scores = scores_list$calibration_metrics

# Get mean values for plotting
cal_score_vals = scores_list$calibration_metrics %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity") %>%
	filter(!taxon_model_rank %in% unconverged$taxon_model_rank)

cal_score_vals_to_plot = cal_score_vals %>% ungroup() %>%
	group_by(fcast_type, pretty_group, model_name, pretty_name, rank)  %>%
	summarize(mean_RSQ=mean(RSQ, na.rm=T),
						sd_RSQ = sd(RSQ, na.rm=T)) %>%
	mutate(ymax = mean_RSQ + 1.96*sd_RSQ,
				 ymin = mean_RSQ - 1.96*sd_RSQ) %>% arrange(pretty_group, pretty_name)


ggplot(cal_score_vals_to_plot,
			 aes(x = pretty_name, y = mean_RSQ, color = pretty_group)) +
	facet_grid(cols=vars(pretty_group), scales="free") +
	theme_minimal(base_size = 18) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("Calibration RSQ") + xlab("Rank") +
	#scale_color_discrete(name=NULL, labels=c("Within-site","New site")) +
	#scale_y_continuous(trans=scales::pseudo_log_trans(base = 2)) +
	geom_pointrange(aes(ymin = ymin, ymax = ymax, color = pretty_group), show.legend = F) + ggtitle(NULL)  +
	theme(plot.margin = margin(1,2,.5,1, "cm")) +
	geom_rug(data = cal_score_vals %>% filter(pretty_group == "Bacteria"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="r", show.legend = F, color = "#F8766D") +
	geom_rug(data = cal_score_vals %>% filter(pretty_group == "Fungi"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="l", show.legend = F, color = "#00BFC4") #+
#geom_rug(data = rank_vals, aes(x = pretty_name, y = RSQ, color = pretty_group), alpha=.5, sides="l", show.legend = F)


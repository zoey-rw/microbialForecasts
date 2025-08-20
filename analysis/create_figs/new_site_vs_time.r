
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
#p_load(forestmangr)
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
#scores_list = readRDS(here("data/summary/scoring_metrics_cv_jan18.rds"))


n_timepoints <- scores_list$calibration_truth_vals %>%
	select(group,siteID,taxon,dateID) %>% distinct() %>%
	group_by(group,siteID,taxon) %>%
	count(name = "n_timepoints")

n_sampling_events <- scores_list$calibration_truth_vals %>%
	select(group,siteID,taxon,dateID,plotID) %>% distinct() %>%
	group_by(group,siteID,taxon) %>%
	count(name = "n_sampling_events")

fcast_success = scores_list$scoring_metrics_site %>% filter(!taxon %in% scores_list$unconverged_list$taxon) %>%  merge(n_timepoints) %>%
	merge(n_sampling_events) %>%
	filter(site_prediction == "New time (observed site)")

ggplot(fcast_success, aes(x = n_timepoints, y = `RSQ.1`)) +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Are predictions better at sites with more sampling dates?", subtitle = "Yes, slightly") +
	theme_minimal() +
	facet_wrap(~pretty_name, scales="free") + geom_smooth(method="lm") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

ggplot(fcast_success, aes(x = n_sampling_events, y = `RSQ.1`)) +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Are predictions better at sites with more sampling dates/plots?",
					subtitle = "no") +
	theme_minimal() +
	facet_wrap(~pretty_name, scales="free") + geom_smooth(method="lm") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))




#keep = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
unconverged = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
converged = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
#keep = keep[keep$median_gbr < 3,]
unconverged = unconverged[unconverged$median_gbr > 7,]
#unconverged = unconverged[unconverged$median_gbr > 3,]
#converged = converged[converged$median_gbr < 4,]
all_cov = converged[converged$model_name=="all_covariates",]

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
	stat_compare_means(data = rank_vals,
										 aes(x = pretty_name, y = RSQ, label = after_stat(p.signif)),
										 method = "t.test", ref.group = "Functional group",
										 label.y.npc = .8, hide.ns = T, size=10, inherit.aes = F) +
	theme(plot.margin = margin(1,2,.5,1, "cm")) +
	geom_rug(data = rank_vals %>% filter(pretty_group == "Bacteria"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="r", show.legend = F, color = "#F8766D") +
	geom_rug(data = rank_vals %>% filter(pretty_group == "Fungi"),
					 aes(x = pretty_name, y = RSQ), alpha=.5, sides="l", show.legend = F, color = "#00BFC4") #+
#geom_rug(data = rank_vals, aes(x = pretty_name, y = RSQ, color = pretty_group), alpha=.5, sides="l", show.legend = F)
newtime_by_rank + 	geom_point(data = rank_vals,
															aes(x = pretty_name, y = RSQ, color = pretty_group), size=3, alpha=.1, show.legend = F)

my_comparisons <- list( c("Genus", "Family"), c("Family", "Order"),
												c("Order", "Class") ,  c("Class","Phylum"),  c("Phylum", "Functional group"))



skill_score_to_plot = scores_list$skill_score_taxon %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity") %>%
	pivot_longer(cols=c(#skill_score, skill_score_random,
											`New time (observed site)`, `New time x site (modeled effect)`),											#`New time x site (random effect)`),
							 names_to = "Prediction type")
rank_vals_newsite = site_time_scores_to_plot %>% filter(!grepl("observed", site_prediction))
rank_means_newsite = rank_vals_newsite %>% ungroup() %>%
	group_by(fcast_type, pretty_group, model_name, pretty_name, rank, site_prediction) %>%
	summarize(mean_RSQ=mean(RSQ, na.rm=T),
						sd_RSQ = sd(RSQ, na.rm=T)) %>%
	mutate(ymax = mean_RSQ + 1.96*sd_RSQ,
				 ymin = mean_RSQ - 1.96*sd_RSQ) %>% arrange(pretty_group, pretty_name)




# Compare across taxonomic ranks and functional groups
site_time_by_rank = ggplot(site_time_scores_to_plot,
											 aes(x = pretty_name, color = site_prediction, y = RSQ)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(width = .1, height=.01, alpha=.01, size=3, position=position_jitterdodge()) +
	facet_grid(rows=vars(pretty_group), scales="free") +
	theme_minimal(base_size = 18) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ") + xlab("Rank") + scale_color_discrete(name=NULL, labels=c("Within-site","New site")) +
	scale_y_continuous(trans=scales::pseudo_log_trans(base = 2))
site_time_by_rank
site_time_by_rank + ggtitle(NULL)  +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .2, show.legend = F, hide.ns = T, size=5)


# Compare f vs b
site_time_f_b = ggplot(site_time_scores_to_plot %>% filter(!taxon_model_rank %in% unconverged$taxon_model_rank),
													 aes(x = pretty_group, y = RSQ)) +
	geom_violin() +
	geom_jitter(width = .1, height=.01, alpha=.5, size=3) +
	facet_wrap(~site_prediction, scales="free", ncol = 1) +
	theme_minimal() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ") + xlab("Domain") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .2, show.legend = F, hide.ns = T, size=5)
site_time_f_b



skill_score_to_plot = scores_list$skill_score_taxon %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity") %>%
	pivot_longer(cols=c(#skill_score, skill_score_random,
											`New time (observed site)`, `New time x site (modeled effect)`),											#`New time x site (random effect)`),
							 names_to = "Prediction type")
#skill_score_to_plot <- skill_score_to_plot %>% filter(!taxon %in% scores_list$unconverged_list$taxon)


# Compare across taxonomic ranks and functional groups
compare_ranks = ggplot(skill_score_to_plot,
			 aes(x = pretty_name, color = `Prediction type`, y = value)) +
	geom_violin() +
	geom_point(width = .1, height=.01, alpha=.5, size=3, position=position_jitterdodge()) +
	facet_wrap(~pretty_group, scales="free", ncol = 1) +
	theme_minimal() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("CRPS") + xlab("Rank")
compare_ranks + ggtitle("Fungal functional group models are most transferable")  +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..), label.x = 1.5, label.y = .2, show.legend = F, hide.ns = T, size=5)

scores_to_plot = scores_list$scoring_metrics_long %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity" &
				 	site_prediction == "New time (observed site)" &
				 	metric %in% c("RSQ","RSQ.1","CRPS"))

bias = scores_list$scoring_metrics_long %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity" &
				 	site_prediction == "New time (observed site)" &
				 	metric %in% c("RSQ","BIAS"))
#scores_to_plot <- scores_to_plot %>% filter(!taxon %in% scores_list$unconverged_list$taxon)


compare_scores = ggplot(scores_to_plot,
												aes(x = fcast_type, color = pretty_group, y = score)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(width = .1, height=.01, alpha=.5, size=3, position=position_jitterdodge()) +
	facet_grid(pretty_group~metric, scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("Score") + xlab("Forecast type")

new_times= compare_scores  +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..), label.x = 1.5, label.y = .2,
										 show.legend = F, hide.ns = F, size=5)


# Compare across taxonomic ranks VS functional groups
# Best so far
new_site <- ggplot(skill_score_to_plot,
			 aes(x = fcast_type, y = skill_score)) +
	geom_violin(aes(color = pretty_group), draw_quantiles = c(.5)) +
	geom_point(aes(color = pretty_group), alpha=.5, size=3, position=position_jitterdodge()) +
	#geom_point(aes(color = pretty_group),alpha=.5, size=3) +
	ggtitle("Predicting to new times and places", "Functional groups outperform taxonomic groups") +
	facet_grid(rows=vars(pretty_group), scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("Decrease in accuracy at new sites") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .2, show.legend = F, #hide.ns = T,
										 size=5) +
	guides(color=guide_legend(title="Domain")) + xlab("Forecast type")

ggpubr::ggarrange(new_times, new_site, common.legend = T)

# Compare fungi vs bacteria
# Compare fungi vs bacteria
ggplot(site_vals_to_plot %>% filter(!grepl("random", site_prediction)),
			 aes(x = site_prediction, y = CRPS_truncated, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(alpha=.1, size=3, position=position_jitterdodge()) +
	ggtitle("Predicting to new times and sites: site-level CRPS") +
	facet_grid(~fcast_type, scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("CRPS") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .12, show.legend = F, #hide.ns = T,
										 size=6) +
	xlab("Forecast type")
# Compare fungi vs bacteria
ggplot(site_vals_to_plot %>% filter(!grepl("random", site_prediction)),
			 aes(x = site_prediction, y = RSQ, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(alpha=.1, size=3, position=position_jitterdodge()) +
	ggtitle("Predicting to new times and sites: site-level RSQ") +
	facet_grid(~fcast_type, scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .12, show.legend = F, #hide.ns = T,
										 size=6) +
	xlab("Forecast type")










# Compare across taxonomic ranks VS functional groups
ggplot(skill_score_to_plot,
			 aes(x = fcast_type, y = value)) +
	geom_violin(aes(color = pretty_group), draw_quantiles = c(.5)) +
	geom_point(aes(color = pretty_group), alpha=.5, size=3, position=position_jitterdodge()) +
	ggtitle("Predicting to new sites and times (group mean)") +
	facet_grid(`Prediction type`~pretty_group, scales="free") +
	theme_minimal() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("CRPS") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .05, show.legend = F, #hide.ns = T,
										 size=6) + xlab("Forecast type")


site_vals_to_plot = scores_list$scoring_metrics_site %>%
	filter(model_name=="all_covariates" & fcast_type !="Diversity")
ggplot(site_vals_to_plot %>% filter(grepl("modeled", site_prediction)),
			 aes(x = fcast_type, y = CRPS_truncated, #color = pretty_group,
			 		color = site_prediction)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(alpha=.5, size=3, position=position_jitterdodge()) +
	ggtitle("Predicting to new times and sites: site-level CRPS") +
	facet_grid(site_prediction~pretty_group, scales="free") +
	theme_minimal() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("CRPS") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .12, show.legend = F, #hide.ns = T,
										 size=6)  + xlab("Forecast type")



f_vs_b_rsq = ggplot(site_vals_to_plot %>% filter(!grepl("random", site_prediction)),
			 aes(x = pretty_group, y = RSQ, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(alpha=.1, size=1, position=position_jitter(), show.legend = F) +
	ggtitle("Predicting to new times and sites: site-level RSQ") +
	facet_grid(~site_prediction, scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .8, show.legend = F, #hide.ns = T,
										 size=6)  + xlab("Domain")


# Fuctional vs tax
ggplot(site_vals_to_plot %>% filter(!grepl("random", site_prediction)),
			 aes(x = fcast_type, y = RSQ.1)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_point(alpha=.1, size=3, position=position_jitter()) +
	ggtitle("Predicting to new times and sites: site-level RSQ.1") +
	facet_grid(pretty_group~site_prediction, scales="free") +
	theme_minimal(base_size = 16) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	ylab("RSQ.1") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x = 1.5, label.y = .12, show.legend = F, #hide.ns = T,
										 size=6) +
	xlab("Forecast type") +
	scale_y_log10()




ecto_fcast <- readRDS(here("data/summary/all_hindcasts.rds"))


ecto_fcast_not_na <- ecto_fcast %>% filter(!is.na(truth) & !is.na(mean))

ecto_fcast_scored_site = ecto_fcast_not_na %>%
	group_by(model_name,group,category,taxon,siteID, new_site,# predicted_site_effect,
					 site_prediction) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = mean,
																sd_predicted = sd))
ecto_fcast_scored_site$RSQ.1 = ifelse(ecto_fcast_scored_site$RSQ.1 < 0, 0, ecto_fcast_scored_site$RSQ.1)

ggplot(ecto_fcast_scored_site, aes(x = site_prediction, y = `RSQ.1`)) +
	geom_violin() +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Predicting spatiotemporally") + theme_minimal() + facet_wrap(~category)

ggplot(ecto_fcast_scored_site, aes(x = site_prediction, y = `CRPS`)) +
	geom_violin() +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Predicting spatiotemporally") + facet_wrap(~category) +theme_minimal()



new_sites_bacteria <- ecto_fcast_not_na %>% ungroup %>%
	filter(site_prediction == "New time x site (random effect)" & group == "16S") %>% select(siteID) %>% unique

new_sites_bacteria_modeled <- ecto_fcast_not_na %>% ungroup %>%
	filter(site_prediction == "New time x site (modeled effect)" & group == "16S") %>% select(siteID) %>% unique

new_sites_fungi <- ecto_fcast_not_na %>% ungroup %>%
	filter(site_prediction == "New time x site (random effect)" & group == "ITS") %>% select(siteID) %>% unique

new_sites_fungi_modeled <- ecto_fcast_not_na %>% ungroup %>%
	filter(site_prediction == "New time x site (random effect)" & group == "ITS") %>% select(siteID) %>% unique

new_sites_bacteria
new_sites_bacteria_modeled
new_sites_fungi
new_sites_fungi_modeled

soap = ecto_fcast[ecto_fcast$siteID=="SOAP",]
soap_oligo = soap[soap$taxon=="oligotroph",]


soap = ecto_fcast_not_na[ecto_fcast_not_na$siteID=="SOAP",]
soap_oligo = soap[soap$taxon=="oligotroph",]
soap_benomyl = soap[soap$taxon=="benomyl",]


soap_ecto[soap_ecto$taxon=="ectomycorrhizal",]


ecto_fcast_scored_site[ecto_fcast_scored_site$siteID=="SOAP",]

n_timepoints <- ecto_fcast %>% filter(!is.na(truth) & fcast_period=="calibration") %>%
	select(group,siteID,taxon,dateID) %>% distinct() %>%
	group_by(group,siteID,taxon) %>%
	count(name = "n_timepoints")

n_sampling_events <- ecto_fcast %>% filter(!is.na(truth) & fcast_period=="calibration") %>%
	select(group,siteID,taxon,dateID,plotID) %>% distinct() %>%
	group_by(group,siteID,taxon) %>%
	count(name = "n_sampling_events")

fcast_success = ecto_fcast_scored_site %>% merge(n_timepoints) %>% merge(n_sampling_events)

ggplot(fcast_success %>% filter(!new_site), aes(x = n_sampling_events, y = `RSQ.1`)) +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Are predictions better at sites with more sampling dates/plots?", subtitle = "Also no") + theme_minimal() +
	facet_wrap(~taxon, scales="free") + geom_smooth(method="lm")

ggplot(fcast_success %>% filter(!new_site), aes(x = n_timepoints, y = `RSQ.1`)) +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Are predictions better at sites with more sampling dates?", subtitle = "No") + theme_minimal() +
	facet_wrap(~taxon, scales="free") + geom_smooth(method="lm")

ecto_fcast %>% filter(!is.na(truth) & fcast_period=="calibration" & taxon=="oligotroph") %>% select(group,siteID,taxon,dateID) %>% distinct()



ggplot(ecto_fcast_scored_site %>% filter(!site_prediction == "New time x site (random effect)"),
			 aes(x = site_prediction, y = `CRPS`)) +
	geom_violin() +
	geom_jitter(aes(color = siteID), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Predicting spatiotemporally") + facet_wrap(~category) +theme_minimal()



ecto_fcast_scored_taxon = ecto_fcast_not_na %>%
	group_by(model_name,group,category,taxon, new_site,# predicted_site_effect,
					 site_prediction) %>%
	summarize(add_scoring_metrics(observed = truth,
																mean_predicted = mean,
																sd_predicted = sd))
ecto_fcast_scored_taxon$RSQ.1 = ifelse(ecto_fcast_scored_taxon$RSQ.1 < 0, 0, ecto_fcast_scored_taxon$RSQ.1)



# Skill score from Scavia 2021
skill_score_taxon <- ecto_fcast_scored_taxon %>%
	pivot_wider(id_cols = c("category","group","model_name","taxon"),
							values_from = "CRPS_truncated", names_from = "site_prediction") %>%
	mutate(skill_score = (1 - (`New time x site (modeled effect)`/`New time (observed site)`)),
				 skill_score_random_eff = (1 - (`New time x site (random effect)`/`New time (observed site)`)),
	)
skill_score_rank <- skill_score_taxon %>%
	group_by(group, category) %>%
	summarize(mean_skill_score = mean(skill_score, na.rm=T),
						mean_skill_score_random_eff = mean(skill_score_random_eff, na.rm=T)) %>% pivot_longer(cols=c(mean_skill_score, mean_skill_score_random_eff))


ggplot(skill_score_rank,
			 aes(x = category)) +
	#geom_violin() +
	geom_jitter(aes(y = value, color = name), width = .1, height=.01, alpha=.5, size=3) +
	ggtitle("Predicting to new  and places") +
	facet_wrap(~group, scales="free") +theme_minimal() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))


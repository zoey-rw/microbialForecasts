source("source.R")

library(ggallin)

# Load scores data first
scores_list <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
converged = scores_list$converged_list
#converged = scores_list$converged_strict_list
# hindcast_in not available, using available data instead
hindcast_data = readRDS(here("data/summary/all_hindcasts.rds"))
hindcast_data[hindcast_data$truth==0,]$truth <- .0001

old_scores_list <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
new_scores_list <- readRDS(here("data/summary/scoring_metrics_plsr.rds"))
new_scores_list2 <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))

old_site_scores = old_scores_list$scoring_metrics_site_long %>%
	filter(model_id %in% converged) %>% mutate(site_prediction_method = "model_average")
new_site_scores = new_scores_list$scoring_metrics_site_long %>% filter(model_id %in% converged) %>%
	mutate(site_prediction_method = "plsr")
site_scores = rbind(old_site_scores, new_site_scores)



#site_scores = new_scores_list$scoring_metrics_site_long #%>%
	#filter(site_prediction != "New time x site (modeled effect)")


fcast_n_obs = hindcast_data %>%
	filter(!is.na(truth) & fcast_period=="hindcast") %>%
	group_by(model_id, siteID) %>%
	tally(name = "n_obs")
# Remove metrics calculated from only <4 observations
site_scores_filt <- merge(site_scores, fcast_n_obs) %>% filter(n_obs > 3)

site_scores_mean = site_scores_filt %>%
	filter(model_id %in% converged) %>%
	group_by(pretty_group, rank_name, taxon,pretty_name, fcast_type, model_id, model_name, site_prediction, metric, site_prediction_method) %>% summarize(score = mean(score, na.rm=T), mean_crps = mean(mean_crps_sample, na.rm=T)) %>%
	distinct()


ggplot(site_scores_mean  %>%
			 	filter(metric %in% c("CRPS_truncated")),
			 aes(x = site_prediction,y = score,
			 		color = site_prediction_method)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), alpha=.1) +
	xlab(NULL) +
	ylab("Hindcast predictability (RSQ 1:1)") +
	xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	facet_wrap(model_name~pretty_group, nrow=2, scales="free_y") +
	scale_x_discrete(labels= model.labs)



# Skill score from Scavia 2021
skill_score_taxon <- site_scores_mean %>%
	filter(metric=="CRPS_truncated") %>%
	pivot_wider(id_cols = c("model_id","fcast_type","pretty_group","model_name","pretty_name","rank_name","taxon","site_prediction_method"),
							values_from = "mean_crps", names_from = "site_prediction") %>%
	mutate(skill_score = (1 - (`New time x site (modeled effect)`/`New time (observed site)`)),
				 skill_score_random = (1 - (`New time x site (random effect)`/`New time (observed site)`)),
	)
skill_score_taxon_long <- skill_score_taxon %>%
	pivot_longer(cols=c("skill_score","skill_score_random"), names_to = "newsite_effect_type", values_to = "skill_score")


ggplot(skill_score_taxon_long %>% filter(site_prediction_method=="plsr"),
			 aes(x = newsite_effect_type,y = skill_score,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), alpha=.1) +
	xlab(NULL) +
	ylab("Skill at new sites (% change in CRPS)") +
	#xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	facet_grid(model_name~pretty_group, scales="free_y") +
	scale_x_discrete(labels= model.labs)




stat_pvalue <- skill_score_taxon %>%
	ungroup() %>%
	#filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	rstatix::tukey_hsd(skill_score ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .4) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

a1 = ggplot(skill_score_taxon,#  %>% filter(fcast_type=="Functional") ,
						aes(x = model_name, y = skill_score,
								color = model_name)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(aes(x = model_name, y = skill_score),
						 position = position_jitterdodge(jitter.height = 0, jitter.width = .2),
						 alpha = .2, size=4, show.legend = F) +
	ylab("Skill at new sites (% change in CRPS)") + xlab("Linear model components") +
	theme_classic(base_size=22) + #ggtitle("Model transferability to new sites") +
	theme(text = element_text(size = 22),
					axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
						angle = 320, vjust=1.5, hjust = -0.03),
					strip.text.y = element_text(size=12,face="bold"), plot.margin = unit(c(1,2,1,1), "cm")) +
	scale_x_discrete(labels= model.labs) +
	geom_hline(yintercept = 0) +
	xlab(NULL) +
	#facet_grid(pretty_name~pretty_group) +
	#facet_grid(rows=vars(pretty_group)) +
	scale_y_continuous(trans = pseudolog10_trans)  +
	ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.adj.signif", bracket.nudge.y = -.3, size=6)#)	#+ ylim(c(-15, 5))

png(here("figures","linear_model_skill_score_bygroup.png"), width = 7, height=7, res = 200, units = "in")
print(a1)

dev.off()



## ESA 2023 FIG
taxon_scores = scores_list$scoring_metrics %>% filter(model_id %in% converged)



stat_pvalue_model <- taxon_scores %>%
	ungroup() %>%
	#group_by(pretty_group) %>%
	rstatix::tukey_hsd(mean_crps_sample ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .7) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
ggplot(taxon_scores,
			 aes(x = model_name,y = mean_crps_sample,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), show.legend = F, alpha=.3) +
	xlab(NULL) +
	ylab("Forecast error (CRPS)") +
	#xlab("Microbial kingdom") +
	xlab(NULL) +
	theme_classic(base_size=22) + theme(text = element_text(size = 22),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1.5, hjust = -0.03),
										 strip.text.y = element_text(size=12,face="bold"), plot.margin = unit(c(1,2,1,1), "cm")) +
	scale_x_discrete(labels= model.labs) +
	#scale_y_continuous(trans = pseudolog10_trans)  +
	scale_y_log10()  +
	ggpubr::stat_pvalue_manual(stat_pvalue_model, label = "p.adj.signif",
														 bracket.nudge.y = -.6,
														 size=7)#)	#+ ylim(c(-15, 5))





stat_pvalue_skill_score_taxon <- skill_score_taxon %>%
	filter(model_id %in% converged) %>%
	ungroup() %>%
	#group_by(pretty_group) %>%
	rstatix::tukey_hsd(skill_score ~ pretty_group) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .4) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
ggplot(skill_score_taxon %>% filter(model_name == "env_cycl"),
			 aes(x = pretty_group,y = skill_score,
			 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), show.legend = F) +
	xlab(NULL) +
	ylab("Skill at new sites (decrease in CRPS)") +
	xlab("Microbial kingdom") +
	theme_classic(base_size=22) +	theme(text = element_text(size = 22),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1.5, hjust = -0.03),
										 strip.text.y = element_text(size=12,face="bold"), plot.margin = unit(c(1,2,1,1), "cm")) +
	scale_x_discrete(labels= model.labs) +
	scale_y_continuous(trans = pseudolog10_trans)  +
	ggpubr::stat_pvalue_manual(stat_pvalue_skill_score_taxon, label = "p.adj.signif",
														 bracket.nudge.y = -.5,
														 size=4)#)	#+ ylim(c(-15, 5))






stat_pvalue_model_functional <- taxon_scores %>%
	ungroup() %>%
	filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	rstatix::tukey_hsd(RSQ ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .7) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
ggplot(taxon_scores %>% filter(fcast_type=="Functional"),
			 aes(x = model_name,y = RSQ,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), show.legend = F, alpha=.3) +
	xlab(NULL) +
	ylab("Forecast RSQ") +
	#xlab("Microbial kingdom") +
	xlab(NULL) +
	theme_classic(base_size=22) + theme(text = element_text(size = 22),
																			axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
																				angle = 320, vjust=1.5, hjust = -0.03),
																			strip.text.y = element_text(size=12,face="bold"), plot.margin = unit(c(1,2,1,1), "cm")) +
	scale_x_discrete(labels= model.labs) +
	#scale_y_continuous(trans = pseudolog10_trans)  +
	ggpubr::stat_pvalue_manual(stat_pvalue_model_functional, label = "p.adj.signif",
														 #bracket.nudge.y = -.6,
														 size=7)#)	#+ ylim(c(-15, 5))





stat_pvalue_fg <- skill_score_taxon %>% filter(model_name == "env_cycl") %>%
	ungroup() %>%
	group_by(pretty_group) %>%
	rstatix::tukey_hsd(skill_score ~ fcast_type) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .01) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

ggplot(skill_score_taxon %>% filter(model_name == "env_cycl"),
			 aes(x = fcast_type,y = skill_score,
			 		color = pretty_group)) +
	facet_grid(~pretty_group) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(
		position = position_jitterdodge(jitter.height = 0, jitter.width = .2), show.legend = F, alpha=.4) +
	xlab(NULL) +
	ylab("Skill at new sites (decrease in CRPS)") +
	xlab("Microbial kingdom") +
	theme_classic(base_size=22) +	theme(text = element_text(size = 22),
																			axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
																				angle = 320, vjust=1.5, hjust = -0.03),
																			strip.text.y = element_text(size=12,face="bold"), plot.margin = unit(c(1,2,1,1), "cm")) +
	scale_x_discrete(labels= model.labs) +
	scale_y_continuous(trans = pseudolog10_trans)  +
	ggpubr::stat_pvalue_manual(stat_pvalue_fg, label = "p.adj.signif", step.group.by = "pretty_group",
														 #bracket.nudge.y = -.5,
														 size=4)#)	#+ ylim(c(-15, 5))


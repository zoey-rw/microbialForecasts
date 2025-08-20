# Compare different models to evaluate the benefit of including environmental predictors vs just seasonality
library(ggallin)
source("source.R")


# Scoring metrics aggregated by taxon/model
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr2.rds")))
converged = scores_list$converged_list
converged = scores_list$converged_strict_list


hindcast_filter <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(model_name != "all_covariates" &
		metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm") &
			site_prediction == "New time (observed site)")


# Raw hindcast crps values
hindcast_data <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))
hindcast_data_df = hindcast_data %>%
	filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period=="hindcast" & new_site==FALSE)


scatter_overall =
	ggplot(hindcast_data_df, aes(x = med, y = truth, group=model_name)) +
	geom_point(aes(color=pretty_group),
						 alpha=.15, show.legend = F) +
	geom_smooth(method = "lm", se = FALSE, color = 1) +
	theme_classic(base_size = 18) +
	facet_grid(pretty_group~model_name, labeller=labeller(model_name = model.labs), scales="free") +
	xlab("Plot forecast") +
	ylab("Plot observation") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), p.accuracy = 0.001,
		label.y = .7, label.x.npc = .2, size=4, label.sep='\n') +
	geom_abline(slope=1, intercept=0, linetype=2)
#	ggtitle("Overall prediction accuracy")
scatter_overall <- tag_facet(scatter_overall, tag_pool = LETTERS, x = .001)
scatter_overall


png(here("figures","model_r2_by_kingdom.png"), width = 1200, height=800)
print(scatter_overall)
dev.off()

# Run "fig_horizon_by_seasonality.r" to generate fig3g
fig3 = ggarrange(plotlist=list(scatter_overall, fig3g), nrow=2, heights=c(2,1), labels=c(NA, "(G)"))

png(here("figures","fig3.png"), width = 1600, height=1600)
print(fig3)
dev.off()

ggscatter(hindcast_data_df, x = "med", y = "truth",
					add = "reg.line", color = "pretty_group"
) + facet_grid(pretty_group~model_name) +
	xlab("Median prediction estimate") + ylab("Observed plot mean") +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
	ggtitle("Overall prediction accuracy")  +
	geom_abline(slope=1, intercept=0, linetype=2)




skill_scores = scores_list$skill_score_taxon %>%
	filter(model_name != "all_covariates") %>%
	filter(model_id %in% converged)


skill_scores_rmse = scores_list$skill_score_taxon_RMSE %>%
	filter(model_name != "all_covariates") %>%
	filter(model_id %in% converged)

site_scores = scores_list$scoring_metrics_site_long %>%
	filter(!siteID %in% "MLBS") %>%
	filter(model_id %in% converged &
				 	site_prediction != "New time x site (modeled effect)") %>%
	filter(metric %in% c("RSQ"))

site_scores_allmetrics = scores_list$scoring_metrics_site_long %>%
	filter(!siteID %in% "MLBS") %>%
	filter(model_id %in% converged &
				 	site_prediction == "New time x site (random effect)")




stat_pvalue <- site_scores_allmetrics  %>%
	filter(metric %in% "CRPS_truncated" & pretty_name %in% c("phylum","functional")) %>%
	group_by(pretty_group) %>%
	rstatix::tukey_hsd(mean_crps_sample ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_xy_position() %>%
	#rstatix::add_y_position(step.increase = .1) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

# View model structure effects on new-site hindcast accuracy
fig2 <- ggplot(site_scores_allmetrics %>%
							 	filter(metric %in% c("CRPS_truncated") & pretty_name %in% c("phylum","functional")),
							 aes(x = model_name, y = mean_crps_sample, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = .2),
						 alpha = .1, size=3, show.legend = F) +
	ylab("Absolute forecast error (CRPS)") +
	xlab("Linear model components") +
	theme_bw(base_size=22) + #ggtitle("Model transferability to new sites") +
	theme(axis.text.x=element_text(angle = 280, vjust=.8, hjust = .2),
				axis.title=element_text(size=22))  +
	scale_x_discrete(labels= model.labs) +
	facet_grid(~pretty_group, scales="free") +
	#facet_grid(cols=vars(pretty_group)) +
	ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.adj.signif", bracket.nudge.y = -.8,
														 size=8,
														 hide.ns = T#, y.position = c(.001,.0005,.0006)
														 ) +
#	scale_y_continuous(trans = pseudolog10_trans) +
scale_y_log10()

fig2

# View model structure effects on new-site hindcast accuracy
a <- ggplot(skill_scores,
			 aes(x = model_name, y = skill_score*100, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(aes(x = model_name, y = skill_score*100),
						 position = position_jitterdodge(jitter.height = 0, jitter.width = .2),
						 alpha = .2, size=4, show.legend = F) +
	ylab("Skill at new sites (% change in CRPS)") + xlab("Linear model components") +
	theme_bw(base_size=18) + #ggtitle("Model transferability to new sites") +
	theme(axis.text.x=element_text(angle = 300, vjust=.5, hjust = .2),
				axis.title=element_text(size=16))  +
	scale_x_discrete(labels= model.labs) +
	geom_hline(yintercept = 0) +
	facet_grid(cols=vars(pretty_group)) +
	stat_compare_means(aes(x = model_name, y = skill_score*100, group=pretty_group,
										 label = paste0("p = ", after_stat(p.signif))),
										 inherit.aes = F, comparisons = list(c('env_cov','env_cycl'),
										 																			c('env_cov','cycl_only'),
										 																			c('env_cycl','cycl_only')), step.increase = .05,
										 hide.ns = F)
a <- tag_facet(a)  +
	theme(plot.margin = unit(c(1,2,1,1), "cm"))


png(here("figures","linear_model_skill_score.png"), width = 7, height=8, res = 200, units = "in")
print(a)

dev.off()

# a1 <- ggplot(skill_scores,
# 						 aes(x = model_name, y = skill_score,
# 						 		shape = model_name, color = pretty_group)) +
# 	geom_violin(draw_quantiles = c(.5), show.legend = F) +
# 	geom_point(aes(x = model_name, y = skill_score),
# 						 position = position_jitterdodge(jitter.height = 0, jitter.width = .2),
# 						 alpha = .2, size=4, show.legend = F) +
# 	ylab("Skill at new sites (% change in CRPS)") + xlab("Linear model components") +
# 	theme_bw(base_size=18) + #ggtitle("Model transferability to new sites") +
# 	theme(axis.text.x=element_text(angle = 320, vjust=.5, hjust = .2),
# 				axis.title=element_text(size=16))  +
# 	scale_x_discrete(labels= model.labs) +
# 	geom_hline(yintercept = 0) +
# 	#facet_grid(pretty_name~pretty_group) +
# 	facet_grid(rows=vars(pretty_group)) +
# 	scale_y_continuous(trans = pseudolog10_trans)  +
# 	stat_compare_means(data = skill_scores,
# 										 aes(x = model_name, y = skill_score, color = pretty_group),
# 										 method = "anova", inherit.aes = F, size=5, label.y.npc = .3, label.x.npc = .7, show.legend = F) +
# 	stat_compare_means(aes(label = after_stat(p.signif)),
# 										 method = "t.test", ref.group = "cycl_only", label.y.npc = .4,
# 										 show.legend = F, hide.ns = T, size=10)
#

stat_pvalue <- skill_scores %>% filter(model_name != "env_cov") %>%
	group_by(pretty_group) %>%
	rstatix::tukey_hsd(skill_score_random ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .4) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

a1 = ggplot(skill_scores %>% filter(model_name != "env_cov"),
			 aes(x = model_name, y = skill_score_random,
			 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_point(aes(x = model_name, y = skill_score),
						 position = position_jitterdodge(jitter.height = 0, jitter.width = .2),
						 alpha = .2, size=4, show.legend = F) +
	ylab("Skill at new sites (% change in CRPS)") + xlab("Linear model components") +
	theme_classic(base_size=22) + #ggtitle("Model transferability to new sites") +
	theme(axis.text.x=element_text(angle = 290, vjust=.5, hjust = .2),
				axis.title=element_text(size=16))  +
	scale_x_discrete(labels= model.labs) +
	geom_hline(yintercept = 0) +
	#facet_grid(pretty_name~pretty_group) +
	#facet_grid(rows=vars(pretty_group)) +
	facet_wrap(~pretty_group, nrow=1, scales="free_y") +
	scale_y_continuous(trans = pseudolog10_trans)  +
	ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.adj.signif", bracket.nudge.y = -.4, size=4)#)	#+ ylim(c(-15, 5))

a1 <- tag_facet(a1, tag_pool = c("C","D"))

png(here("figures","linear_model_skill_score_bygroup.png"), width = 3, height=7, res = 200, units = "in")
print(a1)

dev.off()

# View model structure effects on within-site hindcast accuracy
b <- ggplot(hindcast_filter %>%
			 	filter(metric %in% c("RSQ")),
			 aes(x = model_name,y = score,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	xlab(NULL) +
	ylab("Hindcast predictability (RSQ 1:1)") +
	xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
 facet_grid(#pretty_name
 					 ~pretty_group) + scale_x_discrete(labels= model.labs) +
	stat_compare_means(data = hindcast_filter %>%
										 	filter(metric %in% c("RSQ.1")),
										 aes(x = model_name, y = score),
										 method = "anova", inherit.aes = F, size=5, label.y.npc = .7, label.x.npc = .7, show.legend = F) +
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", ref.group = "cycl_only", label.y.npc = .5,
										 show.legend = F, hide.ns = T, size=10)

#b <- tag_facet(b, size=7)

png(here("figures","linear_model_comparison.png"), width = 800, height=1200)
print(b)

dev.off()

stat_pvalue_hindcast <- hindcast_filter %>%
	filter(metric %in% c("RMSE.norm")) %>%
	group_by(pretty_group) %>%
	rstatix::tukey_hsd(score ~ model_name) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .2) %>%
	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))


converged_all_3 = hindcast_data_df %>%
	distinct(taxon, model_name) %>%
	select(taxon) %>% table
converged_all_3 = names(converged_all_3[converged_all_3 > 2])
ggplot(mean_crps_hindcast %>% filter(taxon %in% converged_all_3),
			 aes(x = model_name,y = mean_crps,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .1, show.legend = F) +

	xlab(NULL) +
	ylab("Hindcast error") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	labs(color = "Domain") + facet_grid(~pretty_group) + scale_x_discrete(labels= model.labs) +
	stat_compare_means(data=mean_crps_hindcast,
										 aes(x = model_name, y = mean_crps),
										 method = "anova", inherit.aes = F, size=5, label.y.npc = .6) +
	scale_y_log10() + # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", ref.group = "cycl_only", label.y.npc = .5,
										 show.legend = F, hide.ns = T, size=10)


# View model structure effects on mean CRPS scores
mean_crps_hindcast = hindcast_data_df %>%
	group_by(model_id, model_name, pretty_group, pretty_name, taxon) %>%
	summarize(mean_crps=mean(crps, na.rm=T))
ggplot(mean_crps_hindcast,
			 aes(x = model_name,y = mean_crps,
			 		color = model_name)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .1, show.legend = F) +

	xlab(NULL) +
	ylab("Hindcast error (CRPS)") +
	xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	labs(color = "Domain") + facet_grid(~pretty_group) + scale_x_discrete(labels= model.labs) +
	stat_compare_means(data=mean_crps_hindcast,
										 aes(x = model_name, y = mean_crps),
										 method = "anova", inherit.aes = F, size=5, label.y.npc = .6) +
	scale_y_log10() + # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", ref.group = "cycl_only", label.y.npc = .5,
										 show.legend = F, hide.ns = T, size=10)



hindcast_data_df %>%
	group_by(model_id, model_name, pretty_group, pretty_name) %>%
	summarize(mean_crps=mean(crps, na.rm=T))
ggplot(mean_crps_hindcast,
			 aes(x = model_name,y = mean_crps,
			 		color = model_name)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .1, show.legend = F) +
	geom_violin(draw_quantiles = c(.5)) +

	xlab(NULL) +
	ylab("Hindcast predictability (RSQ 1:1)") +
	xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	labs(color = "Domain") + facet_grid(pretty_name~pretty_group) + scale_x_discrete(labels= model.labs) +
	stat_compare_means(data=mean_crps_hindcast,
										 aes(x = model_name, y = mean_crps),
										 method = "anova", inherit.aes = F, size=5, label.y.npc = .6) +
	scale_y_log10() + # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", ref.group = "cycl_only", label.y.npc = .5,
										 show.legend = F, hide.ns = T, size=10)

skill_scores_long = skill_scores %>% pivot_longer(cols=c(skill_score, skill_score_random))

# Check average improvement by modeled vs random effects
skill_scores %>% group_by(pretty_group, model_name) %>% summarize(mean_skill_score=mean(skill_score, na.rm=T), mean_skill_score_random = mean(skill_score_random, na.rm=T))

ggplot(skill_scores_long,
			 aes(y = name,x = value)) +
	geom_point() +
#	geom_jitter(width=.2, height = 0, size=4, alpha = .1, show.legend = F) +
	geom_violin(draw_quantiles = c(.5)) +

	# ylab("Hindcast predictability (RSQ 1:1)") +
	# xlab("Taxonomic rank") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	labs(color = "Domain") + facet_grid(model_name~pretty_group) +
	scale_x_discrete(labels= model.labs) +
	stat_compare_means(data=skill_scores_long,
										 aes(x = value, y = name),
										 method = "anova", #inherit.aes = F,
										 size=5, label.y.npc = .6)
	# scale_y_log10() + # Add global p-value
	# stat_compare_means(aes(label = after_stat(p.signif)),
	# 									 method = "t.test", ref.group = "cycl_only", label.y.npc = .5,
	# 									 show.legend = F, hide.ns = T, size=10)

# Just checking which model performed best overall for each group
# Env_cycl has higher crps for some groups, but higher RSQ
cal_metrics = scores_list$calibration_metrics %>%
	filter(model_id %in% converged)
ggplot(cal_metrics,
			 aes(x=pretty_name, y=CRPS, colour = model_name)) +
	geom_boxplot(aes(colour = model_name)) +
	geom_point(position=position_jitterdodge(jitter.width = .1,jitter.height = 0, dodge.width = 1), alpha=.3, size=3) +
	theme_bw(base_size = 20) +
	ggtitle(paste0("In-sample prediction accuracy")) +
	xlab(NULL) + facet_grid(~pretty_group)

ggplot(cal_metrics,
			 aes(x=model_name, y=RSQ, colour = model_name)) +
	geom_boxplot(aes(colour = model_name)) +
	geom_point(position=position_jitterdodge(jitter.width = .1,jitter.height = 0, dodge.width = 1), alpha=.3, size=3) +
	theme_bw(base_size = 20) +
	xlab(NULL) + facet_grid(pretty_name~pretty_group)  +
	stat_compare_means(data = cal_metrics,
										 aes(x = model_name, y = RSQ, group=model_name),
										 method = "anova", inherit.aes = F, size=5, label.y.npc = .6) + # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", ref.group = "cycl_only")





ggscatter(hindcast_data_df, x = "mean", y = "truth",
					add = "reg.line", color = "model_name"
) + facet_grid(fcast_period~model_name) +
	xlab("Mean prediction estimate") + ylab("Observed plot mean") +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) + ggtitle("Overall prediction accuracy")


sum.all <- readRDS(here("data", "summary/predictor_effects.rds"))

df_cal <- sum.all %>%
	filter(beta %in% beta_names & #(workaround since the model names aren't all saving)
				 	model_name %in% c("env_cov","env_cycl") &
				 	!grepl("other", taxon) &
				 	time_period == "2015-11_2018-01" &
				 	!time_period %in% c("2015-11_2020-01")) %>% mutate(rank_only=only_rank)
df_cal$beta <- droplevels(df_cal$beta)

df_cal_fg_tax <- df_cal
# CHeck whether model structures (i.e. inclusion of cycl) lead to different effect sizes
ggplot(df_cal_fg_tax, aes(x = rank_only,y = effSize,
													color = model_name)) +
	geom_violin(position=position_dodge(), draw_quantiles = c(.5), show.legend = F) +
	geom_point(position=position_jitterdodge(jitter.width = .5,jitter.height = 0, dodge.width = 1),size=4, alpha = .5) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(rows = vars(beta), cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme(axis.text.x=element_text(
		angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold"))

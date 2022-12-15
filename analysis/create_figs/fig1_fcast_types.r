# Comparing functional groups vs taxonomy


source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
#p_load(multcomp)
library(ggallin)


scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))

sum.all <- readRDS("./data/summary/all_fcast_effects.rds")
sum.div = readRDS("./data/summary/div_summaries_effects.rds")
sum.div = sum.div %>% mutate(beta = ifelse(beta == "Ectomycorrhizal trees", "Ectomycorrhizal\ntrees", beta),
														 pretty_group = ifelse(taxon=="div_16S", "Bacteria", "Fungi"))

df_refit = rbindlist(list(sum.all, sum.div), fill=T) %>% filter(time_period=="2015-11_2020-01")
df_cal = rbindlist(list(sum.all, sum.div), fill=T) %>% filter(time_period=="2015-11_2018-01")


# Test for differences in effect sizes between fcast types
beta_tukey_group = list()
for (group in c("Fungi", "Bacteria")){
	beta_tukey = list()
	for (beta in c("Ectomycorrhizal\ntrees",
								 "LAI", "pC", "pH", "Temperature", "Moisture")) {
		df_group=df_refit %>% filter(pretty_group==!!group)
		out = df_group %>% filter(beta %in% !!beta) %>%
			#group_by(beta) %>%
			aov(effSize~fcast_type,.) %>%
			emmeans::emmeans(object = ., pairwise ~ "fcast_type", adjust = "tukey") %>% .$emmeans %>%
			multcomp::cld(object = ., Letters = letters) %>% as.data.frame() %>% rename(Letters_Tukey = `.group`) %>%
			#rownames_to_column("beta") %>%
			mutate(pretty_group = !!group, beta = !!beta)
		beta_tukey[[beta]] = out
	}
	beta_tukey_group[[group]] = do.call(rbind, beta_tukey)
}
tukey_refit_beta = do.call(rbind, beta_tukey_group)
tukey_refit_beta$tot = tukey_refit_beta$upper.CL + .1

# Subset to only predictors that have a difference between types
diff_letters <- tukey_refit_beta %>% group_by(beta) %>% #any(Letters_Tukey != "a")  %>%
	mutate(eq = replace(Letters_Tukey, n_distinct(Letters_Tukey)==1, '') ) #%>% filter(eq != "")
df_refit_diff = merge(df_refit, diff_letters, all.x = T) %>% filter(eq != "")

###### # Plot with Tukey, effect sizes across fcast types ----
c <- ggplot(df_refit_diff %>%
							filter(!beta %in% c("sin","cos","rho")),
						aes(x = fcast_type,y = effSize,color = pretty_group)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(rows = vars(beta),
						 cols = vars(pretty_group),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() +
	theme(
		text = element_text(size = 16),
		axis.text.x=element_text(
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18),
		strip.text.y = element_text(size=12)
	) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .1, show.legend = F)

tukey_refit_beta = df_refit %>% filter(!beta %in% c("rho","sin","cos")) %>%
	group_by(pretty_group, beta) %>%
	summarize(tukey(x = fcast_type, y = effSize)) %>%
	rename(fcast_type = x)



# beta_aov = df_refit %>% filter(!beta %in% c("rho","sin","cos")) %>% group_by(beta) %>% aov(effSize~pretty_group+fcast_type,.)
# aov1.1=aov(skill_score~pretty_group+fcast_type,skill_scores)
# beta_aov_tukey <- emmeans::emmeans(object = beta_aov,
# 																		pairwise ~ "fcast_type",
# 																		adjust = "tukey")
# tukey_skill_score <- multcomp::cld(object = mod_means_contr$emmeans,
# 																	 Letters = letters) %>% as.data.frame() %>% rename(Letters_Tukey = `.group`) %>% mutate(tot = 1)


tukey_beta_diff = merge(tukey_refit_beta,
												diff_letters[,c("fcast_type","pretty_group","beta","eq")], all.x = T) %>% filter(eq != "")

fcast_type_beta_plot <- c +
	geom_text(data = tukey_beta_diff,

						aes(x = fcast_type, y = tot, label = Letters_Tukey), show.legend = F, color = 1, size =6)
######
fcast_type_beta_plot



fcast_info_simple <- scores_list$scoring_metrics_site_lon %>%
	select(fcast_type, pretty_group, model_name, pretty_name, taxon) %>% distinct()


cal_rsq_filter <- scores_list$calibration_metrics %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1))  %>%
	filter(model_name == "all_covariates" & RSQ > .1)
saveRDS(cal_rsq_filter, here("data/summary/tax_filter_pass.rds"))


cal_rsq <- scores_list$calibration_metrics %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS, CRPS_truncated, CRPS_penalty, RSQ, RSQ.1,
												RMSE.norm,residual_variance, predictive_variance, total_PL),
																																		 names_to = "metric") %>%

	filter(model_name == "all_covariates") %>%
	distinct(.keep_all = T) %>% merge(fcast_info_simple, all.x=T, all.y=F)

hindcast_rsq <- scores_list$scoring_metrics_site_lon %>%
	mutate(value = ifelse(metric == "RSQ.1" & value < 0, 0, value)) %>%
	filter(model_name == "all_covariates" &
																																	metric %in% c("CRPS_penalty","RSQ","RSQ.1","CRPS","residual_variance", "predictive_variance", "total_PL","CRPS_truncated")) %>%
	distinct() #%>% filter(pretty_group == "Bacteria")


hindcast_rsq <- scores_list$scoring_metrics %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	filter(taxon %in% cal_rsq_filter$taxon | pretty_name == "Diversity") %>%
	filter(model_name == "all_covariates") %>%
	distinct()  %>%
	merge(fcast_info_simple, all.x=T, all.y=F)

hindcast_wide_site <- scores_list$scoring_metrics_site %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	filter(taxon %in% cal_rsq_filter$taxon | pretty_name == "Diversity") %>%
	filter(model_name == "all_covariates") %>%
	distinct()

calibration_wide_site <- scores_list$calibration_metrics_site %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity") %>%
	merge(fcast_info_simple, all.x=T, all.y=F)

cal_site_rsq_filter <- calibration_wide_site %>%
	filter(model_name == "all_covariates" & RSQ > .1) %>% mutate(filt = paste(taxon, siteID))


hindcast_site_filtered <- hindcast_wide_site %>%
	filter(paste(taxon, siteID) %in% cal_site_rsq_filter$filt | pretty_name == "Diversity")


tukey_fcast_CRPS = hindcast_site_filtered %>% group_by(pretty_group) %>% summarize(tukey(x = fcast_type, y = RSQ.1)) %>% rename(fcast_type = x)
tukey_fcast_CRPS = hindcast_wide_site %>% group_by(pretty_group) %>% summarize(tukey(x = fcast_type, y = RSQ.1)) %>% rename(fcast_type = x)
tukey_fcast_CRPS = hindcast_rsq %>% group_by(pretty_group) %>% summarize(tukey(x = fcast_type, y = RSQ)) %>% rename(fcast_type = x)
tukey_cal_CRPS = calibration_wide_site %>% group_by(pretty_group) %>% summarize(tukey(x = fcast_type, y = CRPS)) %>% rename(fcast_type = x)



# Test for differences in effect sizes between fcast types
score_tukey_group = list()
for (group in c("Fungi", "Bacteria")){
		df_group=hindcast_rsq %>% filter(pretty_group==!!group)
		out = df_group %>%
			#group_by(beta) %>%
			aov(RSQ.1~fcast_type,.) %>%
			emmeans::emmeans(object = ., pairwise ~ "fcast_type", adjust = "tukey") %>% .$emmeans %>%
			multcomp::cld(object = ., Letters = letters) %>% as.data.frame() %>%
			rename(Letters_Tukey = `.group`) %>%
			#rownames_to_column("beta") %>%
			mutate(pretty_group = !!group)
		score_tukey_group[[group]] = out
}
score_tukey_group
tukey_score = do.call(rbind, score_tukey_group)
tukey_score$tot = tukey_score$upper.CL + .3

# hindcast_taxon <- scores_list$scoring_metrics %>%
# 	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1)) %>%
# 	filter(taxon %in% cal_rsq_filter$taxon | pretty_name == "Diversity") %>%
# 	filter(model_name == "all_covariates") %>%
# 	distinct()  %>% merge(fcast_info_simple, all.x=T, all.y=F)

# Test for differences in hindcast accuracy

b <-
	ggplot(hindcast_rsq,
				 aes(x = fcast_type, y = RSQ.1,
				 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .4), alpha=.1, show.legend = F) +
	facet_grid(cols=vars(pretty_group), drop = T, scales="free") +
	# geom_jitter(aes(x = metric, y = value), width=.1,
	# 						height = 0, alpha = .8, size=4) +
	ylab("RSQ 1:1 (hindcast accuracy)") + xlab(NULL) +
	theme_bw(base_size=16) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16), legend.position = c(.9,1.1)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm"))
b

fcast_type_crps_plot <- b +
	geom_text(data = tukey_score,

						aes(x = fcast_type, y = tot, label = Letters_Tukey), show.legend = F, color = 1, size =6)
fcast_type_crps_plot

# ggplot(calibration_wide,
# 			 aes(x = RSQ.1, y = CRPS_truncated,
# 			 		color = pretty_group)) +
# 	geom_point(size = 4,  alpha=.1) + geom_abline(slope=1)  +
# 	theme_bw(base_size=18)

# ggplot(hindcast_wide,
# 			 aes(x = RSQ.1, y = CRPS_truncated,
# 			 		color = pretty_group)) +
# 	geom_point(size = 4,  alpha=.1) + geom_abline(slope=1)  +
# 	theme_bw(base_size=18)


# ggplot(cal_rsq %>% filter(metric %in% c("residual_variance", "predictive_variance", "total_PL")),
# 			 aes(x = pretty_name, y = value,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = 1), alpha=.1, show.legend = F) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	facet_grid(pretty_group ~ metric, drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Metric scores") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	ggtitle("calibration: Predictive loss increases (predictability worsens) at broader ranks")  +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))  + scale_y_log10()
#
# ggplot(cal_rsq %>% filter(metric %in% c("RSQ","RSQ.1","CRPS_truncated")),
# 			 aes(x = pretty_name, y = value,
# 			 		color = pretty_name)) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = 1), alpha=.1, show.legend = F) +
# 	facet_grid( metric ~ pretty_group, drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Metric scores") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	ggtitle("calibration: Rsq predictability increases at broader ranks")  +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))

skill_scores = scores_list$skill_score_taxon %>% filter(model_name == "all_covariates") %>%
	filter(taxon %in% cal_rsq_filter$taxon | pretty_name == "Diversity")
skill_scores$fcast_type = ifelse(skill_scores$fcast_type=="Diversity", "Evenness", skill_scores$fcast_type)

# View skill scores by forecast type
a <- ggplot(skill_scores,
			 aes(x = fcast_type, y = skill_score,
			 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = fcast_type, y = skill_score),
						 position = position_jitterdodge(jitter.height = 0, jitter.width = .4), alpha = .2, size=4) +
	ylab("Skill at new sites (% change in CRPS)") + xlab(NULL) +
	theme_bw(base_size=18) + #ggtitle("Model transferability to new sites") +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16))  +
	guides(color=guide_legend(title="Domain"),
				 shape=guide_legend(title="Forecast type")) +
	#facet_grid(cols=vars(pretty_group)) +
	scale_y_continuous(trans = pseudolog10_trans)

aov1.1=aov(skill_score~pretty_group+fcast_type,skill_scores)
mod_means_contr <- emmeans::emmeans(object = aov1.1,
																		pairwise ~ "fcast_type",
																		adjust = "tukey")
tukey_skill_score <- multcomp::cld(object = mod_means_contr$emmeans,
													 Letters = letters) %>% as.data.frame() %>% rename(Letters_Tukey = `.group`) %>% mutate(tot = 1)

skill_score_plot <- a +
	geom_text(data = tukey_skill_score,

						aes(x = fcast_type, y = tot, label = Letters_Tukey), show.legend = F, color = 1, size =6)
skill_score_plot

# ggplot(hindcast_rsq %>% filter(metric %in% c("residual_variance", "predictive_variance", "total_PL")),
# 			 aes(x = pretty_name, y = value,
# 			 		color = pretty_name)) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = 1), alpha=.1, show.legend = F) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	facet_grid(pretty_group ~ metric, drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Metric scores") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	ggtitle("hindcasts: Predictive loss increases (predictability worsens) at broader ranks")  +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm")) + scale_y_log10()
#
#
# ggplot(hindcast_rsq %>% filter(metric %in% c("RSQ","RSQ.1","CRPS_truncated")),
# 			 aes(x = pretty_name, y = value,
# 			 		color = pretty_name)) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = 1), alpha=.1, show.legend = F) +
# 	facet_grid( metric ~ pretty_group, drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Metric scores") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	ggtitle("hindcasts: Rsq predictability DECREASES at broader ranks")  +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))


fcast_type_beta_plot
fcast_type_crps_plot
skill_score_plot

library(ggpubr)
ggarrange(fcast_type_beta_plot, fcast_type_crps_plot, skill_score_plot, common.legend = T, nrow=1)

# Including C from seasonalEffSize script
ggarrange(fcast_type_beta_plot, fcast_type_crps_plot, skill_score_plot, amplitude_plot,
					common.legend = T, labels = c("A","B","C","D"))

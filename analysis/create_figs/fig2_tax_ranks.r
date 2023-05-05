source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(ggallin)

scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))

sum.all <- readRDS(here("data", paste0("summary/predictor_effects.rds")))
hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds"))

converged <- scores_list$converged_strict_list
converged <- scores_list$converged_list

hindcast_filter <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(model_name == "env_cycl" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm") &
				 	site_prediction == "New time (observed site)")
hindcast_filter$score = ifelse(hindcast_filter$metric=="RMSE.norm" & hindcast_filter$score > 5, 5, hindcast_filter$score)


calibration_filter <- scores_list$calibration_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(model_name == "env_cycl" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm"))

hindcast_crps_mean =  hindcast_data %>%
	filter(model_id %in% converged &
				 	fcast_period=="hindcast" & new_site==FALSE & model_name == "env_cycl") %>%
	group_by(model_name,pretty_group,pretty_name, rank_only, model_id) %>%
	summarize(score = mean(crps, na.rm=T), metric="mean_crps")

tukey_rmse_rank = hindcast_filter %>%
	filter(metric %in% c(#"CRPS_truncated","RSQ",
											 "RMSE.norm")) %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="RMSE.norm")

tukey_crps_rank = hindcast_crps_mean %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="mean_crps")


plotting_df_rank_scores = rbind(hindcast_crps_mean, hindcast_filter %>%
																	filter(metric %in% c("RMSE.norm")))
plotting_tukey =rbind(tukey_crps_rank, tukey_rmse_rank)


metric.labs <- c("Relative forecast error (nRMSE)", "Absolute forecast error (CRPS)")
names(metric.labs) <- c("RMSE.norm", "mean_crps")

a <- ggplot(plotting_df_rank_scores,
			 aes(x = pretty_name, y = score,
			 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .2), alpha=.2, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	facet_wrap(metric~pretty_group, drop = T, scales="free_y",
						 labeller = labeller(metric = metric.labs)) +
	ylab("Forecast error") + xlab(NULL) +
	theme_bw(base_size=20) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	geom_text(data = plotting_tukey,
						aes(x = pretty_name, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6) +
	scale_y_continuous(trans = pseudolog10_trans)



a <- tag_facet(a, size=7, tag_pool = LETTERS)

a



b <- ggplot(plotting_df_rank_scores,
						aes(x = pretty_group, y = score,
								color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitter(width = .2), alpha=.2, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	facet_wrap(nrow = 2, ~metric, drop = T, scales="free_y",
						 labeller = labeller(metric = metric.labs)) +
	ylab(NULL) + xlab(NULL) +
	theme_bw(base_size=20) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	scale_y_continuous(trans = pseudolog10_trans) +
	stat_compare_means(show.legend = F, label.x = 1.5,
										 label.y.npc = .6, size=6, aes(label = paste0("p = ", after_stat(p.format))))
b <- tag_facet(b, size=7, tag_pool = c("E","F"))


fig2 <- ggarrange(a, NULL, b, widths = c(2,-.1, 1), nrow=1)



png(here("figures","rank_scores.png"), width = 1600, height=1000)
print(fig2)

dev.off()
# Same thing but split onto two plots

# rank_rmse_plot <- ggplot(hindcast_filter %>%
# 												 	filter(metric %in% c("RMSE.norm")),
# 												 aes(x = pretty_name, y = score,
# 												 		color = pretty_group)) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
# 	facet_grid(rows=vars(pretty_group), drop = T, scales="free") +
# 	# geom_jitter(aes(x = metric, y = value), width=.1,
# 	# 						height = 0, alpha = .8, size=4) +
# 	ylab("Relative error (NRMSE)") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm")) +
# 	scale_y_continuous(trans = pseudolog10_trans)
# #+ ylim(c(0,5))
# rank_rmse_plot <- rank_rmse_plot +
# 	geom_text(data = tukey_rmse_rank,
# 						aes(x = pretty_name, y = tot, label = Letters_Tukey),
# 						show.legend = F, color = 1, size =6)
# rank_rmse_plot
#
#
#
# rank_crps_plot <- ggplot(hindcast_crps_mean,
# 												 aes(x = pretty_name, y = mean_crps,
# 												 		color = pretty_group)) +
# 	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
# 	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
# 	facet_grid(rows=vars(pretty_group), drop = T, scales="free") +
# 	ylab("Absolute error (CRPS)") + xlab(NULL) +
# 	theme_bw(base_size=18) +
# 	#scale_color_manual(values = c(1,2))	+
# 	theme(text = element_text(size = 22),
# 				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
# 				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
# 	#guides(color=guide_legend(title=NULL)) +
# 	geom_hline(yintercept = 0, linetype=2) +
# 	theme(plot.margin = margin(1,2,1,1, "cm"))
# rank_crps_plot <- rank_crps_plot +
# 	geom_text(data = tukey_crps_rank,
# 						aes(x = pretty_name, y = tot, label = Letters_Tukey),
# 						show.legend = F, color = 1, size =6)
# rank_crps_plot

#ggarrange(rank_rmse_plot, NULL, rank_crps_plot,  widths = c(1, -0.15, 1), nrow=1)

source("source.R")
library(ggallin)
library(deeptime) # devtools::install_github("willgearty/deeptime")
library(tagger) #devtools::install_github("eliocamp/tagger")



#scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr.rds")))

scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr2.rds")))
sum.all <- readRDS(here("data", paste0("summary/predictor_effects.rds")))
# Check if Parquet file exists, otherwise use RDS
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(parquet_file)) {
  cat("Using Parquet file for memory efficiency...\n")
  hindcast_in <- arrow::read_parquet(parquet_file)
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  hindcast_in <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}

converged <- scores_list$converged_strict_list
#converged <- readRDS("data/summary/stricter_converged_taxa_list.rds")
#converged <- scores_list$converged_list

hindcast_filter <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(#model_name == "env_cycl" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm","RMSE.iqr") &
				 	site_prediction == "New time (observed site)")
# hindcast_filter$score = ifelse(hindcast_filter$metric=="RMSE.norm" &
# 															 	hindcast_filter$score > 5, 5, hindcast_filter$score)
# 
# hindcast_filter$score = ifelse(hindcast_filter$metric=="RMSE.iqr" &
# 															 	hindcast_filter$score > 5, 5, hindcast_filter$score)



# scores_list$scoring_metrics_long %>%
# 	filter(model_id %in% converged) %>%
# 	filter(model_name == "env_cycl" &
# 				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm","RMSE.iqr") &
# 				 	site_prediction == "New time (observed site)") %>% ungroup %>%
# 	filter(metric %in% c("RMSE.iqr")) %>%
# 	group_by(pretty_group) %>%
# 	summarize(tukey(x = pretty_name, y = score))


calibration_filter <- scores_list$calibration_metrics_long %>%
	filter(model_id %in% converged) %>%
	filter(model_name == "cycl_only" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1","RMSE.norm"))

# RMSE calculated and normalized per group
rmse_values = hindcast_filter %>% ungroup %>%
	filter(metric %in% c(#"CRPS_truncated","RSQ",
		"RMSE.norm"))

score_values = hindcast_filter %>% ungroup %>%
	filter(metric %in% c("CRPS_truncated","RSQ","RSQ.1","MAE","RMSE.norm"))

rsq1_values = hindcast_filter %>% ungroup %>%
	filter(metric %in% c("RSQ.1"))


rsq_values = hindcast_filter %>% ungroup %>%
	filter(metric %in% c("RSQ"))


crps_values = hindcast_filter %>% ungroup %>%
	filter(metric %in% c("CRPS_truncated")) 
# CRPS calculated per-prediction
# crps_values =  hindcast_in %>%
# 	filter(model_id %in% converged &
# 				 	fcast_period=="hindcast" & new_site==FALSE &
# 				 	model_name == "cycl_only" & !is.na(crps)) %>%
# 	group_by(model_name,pretty_group,pretty_name, rank_only, model_id) %>%
# 	summarize(score = mean(crps, na.rm=T), metric="mean_crps") %>% ungroup()

# Split by f/b
tukey_rmse_rank = rmse_values %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="RMSE.norm")


tukey_rsq_rank = rsq_values %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="RSQ")


tukey_rmse_rank = rmse_values %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="RMSE.norm")

tukey_crps_rank = crps_values %>%
	group_by(pretty_group) %>%
	summarize(tukey(x = pretty_name, y = score)) %>%
	rename(pretty_name = x) %>%
	mutate(metric="mean_crps")

plotting_df_rank_scores = plyr::rbind.fill(rmse_values, crps_values)
plotting_tukey =rbind(tukey_crps_rank, tukey_rmse_rank)



metric.labs <- c("Relative forecast error (nRMSE)", "Absolute forecast error (CRPS)")
names(metric.labs) <- c("RMSE.norm", "mean_crps")

crps_rank_plot <- ggplot(plotting_df_rank_scores %>% filter(metric=="mean_crps"),
						aes(x = pretty_name, y = score,
								color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .2), alpha=.2, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	facet_wrap(~pretty_group, drop = T, scales="free_y",
						 labeller = labeller(metric = metric.labs)) +
	ylab("Absolute forecast error (CRPS)") + xlab(NULL) +
	theme_classic(base_size=20) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16), legend.position = c(.9,1.1)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	geom_text(data = plotting_tukey %>% filter(metric=="mean_crps"),
						aes(x = pretty_name, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6) +
	scale_y_continuous(trans = pseudolog10_trans)

rmse_rank_plot <- ggplot(plotting_df_rank_scores %>% filter(metric=="RMSE.norm"),
										aes(x = pretty_name, y = score)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, aes(
		color = pretty_group),
						 position = position_jitter(width = .2),
						 #position = position_jitterdodge(jitter.width = .2),
						 alpha=.5, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	facet_wrap(~pretty_group, drop = T, scales="free_y", nrow = 2) +
	# 					 labeller = labeller(metric = metric.labs)) +
	ylab("Relative forecast error (nRMSE)") + xlab(NULL) +
	theme_classic(base_size=16) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16), legend.position = c(.9,1.1)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	geom_text(data = tukey_rmse_rank,
						aes(x = pretty_name, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6) +
	scale_y_continuous(trans = pseudolog10_trans)


crps_fb_plot <- ggplot(plotting_df_rank_scores %>% filter(metric=="mean_crps"),
						aes(x = pretty_group, y = score,
								color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitter(width = .1), alpha=.6, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	# facet_wrap(nrow = 2, ~metric, drop = T, scales="free_y",
	# 					 labeller = labeller(metric = metric.labs)) +
  xlab(NULL) +
	theme_classic(base_size=16) +
	ylab("Absolute forecast error (CRPS)") + xlab(NULL) +
		#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	scale_y_continuous(trans = pseudolog10_trans) +
	stat_compare_means(show.legend = F, label.x = 1.5,
										 label.y.npc = .6, size=6, aes(label = paste0("p = ", after_stat(p.format))))

rmse_fb_plot <- ggplot(plotting_df_rank_scores %>% filter(metric=="RMSE.norm"),
											 aes(x = pretty_group, y = score,
											 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitter(width = .1), alpha=.6, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	# facet_wrap(nrow = 2, ~metric, drop = T, scales="free_y",
	# 					 labeller = labeller(metric = metric.labs)) +
	ylab("Relative forecast error (nRMSE)") + xlab(NULL) +
	theme_classic(base_size=20) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	scale_y_continuous(trans = pseudolog10_trans) +
	stat_compare_means(show.legend = F, label.x = 1.33,
										 label.y.npc = .8, size=5,
										 aes(label = paste0("p = ", after_stat(p.format))))


newsite_fb_plot <-
	ggplot(scores_list$skill_score_taxon 	%>% filter(model_id %in% converged #& model_name == "cycl_only"
	),# %>% filter(model_name=="cycl_only"),
											 aes(x = pretty_group, y = skill_score * 100,
											 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitter(width = .1), alpha=.6, show.legend = F) +
	#facet_grid(rows=vars(pretty_group), cols=vars(metric), drop = T, scales="free") +
	# facet_wrap(nrow = 2, ~model_name, drop = T, scales="free_y",
	# 					 labeller = labeller(metric = metric.labs)) +
		ylab("Transferability to new sites\n(% change in CRPS)") + xlab(NULL) +
		theme_classic(base_size=16) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 16),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=16), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	#scale_y_continuous(trans = pseudolog10_trans) +
	stat_compare_means(show.legend = F, label.x = 1.33,
										 label.y.npc = .6, size=5,
										 aes(label = paste0("p = ", after_stat(p.format))))


fcast_horizon_in = readRDS(here("data", "summary/fcast_horizon_clean.rds"))
fit_combined = fcast_horizon_in[[1]]
fcast_horizon_x_RSQ = fcast_horizon_in[[2]]
fcast_horizon_x_RMSE = fcast_horizon_in[[3]]
fcast_horizon_x_CRPS = fcast_horizon_in[[4]]


horizon_plot_f = ggplot(fcast_horizon_x_CRPS %>% 	filter(model_id %in% converged #& model_name == "cycl_only"
),# %>% filter(model_id %in% converged & model_name=="cycl_only"),
											 aes(x = pretty_group, y = months_since_obs, color=pretty_group)) +
	coord_flip() +
	geom_violin(draw_quantiles = c(0.5), show.legend = F) +
	geom_point(position=position_jitterdodge(jitter.width = .5, jitter.height = .5), alpha=.3, show.legend = F, size=2) +
	ylab("Forecast horizon\n(months since last observation)") + xlab(NULL) +
	stat_compare_means(show.legend = F, label.x = 1.5,
										 label.y.npc = .6, size=5, aes(label = paste0("p = ", after_stat(p.format)))) +
	theme_classic(base_size=16)  + theme(plot.margin = margin(0.1, 0.1, .1, 0.1, "cm"),
																			 axis.text.x = element_text(vjust = -1.6))  +
	tag_facets(tag_pool = c("f")) +
	theme(tagger.panel.tag.text = element_text(size = 16),
				tagger.panel.tag.background = element_rect(fill = "white", color = "white"))
 

panel_ac <- ggarrange2(
	# crps_rank_plot + rremove("x.text") +
	# 	theme(plot.margin = margin(0.1, 0.15, .1, 0.1, "cm")) +
	# 	tag_facets(tag_pool =	c("a", "b")) +
	# 	theme(tagger.panel.tag.text = element_text(size = 20),
	# 		tagger.panel.tag.background = element_rect(fill = "white", color = "white")),
	rmse_rank_plot + theme(plot.margin = margin(0.1, 0.15, .2, 0.1, "cm"),
												 panel.spacing.x = unit(2, "lines")) +
		tag_facets(tag_pool = c("a", "c")) +
		theme(tagger.panel.tag.text = element_text(size = 16),
			tagger.panel.tag.background = element_rect(fill = "white", color = "white")))

panel_bd <- ggarrange2(
	rmse_fb_plot +
		theme(plot.margin = margin(0.5, .1, .2, 0.1, "cm"),
					axis.text.y = element_text(vjust = -1.8)) + tag_facets(tag_pool = c("b")) +
		theme(tagger.panel.tag.text = element_text(size = 16),
			tagger.panel.tag.background = element_rect(fill = "white", color = "white")),
	newsite_fb_plot + theme(plot.margin = margin(0.1, 0.1, .2, 0.1, "cm"))  +
		tag_facets(tag_pool = c("d")) +
		theme(tagger.panel.tag.text = element_text(size = 16),
			tagger.panel.tag.background = element_rect(fill = "white", color = "white")),
	nrow = 2)


panel_ad <-
	ggarrange(
		panel_ac,panel_bd, nrow=1) 
#fig2 <- ggarrange2(panel_ad, panel_eg, nrow=1)
fig2 <- ggarrange(panel_ad, horizon_plot_f, nrow = 2, 		heights = c(3, 1))
fig2

png(here("figures","rank_scores.png"), width = 1200, height=1800)
print(fig2)
dev.off()




panel_a = rmse_rank_plot +
	theme(plot.margin = margin(0.2, 1, .2, 0.1, "cm"),
				axis.text.y = element_text(vjust = -1.8)) + tag_facets(tag_pool = c("(A")) +
	theme(tagger.panel.tag.text = element_text(size = 20),
				tagger.panel.tag.background = element_rect(fill = "white", color = "white"))

panel_abc <-
	ggarrange2(
		panel_a,rmse_fb_plot +
			theme(plot.margin = margin(0.5, .1, .2, 0.1, "cm"),
						axis.text.y = element_text(vjust = -1.8)) + tag_facets(tag_pool = c("(B")) +
			theme(tagger.panel.tag.text = element_text(size = 20),
						tagger.panel.tag.background = element_rect(fill = "white", color = "white")), nrow = 1, widths=c(2,1))

panel_de <- ggarrange2(horizon_plot_f + theme(plot.margin = margin(0.1, 1, .1, 0.1, "cm"),
																	axis.text.x = element_text(vjust = -1.6))  +
									 		tag_facets(tag_pool = c("(C")) +
									 		theme(tagger.panel.tag.text = element_text(size = 20),
									 					tagger.panel.tag.background = element_rect(fill = "white", color = "white")),
	newsite_fb_plot + theme(plot.margin = margin(0.1, 0.1, .2, 0.1, "cm"))  +
		tag_facets(tag_pool = c("(D")) +
		theme(tagger.panel.tag.text = element_text(size = 20),
					tagger.panel.tag.background = element_rect(fill = "white", color = "white")), nrow=1, widths=c(2,1))

fig2_v2 <-
	ggarrange2(
		panel_abc,
		panel_de,
		nrow = 2) #+
#	theme(plot.margin = unit(c(2,2,1,1), "cm"))
#fig2 <- ggarrange2(panel_ad, panel_eg, nrow=1)
fig2_v2

png(here("figures","figure_2ad.png"), width = 1200, height=1000)
print(fig2_v2)
dev.off()



metric.labs = c(RMSE.norm = "Relative forecast\nerror (nRMSE)",
								mean_crps = "Absolute forecast\nerror (CRPS)",
								CRPS_truncated = "Absolute forecast\nerror (CRPS)",
								RSQ="R-squared",
								RSQ.1 ="R-squared\nrelative to 1:1 line"
)


# For supplement
rank_models = ggplot(scores_list$scoring_metrics_long %>%
			 	filter(model_id %in% converged) %>%
			 	filter(metric %in% c("CRPS_truncated","RMSE.norm","RSQ","RSQ.1") &
			 				 	site_prediction == "New time (observed site)"),
			 aes(x = pretty_name, y = score,
			 		color = model_name)) +

	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .2), alpha=.3) +

	facet_grid(metric~pretty_group, drop = T, scales="free_y",
						 labeller = labeller(metric = metric.labs, model_name=model.labs)) +
	ylab("Forecast metric") + xlab(NULL) +
	theme_classic(base_size=20) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm")) +
	scale_y_continuous(trans = pseudolog10_trans) + labs(color="Model predictors") +
	scale_color_discrete(labels=model.labs) + theme(legend.position = "right")  +
	tag_facets() +
	theme(tagger.panel.tag.text = element_text(size = 20),
				tagger.panel.tag.background = element_rect(fill = "white", color = "white"))

rank_models
png(here("figures","supp_rank_model_scores.png"), width = 1000, height=2000)
print(rank_models)

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
# 	theme_classic(base_size=18) +
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
# 	theme_classic(base_size=18) +
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

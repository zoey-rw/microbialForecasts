
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))
scores_list = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/scoring_metrics_cv.rds")

sum.all <- readRDS("./data/summary/all_fcast_effects.rds")

hindcast_filter <- scores_list$scoring_metrics_long %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	filter(!taxon_model_rank %in% scores_list$unconverged$taxon_model_rank) %>%
	filter(fcast_type != "Diversity" &
				 	model_name == "all_covariates" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1") &
				 	site_prediction == "New time (observed site)")
calibration_filter <- scores_list$calibration_metrics_long %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	filter(!taxon_model_rank %in% scores_list$unconverged$taxon_model_rank) %>%
	filter(fcast_type != "Diversity" &
				 	model_name == "all_covariates" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1"))

# Taxonomic rank effect sizes
ggplot(df_refit %>% filter(!beta %in% c("sin","cos")),
			 aes(x = only_rank,y = effSize,
			 		color = pretty_group)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .4) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(rows = vars(beta), cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw(base_size = 18) + theme(axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		strip.text.y = element_text(size=12,face="bold")
	) +
	geom_smooth(method = "lm", aes(x = as.numeric(only_rank), y = effSize), show.legend = F) +
	stat_cor(aes(x = as.numeric(only_rank), y = effSize), color = 1, size =5, p.accuracy = .001)


all_cov_vals_scores_filter = all_cov_vals_scores %>% filter(taxon %in% pass_filter$taxon)

# Taxonomic rank seasonal amplitude sizes
ggplot(all_cov_vals_scores_filter %>% filter(fcast_type=="Taxonomic"),
																			aes(x = only_rank,y = amplitude,
													 color = pretty_group)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	#ggtitle("Residual seasonal amplitude") +
	xlab(NULL) +
	ylab("Residual seasonal amplitude") +
	facet_grid(cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	geom_smooth(method = "lm", aes(x = as.numeric(only_rank), y = amplitude), show.legend = F) +
	stat_cor(aes(x = as.numeric(only_rank), y = amplitude), color = 1, size =5, p.accuracy = .001)






ggplot(hindcast_filter %>% filter(pretty_name !="Functional group") %>%
			 	filter(metric %in% c("RSQ.1")),
			 aes(x = pretty_name,y = score,
			 		color = pretty_group)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	xlab(NULL) +
	ylab("Hindcast predictability (RSQ 1:1)") + xlab("Taxonomic rank") +

	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	geom_smooth(method = "lm", aes(x = as.numeric(pretty_name), y = score, color = pretty_group), show.legend = F) +
	stat_cor(aes(x = as.numeric(pretty_name), y = score, color = pretty_group), size =5, p.accuracy = .001) + labs(color = "Domain")




ggplot(scores_list$calibration_metrics_long %>%
			 	filter(metric %in% c("RSQ.1")),
			 aes(x = pretty_name,y = score,
			 		color = pretty_group)) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .4, show.legend = F) +
	geom_violin(draw_quantiles = c(.5)) +
	xlab(NULL) +
	ylab("Hindcast predictability (RSQ 1:1)") + xlab("Taxonomic rank") +

	theme_bw() + theme(text = element_text(size = 16),
										 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
										 	angle = 320, vjust=1, hjust = -0.05),
										 strip.text.y = element_text(size=12,face="bold")) +
	labs(color = "Domain")




rank_score_plot <- ggplot(hindcast_filter %>%
														filter(metric %in% c("RSQ.1")),
													aes(x = pretty_name, y = score,
															color = pretty_name)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
	facet_grid(rows=vars(pretty_group), drop = T, scales="free") +
	# geom_jitter(aes(x = metric, y = value), width=.1,
	# 						height = 0, alpha = .8, size=4) +
	ylab("Predictability (RSQ 1:1)") + xlab(NULL) +
	theme_bw(base_size=18) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm"))  +
	geom_smooth(method = "lm", aes(x = as.numeric(pretty_name), y = value), show.legend = F)


rank_score_plot <- rank_score_plot +
	geom_text(data = tukey_hindcast_rank,
						aes(x = pretty_name, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6)
rank_score_plot


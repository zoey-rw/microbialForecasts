
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))

cal_rsq <- scores_list$calibration_metrics_long %>% filter(model_name == "all_covariates" &
																													 	pretty_name != "Diversity" &
																															metric %in% c("RSQ","RSQ.1","CRPS","residual_variance", "predictive_variance", "total_PL","CRPS_truncated")) %>%
	mutate(value = ifelse(metric == "RSQ.1" &  value < 0, 0, value)) %>%
	distinct()

hist(cal_rsq[cal_rsq$metric=="RSQ",]$value, breaks=50)

# Check for natural break in predictability: looks like .25
p<-ggplot(cal_rsq %>% filter(metric == "RSQ"),
					aes(x=value)) +
	theme_bw(base_size = 18) +
	geom_histogram() + xlab("Calibration R-Squared") + ylab("Frequency") +
	geom_vline(xintercept=.25, linetype =2, color = 2) +
	ggtitle("Natural break in calibration predictability")
p

pass_filter <- cal_rsq %>% filter(metric == "RSQ" & value > .25) %>%
	select(pretty_group,model_name,pretty_name,taxon)
saveRDS(pass_filter, here("data/summary/tax_filter_pass.rds"))


hindcast_site = scores_list$scoring_metrics_site  %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	select(fcast_type, pretty_group, model_name, pretty_name, taxon, siteID,
				 CRPS_truncated, CRPS_penalty, RSQ, `RSQ.1`) %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity") %>%
	distinct() %>%
	pivot_longer(7:10, names_to = "metric")


hindcast_simple <- scores_list$scoring_metrics %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	select(pretty_group, model_name, pretty_name, taxon, CRPS_truncated, CRPS_penalty, RSQ, `RSQ.1`) %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity") %>%
	distinct() %>%
	pivot_longer(5:8, names_to = "metric")

hindcast_filter = hindcast_simple %>% merge(pass_filter, all.y=T)
hindcast_site_filter = hindcast_site %>% merge(pass_filter, all.y=T)


# Hindcast by taxon
hindcast_site_filter
# Hindcast by taxon x site
hindcast_filter

rsq_hindcast <- hindcast_site_filter %>% filter(metric == "RSQ.1" & pretty_group=="Fungi")
summary(lm(value~as.numeric(pretty_name), rsq_hindcast))
plot(value~as.numeric(pretty_name), rsq_hindcast)

# Test for differences in predictability between ranks
hindcast_tukey_group = list()
for (group in c("Fungi", "Bacteria")){
	df_group=hindcast_filter %>% filter(pretty_group==!!group)
	out = df_group %>%
		filter(metric == "RSQ.1") %>%
		#group_by(metric) %>%
		aov(value~pretty_name,.) %>%
		emmeans::emmeans(object = ., pairwise ~ "pretty_name", adjust = "tukey") %>% .$emmeans %>%
		multcomp::cld(object = ., Letters = letters) %>% as.data.frame() %>%
		rename(Letters_Tukey = `.group`) %>%
		#rownames_to_column("beta") %>%
		mutate(pretty_group = !!group)
	hindcast_tukey_group[[group]] = out
}
tukey_hindcast_rank = do.call(rbind, hindcast_tukey_group)
tukey_hindcast_rank$tot = tukey_hindcast_rank$upper.CL + .2
tukey_hindcast_rank





rank_score_plot <- ggplot(hindcast_filter %>%
														filter(metric %in% c(#"CRPS", "RSQ",
																								 "RSQ.1")),
			 aes(x = pretty_name, y = value,
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
	theme(plot.margin = margin(1,2,1,1, "cm"))

rank_score_plot <- rank_score_plot +
	geom_text(data = tukey_hindcast_rank,
						aes(x = pretty_name, y = tot, label = Letters_Tukey),
						show.legend = F, color = 1, size =6)
rank_score_plot




# Calibration scores
ggplot(cal_rsq %>%
			 	filter(metric %in% c("RSQ","RSQ.1")),
			 aes(x = pretty_name, y = value,
			 		color = pretty_name)) +
	geom_violin(draw_quantiles = c(0.5), show.legend=F) +
	geom_point(size = 4, position = position_jitterdodge(jitter.width = .5), alpha=.2, show.legend = F) +
	facet_grid(metric~pretty_group, drop = T, scales="free") +
	# geom_jitter(aes(x = metric, y = value), width=.1,
	# 						height = 0, alpha = .8, size=4) +
	ylab("Predictability") + xlab(NULL) +
	ggtitle("Predictability in calibration time period") +
	theme_bw(base_size=18) +
	#scale_color_manual(values = c(1,2))	+
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) +
	#guides(color=guide_legend(title=NULL)) +
	geom_hline(yintercept = 0, linetype=2) +
	theme(plot.margin = margin(1,2,1,1, "cm"))


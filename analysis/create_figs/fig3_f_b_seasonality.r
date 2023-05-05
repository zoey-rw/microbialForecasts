source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

library("cowplot")
library("ggExtra")

seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))


cycl_vals_scores = seasonal_amplitude_in[[6]]

scores_list = readRDS(here("data", "summary/scoring_metrics_cv.rds"))
converged = scores_list$converged_list
#converged = scores_list$converged_strict_list
hindcast_rsq <- scores_list$scoring_metrics %>%
	filter(model_id %in% converged & model_name == "env_cycl" &
	grepl("observed", site_prediction)) %>% ungroup %>%
	distinct() %>% select(-c(model_id, model_name))
hindcast_rsq$RMSE.norm = ifelse(hindcast_rsq$RMSE.norm > 5, 5, hindcast_rsq$RMSE.norm)

cycl_vals_scores <- merge(cycl_vals_scores, hindcast_rsq, all=F)

to_plot_seasonality <- cycl_vals_scores %>% filter(model_name=="cycl_only")
to_plot_seasonality_resid <- cycl_vals_scores %>% filter(model_name=="env_cycl")


hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds"))
calibration_abun  = hindcast_data %>%
	filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period!="hindcast")
mean_group_abun <- calibration_abun %>% group_by(model_id) %>% 
	summarize(mean_abun = mean(truth, na.rm=T),
						median_abun = median(truth, na.rm=T))
to_plot_seasonality_abun = merge(to_plot_seasonality, mean_group_abun)


ggplot(to_plot_seasonality_abun,
			 aes(x = mean_abun, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), label.x.npc = .1, label.y.npc = .8, size=6) +
	labs(color = "Kingdom")  +
	theme_bw(base_size = 18) +	ylab("Seasonal amplitude") + xlab("Mean abundance across samples") 


f1 <- ggplot(to_plot_seasonality,
						aes(x = RMSE.norm, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Seasonal amplitude") + xlab("Forecast error (nRMSE)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), label.x.npc = .1, label.y.npc = .8, size=6) +
	labs(color = "Kingdom")
f1


xplot <- ggdensity(to_plot_seasonality, "RMSE.norm", fill = "pretty_group", ) + clean_theme() + rremove("legend")
yplot <- ggdensity(to_plot_seasonality, "amplitude", fill = "pretty_group") +
#, color="pretty_group", add = "mean") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(f1)$grobs
legend_grob <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

sp <- f1 + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

# alternate plotting approach?
#xplot <- ggplot(to_plot_seasonality) + geom_density(aes(amplitude, fill = pretty_group)) + rotate() + theme_minimal() + rremove("legend")


fig3 <- plot_grid(xplot, legend_grob,
					sp, yplot, ncol = 2, align = "hv",
					rel_widths = c(2, 1), rel_heights = c(1, 2))
fig3 <- f1

png(here("figures","seasonality_error.png"), width = 800, height=800)
print(fig3)
dev.off()



g1 <- ggplot(to_plot_seasonality_resid,
						 aes(x = RMSE.norm, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Residual seasonal amplitude") + xlab("Forecast error (nRMSE)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Kingdom")

xplot <- ggdensity(to_plot_seasonality_resid, "RMSE.norm", fill = "pretty_group") + clean_theme() + rremove("legend")
yplot <- ggdensity(to_plot_seasonality_resid, "amplitude", fill = "pretty_group") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(g1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

sp <- g1 + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

plot_grid(xplot, legend, sp, yplot, ncol = 2, align = "hv",
					rel_widths = c(2, 1), rel_heights = c(1, 2))

#
# ggarrange(xplot, sp,
# 	yplot,
# 	nrow=2,ncol=2,
# 	common.legend = TRUE,
# )



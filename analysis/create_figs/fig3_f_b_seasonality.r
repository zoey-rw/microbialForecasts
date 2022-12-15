library(cowplot)
library("ggExtra")

to_plot_seasonality <- cycl_vals_scores %>% filter(model_name=="cycl_only" & !grepl("effect", site_prediction))
to_plot_seasonality_resid <- all_cov_vals_scores %>% filter(model_name=="all_covariates" & !grepl("effect", site_prediction))

# from seasonalityAccuracy.r
f <- ggplot(to_plot_seasonality,
			 aes(x = CRPS_truncated, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Seasonal amplitude") + xlab("Predictability (CRPS)") +
	stat_cor() + labs(color = "Domain")


f1 <- ggplot(to_plot_seasonality,
						aes(x = RSQ.1, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Seasonal amplitude") + xlab("Predictability (RSQ 1:1)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")

xplot <- ggdensity(all_cov_vals_scores, "RSQ.1", fill = "pretty_group") + clean_theme() + rremove("legend")
yplot <- ggdensity(all_cov_vals_scores, "amplitude", fill = "pretty_group") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(g)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

sp <- f + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

plot_grid(xplot, legend, sp, yplot, ncol = 2, align = "hv",
					rel_widths = c(2, 1), rel_heights = c(1, 2))





g1 <- ggplot(to_plot_seasonality_resid,
						 aes(x = RSQ.1, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Residual seasonal amplitude") + xlab("Predictability (RSQ 1:1)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")
g <-  ggplot(to_plot_seasonality_resid,
						 aes(x = CRPS_truncated, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Residual seasonal amplitude") + xlab("Predictability (CRPS)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")
	stat_cor() + labs(color = "Domain")

xplot <- ggdensity(to_plot_seasonality_resid, "RSQ.1", fill = "pretty_group") + clean_theme() + rremove("legend")
yplot <- ggdensity(to_plot_seasonality_resid, "amplitude", fill = "pretty_group") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(g)$grobs
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



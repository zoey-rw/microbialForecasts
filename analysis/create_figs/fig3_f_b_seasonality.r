source("source.R")

library("cowplot")
library("ggExtra")

seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))


cycl_vals_scores = seasonal_amplitude_in[[6]]

scores_list = readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list
#converged = scores_list$converged_strict_list
hindcast_rsq <- scores_list$scoring_metrics %>%
	filter(model_id %in% converged & model_name == "env_cycl" &
	grepl("observed", site_prediction)) %>% ungroup %>%
	distinct() %>% select(-c(model_id, model_name))
#hindcast_rsq$RMSE.norm = ifelse(hindcast_rsq$RMSE.norm > 5, 5, hindcast_rsq$RMSE.norm)

cycl_vals_scores <- merge(cycl_vals_scores, hindcast_rsq, all=F)

to_plot_seasonality <- cycl_vals_scores %>% filter(model_name=="cycl_only")
to_plot_seasonality_resid <- cycl_vals_scores %>% filter(model_name=="env_cycl")


# Check if Parquet file exists, otherwise use RDS
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(parquet_file)) {
  cat("Using Parquet file for memory efficiency...\n")
  hindcast_data <- arrow::read_parquet(parquet_file)
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  hindcast_data <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}
calibration_abun  = hindcast_data %>%
	filter(model_id %in% converged) %>%
	filter(!is.na(truth) & fcast_period!="hindcast")
mean_group_abun <- calibration_abun %>% group_by(model_id) %>%
	summarize(mean_abun = mean(truth, na.rm=T),
						median_abun = median(truth, na.rm=T))
to_plot_seasonality_abun = merge(to_plot_seasonality, mean_group_abun)
to_plot_seasonality_abun$nRMSE = to_plot_seasonality_abun$RMSE.norm
to_plot_long = to_plot_seasonality_abun %>% pivot_longer(cols=c(mean_abun, nRMSE), names_to = "metric")

metric.labs <- c("Relative forecast error (nRMSE)", "Mean abundance across samples")
names(metric.labs) <- c("nRMSE", "mean_abun")


f <- ggplot(to_plot_seasonality_abun,
			 aes(x = mean_abun, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01, show.legend = F) +
	scale_x_sqrt() +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F, show.legend = F)  +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), label.x.npc = .1, label.y.npc = .8, size=6, show.legend = F) +
	labs(color = "Kingdom")  +
	theme_bw(base_size = 18) +
	ylab("Seasonal amplitude") + xlab("Mean abundance across samples")
#f <- tag_facet(f, tag_pool = "A")


f1 <- ggplot(to_plot_seasonality_abun,
						aes(x = mean_abun, y = RMSE.norm, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	xlab("Mean abundance across samples") +	ylab("Forecast error (nRMSE)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")),
					 label.x.npc = .1, label.y.npc = .8, size=6) +
	labs(color = "Kingdom")
#f1 <- tag_facet(f1, tag_pool = "B", x = .01)


g1 <- ggplot(to_plot_seasonality_abun,
						 aes(x = RMSE.norm, y = amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01, show.legend = F) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F, show.legend = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	scale_y_sqrt() +
	ylab("Seasonal amplitude") + xlab("Forecast error (nRMSE)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), show.legend = F, size=6) +
	labs(color = "Kingdom")
# alternate plotting approach?
#xplot <- ggplot(to_plot_seasonality) + geom_density(aes(amplitude, fill = pretty_group)) + rotate() + theme_minimal() + rremove("legend")

library(deeptime)
# Skip early composite save; legend is built later and an alternate composite is saved


xplot <- ggdensity(to_plot_seasonality, "RMSE.norm", fill = "pretty_group", ) + clean_theme() + rremove("legend")
yplot <- ggdensity(to_plot_seasonality, "amplitude", fill = "pretty_group") +
	#, color="pretty_group", add = "mean") +
	rotate() + clean_theme() + rremove("legend")

# Create legend_grob from f1 plot first
grobs <- ggplotGrob(f1)$grobs
legend_indices <- which(sapply(grobs, function(x) x$name) == "guide-box")
if (length(legend_indices) > 0) {
  legend_grob <- grobs[[legend_indices[1]]]
} else {
  # No legend found, create a simple placeholder or skip legend
  cat("Warning: No legend found in plot, skipping legend extraction\n")
  legend_grob <- NULL
}

sp <- f1 + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")


fig3 <- plot_grid(xplot, legend,
									sp, yplot, ncol = 2, align = "hv",
									rel_widths = c(2, 1), rel_heights = c(1, 2))



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


no_kingdom_diff <- ggplot(to_plot_long,
													 aes(x = value, y = amplitude)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	facet_grid(~metric, scales="free", labeller = labeller(metric = metric.labs)) +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), label.x.npc = .3, label.y.npc = .9, size=6) +
	labs(color = "Kingdom")  +
	theme_bw(base_size = 18) +
	ylab("Seasonal amplitude") + scale_y_log10() + xlab(NULL)



no_kingdom_diff <- tag_facet(no_kingdom_diff, tag_pool = c("A","B"))


png(here("figures","seasonality_error2.png"), width = 800, height=800)
print(no_kingdom_diff)
dev.off()

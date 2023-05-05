source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(cowplot)
library(MuMIn)
options(na.action = "na.fail")

# Read in forecast scores
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
unconverged = scores_list$unconverged_list
converged = scores_list$converged_list
converged_strict = scores_list$converged_strict_list

# Read in Moran's I (spatial autocorrelation)
moran.stat_all_rank =
	readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/moran_stat.rds") %>%
	filter(taxon != "other")

# Read in model estimates for rho (temporal autocorrelation) and core_sd (spatial variation)
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(model_name != "all_covariates" & model_id %in% converged_strict) %>%
	select(-pretty_name)
core_sd = rho_core_in  %>% filter(rowname == "core_sd")
rho = rho_core_in  %>% filter(rowname == "rho")

# Merge the data that has been read in so far
model_scores_vals = merge(scores_list$scoring_metrics, rho_core_in)  %>%
	filter(model_id %in% converged_strict &
		grepl("observed", site_prediction))
model_scores_vals <- merge(model_scores_vals, moran.stat_all_rank, all.x = T)


# Read in model estimates of seasonal amplitude
seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_only_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="cycl_only") %>%
	mutate(cycl_amplitude=amplitude, cycl_max = max)
env_cycl_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="env_cycl")  %>%
	mutate(env_cycl_amplitude=amplitude, env_cycl_max = max)

model_scores_vals <- merge(model_scores_vals, cycl_only_vals_scores %>%
													 	select(taxon, time_period, cycl_amplitude, cycl_max))
model_scores_vals <- merge(model_scores_vals, env_cycl_vals_scores %>%
													 	select(taxon, time_period, env_cycl_amplitude, env_cycl_max))
model_scores_vals_wide = model_scores_vals %>% select(-c( SD, `Naive SE`, `Time-series SE`)) %>%
	pivot_wider(names_from = "rowname", values_from = "Mean") %>%
	filter(model_name=="env_cycl")

# Read in model estimates of environmental predictors
beta_in <- readRDS(here("data", "summary/predictor_effects.rds")) %>%
	filter(time_period=="2015-11_2018-01" & model_name=="env_cycl" &
				 	model_id %in% converged_strict) %>%
	pivot_wider(id_cols = c("model_name","pretty_group","taxon","time_period"),
							names_from = "beta", values_from = "effSize")
beta_wide <- beta_in[complete.cases(beta_in),]
model_scores_vals_wide_betas <- merge(model_scores_vals_wide, beta_wide)


# Check correlation between moran's and seasonality (negative)
cycl_vals = cycl_only_vals_scores %>% filter(model_id %in% converged_strict &
																						 	time_period == "2015-11_2018-01")
plotting_df = merge(cycl_vals, moran.stat_all_rank)
ggplot(model_scores_vals_wide_betas,
			 aes(x = mean_morans, y = cycl_amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Seasonal amplitude") + xlab("Variation among cores") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")

testing_df = merge(core_sd, moran.stat_all_rank) %>% filter(model_name=="cycl_only" &
																														model_id %in% converged &
																															time_period == "2015-11_2018-01")

a <- ggplot(model_scores_vals_wide,
						 aes(x = mean_morans, y = cycl_amplitude, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	ylab("Seasonal amplitude") + xlab("Variation among cores") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")

xplot <- ggdensity(model_scores_vals_wide, "core_sd", fill = "pretty_group") + clean_theme() + rremove("legend")
yplot <- ggdensity(model_scores_vals_wide, "cycl_amplitude", fill = "pretty_group") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(a)$grobs
legend_grob <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

sp <- a + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

plot_grid(xplot, legend_grob,
					sp, yplot, ncol = 2, align = "hv",
					rel_widths = c(2, 1), rel_heights = c(1, 2))

summary(lm(RSQ.1 ~ mean_morans * cycl_amplitude * pretty_group, model_scores_vals_wide))




b <- ggplot(model_scores_vals_wide,
						aes(x = mean_morans, y = RSQ.1, color = pretty_group)) +
	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
	geom_smooth(method="lm",
							linewidth=2, alpha = .2, na.rm = T, se = F)  +
	theme_bw(base_size = 18) +
	#facet_grid(metric ~pretty_group, scales="free") +
	scale_x_sqrt() +
	xlab("Seasonal amplitude") + ylab("Predictability") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
	labs(color = "Domain")

xplot <- ggdensity(model_scores_vals_wide, "core_sd", fill = "pretty_group") + clean_theme() + rremove("legend")
yplot <- ggdensity(model_scores_vals_wide, "cycl_amplitude", fill = "pretty_group") +
	rotate() + clean_theme() + rremove("legend")

grobs <- ggplotGrob(a)$grobs
legend_grob <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

sp <- b + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

plot_grid(xplot, legend_grob,
					sp, yplot, ncol = 2, align = "hv",
					rel_widths = c(2, 1), rel_heights = c(1, 2))

summary(lm(RSQ.1 ~ mean_morans * pretty_group, model_scores_vals_wide))

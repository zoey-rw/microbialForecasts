# View predictability by beta effects
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

predict_scores <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
cv_metric_scaled <- predict_scores$cv_metric_scaled


beta_effects <- readRDS(here("data/summary/all_fcast_effects.rds"))
beta_effects_min <- beta_effects %>% filter(taxon != "other") %>%
	select(siteID, model_name, group, rank_only, fcast_type, pretty_group, rank_name,
				 Mean, SD, beta, taxon, significant, effSize,time_period)
score_beta <- merge(beta_effects_min, cv_metric_scaled, allow.cartesian=TRUE)

ggplot(score_beta %>%
			 	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
			 				 	cv_type == "mean_per_plot_cv" & metric == "CRPS"),
			 aes(x = effSize, y = metric_scale,
			 		shape = pretty_group,
			 		color = pretty_name)) +
	geom_point(size = 4, alpha=.5, show.legend = T) +
	facet_wrap(~beta, drop = T, scales="free") + theme_bw(base_size = 18) +
	xlab("Parameter effect size") +
	ylab("CRPS (scaled by rank)") + ggtitle("Predictability ~ parameter effects: no clear trends") +
	guides(color=guide_legend(title="Rank"),
				 shape=guide_legend(title="Domain"))


# Model predictability as a function of env. sensitivities
library(MuMIn)
options(na.action = "na.fail")
to_model <- score_beta %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
				 	cv_type == "mean_per_plot_cv" & metric == "CRPS") %>% distinct()
to_model_wide <- to_model %>% pivot_wider(id_cols = c("model_name","pretty_group","pretty_name","taxon",
																											"metric_scale","CV_scale","time_period"),
																					names_from = "beta", values_from = "effSize")
to_model_wide <- to_model_wide[complete.cases(to_model_wide),]

to_model_wide <- to_model_wide %>% filter(time_period=="2015-11_2018-01")
predictability_lm <-lm(metric_scale ~ Temperature +
											 	Moisture + pH + pC +
											 	`Ectomycorrhizal\ntrees` + LAI +
											 	sin +
											 	cos,
											 data = to_model_wide)
models <- lapply(dredge(predictability_lm, evaluate = FALSE, m.lim = c(1,5)), eval)
ma <- model.avg(models,  subset = delta <= 2)
avg_results <- summary(ma)
avg <- as.data.frame(avg_results$coefmat.subset)

imp1 <- cbind.data.frame(predictor = names(ma$sw), values = ma$sw)
attr(imp1$values, "n.models") <- NULL
imp <- imp1 %>% mutate(taxon_rank = taxon_rank, taxon = !!s)
rownames(imp) <- NULL


to_model_wide <- to_model_wide %>% group_by(taxon) %>%
	mutate(temporal_env_sensitivity = sum(Temperature + Moisture + sin + cos + LAI + `Ectomycorrhizal\ntrees`),
				 spatial_env_sensitivity = sum(pH + pC),
				 env_sensitivity = temporal_env_sensitivity + spatial_env_sensitivity) %>%
	group_by(pretty_group, pretty_name) %>%
	mutate(env_sensitivity_scaled = scale(env_sensitivity)[,1],
				 temporal_env_sensitivity_scaled = scale(temporal_env_sensitivity)[,1],
				 quadrant = ifelse(env_sensitivity_scaled < 0 & metric_scale < 0, 4,
				 									 ifelse(env_sensitivity_scaled < 0 & metric_scale > 0, 3,
				 									 			 ifelse(env_sensitivity_scaled > 0 & metric_scale > 0, 2,
				 									 			 			 1))))


to_model_wide
master_plot <- ggplot(to_model_wide,
			 aes(x = metric_scale, y = env_sensitivity_scaled,
			 		shape = pretty_group,
			 		color = pretty_name)) +
	geom_point(#aes(size = CV_scale),
		size = 3,
						 alpha=.5, show.legend = T) +
	#facet_wrap(~beta, drop = T, scales="free") + theme_bw(base_size = 18) +
	ylab("Environmental sensitivity (scaled by rank)") +
	xlab("Predictability (scaled by rank)") +
	ggtitle("Predictability and parameter sensitivity") +
	guides(color=guide_legend(title="Rank"),
				 shape=guide_legend(title="Domain")) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal() + ylim(c(-3,3)) + xlim(c(-3.5,3.5))



for_examples <- to_model_wide %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == max(abs(env_sensitivity_scaled)) | metric_scale == max(abs(metric_scale))) %>% as.data.frame()

for_examples2 <- to_model_wide %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == min(env_sensitivity_scaled) | metric_scale == min(metric_scale)) %>% as.data.frame()

for_examples3 <- to_model_wide %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == min(env_sensitivity_scaled)) %>% as.data.frame()

# Sensitive, but unpredictable (Q1)
plot_model(hindcasts, taxon = "denitrification", siteID = "HARV", site_plots = "overlap")
to_model_wide[to_model_wide$env_sensitivity_scaled > 1 & to_model_wide$temporal_env_sensitivity_scaled > 1 & to_model_wide$quadrant == 1,] %>% as.data.frame()

q1_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "animal_pathogen", siteID = "CPER")

# Sensitive, and predictable (Q2)
to_model_wide[to_model_wide$env_sensitivity_scaled > 1 & to_model_wide$temporal_env_sensitivity_scaled > 1 & to_model_wide$quadrant == 2,] %>% as.data.frame()
plot_model(hindcasts, taxon = "copiotroph", siteID = "HARV", site_plots = "overlap")
q2_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "mycobacterium", siteID = "CPER")


# Insensitive, but predictable (Q3)
plot_model(hindcasts, taxon = "rhizobiales", siteID = "HARV", site_plots = "overlap")
plot_model(hindcasts, taxon = "actinobacteria", siteID = "HARV", site_plots = "overlap")
plot_model(hindcasts, taxon = "russula", siteID = "HARV", site_plots = "overlap")
to_model_wide[to_model_wide$env_sensitivity_scaled < -1 & to_model_wide$temporal_env_sensitivity_scaled < -1 & to_model_wide$quadrant == 3,] %>% as.data.frame()

q3_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "rhizobiales", siteID = "CPER")


# Insensitive, and unpredictable (Q4)
plot_model(hindcasts, taxon = "sphingomonas", siteID = "HARV", site_plots = "overlap")
plot_model(hindcasts, taxon = "tremellales", siteID = "HARV", site_plots = "overlap")
plot_model(hindcasts, taxon = "steroidobacter", siteID = "HARV", site_plots = "overlap")
plot_model(hindcasts, taxon = "rozellomycota", siteID = "HARV", site_plots = "overlap")
plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "rozellomycota", siteID = "CPER")
to_model_wide[to_model_wide$env_sensitivity_scaled < -1 & to_model_wide$temporal_env_sensitivity_scaled < -1 & to_model_wide$quadrant == 4,] %>% as.data.frame()
q4_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "sphingomonas", siteID = "CPER")


master_out <- master_plot +
	inset_element(q1_inset, 0, 0.6, 0.4, 1, align_to = 'plot') +
	inset_element(q2_inset, 0.6, 0.6, 1, 1, align_to = 'plot') +
	inset_element(q3_inset, 0.7, 0, 1, .7, align_to = 'plot') +
	inset_element(q4_inset, 0, 0, .7, .7, align_to = 'plot')
	# inset_element(
	# 	q1_inset,
	# 	left = 0.5,
	# 	bottom = 0.5,
	# 	right = 1,
	# 	top=1,
	# 	align_to = "full"
	# )
master_out

patch1 <- q1_inset / q4_inset
patch2 <- q2_inset / q3_inset
patch1 | master_plot | patch2

master_plot

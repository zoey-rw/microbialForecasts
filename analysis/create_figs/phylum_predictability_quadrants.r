# View predictability by beta effects
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
library(ggrepel)

predict_scores <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
cv_metric_scaled <- predict_scores$cv_metric_scaled


tax_scores <- predict_scores$scoring_metrics_cv %>% filter(model_name == "all_covariates") %>%
	filter(metric %in% c("CRPS_truncated","RSQ","RSQ.1"))
	# group_by(metric) %>%
	# mutate(CRPS_scale = scale(CRPS_tr)[,1],
	# 			 predictability = -CRPS_scale)

ggplot(tax_scores, aes(x = score, y = per_site_cv, color = pretty_group)) + geom_point() + geom_smooth() + facet_grid(~metric)

ggplot(tax_scores, aes(x = pretty_group, y = per_site_cv, color = pretty_group)) + geom_point()


beta_effects <- readRDS(here("data/summary/all_fcast_effects.rds"))
beta_effects_min <- beta_effects %>% filter(!grepl("other", taxon)) %>%
	select(siteID, model_name, group, rank_only, pretty_group, rank_name,
				 Mean, SD, beta, taxon, significant, effSize,time_period)
score_beta <- merge(beta_effects_min, tax_scores, allow.cartesian=TRUE)

to_model <- score_beta %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity") %>% distinct()
to_model_wide <- to_model %>% pivot_wider(id_cols = c("model_name","pretty_group","pretty_name","taxon",
																											"time_period","predictability"),
																					names_from = "beta", values_from = "effSize")

to_model_wide <- to_model_wide %>% filter(time_period=="2015-11_2018-01")
to_model_phylum <- to_model_wide %>% group_by(taxon) %>%
	filter(pretty_name == "Phylum") %>%
	mutate(temporal_env_sensitivity = sum(Temperature + Moisture + sin + cos + LAI + `Ectomycorrhizal\ntrees`),
				 spatial_env_sensitivity = sum(pH + pC),
				 env_sensitivity = temporal_env_sensitivity + spatial_env_sensitivity) %>%
	group_by(pretty_group, pretty_name) %>%
	mutate(env_sensitivity_scaled = scale(env_sensitivity)[,1],
				 temporal_env_sensitivity_scaled = scale(temporal_env_sensitivity)[,1],
				 quadrant = ifelse(env_sensitivity_scaled < 0 & predictability < 0, 4,
				 									ifelse(env_sensitivity_scaled < 0 & predictability > 0, 3,
				 												 ifelse(env_sensitivity_scaled > 0 & predictability > 0, 2,
				 												 			 1))))


to_model_wide
master_plot <- ggplot(to_model_phylum,
											aes(x = predictability, y = env_sensitivity_scaled,
													shape = pretty_group,
													color = pretty_group)) +
	geom_point(#aes(size = CV_scale),
		size = 3,
		alpha=.5, show.legend = F) +
	#facet_wrap(~beta, drop = T, scales="free") +
	theme_bw(base_size = 22) +
	ylab("Environmental sensitivity") +
	xlab("Predictability") +
	ggtitle("Predictability and parameter sensitivity") +
	guides(color=guide_legend(title="Domain")) +
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal() +
	ylim(c(-2,2)) +
	xlim(c(-3,3)) +
	geom_label_repel(
									 aes(label = taxon), size = 8,
									 nudge_x = .1,
									 na.rm = TRUE) + theme(axis.title = element_text(size=26))  +
	annotate("text", x = 2, y = 2, label = "More predictable, \nmore sensitive to environment", size = 6) +
	annotate("text", x = -2, y = -2, label = "Less predictable, \nless sensitive to environment", size = 6)
master_plot

to_model_phylum[to_model_phylum$quadrant==1,] %>% as.data.frame()
to_model_phylum[to_model_phylum$quadrant==2,] %>% as.data.frame()
to_model_phylum[to_model_phylum$quadrant==3,] %>% as.data.frame()
to_model_phylum[to_model_phylum$quadrant==4,] %>% as.data.frame()

for_examples <- to_model_phylum %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == max(abs(env_sensitivity_scaled)) | predictability == max(abs(predictability))) %>% as.data.frame()

for_examples2 <- to_model_phylum %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == min(env_sensitivity_scaled) | predictability == min(predictability)) %>% as.data.frame()

for_examples3 <- to_model_phylum %>% group_by(quadrant) %>%
	filter(env_sensitivity_scaled == min(env_sensitivity_scaled)) %>% as.data.frame()

# Sensitive, but unpredictable (Q1)
to_model_phylum[to_model_phylum$env_sensitivity_scaled > 1 & to_model_phylum$temporal_env_sensitivity_scaled > 1 & to_model_phylum$quadrant == 1,] %>% as.data.frame()

#q1_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "acidobacteriota", siteID = "CPER")
q1_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "verrucomicrobiota", siteID = "CPER")
q1_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "ascomycota", siteID = "CPER")

# Sensitive, and predictable (Q2)
to_model_wide[to_model_wide$env_sensitivity_scaled > 1 & to_model_wide$temporal_env_sensitivity_scaled > 1 & to_model_wide$quadrant == 2,] %>% as.data.frame()
plot_model(hindcasts, taxon = "copiotroph", siteID = "HARV", site_plots = "overlap")
q2_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "mortierellomycota", siteID = "CPER")
q2_inset + ggtitle("Mortierellomycota at Harvard Forest plots")

# Insensitive, but predictable (Q3)
to_model_wide[to_model_wide$env_sensitivity_scaled < -1 & to_model_wide$temporal_env_sensitivity_scaled < -1 & to_model_wide$quadrant == 3,] %>% as.data.frame()
#q3_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "basidiomycota", siteID = "CPER")
q3_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "chloroflexi", siteID = "CPER")


# Insensitive, and unpredictable (Q4)
plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "rozellomycota", siteID = "CPER")
to_model_wide[to_model_wide$env_sensitivity_scaled < -1 & to_model_wide$temporal_env_sensitivity_scaled < -1 & to_model_wide$quadrant == 4,] %>% as.data.frame()
#q4_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "rozellomycota", siteID = "CPER")
q4_inset <- plot_summary_inset(hindcasts, beta_df = beta_df, taxon = "actinobacteriota", siteID = "CPER")
q4_inset <- q4_inset + ylim(c(0,1))

patch1 <- q1_inset / q4_inset
patch2 <- q2_inset / q3_inset
plot_out <- patch1 | master_plot | patch2
plot_out

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(ggallin)


weak_keep_list <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
strict_keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))

converged = strict_keep_list

rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(model_name != "all_covariates" & model_id %in% converged) %>%
	select(-pretty_name)

core_sd = rho_core_in  %>%
	filter(rowname == "core_sd")
rho = rho_core_in  %>%
	filter(rowname == "rho")

rho$hi <- rho$Mean + rho$SD*1.96
rho$lo <- rho$Mean - rho$SD*1.96
rho$significant <- microbialForecast:::is_significant(rho$lo, rho$hi)



ggplot(core_sd %>% filter(model_name == "cycl_only")) +
	geom_boxplot(aes(x = only_rank,
									 y = Mean, color = pretty_group),
							 position=position_dodge()) +
	geom_point(aes(x = rank_only,
								 y = Mean, color = pretty_group),
						 position=position_jitterdodge(), alpha=.3) +
	facet_grid(~pretty_group) + theme_minimal(base_size = 18) + ggtitle("Variation among soil cores")

# Supplementary figure
ggplot(rho %>% filter(model_name == "cycl_only"), 
			 aes(x = only_rank,
																											y = Mean, color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F, color = 1) +
	# geom_boxplot(show.legend = F) +
	# 						 position=position_dodge()) +
	geom_point(aes(shape=as.factor(significant)), #show.legend = F, 
						 size=3,
						 alpha=.3, position=position_jitter(height=0)) +
	 facet_grid( ~ pretty_group,  scales="free") +
	# 					 labeller = labeller(model_name = model.labs)) + 
						 #rows=vars(pretty_group), scales="free") + 
	xlab("Rank") +
	theme_minimal(base_size = 18)  +
	ggtitle("Temporal autocorrelation effect size")  + 
	scale_y_continuous(trans=ssqrt_trans) + 
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05))  + 
	scale_shape_manual(values = c(21, 16), name = NULL,																																						 labels = c("Not significant","Significant")) + guides(color="none")


scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
unconverged = scores_list$unconverged_list
converged = scores_list$converged_list
converged_strict = scores_list$converged_strict_list

model_scores_vals = merge(scores_list$scoring_metrics, rho_core_in)  %>%
	filter(#model_id %in% converged &
				 	grepl("observed", site_prediction))


# Get mean values for plotting
model_scores_vals_env_cycl = model_scores_vals %>%
	filter(model_name=="env_cycl")



ggplot(model_scores_vals_env_cycl %>% filter(rowname=="core_sd"),
			 aes(x = RSQ,
															y = Mean)) +
	geom_point(aes( color = pretty_group), alpha=.5) +
	geom_smooth() #+
	#facet_grid(~pretty_group, scales="free")

model_scores_vals_wide = model_scores_vals_env_cycl %>% select(-c( SD, `Naive SE`, `Time-series SE`)) %>%
	pivot_wider(names_from = "rowname", values_from = "Mean")
summary(lm(RSQ ~ core_sd * rho, model_scores_vals_wide))





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
	pivot_wider(names_from = "rowname", values_from = "Mean")


summary(lm(RSQ ~ core_sd * cycl_amplitude, model_scores_vals_wide %>% filter(model_name=="env_cycl")))
summary(lm(RSQ.1 ~ core_sd * cycl_amplitude, model_scores_vals_wide %>% filter(model_name=="env_cycl")))


ggplot(model_scores_vals_wide %>% filter(model_name=="cycl_only"),
			 aes(x = core_sd,
			 		y = cycl_amplitude,  color = pretty_group)) +
	geom_point(aes( color = pretty_group), alpha=.5) +
	geom_smooth()

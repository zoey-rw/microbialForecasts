source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(cowplot)
library(MuMIn)
options(na.action = "na.fail")

# Read in forecast scores
scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
unconverged = scores_list$unconverged_list
converged = scores_list$converged_list
converged  = scores_list$converged_strict_list

# Read in Moran's I (spatial autocorrelation)
moran.stat_all_rank =
	readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/moran_stat.rds") %>%
	filter(taxon != "other")

# Read in model estimates for rho (temporal autocorrelation) and core_sd (spatial variation)
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(model_name != "all_covariates" & model_id %in% converged) %>%
	select(-pretty_name)
core_sd = rho_core_in  %>% filter(rowname == "core_sd")
rho = rho_core_in  %>% filter(rowname == "rho")

# Merge the data that has been read in so far
model_scores_vals = merge(scores_list$scoring_metrics, rho_core_in)  %>%
	filter(model_id %in% converged &
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
				 	model_id %in% converged) %>%
	pivot_wider(id_cols = c("model_name","pretty_group","taxon","time_period"),
							names_from = "beta", values_from = "effSize")
beta_wide <- beta_in[complete.cases(beta_in),]
model_scores_vals_wide_betas <- merge(model_scores_vals_wide, beta_wide)




# Model predictability as a function of env. sensitivities
to_model <- model_scores_vals_wide_betas %>% mutate(`temporal autocorrelation` = abs(rho),
																										seasonality = cycl_amplitude,
																										`variation among soil cores` = core_sd,
																										`Spatial autocorrelation` = mean_morans)
to_model_bacteria =to_model %>% filter(pretty_group=="Bacteria")
predictability_lm_bac <-lm(RSQ ~ Temperature +
													 	Moisture + pH + pC +
													 	`Ectomycorrhizal\ntrees` + LAI + 	`Spatial autocorrelation` +
													 	`temporal autocorrelation` + seasonality,
													 data = to_model_bacteria)
models_bac <- lapply(dredge(predictability_lm_bac, evaluate = FALSE, m.lim = c(1,9)), eval)
ma_bac <- model.avg(models_bac,  subset = delta <= 2)

plot(ma_bac, full = F, main = "Effects on bacterial predictability (RSQ)", intercept = F)


ma_bac_avg_results <- summary(ma_bac)
avg_bac <- as.data.frame(ma_bac_avg_results$coefmat.subset)
bac_importance <- cbind.data.frame(predictor = names(ma_bac_avg_results$sw), values = ma_bac_avg_results$sw) %>% mutate(Kingdom = "Bacteria")


to_model_fungi =to_model %>% filter(pretty_group=="Fungi")
predictability_lm_fun <-lm(RSQ ~ Temperature +
													 	Moisture + pH + pC +
													 	`Ectomycorrhizal\ntrees` + LAI + 	`Spatial autocorrelation` +
													 	`temporal autocorrelation` + seasonality,
													 data = to_model_fungi)
models_fun <- lapply(dredge(predictability_lm_fun, evaluate = FALSE, m.lim = c(1,9)), eval)
ma_fun <- model.avg(models_fun,  subset = delta <= 2)
plot(ma_fun, full = F, main = "Effects on fungal predictability (RSQ)", intercept = F)

ma_fun_avg_results <- summary(ma_fun)
avg_fun <- as.data.frame(ma_fun_avg_results$coefmat.subset)
fun_importance <- cbind.data.frame(predictor = names(ma_fun_avg_results$sw),
																	 values = ma_fun_avg_results$sw) %>% mutate(Kingdom = "Fungi")



all_importance <- rbind(fun_importance, bac_importance)
importance_f_b = ggplot(all_importance,
												aes(x = reorder(predictor, -values), y = values, color = Kingdom)) +
	geom_point(position = position_jitterdodge(jitter.width = .1, jitter.height = 0),
						 size=3, alpha=.8, show.legend = F)  +
	theme_minimal(base_size = 18) +
	geom_hline(yintercept = 0) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	xlab("Variables explaining predictability across groups") +
	ylab("Variable importance") +
	scale_color_manual(values = c(#"Fungi" = "#00BA38", # green
		"Fungi" = "#00BFC4",
		"Bacteria" = "#F8766D",
		"Overall" = "black"))




ph_pred <- ggplot(to_model, aes(x = pH,
																y = RSQ, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("pH (absolute effect)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6)

patchiness_pred <- ggplot(to_model, aes(x = `Spatial autocorrelation`,
																				y = RSQ, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Spatial autocorrelation") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6)

multiplot = ggarrange(importance_f_b, ggarrange(ph_pred, patchiness_pred, ncol=1,common.legend = T, labels = c("B","C")), ncol=2, common.legend = T, labels = c("A", NULL
))
multiplot

png(here("figures","explain_predictability.png"), width = 1200, height=800)
print(multiplot)
dev.off()









ph_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = pH,
																												 y = mean_crps_sample, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("pH (absolute effect)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")




pc_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = pC,
																												 y = mean_crps_sample, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("pH (absolute effect)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")



core_sd_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = core_sd,
																															y = mean_crps_sample, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Variation among soil cores") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")


ggarrange(ph_pred_crps, pc_pred_crps, core_sd_pred_crps, common.legend = T)


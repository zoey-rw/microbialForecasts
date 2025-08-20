source("source.R")
library(cowplot)
library(MuMIn)
options(na.action = "na.fail")

# Read in forecast scores
scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list
converged  = scores_list$converged_strict_list

# Read in Moran's I (spatial autocorrelation)
# moran.stat_all_rank data not available, skipping this analysis
# moran.stat_all_rank =
# 	readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/moran_stat.rds") %>%
# 	filter(taxon != "other")
moran.stat_all_rank = data.frame() # Empty dataframe to avoid errors

# Read in model estimates for rho (temporal autocorrelation) and core_sd (spatial variation)
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(time_period=="2015-11_2018-01") %>%
	filter(model_name != "all_covariates" & model_id %in% converged) %>%
	select(model_id, taxon, rowname, Mean) %>%
	pivot_wider(names_from = "rowname", values_from = "Mean")


# Merge the data that has been read in so far
scores_df = scores_list$scoring_metrics %>%
	filter(model_id %in% converged &
				 	site_prediction != "New time x site (random effect)")

scores_site_df = scores_list$scoring_metrics_site %>%
	filter(model_id %in% converged &
				 	site_prediction != "New time x site (random effect)")



# Read in model estimates of seasonal amplitude
seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
seas_amplitude_long <- seasonal_amplitude_in[[6]] %>% filter(time_period=="2015-11_2018-01")  %>%
	select(model_id, taxon, model_name, amplitude)

# Read in model estimates of environmental predictors
beta_in <- readRDS(here("data", "summary/predictor_effects.rds")) %>%
	filter(time_period=="2015-11_2018-01" & model_id %in% converged) %>%
	select(model_id, beta, effSize) %>%
	pivot_wider(id_cols = c("model_id"),
							names_from = "beta", values_from = "effSize")
# beta_wide <- beta_in[complete.cases(beta_in),]


master_df <- merge(scores_df, rho_core_in)


#master_df <- merge(scores_site_df, rho_core_in)
master_df <- merge(master_df, moran.stat_all_rank, all.x = T)
master_df <- merge(master_df, seas_amplitude_long, all.x = T)
master_df <- merge(master_df, beta_in, all.x = T)



# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))
site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups
site_descr$latitude_category = ifelse(site_descr$latitude > 44, "High-latitude",
																			ifelse(site_descr$latitude < 31, "Low-latitude",
																						 "Mid-latitude"))
# This factor should be ordered
site_descr$latitude_category =
	factor(site_descr$latitude_category, ordered = T, levels = c("Low-latitude",  "Mid-latitude","High-latitude"))


#master_df <- merge(master_df, site_descr, all.x = T)

# Model predictability as a function of env. sensitivities
to_model <- master_df %>% mutate(#`temporal autocorrelation` = abs(rho),
																 `temporal memory` = rho,
																										seasonality = amplitude,
																										`variation among soil cores` = core_sd
																										# `spatial autocorrelation` = mean_morans  # Commented out - data not available
																										) %>%
	filter(!is.na(RSQ))

to_model$seasonal_predictors = ifelse(to_model$model_name %in% c("env_cycl","cycl_only"), 1, 0)
to_model$environmental_predictors = ifelse(to_model$model_name %in% c("env_cycl","env_cov"), 1, 0)
to_model$new_site = ifelse(to_model$site_prediction =="New time x site (random effect)", 1, 0)
to_model$Fungi = ifelse(to_model$pretty_group =="Fungi", 1, 0)
to_model$Taxonomic = ifelse(to_model$fcast_type =="Taxonomic", 1, 0)


to_model_scaled = to_model %>% group_by(pretty_group, model_name) %>%
	filter(model_name == "env_cycl" & site_prediction == "New time (observed site)") %>%
	mutate(#latitude = scale(latitude),
																			RMSE.norm = scale(RMSE.norm),
																			mean_crps_sample = scale(mean_crps_sample),
																			RSQ = scale(RSQ),
																			`seasonal amplitude` = scale(seasonality),
																			`temporal memory` = scale(`temporal memory`),
																			# mean_morans = scale(mean_morans),  # Commented out - data not available
																			amplitude = scale(amplitude),
																			temperature = scale(Temperature),
																			moisture = scale(Moisture),
																			pH = scale(pH),
																			`percent carbon` = scale(pC),
																		`ectomycorrhizal trees` = scale(`Ectomycorrhizal\ntrees`),
LAI = scale(LAI)
)
env_df = to_model_scaled #%>% filter(model_name == "env_cycl" & site_prediction == "New time (observed site)")


library(sjlabelled)
set_label(env_df$Taxonomic) <- "Group type: Taxonomic"
set_label(env_df$Fungi) <- "Microbial kingdom: Fungi"
set_label(env_df$seasonal_predictors) <- "Seasonal predictors"
set_label(env_df$environmental_predictors) <- "Environmental predictors"
set_label(env_df$new_site) <- "New location"
set_label(env_df$RMSE.norm) <- "Effects on error (nRMSE)"
set_label(env_df$amplitude) <- "Seasonal amplitude"
set_label(env_df$`ectomycorrhizal trees`) <- "Ectomycorrhizal trees"
set_label(env_df$LAI) <- "Leaf area index"
set_label(env_df$pC) <- "percent carbon"


#
# to_model =to_model %>% filter(site_prediction != "New time x site (modeled effect) " &
# 																!is.na(RSQ)) #&
# 															#	pretty_name %in% c("phylum","functional")
# 															)
	#pretty_group=="Bacteria")

# Check if we have enough data for modeling
if (nrow(env_df) > 0 && sum(!is.na(env_df$RMSE.norm)) > 0) {
  overall_lm <-lm(RMSE.norm ~
  							temperature +
  							moisture + pH + `percent carbon` +
  							`ectomycorrhizal trees` + LAI +
  							`temporal memory` + `seasonal amplitude` ,
  						data = env_df)
  
  overall_plot <- sjPlot::plot_models(model = list(overall_lm),
  																		vline.color = "black", show.legend = F, colors = "gs") +
  	theme_bw(base_size = 20) + ylab("Parameter effect on forecast error (nRMSE)") +
  	ylim(c(-1.5,1.5)) +
  	theme(plot.margin = unit(c(1,1,1,0), "cm")) +
  	annotate("rect", xmin = 8.6, xmax = 9.7, ymin = -1.5, ymax = 1.5,
  					 color="black", fill = "lightgrey")+
  	annotate(label = "Environmental predictors\n& seasonality", x= 9.2, y=0, geom="text", size=5.5) +  coord_flip(expand = FALSE)
  print(overall_plot)
} else {
  cat("Insufficient data for overall modeling\n")
  overall_plot <- NULL
}



# create fig4a using analysis/misc/view_insensitive_taxa.r
fig4 <- ggarrange(fig4a + scale_x_discrete(guide = guide_axis(angle = -60)),
overall_plot, nrow=1, labels = c("(A)","(B)"), heights = c(1.5,1))
fig4

png(here("figures","seasonality_predictors.png"), width = 1200, height=800)
print(fig4)
# dev.off()


# Check if we have enough data for who_what_where modeling
if (nrow(to_model_scaled) > 0 && sum(!is.na(to_model_scaled$RMSE.norm)) > 0) {
  # Check if required columns exist for who_what_where modeling
required_cols <- c("Taxonomic", "Fungi", "seasonal_predictors", "environmental_predictors", 
                   "new_site", "latitude", "spatial autocorrelation")
missing_cols <- required_cols[!required_cols %in% names(to_model_scaled)]

if (length(missing_cols) == 0) {
  who_what_where_lm <-lm(RMSE.norm ~
  											 	Taxonomic +
  											 	Fungi +
  											 	seasonal_predictors +
  											 	environmental_predictors +
  											 	new_site +
  											 	latitude +
  											 `spatial autocorrelation`,
  											 data = to_model_scaled)
} else {
  cat("Missing columns for who_what_where modeling:", paste(missing_cols, collapse = ", "), "\n")
  who_what_where_lm <- NULL
}
  
  who_what_where_plot <- sjPlot::plot_model(who_what_where_lm,  intercept = F,
  																					sort.est = T,  vline.color = "black", color= c("black", "black")) +
  	theme_bw(base_size = 18) +
  	ggtitle("A priori effects on forecast error (nRMSE)") + ylim(-1.5, 1)
  print(who_what_where_plot)
} else {
  cat("Insufficient data for who_what_where modeling\n")
  who_what_where_plot <- NULL
}

# Check if we have enough data for bacterial modeling
bac_data <- env_df %>% filter(pretty_group=="Bacteria")
if (nrow(bac_data) > 0 && sum(!is.na(bac_data$RMSE.norm)) > 0) {
  bac_lm <-lm(RMSE.norm ~
  							temperature +
  							moisture + pH + `percent carbon` +
  							`ectomycorrhizal trees` + LAI +
  							`temporal memory` + `seasonal amplitude`,
  						data = bac_data)
} else {
  cat("Insufficient data for bacterial modeling\n")
  bac_lm <- NULL
}

# Check if we have enough data for fungal modeling
fun_data <- env_df %>% filter(pretty_group=="Fungi")
if (nrow(fun_data) > 0 && sum(!is.na(fun_data$RMSE.norm)) > 0) {
  fun_lm <-lm(RMSE.norm ~
  							temperature +
  							moisture + pH + `percent carbon` +
  							`ectomycorrhizal trees` + LAI +
  							`temporal memory` + `seasonal amplitude`,
  						data = fun_data)
} else {
  cat("Insufficient data for fungal modeling\n")
  fun_lm <- NULL
}

# Only create param_plot if both models succeeded
if (!is.null(bac_lm) && !is.null(fun_lm)) {
  param_plot <- sjPlot::plot_models(model = list(bac_lm,fun_lm),
  																	m.labels = c("Bacteria","Fungi"),
  										title = "Model parameter effects on forecast error (nRMSE)",
  										legend.title = "Microbial kingdom", vline.color = "black") +
  	theme_bw(base_size = 18) +
  	scale_color_manual(values = c("#F8766D", "#00BFC4"))
  print(param_plot)
} else {
  cat("Cannot create param_plot due to modeling failures\n")
  param_plot <- NULL
}

# Only combine plots if they exist
if (!is.null(who_what_where_plot) && !is.null(param_plot)) {
  ggarrange(who_what_where_plot, param_plot)
} else {
  cat("Cannot combine plots due to missing components\n")
}

# Check if we have enough data for pH plot
if (nrow(env_df) > 0 && sum(!is.na(env_df$pH)) > 0 && sum(!is.na(env_df$RMSE.norm)) > 0) {
  pH_pred <- ggplot(env_df, aes(x = pH,
  																				y = RMSE.norm, color = pretty_group)) +
  	geom_point(alpha=.4, size=3, position=position_jitter(height=0)) +
  	geom_smooth(method="lm") +
  	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
  	xlab("Absolute pH effect") +
  	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
  	ylab("Forecast error (nRMSE) of group\nat observed sites")
  print(pH_pred)
} else {
  cat("Insufficient data for pH plot\n")
  pH_pred <- NULL
}

# Check if spatial autocorrelation column exists
if ("spatial autocorrelation" %in% names(env_df)) {
  patchiness_pred <- ggplot(env_df, aes(x = `spatial autocorrelation`,
  																			y = RMSE.norm, color = pretty_group)) +
  	geom_point(alpha=.4, size=3, position=position_jitter(height=0)) +
  	geom_smooth(method="lm") +
  	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
  	xlab("Spatial autocorrelation (Moran's I)") +
  	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
  	ylab("Forecast error (nRMSE) of group\nat observed sites")
  print(patchiness_pred)
} else {
  cat("Spatial autocorrelation column not found - skipping patchiness plot\n")
  patchiness_pred <- NULL
}

# Check if latitude column exists
if ("latitude" %in% names(env_df)) {
  env_df$orig_latitude = (env_df$latitude * attr(env_df$latitude, "scaled:scale")) + attr(env_df$latitude, "scaled:center")
  
  latitude_pred <- ggplot(env_df, aes(x = orig_latitude,
  																			y = RMSE.norm, color = pretty_group)) +
  	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
  	geom_smooth(method="lm") +
  	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
  	xlab("Latitude") +
  	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) + ylab("Forecast error (nRMSE)")
  print(latitude_pred)
} else {
  cat("Latitude column not found - skipping latitude plot\n")
  latitude_pred <- NULL
}

# Only create plot combinations if all components exist
if (!is.null(who_what_where_plot) && !is.null(param_plot)) {
  pt1 = ggarrange(who_what_where_plot, param_plot, ncol=2, common.legend = T, labels = c("A","B"), legend = "none")
} else {
  cat("Cannot create pt1 due to missing components\n")
  pt1 <- NULL
}

if (!is.null(patchiness_pred) && !is.null(latitude_pred) && !is.null(pH_pred)) {
  pt2 = ggarrange(patchiness_pred, latitude_pred + ylab(NULL), pH_pred + ylab(NULL), nrow=1, common.legend = T, labels = c("C","D","E"))
} else {
  cat("Cannot create pt2 due to missing components\n")
  pt2 <- NULL
}

if (!is.null(pt1) && !is.null(pt2)) {
  ggarrange(pt1, pt2, ncol=1)
} else {
  cat("Cannot create final plot combination due to missing components\n")
}




r2_lat_gradients = c("cortinarius","aspergillus")
rmse_lat_gradients = c("agaricales","aspergillus","mortierellomycota","firmicutes_a","burkholderiaceae")
ggplot(env_df %>% filter(taxon %in% rmse_lat_gradients), aes(x = orig_latitude,
																														 y = RMSE.norm)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_wrap(~taxon, scales = "free")



ggplot(env_df, aes(x = deciduous,
									 y = RMSE.norm, color = pretty_group)) +
	geom_point(alpha=.4, size=3, position=position_jitter(height=0, width = .3)) +
	geom_boxplot() +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("Latitude") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6, label.x.npc = .5) +
	facet_wrap(~pretty_group, scales = "free")



categorical_predictors = no_env_df %>%
	select(c("model_name", "pretty_group", "pretty_name",
					 "fcast_type","site_prediction"))
upset(fromList(categorical_predictors), order.by = "freq")



ma_bac_avg_results <- summary(ma_bac)
avg_bac <- as.data.frame(ma_bac_avg_results$coefmat.subset)
bac_importance <- cbind.data.frame(predictor = names(ma_bac_avg_results$sw), values = ma_bac_avg_results$sw) %>% mutate(Kingdom = "Bacteria")


to_model_fungi =to_model %>% filter(pretty_group=="Fungi")
# Check if required columns exist for fungal predictability modeling
required_cols_fun <- c("Temperature", "Moisture", "pH", "pC", "Ectomycorrhizal\ntrees", 
                       "LAI", "Spatial autocorrelation", "temporal autocorrelation", "seasonality")
missing_cols_fun <- required_cols_fun[!required_cols_fun %in% names(to_model_fungi)]

if (length(missing_cols_fun) == 0) {
  predictability_lm_fun <-lm(RSQ ~ Temperature +
  													 	Moisture + pH + pC +
  													 	`Ectomycorrhizal\ntrees` + LAI + 	`Spatial autocorrelation` +
  													 	`temporal autocorrelation` + seasonality,
  													 data = to_model_fungi)
} else {
  cat("Missing columns for fungal predictability modeling:", paste(missing_cols_fun, collapse = ", "), "\n")
  predictability_lm_fun <- NULL
}
# Only proceed if the model was created successfully
if (!is.null(predictability_lm_fun)) {
  tryCatch({
    models_fun <- lapply(dredge(predictability_lm_fun, evaluate = FALSE, m.lim = c(1,9)), eval)
    ma_fun <- model.avg(models_fun,  subset = delta <= 2)
    plot(ma_fun, full = F, main = "Effects on fungal predictability (RSQ)", intercept = F)
  }, error = function(e) {
    cat("Error in dredge/model averaging:", e$message, "\n")
    ma_fun <- NULL
  })
} else {
  cat("Cannot perform dredge analysis - model not available\n")
  ma_fun <- NULL
}

# Only proceed if ma_fun was created successfully
if (!is.null(ma_fun)) {
  tryCatch({
    ma_fun_avg_results <- summary(ma_fun)
    avg_fun <- as.data.frame(ma_fun_avg_results$coefmat.subset)
    fun_importance <- cbind.data.frame(predictor = names(ma_fun_avg_results$sw),
  																	 values = ma_fun_avg_results$sw) %>% mutate(Kingdom = "Fungi")
  }, error = function(e) {
    cat("Error in fungal model summary:", e$message, "\n")
    fun_importance <- data.frame() # Empty dataframe to avoid errors
  })
} else {
  cat("Cannot create fungal importance - model not available\n")
  fun_importance <- data.frame() # Empty dataframe to avoid errors
}



# Only create importance plot if we have data
if (nrow(fun_importance) > 0 || nrow(bac_importance) > 0) {
  all_importance <- rbind(fun_importance, bac_importance)
  
  if (nrow(all_importance) > 0) {
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
    print(importance_f_b)
  } else {
    cat("No importance data available for plotting\n")
    importance_f_b <- NULL
  }
} else {
  cat("No importance data available for plotting\n")
  importance_f_b <- NULL
}




ph_pred <- ggplot(to_model, aes(x = pH,
																y = RSQ, color = pretty_group)) +
	geom_point(alpha=.5, size=3) +
	geom_smooth(method="lm") +
	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
	xlab("pH (absolute effect)") +
	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6)

# Check if Spatial autocorrelation column exists in to_model
if ("Spatial autocorrelation" %in% names(to_model)) {
  patchiness_pred <- ggplot(to_model, aes(x = `Spatial autocorrelation`,
  																				y = RSQ, color = pretty_group)) +
  	geom_point(alpha=.5, size=3) +
  	geom_smooth(method="lm") +
  	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
  	xlab("Spatial autocorrelation") +
  	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~")), size=6)
} else {
  cat("Spatial autocorrelation column not found in to_model - skipping plot\n")
  patchiness_pred <- NULL
}



multiplot = ggarrange(importance_f_b, ggarrange(ph_pred, patchiness_pred, ncol=1,common.legend = T, labels = c("B","C")), ncol=2, common.legend = T, labels = c("A", NULL
))
multiplot

png(here("figures","explain_predictability.png"), width = 1200, height=800)
print(multiplot)
dev.off()









# model_scores_vals_wide_betas not defined - skipping CRPS plots
cat("model_scores_vals_wide_betas not defined - skipping CRPS plots\n")

# ph_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = pH,
# 																												 y = mean_crps_sample, color = pretty_group)) +
# 	geom_point(alpha=.5, size=3) +
# 	geom_smooth(method="lm") +
# 	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
# 	xlab("pH (absolute effect)") +
# 	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")

# pc_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = pC,
# 																												 y = mean_crps_sample, color = pretty_group)) +
# 	geom_point(alpha=.5, size=3) +
# 	geom_smooth(method="lm") +
# 	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
# 	xlab("pH (absolute effect)") +
# 	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")

# core_sd_pred_crps <- ggplot(model_scores_vals_wide_betas, aes(x = core_sd,
# 																															y = mean_crps_sample, color = pretty_group)) +
# 	geom_point(alpha=.5, size=3) +
# 	geom_smooth(method="lm") +
# 	theme_minimal(base_size = 18) + labs(color = "Kingdom") +
# 	xlab("Variation among soil cores") +
# 	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) + ylab("Forecast error (CRPS)")

# ggarrange(ph_pred_crps, pc_pred_crps, core_sd_pred_crps, common.legend = T)


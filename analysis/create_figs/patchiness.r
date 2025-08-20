source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
library(MuMIn)
options(na.action = "na.fail")

# Read in moran's I stats
morans <- readRDS(here("data/clean", "moran_stat_MRM.rds"))

morans <-  morans %>%  mutate_at(c("Int", "intraAnnual",
																	 "interAnnual", "space", "overall_Rsq"), as.numeric) %>%
	filter(!grepl("species", rank))
morans <- merge(morans, to_model[,c("rank","pretty_group","taxon","fcast_type","CRPS_truncated","RSQ.1","unexplained seasonality")])

morans_long = morans %>% pivot_longer(cols = c("Int", "intraAnnual",
																							 "interAnnual", "space", "overall_Rsq","CRPS_truncated","unexplained seasonality"), names_to = "param")

to_plot = morans_long %>% filter(!is.na(pretty_group) &
																 	param %in% c("unexplained seasonality","interAnnual","space","CRPS_truncated")) %>%
																 	mutate(param = recode(param, "interAnnual"= "temporal autocorrelation (days apart)",
																 																																		"intraAnnual"="temporal autocorrelation (intra-annual)",
																 																																		"space"="spatial autocorrelation (sample distance)", "CRPS_truncated" = "Predictability (CRPS)"))
ggplot(to_plot,
			 aes(x = pretty_group, y = value,
											 color = pretty_group)) +
	geom_boxplot(show.legend = F) +
	geom_point(position = position_jitterdodge(), alpha = .3, show.legend = F)  +
	xlab("Kingdom") +
	theme_minimal(base_size = 18) +
	stat_compare_means(na.rm = T, show.legend = F) +
	geom_hline(yintercept = 0) + facet_wrap(~param, scales="free") +
	ggtitle("Bacteria have stronger spatio-temporal autocorrelation and more predictability than fungi")

ggplot(morans %>% filter(!is.na(pretty_group)),
			 aes(x = pretty_group, y = space,
																										 color = pretty_group)) +
	geom_boxplot(show.legend = F, outlier.alpha = .1) +
	geom_jitter(alpha = .5, show.legend = F)  +
	xlab("Kingdom") +
	theme_minimal(base_size = 20) +
	stat_compare_means(na.rm = T, show.legend = F) +
	geom_hline(yintercept = 0) +
	ylab("Spatial autocorrelation") #+ coord_flip()

# spatial vs temporal variability
fg_summary_beta <- readRDS(here("data", paste0("summary/beta_fg_summaries_20151101_20200101.rds")))
fg_obs <- fg_summary_beta$plot_est %>% filter(!is.na(truth))
fg_obs$rank <- fg_obs$taxon

tax_obs <- readRDS(here("data", "clean/truth_plot_long_tax.rds")) %>% filter(!is.na(truth)) %>% mutate(truth=as.numeric(truth))
tax_obs$pretty_group <- ifelse(grepl("fun", tax_obs$rank), "Fungi", "Bacteria")
obs <- rbindlist(list(fg_obs, tax_obs), fill=T)

# Calculate within-site "patchiness:" CV for a given month, averaged at each site
patchiness <- obs %>% group_by(pretty_group,rank,taxon,siteID,dateID) %>%
	summarize(within_site_date_patchy = calc_cv(truth)) %>%
	ungroup %>%
	group_by(pretty_group,rank,taxon) %>%
	summarize(within_site_patchy = mean(within_site_date_patchy, na.rm=T))

# Calculate across-site "patchiness:" CV for a given month, across all sites
across_site_patchiness <- obs %>% group_by(pretty_group,rank,taxon,dateID) %>%
	summarize(across_site_date_patchy = calc_cv(truth)) %>%
	ungroup %>%
	group_by(pretty_group,rank,taxon) %>%
	summarize(across_site_patchy = mean(across_site_date_patchy, na.rm=T))


morans_patchy <- merge(as.data.frame(morans), across_site_patchiness)

a <- ggplot(patchiness) +
	geom_jitter(aes(x = pretty_group, y = within_site_patchy,
									color = pretty_group))  + xlab("Kingdom") + theme_minimal(base_size = 18)

b <- ggplot(across_site_patchiness) +
	geom_jitter(aes(x = pretty_group, y = across_site_patchy,
									color = pretty_group)) + xlab("Kingdom")  + theme_minimal(base_size = 18)

ggpubr::ggarrange(a,b, common.legend = T) + theme_minimal()


# Predictability scores
predict_scores <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
scores <- predict_scores$scoring_metrics %>%
	filter(model_name == "all_covariates" & pretty_name != "Diversity") %>% distinct()
scores$RSQ.1 <- ifelse(scores$RSQ.1 > 0, scores$RSQ.1, 0)
# scores <- predict_scores$cv_metric_scaled %>%
# 	filter(model_name == "all_covariates" & pretty_name != "Diversity" &
# 				 	cv_type == "mean_per_plot_cv" & metric == "RSQ.1") %>% distinct()

# Beta estimates
beta_effects <- readRDS(here("data/summary/all_fcast_effects.rds"))
beta_effects_min <- beta_effects %>% filter(taxon != "other" & !is.na(model_name)) %>%
	select(siteID, model_name, group, rank_only, fcast_type, pretty_group, rank_name,
				 Mean, SD, beta, taxon, significant, effSize,time_period) %>% distinct()
beta_effects_wide <- beta_effects_min %>% pivot_wider(id_cols = c("model_name","rank_name",
																											"pretty_group","taxon","rank_only",
																											"time_period"),
																					names_from = "beta", values_from = "effSize")
beta_effects_wide <- beta_effects_wide %>% filter(time_period=="2015-11_2018-01" & model_name == "all_covariates")
# Read in seasonality values
seasonality <- readRDS(here("data/summary/seasonalAmplitude.rds"))
all_cov_vals <- seasonality[[1]] %>% select(taxon, fcast_type,
																						pretty_name, pretty_group,
																						resid_max = max, residual_amplitude = amplitude)
cycl_vals <- seasonality[[2]] %>% select(taxon, fcast_type,
																				 pretty_name, pretty_group, max, amplitude)

# Merge data types
to_model <- list(cycl_vals, all_cov_vals, beta_effects_wide, scores, patchiness, across_site_patchiness) %>%
	reduce(full_join) #%>% filter(fcast_type=="Taxonomic")
to_model <- merge(to_model, morans, all=T)


#to_model[,c("RSQ.1","amplitude","space")]


# so dumb
to_model$fg_cats <- NA
to_model$group <- NA
to_model[is.na(to_model$pretty_group),]$fg_cats = assign_fg_categories(to_model[is.na(to_model$pretty_group),]$rank)
to_model[is.na(to_model$pretty_group),]$group = assign_fg_kingdoms(to_model[is.na(to_model$pretty_group),]$fg_cats)
to_model[is.na(to_model$pretty_group),]$pretty_group <- ifelse(to_model[is.na(to_model$pretty_group),]$group=="16S", "Bacteria", "Fungi")


ggplot(to_model %>% filter(!is.na(pretty_group)),
			 aes(x = space, y = mean_morans, color = pretty_group)) +
	geom_point() +
	theme_minimal(base_size = 18) +
	stat_smooth(method = "lm") +
	ggtitle("Moran's and MRM space coefficient agree..?")


library(ggpubr)
ggplot(to_model %>% filter(!is.na(pretty_group)),aes(x = pretty_group, y = space,
																										 color = pretty_group)) +
	geom_boxplot() +
	geom_jitter()  +
	xlab("Kingdom") +
	theme_minimal(base_size = 18) +
	stat_compare_means(na.rm = T) +
	geom_hline(yintercept = 0)

ggplot(to_model %>% filter(!is.na(pretty_group)),aes(x = pretty_group, y = interAnnual,
																										 color = pretty_group)) +
	geom_boxplot() +
	geom_jitter()  +
	xlab("Kingdom") +
	theme_minimal(base_size = 18) +
	stat_compare_means(na.rm = T) +
	geom_hline(yintercept = 0)

to_model <- to_model %>% select(c("taxon", "rank", "fcast_type", "pretty_group",
																	"seasonality"="amplitude",
																	"unexplained seasonality"="residual_amplitude",
																	"Temperature", "Moisture",
																	"pH", "pC", "Ectomycorrhizal\ntrees", "LAI","CRPS", "RSQ", "RSQ.1",  "CRPS_truncated","mean_morans", "intraAnnual",
																	"interAnnual", "space")) %>% drop_na

#to_model <-  to_model %>% mutate_at(c("mean_morans","intraAnnual",
#																	 "interAnnual", "space"), function(x) scale(abs(x)))

# MEAN/CENTER SCALE EVERYTHING?
to_model <-  to_model %>% mutate_at(c("seasonality", "unexplained seasonality", "Temperature", "Moisture",
																			"pH", "pC", "Ectomycorrhizal\ntrees", "LAI","CRPS","mean_morans", "intraAnnual",
																			"interAnnual", "space"), function(x) scale(x))


# Model predictability as a function of everything!
predictability_lm <-lm(CRPS_truncated ~ Temperature +
											 	Moisture + pH + pC +
											 	`Ectomycorrhizal\ntrees` + LAI +
											 	`unexplained seasonality` +
											 	seasonality +
											 	interAnnual + intraAnnual + space,
											 data = to_model)
#models <- lapply(dredge(predictability_lm, evaluate = FALSE, m.lim = c(1,5)), eval)
models <- lapply(dredge(predictability_lm, evaluate = FALSE), eval)
ma <- model.avg(models, subset = delta <= 2)
plot(ma, full = F, main = "Effects on predictability (CRPS)", intercept = F)
#plot(ma, full = F, main = "Effects on predictability (RSQ)")
#ggpubr::ggarrange(c,d)
avg_results <- summary(ma)
avg <- as.data.frame(avg_results$coefmat.subset)
importance <- cbind.data.frame(predictor = names(ma$sw), values = ma$sw) %>% mutate(Kingdom = "Overall")

importance = importance %>% mutate(predictor = recode(predictor, "interAnnual"= "temporal autocorrelation (days apart)",
																 "intraAnnual"="temporal autocorrelation (intra-annual)",
																	"space"="spatial autocorrelation (sample distance)" ))
ggplot(importance,
			 aes(x = reorder(predictor, -values), y = values)) +
	geom_point(size=3, alpha=.8)  +
	theme_minimal(base_size = 18) +
	geom_hline(yintercept = 0) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	xlab("Variables explaining predictability across 141 groups") +
	ylab("Variable importance")



bac <- to_model[to_model$pretty_group=="Bacteria",]
predictability_lm <-lm(CRPS_truncated ~ Temperature +
											 	Moisture + pH + pC +
											 	`Ectomycorrhizal\ntrees` + LAI +
											 	`unexplained seasonality` +
											 	seasonality +
											 	interAnnual + intraAnnual + space,
											 ,
											 data = bac)
models <- lapply(dredge(predictability_lm, evaluate = FALSE)#, m.lim = c(1,5))
								 , eval)
ma <- model.avg(models, subset = delta <= 2)
plot(ma, full = F, main = "Effects on predictability (CRPS), bacteria", intercept = F)
avg_results <- summary(ma)
bac_avg <- as.data.frame(avg_results$coefmat.subset)
bac_importance <- cbind.data.frame(predictor = names(ma$sw), values = ma$sw) %>% mutate(Kingdom = "Bacteria")



fun <- to_model[to_model$pretty_group=="Fungi",]
predictability_lm <-lm(CRPS_truncated ~ Temperature +
											 	Moisture + pH + pC +
											 	`Ectomycorrhizal\ntrees` + LAI +
											 	#residual_amplitude +
											 	`unexplained seasonality` +
											 	seasonality +
											 	interAnnual + intraAnnual + space,
											 data = fun)
models <- lapply(dredge(predictability_lm, evaluate = FALSE), #, m.lim = c(1,5)),
								 eval)
ma <- model.avg(models, subset = delta <= 2)
plot(ma, full = F, main = "Effects on predictability (CRPS), fungi", intercept = F)
avg_results <- summary(ma)
fun_avg <- as.data.frame(avg_results$coefmat.subset)
fun_importance <- cbind.data.frame(predictor = names(ma$sw), values = ma$sw) %>% mutate(Kingdom = "Fungi")



all_importance <- rbind(importance, fun_importance, bac_importance)
ggplot(all_importance,
			 aes(x = reorder(predictor, -values), y = values, color = Kingdom)) +
	geom_point(position = position_jitterdodge(jitter.width = .1, jitter.height = 0),
						 size=3, alpha=.8)  +
	theme_minimal(base_size = 18) +
	geom_hline(yintercept = 0) +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	xlab("Variables explaining predictability across groups") +
	ylab("Variable importance") +
	scale_color_manual(values = c(#"Fungi" = "#00BA38", # green
																"Fungi" = "#00BFC4",
																"Bacteria" = "#F8766D",
																"Overall" = "black"))


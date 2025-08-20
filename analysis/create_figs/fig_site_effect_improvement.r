source("source.R")
library(ggallin)
old_scores_list <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
new_scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))

# Check if Parquet files exist, otherwise use RDS
old_parquet <- here("data/summary/parquet/all_hindcasts.parquet")
old_rds <- here("data/summary/all_hindcasts.rds")
new_parquet <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
new_rds <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(old_parquet)) {
  cat("Using old hindcast Parquet file...\n")
  old_hindcast <- arrow::read_parquet(old_parquet)
} else if (file.exists(old_rds)) {
  cat("Old hindcast Parquet not found, using RDS...\n")
  old_hindcast <- readRDS(old_rds)
} else {
  stop("Old hindcast file not found!")
}

if (file.exists(new_parquet)) {
  cat("Using new hindcast Parquet file...\n")
  new_hindcast <- arrow::read_parquet(new_parquet)
} else if (file.exists(new_rds)) {
  cat("New hindcast Parquet not found, using RDS...\n")
  new_hindcast <- readRDS(new_rds)
} else {
  stop("New hindcast file not found!")
}




# Check if Parquet files exist, otherwise use RDS
old_parquet <- here("data/summary/parquet/all_hindcasts.parquet")
old_rds <- here("data/summary/all_hindcasts.rds")
new_parquet <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
new_rds <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(old_parquet)) {
  cat("Using old hindcast Parquet file...\n")
  old_hindcast <- arrow::read_parquet(old_parquet)
} else if (file.exists(old_rds)) {
  cat("Old hindcast Parquet not found, using RDS...\n")
  old_hindcast <- readRDS(old_rds)
} else {
  stop("Old hindcast file not found!")
}

if (file.exists(new_parquet)) {
  cat("Using new hindcast Parquet file...\n")
  new_hindcast <- arrow::read_parquet(new_parquet)
} else if (file.exists(new_rds)) {
  cat("New hindcast Parquet not found, using RDS...\n")
  new_hindcast <- readRDS(new_rds)
} else {
  stop("New hindcast file not found!")
}


converged <- new_scores_list$converged_list

ggplot(new_hindcast %>%
			 	filter(plotID %in% c("DEJU_001") &
			 		model_id %in% c("env_cycl_copiotroph_20151101_20180101"))
) +
	#facet_grid(rows=vars(species), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, na.rm = T) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue", na.rm = T) +
	geom_ribbon(aes(x = dates, ymin = lo_25, ymax = hi_75),fill="red", alpha=0.6, na.rm = T) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(model_id+plotID~predicted_site_effect) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')

# These variables are not defined - skipping these plots
cat("hindcast_plsr and hindcast not defined - skipping related plots\n")



old_hindcast_not_na <- old_hindcast %>% filter(!is.na(crps)) %>%
	filter(new_site==T) %>%
	group_by(new_site, pretty_group,model_name,predicted_site_effect) %>% summarize(mean(crps))
new_hindcast_not_na <- new_hindcast %>% filter(!is.na(crps)) %>%
	filter(new_site==T) %>%
	group_by(pretty_group,model_name,predicted_site_effect) %>% summarize(mean(crps))
# new_hindcast2 is not defined - skipping this calculation
cat("new_hindcast2 not defined - skipping this calculation\n")

# This plot will be created after stat_pvalue is defined
# ggplot(new_hindcast_not_na) + geom_point(aes(x = predicted_site_effect, y = `mean(crps)`, color = pretty_group)) + facet_grid(pretty_group~model_name) + stat_pvalue_manual(stat_pvalue, y.position = .023)

new_hindcast <- new_hindcast %>% mutate(site_effect = ifelse(predicted_site_effect==T,"Modeled (PLSR)","Random"))

# Define stat_pvalue for the plots
stat_pvalue <- new_hindcast %>% filter(!is.na(crps)) %>%
	filter(new_site==T & !is.na(site_effect))  %>%
	ungroup() %>%
	group_by(pretty_group,model_name) %>%
	compare_means(crps ~ site_effect, data =.,
								group.by = c("pretty_group","model_name"))

ggplot(new_hindcast %>% filter(!is.na(crps)) %>%
			 	filter(new_site==T & !is.na(site_effect)),
			 aes(x = site_effect, y = crps, color = pretty_group)) +
	#geom_point(aes(x = `site effect`, y = crps, color = pretty_group), position = position_jitter(), alpha=.1) +
	geom_boxplot(outlier.shape = NA) +
	#	geom_point(aes(x = predicted_site_effect, y = `mean(crps)`, color = pretty_group)) +
	facet_grid(pretty_group~model_name) +
	scale_y_continuous(trans = pseudolog10_trans) +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 show.legend = F, hide.ns = F, size=4)




ggplot(new_hindcast %>% filter(!is.na(crps)) %>%
			 	filter(new_site==T & !is.na(site_effect)), aes(x = site_effect, y = crps, color = pretty_group)) +
	#geom_point(aes(x = `site effect`, y = crps, color = pretty_group), position = position_jitter(), alpha=.1) +
	#	geom_point(aes(x = predicted_site_effect, y = `mean(crps)`, color = pretty_group)) +
	facet_grid(pretty_group~model_name) +
	scale_y_continuous(trans = pseudolog10_trans,limits = c(0,.035)) +
	#geom_violin(draw_quantiles = c(.5), show.legend = F) +

	geom_boxplot(outlier.shape = NA) +
	stat_pvalue_manual(stat_pvalue, y.position = .005, label="p.signif") + theme_bw(base_size = 22)
	#coord_cartesian(ylim = c(0, .1))


stat_pvalue <- new_hindcast %>% filter(!is.na(crps)) %>%
	filter(new_site==T & !is.na(site_effect))  %>%
	ungroup() %>%
	group_by(pretty_group,model_name) %>%
	#mutate(modeled_site_effect = as.factor(predicted_site_effect)) %>%
	#filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	compare_means(crps ~ site_effect, data =.,
								group.by = c("pretty_group","model_name"))
#	rstatix::tukey_hsd(crps ~ modeled_site_effect)

new_hindcast_wide <- new_hindcast %>% pivot_wider(names_from = site_effect, values_from = crps) %>%
	mutate(site_eff_improvement = Random - `Modeled (PLSR)`)
ggplot(new_hindcast_wide,
			 aes(x = pretty_group, y = site_eff_improvement, group = model_id)) +
	geom_line(alpha=.2) +
	geom_point(size = 2,position = position_jitter(width = .05), alpha=.5) +
	facet_wrap( ~ pretty_group, switch = "x") +
	scale_x_discrete("") +
	theme_minimal() +
	scale_y_continuous(trans = pseudolog10_trans)

# scores_list is not defined - using new_scores_list instead
newsite_scores_taxon = new_scores_list$scoring_metrics_cv %>%
	filter(model_id %in% converged) %>%
	filter(site_prediction != "New time (observed site)" & model_name=="env_cycl") %>%
	select(-c(metric, score)) %>% distinct()

newsite_scores_taxon_wide <- newsite_scores_taxon %>%
	pivot_wider(names_from = site_prediction, values_from = mean_crps_sample) %>%
	mutate(site_eff_improvement = `New time x site (modeled effect)` - `New time x site (random effect)`)

# THIS IS THE APPROACH
# scores_list is not defined - using new_scores_list instead
skill_scores = new_scores_list$skill_score_taxon %>%
	filter(model_name=="env_cycl") %>%

	filter(model_id %in% converged) %>% pivot_longer(cols = c(skill_score,skill_score_random))
ggplot(skill_scores,
			 aes(x = name, y = value, group = model_id)) +
	geom_line(alpha=.2) +
	geom_point(size = 2, aes(color = name), position = position_jitter(width = .05), alpha=.5) +
	facet_wrap( ~ pretty_group, switch = "x") +
	scale_x_discrete("") +
	theme_minimal() +
	theme(legend.position = "top",
				panel.grid = element_blank(),
				axis.text.x = element_blank(),
				axis.line.y = element_line(size = .5)) +
	scale_y_continuous(trans = pseudolog10_trans)#+
#	scale_y_log10()

stat_pvalue <- skill_scores %>% filter(!is.na(mean_crps_sample)) %>%

	#filter(new_site==T & !is.na(site_effect))  %>%
	ungroup() %>%
	group_by(pretty_group,model_name) %>%
	#mutate(modeled_site_effect = as.factor(predicted_site_effect)) %>%
	#filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	compare_means(value ~ name, data =.,
								group.by = c("pretty_group","model_name"))

stat_pvalue <- newsite_scores_taxon %>% filter(!is.na(mean_crps_sample)) %>%
	filter(site_prediction != "New time (observed site)" & model_name=="env_cycl") %>%

	#filter(new_site==T & !is.na(site_effect))  %>%
	ungroup() %>%
	group_by(pretty_group,model_name) %>%
	#mutate(modeled_site_effect = as.factor(predicted_site_effect)) %>%
	#filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	compare_means(mean_crps_sample ~ site_prediction, data =.,
								group.by = c("pretty_group","model_name"))
ggplot(newsite_scores_taxon,
			 aes(x = site_prediction, y = mean_crps_sample, group = model_id)) +
	geom_line(alpha=.2) +
	geom_point(size = 2, aes(color = site_prediction), position = position_jitter(width = .05), alpha=.5) +
	facet_wrap( ~ pretty_group, switch = "x") +
	scale_x_discrete("") +
	theme_minimal() +
	theme(legend.position = "top",
				panel.grid = element_blank(),
				axis.text.x = element_blank(),
				axis.line.y = element_line(size = .5)) +
	scale_y_log10()

stat_pvalue_univ <- new_hindcast %>% filter(!is.na(crps)) %>%
	filter(new_site==T & !is.na(site_effect))  %>%
	ungroup() %>%
	group_by(pretty_group,model_name) %>%
	#mutate(modeled_site_effect = as.factor(predicted_site_effect)) %>%
	#filter(fcast_type=="Functional") %>%
	#group_by(pretty_group) %>%
	compare_means(crps ~ site_effect, data =.)
#	rstatix::tukey_hsd(crps ~ modeled_site_effect)


ggplot(new_hindcast %>% filter(!is.na(crps)) %>%
			 	filter(new_site==T & !is.na(site_effect)), aes(x = site_effect, y = crps)) +
	#geom_point(position = position_jitter(), alpha=.1) +
	#	geom_point(aes(x = predicted_site_effect, y = `mean(crps)`, color = pretty_group)) +
	#geom_violin(draw_quantiles = c(.5), show.legend = F) +

	geom_boxplot(outlier.shape = NA) +
	stat_pvalue_manual(stat_pvalue_univ, y.position = .005, label="p.signif") + theme_bw(base_size = 22) +
	scale_y_continuous(trans = pseudolog10_trans,limits = c(0,.035))

hindcast_crps_wide = new_hindcast %>%
	filter(new_site==T & !is.na(site_effect)) %>%
	select(pretty_group, model_name, model_id, plotID, timepoint, site_effect, crps) %>%
	pivot_wider(names_from = site_effect, values_from = crps)


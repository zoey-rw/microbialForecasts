# Visualize effect size estimates (beta covariates) from all model
source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr)
library(rstatix)
library(ggh4x)

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

sum.all <- readRDS(here("data/summary/predictor_effects.rds"))

seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_only_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="cycl_only") %>%
	mutate(cycl_amplitude=amplitude) %>%
	select(-c(sin,cos,max,amplitude_orig,max)) %>%
	pivot_longer(#id_cols = c("model_name","pretty_group","taxon","time_period"),
		cols=cycl_amplitude, values_to = "effSize", names_to="beta")
env_cycl_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="env_cycl")  %>%
	mutate(residual_amplitude=amplitude) %>%
	select(-c(sin,cos,max,amplitude_orig,max)) %>%
	pivot_longer(#id_cols = c("model_name","pretty_group","taxon","time_period"),
							 cols=residual_amplitude, values_to = "effSize", names_to="beta")



df_cal_fg_tax <- sum.all %>% filter(time_period == "2015-11_2018-01") %>%
	filter(model_name == "env_cycl" & model_id %in% converged)  %>%
	filter(!beta %in% c("sin","cos"))

df_cal_fg_tax <- rbindlist(list(df_cal_fg_tax, env_cycl_vals_scores, cycl_only_vals_scores), fill=T)


df_cal_fg_tax$beta_pretty <- recode(df_cal_fg_tax$beta, "residual_amplitude" = "seasonality",
																		"cycl_amplitude" = "seasonality",
																		"pC" = "percent carbon",
																		"LAI" = "Leaf area index")


df_cal_fg_tax$beta_pretty <- factor(df_cal_fg_tax$beta_pretty,
																	 levels=c("seasonality", "Ectomycorrhizal\ntrees", "Leaf area index",
																	 				 "percent carbon", "pH", "Temperature", "Moisture"))

beta.labs = c("sin", "cos",
							"Ectomycorrhizal\ntrees", "Leaf area index", "percent carbon",
							"pH", "Temperature", "Moisture", "rho", "seasonality")
names(beta.labs) = beta.labs
group.labs = c("Functional", "Taxonomic")
names(group.labs) = group.labs
label_fn2 <- as_labeller(c(model.labs, beta.labs, group.labs))


stat.test <- df_cal_fg_tax %>%
	group_by(fcast_type,model_name,beta_pretty) %>%
	t_test(effSize ~ pretty_group) %>%
	adjust_pvalue(method = "fdr")  %>%
	add_significance() %>%
		add_xy_position(x = "pretty_group")


###### # Plot with Tukey, fungi vs bacteria effect sizes -----
b_vs_f_fcast_type_plot <- ggplot(data=df_cal_fg_tax,
																 aes(x = pretty_group,
																 		color = pretty_group, y = effSize)) +

	geom_violin(draw_quantiles = c(.5), show.legend = F, color = 1) +
	geom_jitter(aes(color = pretty_group),
							size = 5, height = 0, width=.1, alpha = .2,
		shape = 16, show.legend = F) +
	xlab("Kingdom")+
	ylab("Absolute effect size") +
	facet_nested(fcast_type~model_name + beta_pretty,
							 labeller = label_fn2) +
	theme_classic(base_size = 18) +
	theme(
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05)) +
		scale_y_log10() +
		stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, size=5)
#+ stat_compare_means(method="t.test", hide.ns = T, show.legend = F)  +
	# +coord_flip()
# b_vs_f_fcast_type_plot <- b_vs_f_fcast_type_plot  + stat_compare_means(
# 	aes(label = after_stat(p.adjust)),
# 	method = "t.test",
# 	label.y.npc = .8,
# 	label.x=1.5,
# 	hide.ns = T,
# 	size=10, show.legend = F)

b_vs_f_fcast_type_plot

png(here("figures","effsize_f_b.png"), width = 1600, height=800)
print(b_vs_f_fcast_type_plot)
dev.off()


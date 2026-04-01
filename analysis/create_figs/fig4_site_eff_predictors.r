source("source.R")
library(tidytext)


# Read in the dredged results from 05_predictSiteEffects.r
site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

# Extract the components from the saved list structure
# site_eff_dredged_in is a list with 4 elements:
# 1: dredged_predictor_importance
# 2: pred_sites  
# 3: pred_sites_plsr
# 4: plsr_model_importance

# Get the predictor importance data (element 1)
obs_sites <- site_eff_dredged_in[[1]] %>% 
  filter(model_id %in% converged) %>%
  as.data.frame()

# Define pred_sites from the site effects data (element 2)
pred_sites <- site_eff_dredged_in[[2]] %>% 
  filter(model_id %in% converged) %>%
  mutate(fcast_type = ifelse(grepl("function", rank_only), "Functional group",
                             ifelse(grepl("diversity", rank_only), "Diversity", "Taxonomic group")),
         Mean = pred)  # Use pred as Mean for the plots


obs_sites$fcast_type = ifelse(grepl("function", obs_sites$rank_only), "Functional group",
															 ifelse(grepl("diversity", obs_sites$rank_only), "Diversity", "Taxonomic group"))

# # Add zeros since unimportant predictors were dropped by dredge
# obs_sites$values <- as.numeric(obs_sites$values)
# obs_sites <- obs_sites %>% group_by(pretty_group, fcast_type, model_name, rank_only) %>% complete(taxon, predictor, fill = list(values = 0))

# Define predictor category and name mappings first
pred_cat_key = list("nitrogenTot" = "macronutrient",
										"totalP" = "macronutrient",
										"sulfurTot" = "macronutrient",
										"estimatedOC" = "macronutrient",
										"pMjelm"= "macronutrient",
										"feMjelm" = "metal",
										"alOxalate" = "metal",
										"cecdNh4" = "cations",
										"mgNh4d" = "cations",
										"naNh4d" = "cations",
										"mnMjelm" = "cations",
										"kNh4d"  = "cations",
										"caNh4d" = "cations",
										"MAT"	     = "climate",
										"MAP"      = "climate",
										"latitude_scaled"      = "climate",
										"so4Satx"= "micronutrient",
										"siMjelm"  = "micronutrient")


pred_name_key = list("no2Satx" = "nitrite",
										 "totalP" = "total phosphorus",
										 "nitrogenTot" = "total nitrogen",
										 "estimatedOC" = "organic carbon",
										 "no3Satx" = "nitrate",
										 "pMjelm"= "phosphorus (Mehlich)",
										 "alOxalate"= "aluminum",
										 "feMjelm" = "iron (Mehlich)",
										 "alKcl" = "aluminum",
										 "cecdNh4" = "cation exchange",
										 "mgNh4d" = "magnesium",
										 "naNh4d" = "sodium",
										 "mnMjelm" = "manganese (Mehlich)",
										 "kNh4d"  = "potassium",
										 "caNh4d" = "calcium",
										 "MAT"	     = "mean annual temperature",
										 "MAP"      = "mean annual precipitation",
										 "latitude_scaled" = "latitude",
										 "so4Satx"= "sulfate",
										 "sulfurTot"= "sulfur",
										 "siMjelm"  = "silicon (Mehlich)")

# Now apply the recoding
# The predictor importance data already has 'values' column from the dredged results
obs_sites$predictor_category = recode(obs_sites$predictor, !!!pred_cat_key)
obs_sites$predictor = recode(obs_sites$predictor, !!!pred_name_key)



all.out <- obs_sites %>%
	group_by(predictor, model_name) %>%
	mutate(importance = round(mean(values), 2)) %>%
	ungroup() %>%
	mutate(predictor = factor(predictor),
				 predictor = fct_reorder(predictor, importance))

group_vals = obs_sites %>%
	group_by(pretty_group, model_name, predictor_category) %>%
	mutate(mean_importance = mean(values, na.rm = TRUE),
				 se_importance = sd(values, na.rm = TRUE) / sqrt(sum(!is.na(values)))) %>%
	mutate(ymax = mean_importance + 1.96 * se_importance,
				 ymin = mean_importance - 1.96 * se_importance) %>%
	arrange(pretty_group, predictor_category) %>% ungroup


pred_vals = obs_sites %>%
	group_by(pretty_group, model_name, predictor) %>%
	mutate(mean_importance = mean(values, na.rm = TRUE),
				 se_importance = sd(values, na.rm = TRUE) / sqrt(sum(!is.na(values)))) %>%
	mutate(ymax = mean_importance + 1.96 * se_importance,
				 ymin = mean_importance - 1.96 * se_importance) %>%
	arrange(pretty_group, predictor) %>% ungroup

# -- Shared styling --
kingdom_colors <- c("Bacteria" = "#E07A5F", "Fungi" = "#3D85C6")
shared_theme <- theme_bw(base_size = 13) +
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.ticks = element_line(color = "grey60", linewidth = 0.3),
		legend.position = "top",
		legend.title = element_text(face = "bold", size = 11),
		legend.text = element_text(size = 10),
		plot.tag = element_text(size = 14, face = "bold")
	)

# Capitalize predictor names for display
str_title_case <- function(x) {
	sapply(x, function(s) paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = ""))
}

# -- Panel A: Individual predictor importance --
pred_vals_plot <- pred_vals %>%
	filter(model_name %in% c("env_cycl")) %>%
	mutate(predictor = str_title_case(predictor))

overall_importance <- ggplot(pred_vals_plot,
	aes(x = reorder(predictor, -values), y = values, color = pretty_group)) +
	geom_pointrange(aes(y = mean_importance, ymin = ymin, ymax = ymax),
		position = position_dodge(width = 0.5), size = 0.4, fatten = 2.5) +
	stat_compare_means(method = "t.test",
		aes(y = values, label = after_stat(p.signif)),
		label.y = 0.68, show.legend = FALSE, hide.ns = TRUE, size = 3.5) +
	scale_color_manual(values = kingdom_colors, name = NULL) +
	coord_cartesian(ylim = c(0, 0.75)) +
	labs(x = "Predictor", y = "Importance for\nexplaining site effects", tag = "A") +
	shared_theme +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
		legend.position = "none"
	)

# -- Panel B: Predictor category importance --
group_vals_plot <- group_vals %>%
	filter(model_name %in% c("env_cycl")) %>%
	mutate(predictor_category = str_title_case(predictor_category))

f_b_category <- ggplot(group_vals_plot,
	aes(x = reorder(predictor_category, -mean_importance),
		y = mean_importance, color = pretty_group)) +
	geom_pointrange(aes(ymin = ymin, ymax = ymax),
		position = position_dodge(width = 0.5), size = 0.5, fatten = 3) +
	stat_compare_means(method = "t.test",
		aes(y = values, label = after_stat(p.signif)),
		label.y = 0.68, show.legend = FALSE, hide.ns = TRUE, size = 3.5) +
	scale_color_manual(values = kingdom_colors, name = NULL) +
	coord_cartesian(ylim = c(0, 0.75)) +
	labs(x = "Predictor category", y = NULL, tag = "B") +
	shared_theme +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
	)

# -- Composite figure --
supp_fig <- ggarrange(overall_importance, f_b_category,
	widths = c(1.6, 1), align = "h",
	common.legend = TRUE, legend = "top")

ggsave(here("figures", "figS12_plsr_importance.png"), supp_fig,
	width = 11, height = 5.5, dpi = 300, bg = "white")

ggsave(here("figures", "site_effect_f_b_category.png"), f_b_category,
	width = 5, height = 5.5, dpi = 300, bg = "white")

# ── Panel C: ECDF of variogram p-values (env_cycl, before/after PLSR) ────────
variograms <- readRDS(here("data/summary/site_effect_variograms.rds"))
sig_results <- variograms[[1]] %>%
  mutate(model_name = ifelse(is.na(model_name),
    case_when(grepl("^env_cycl", model_id)  ~ "env_cycl",
              grepl("^env_cov", model_id)   ~ "env_cov",
              grepl("^cycl_only", model_id) ~ "cycl_only",
              TRUE ~ NA_character_), model_name))

# Use same converged list as rest of script
sig_var_long <- sig_results %>%
  filter(model_id %in% converged, model_name == "env_cycl") %>%
  pivot_longer(cols = c("site effect", "site effect residuals"),
               names_to = "stage", values_to = "pval") %>%
  mutate(pval = as.numeric(pval),
         stage_label = factor(
           ifelse(stage == "site effect",
                  "Raw site effects", "After PLSR modeling"),
           levels = c("Raw site effects", "After PLSR modeling")))

# Compute % significant for annotation
pct_sig_var <- sig_var_long %>%
  group_by(stage_label) %>%
  summarise(pct = round(100 * mean(pval < 0.05, na.rm = TRUE)), .groups = "drop")

ecdf_plot <- ggplot(sig_var_long, aes(x = pval, color = stage_label)) +
  stat_ecdf(linewidth = 0.8, pad = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dotted", linewidth = 0.6) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey70") +
  annotate("text", x = 0.72, y = 0.42,
           label = "Uniform reference\n(no autocorrelation)",
           color = "grey50", size = 3.2, fontface = "italic") +
  annotate("text", x = 0.72, y = 0.28,
           label = paste0("Significant (p<0.05):  ",
                          pct_sig_var$pct[pct_sig_var$stage_label == "Raw site effects"], "% raw  vs  ",
                          pct_sig_var$pct[pct_sig_var$stage_label == "After PLSR modeling"], "% after PLSR"),
           size = 3.2) +
  scale_color_manual(values = c("Raw site effects"    = "#D55E00",
                                "After PLSR modeling" = "#0072B2"),
                     name = NULL) +
  labs(x = "Variogram p-value",
       y = "Cumulative proportion of taxa",
       tag = "C") +
  shared_theme +
  theme(legend.position = c(0.35, 0.92),
        legend.background = element_blank(),
        legend.key.width = unit(1.2, "cm"))

# ── 3-panel composite (2 rows) ────────────────────────────────────────────────
row1 <- ggarrange(overall_importance, f_b_category,
                  widths = c(1.6, 1), align = "h",
                  common.legend = TRUE, legend = "top") +
  theme(plot.margin = margin(t = 0, b = 0))

ecdf_plot_tight <- ecdf_plot +
  theme(plot.margin = margin(t = 0, r = 5, b = 5, l = 5))

fig4_composite <- ggarrange(row1, ecdf_plot_tight,
                            nrow = 2, heights = c(1, 0.85))

ggsave(here("figures", "fig4_predictor_sets_accuracy.png"), fig4_composite,
       width = 11, height = 8, dpi = 300, bg = "white")
ggsave(here("figures", "fig4_predictor_sets_accuracy.pdf"), fig4_composite,
       width = 11, height = 8, bg = "white")
cat("Saved: figures/fig4_predictor_sets_accuracy.png\n")


ggplot(all.out,#  %>% filter(rank_only != "functional"),
			 aes(x = pretty_group,
			 		y = values,
			 		color = fcast_type)) +
	geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = pretty_group,
								 y = values,
								 color = fcast_type),
						 size=3, alpha = .2,
						 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = 0)) +
	facet_wrap(predictor~pretty_group, scales="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_blank(),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ylab("Variable importance") +
	ggtitle("Variables best explaining site random effects") +
	scale_color_discrete("Domain") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x= 1.5, label.y.npc = .75,
										 show.legend = F, hide.ns = T, size=8)

##### Misc checks #####



# FACET BY PREDICTOR, compare by MICROBIAL DOMAIN
all_vars_domain <- ggplot(all.out %>% filter(model_name == "env_cycl"),
													#  %>% filter(rank_only != "functional"),
													aes(x = pretty_group,
															y = values,
															color = pretty_group)) +
	#geom_violin(draw_quantiles = c(.5))+
	geom_boxplot(draw_quantiles = c(.5))+
	# geom_point(aes(x = pretty_group,
	# 							 y = values,
	# 							 color = pretty_group),
	# 					 size=3, alpha = .1,
	# 					 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = 0)) +
	facet_wrap(predictor_category~predictor, scales="free") +
	theme_pubclean() +
	theme(
		text = element_text(size = 16),
		axis.text.x=element_blank(),
		#axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ylab("Variable importance") +
	ggtitle("Variables best explaining site random effects") +
	scale_color_discrete(name = NULL) +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..),
										 label.x= 1.5,
										 #label.y = .75,
										 label.y.npc = .75,
										 show.legend = F, hide.ns = T, size=8) +
	scale_y_continuous(trans = pseudolog10_trans)
all_vars_domain

png(here("figures","site_effect_predictors.png"), width = 800, height=1000)
print(all_vars_domain)
dev.off()



# FACET BY PREDICTOR across taxonomic ranks
ggplot(all.out) + geom_point(aes(x = rank_only,
																 y = values,
																 color = pretty_group),
														 size=3, alpha = .2,
														 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(~predictor, scale="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ggtitle("Variables best explaining site random effects, across taxonomic ranks")



### Overall predictability of site effects ###
ggscatter(pred_sites %>% filter(model_name == "env_cycl") %>% filter(fcast_type != "Diversity"),
					x = "Mean", y = "pred",color = "pretty_group",
					add = "reg.line",                                 # Add regression line
					conf.int = TRUE,                                  # Add confidence interval
					add.params = list(color = "pretty_group", alpha=.3))+
	stat_cor(method = "pearson", label.x = -2,  p.digits = 2) +
	geom_abline(slope = 1, intercept = 0, linetype=2) +
	xlab("Observed") + ylab("Predicted") + theme(
		text = element_text(size = 16)) +
	#ggtitle("Fungal functional group site effects are \nless predictable from soil chemistry and climate") +
	facet_grid(pretty_group~fcast_type)

### Visualize absolute size of site effects ###
ggplot(data=pred_sites,
			 aes(x = rank_only,y = abs(Mean)))+
	geom_jitter(aes(color = pretty_group), size = 4, width=.2) +
	labs(col = "Site", title = "Absolute site effect size") +
	xlab("Rank")+ 	facet_grid(#rows = vars(only_rank),
		rows = vars(pretty_group), drop = T,
		scales = "free", space = "free_x") +
	ylab(NULL)+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	)


source("source.R")
library(tidytext)


site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged_env_cycl.rds"))

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))



# Read in predicted site effects
unobs_sites_env_cov <- readRDS(here("data/summary/site_effects_dredged_env_cov.rds"))[[1]]
unobs_sites_cycl_only <- readRDS(here("data/summary/site_effects_dredged_cycl_only.rds"))[[1]]
unobs_sites_env_cycl <- readRDS(here("data/summary/site_effects_dredged_env_cycl.rds"))[[1]]
obs_sites <- list(unobs_sites_env_cov, unobs_sites_cycl_only, unobs_sites_env_cycl) %>%
	rbindlist() %>% filter(model_id %in% converged) %>% as.data.frame()

obs_sites <- readRDS(here("data/summary/site_effects_dredged.rds"))[[4]]


obs_sites$fcast_type = ifelse(grepl("function", obs_sites$rank_only), "Functional group",
															 ifelse(grepl("diversity", obs_sites$taxon_rank), "Diversity", "Taxonomic group"))

# # Add zeros since unimportant predictors were dropped by dredge
# obs_sites$values <- as.numeric(obs_sites$values)
# obs_sites <- obs_sites %>% group_by(pretty_group, fcast_type, model_name, taxon_rank, rank_only) %>% complete(taxon, predictor, fill = list(values = 0))

# Define predictor category and name mappings first
pred_cat_key = list("nitrogenTot" = "macronutrient",
										"totalP" = "macronutrient",
										"sulfurTot" = "macronutrient",
										"estimatedOC" = "macronutrient",
										"pOxalate"= "macronutrient",
										"feOxalate" = "metal",
										"alOxalate" = "metal",
										"cecdNh4" = "cations",
										"mgNh4d" = "cations",
										"naNh4d" = "cations",
										"mnOxalate" = "cations",
										"kNh4d"  = "cations",
										"caNh4d" = "cations",
										"MAT"	     = "climate",
										"MAP"      = "climate",
										"latitude_scaled"      = "climate",
										"so4Satx"= "micronutrient",
										"siOxalate"  = "micronutrient")


pred_name_key = list("no2Satx" = "nitrite",
										 "totalP" = "total phosphorus",
										 "nitrogenTot" = "total nitrogen",
										 "estimatedOC" = "organic carbon",
										 "no3Satx" = "nitrate",
										 "pOxalate"= "phosphate",
										 "alOxalate"= "aluminum",
										 "feOxalate" = "iron",
										 "alKcl" = "aluminum",
										 "cecdNh4" = "cation exchange",
										 "mgNh4d" = "magnesium",
										 "naNh4d" = "sodium",
										 "mnOxalate" = "manganese",
										 "kNh4d"  = "potassium",
										 "caNh4d" = "calcium",
										 "MAT"	     = "mean annual temperature",
										 "MAP"      = "mean annual precipitation",
										 "latitude_scaled" = "latitude",
										 "so4Satx"= "sulfate",
										 "sulfurTot"= "sulfur",
										 "siOxalate"  = "silicon")

# Now apply the recoding
obs_sites$values = obs_sites$importance
obs_sites$predictor_category = recode(obs_sites$predictor, !!!pred_cat_key)
obs_sites$predictor = recode(obs_sites$predictor, !!!pred_name_key)



all.out <- obs_sites %>%
	group_by(predictor, model_name) %>%
	mutate(importance = round(mean(values), 2)) %>%
	ungroup() %>%
	mutate(predictor = factor(predictor),
				 predictor = fct_reorder(predictor, importance))

group_vals = obs_sites %>%
	group_by(pretty_group,model_name, predictor_category) %>%
	mutate(mean_importance=mean(values, na.rm=T),
				 sd_importance = sd(values, na.rm=T)) %>%
	mutate(ymax = mean_importance + 1.96*sd_importance,
				 ymin = mean_importance - 1.96*sd_importance) %>%
	arrange(pretty_group, predictor_category) %>% ungroup


pred_vals = obs_sites %>%
	group_by(pretty_group,model_name,predictor) %>%
	mutate(mean_importance=mean(values, na.rm=T),
				 sd_importance = sd(values, na.rm=T)) %>%
	mutate(ymax = mean_importance + 1.96*sd_importance,
				 ymin = mean_importance - 1.96*sd_importance) %>%
	arrange(pretty_group, predictor) %>% ungroup

# FACET BY PREDICTOR TYPE, compare by MICROBIAL DOMAIN
f_b_category =	ggplot(group_vals %>% filter(model_name %in% c("env_cycl")), #,"env_cov")),
											# %>% filter(rank_only != "functional"),
				 aes(x = reorder(predictor_category, -mean_importance),
				 		y = mean_importance,
				 		color = pretty_group)) +
		geom_pointrange(aes(ymin = ymin, ymax = ymax,
												color = pretty_group), position = position_dodge(width = .5)
										#show.legend = F
										) + ggtitle(NULL)  +
		# geom_point(data = obs_sites, aes(x = predictor_category,
		# 																 y = values,
		# 																 color = pretty_group), size=3, alpha = .1, position = position_jitterdodge(dodge.width = .5, jitter.width = .1, jitter.height = .1)) +
		#facet_wrap(~pretty_group, scales="free") +
		theme_bw() + theme(
			text = element_text(size = 20),
			axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
			axis.title=element_text(size=18)) + xlab("Predictor category") +
		ylab("Importance for \nexplaining site effects") +
		scale_color_discrete("Kingdom") +
		stat_compare_means(method = "t.test",
											 aes(y = values, label = ..p.signif..),
											  label.y = 1,
											 show.legend = F, hide.ns = T, size=8) + ylab(NULL)
#	facet_grid(rows=vars(model_name), labeller = labeller(model_name = model.labs))
	#	scale_y_sqrt()
f_b_category

png(here("figures","site_effect_f_b_category.png"), width = 800, height=1000)
print(f_b_category)
dev.off()


overall_importance =	ggplot(group_vals %>% filter(model_name %in% c("env_cycl")),#  %>% filter(rank_only != "functional"),
											aes(x = reorder(predictor, -values),
													y = values, color=pretty_group)) +
	#geom_point(position = position_jitter(width = .1, height=0.01), show.legend = F, alpha=.01) +
	ggtitle(NULL)  +
	#geom_violin(draw_quantiles = (.5)) +
	#facet_grid(rows=vars(pretty_group), scales="free") +
	theme_bw() + theme(
		text = element_text(size = 20),
		axis.text.x=element_text(angle = 60, hjust = 1, vjust = 1),
		axis.title=element_text(size=18)) + xlab("Predictor") +
	ylab("Importance for \nexplaining site effects") +
	stat_compare_means(method = "t.test",
										 aes(y = values, label = ..p.signif..),
										 label.y = 1.5,
										 show.legend = F, hide.ns = T, size=8)
overall_importance <- overall_importance %>% ggadd(c("mean_sd"), #color = "pretty_group",
																					 show.legend = F) + ylim(0.3,2.)
overall_importance

supp_fig <- ggarrange(overall_importance, f_b_category, labels = c("A","B"))

png(here("figures","site_effect_predictor_importance.png"), width = 1200, height=1000)
print(supp_fig)
dev.off()


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
	scale_color_discrete("Kingdom") +
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


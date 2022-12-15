source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
library(tidytext)


site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))

all.out <- site_eff_dredged_in[[1]]
pred_sites <- site_eff_dredged_in[[2]]

all.out <- all.out %>%
	group_by(predictor) %>%
	mutate(importance = round(mean(values), 2)) %>%
	ungroup() %>%
	#group_by(pretty_group, only_rank) %>%
	mutate(predictor = factor(predictor),
				 predictor = fct_reorder(predictor, importance))

# all.out$predictor_category = recode(all.out$predictor, "siOxalate"  = "micronutrient",
# 																		"no2Satx" = "macronutrient",
# 																		"totalP" = "macronutrient",
# 																		"no3Satx" = "macronutrient",
# 																		"cecdNh4" = "cations",
# 																		"pOxalate"= "macronutrient",
# 																		"so4Satx"= "macronutrient",
# 																		"naNh4d" = "micronutrient",
# 																		"mnOxalate" = "micronutrient",
# 																		"MAT"	     = "climate",
# 																		"feOxalate" = "micronutrient",
# 																		"MAP"      = "climate",
# 																		"kNh4d"  = "macronutrient",
# 																		"caNh4d" = "macronutrient")

all.out$predictor_category = recode(all.out$predictor, "siOxalate"  = "micronutrient",
			 "no2Satx" = "macronutrient",
			 "totalP" = "macronutrient",
			 "no3Satx" = "macronutrient",
			 "cecdNh4" = "cations",
			 "pOxalate"= "macronutrient",
			 "so4Satx"= "micronutrient",
			 "mgNh4d" = "micronutrient",
			 "naNh4d" = "salt",
			 "mnOxalate" = "salt",
			 "MAT"	     = "climate",
			 "feOxalate" = "metal",
			 "alKcl" = "metal",
			 "MAP"      = "climate",
			 "kNh4d"  = "salt",
			 "caNh4d" = "salt")

# FACET BY PREDICTOR TYPE, compare by MICROBIAL DOMAIN
	ggplot(all.out,#  %>% filter(rank_only != "functional"),
				 aes(x = pretty_group,
				 		y = values,
				 		color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = pretty_group,
								 y = values,
								 color = pretty_group),
						 size=3, alpha = .2,
						 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(~predictor_category, scales="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_blank(),
		#axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ylab("Variable importance") +
	ggtitle("Variables best explaining site random effects") +
	scale_color_discrete("Domain") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..), label.x = 1.5, label.y = .75, show.legend = F, hide.ns = T, size=5)



# FACET BY PREDICTOR, compare by MICROBIAL DOMAIN
ggplot(all.out,#  %>% filter(rank_only != "functional"),
			 aes(x = pretty_group,
			 		y = values,
																																	color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5))+
	geom_point(aes(x = pretty_group,
								 y = values,
								 color = pretty_group),
						 size=3, alpha = .2,
						 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(predictor_category~predictor, scales="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_blank(),
		#axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ylab("Variable importance") +
	ggtitle("Variables best explaining site random effects") +
	scale_color_discrete("Domain") +
	stat_compare_means(method = "t.test", aes(label = ..p.signif..), label.x = 1.5, label.y = .75, show.legend = F, hide.ns = T, size=5)



# FACET BY PREDICTOR
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


ggplot(all.out %>% filter(rank_only != "diversity"), aes(x = rank_only,
										y = values,
										color = pretty_group)) + geom_point(
											size=3, alpha = .2,
											position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(~predictor_category, scale="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ggtitle("Variables best explaining site random effects, across taxonomic ranks") +
	geom_smooth(aes(x = as.numeric(rank_only)), span=.7)


pred_sites$fcast_type = ifelse(grepl("function", pred_sites$rank_only), "Functional group",
															 ifelse(grepl("diversity", pred_sites$taxon_rank), "Diversity", "Taxonomic group"))
### Overall predictability of site effects ###
ggscatter(pred_sites %>% filter(fcast_type != "Diversity"), x = "Mean", y = "pred",color = "pretty_group",
					add = "reg.line",                                 # Add regression line
					conf.int = TRUE,                                  # Add confidence interval
					add.params = list(color = "pretty_group"))+
	stat_cor(method = "pearson", label.x = -2,  p.digits = 2) +
	geom_abline(slope = 1, intercept = 0, linetype=2) +
	xlab("Observed") + ylab("Predicted") + theme(
		text = element_text(size = 16)) +
	ggtitle("Fungal functional group site effects are \nless predictable from soil chemistry and climate") + facet_grid(pretty_group~fcast_type)



#####
# older


# FACET BY PREDICTOR
all.out$rank_only  <- ordered(all.out$rank_only , levels = c("genus",
																								 "family",
																								 "order",
																								 "class",
																								 "phylum", "functional", "diversity"))
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

# FACET BY RANK
ggplot(all.out) + geom_point(aes(x = reorder_within(predictor, values, rank_only),
																 y = values,
																 color = pretty_group),
														 size=3, alpha = .2,
														 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(~rank_only, scale="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ggtitle("Variables best explaining site random effects")


all.out_filter = all.out %>% filter(taxon %in% pass_filter$taxon)


# FACET BY PREDICTOR and TYPE
ggplot(all.out  %>% filter(rank_only != "functional"), aes(x = pretty_group,
																																	y = values,
																																	color = pretty_group)) +
	geom_point(aes(x = rank_only,
																			y = values,
																			color = pretty_group),
																	size=3, alpha = .2,
																	position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
			 	facet_wrap(predictor_category~predictor, scale="free") +
			 	theme_bw() + theme(
			 		text = element_text(size = 16),
			 		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
			 		axis.title=element_text(size=22,face="bold")) + xlab(NULL) +
	ylab("Variable importance") +
			 	ggtitle("Variables best explaining site random effects")





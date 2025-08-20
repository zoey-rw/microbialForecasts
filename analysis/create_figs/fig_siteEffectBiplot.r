source("source.R")
pacman::p_load(ggforce, ggrepel)

#https://stackoverflow.com/questions/39137287/plotting-partial-least-squares-regression-plsr-biplot-with-ggplot2

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


read_in = readRDS(here("data/summary/site_effects_dredged.rds"))
component_df = read_in[[4]][[660]]$plsr_scores
# 5 is ok, 9 is interesting, 11 is clearest so far, maybe 45, 46, 47
# 660
# my fav is 643
component_df$predictor_category = recode(rownames(component_df), !!!pred_cat_key)
rownames(component_df) = recode(rownames(component_df), !!!pred_name_key)
p1 <- ggplot(data=component_df, aes(x = `Comp 1`,y = `Comp 2`))+
	ylab("Component 2")+
	xlab("Component 1")+
	ggtitle("Site effects, partial least squares regression", "model_id: env_cycl_saprotroph_20151101_20180101")+
	theme_bw(base_size = 18) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.grid.major = element_blank(),
										 panel.grid.minor = element_blank(),
										 axis.line = element_line(colour = "black"))+
	#geom_text(aes(label=rownames(component_df), color=predictor_category))+
	geom_text_repel(aes(label=rownames(component_df), color=predictor_category), max.overlaps = 20, size=5)+
	coord_fixed(ylim=c(-1, 1),xlim=c(-1,1))+
	theme(axis.ticks = element_line(colour = "red")) +
	geom_circle(aes(x0=0, y0=0, r=sqrt(1/2)), inherit.aes=FALSE) +
	geom_circle(aes(x0=0, y0=0, r=1), inherit.aes=FALSE)  + labs(color = "Predictor type")
p2 <- ggplot(data=component_df, aes(x = `Comp 3`,y = `Comp 4`))+
	ylab("Component 4")+
	xlab("Component 3")+
	ggtitle("","")+
	#ggtitle("")+
	theme_bw(base_size = 18) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
										 panel.grid.major = element_blank(),
										 panel.grid.minor = element_blank(),
										 axis.line = element_line(colour = "black"))+
#	geom_text(aes(label=rownames(component_df), color=predictor_category))+
	geom_text_repel(aes(label=rownames(component_df), color=predictor_category), max.overlaps = 20, size=5)+

	coord_fixed(ylim=c(-1, 1),xlim=c(-1,1))+
	theme(axis.ticks = element_line(colour = "red")) +
	geom_circle(aes(x0=0, y0=0, r=sqrt(1/2)), inherit.aes=FALSE) +
	geom_circle(aes(x0=0, y0=0, r=1), inherit.aes=FALSE) + labs(color = "Predictor type")
ggarrange(p1, p2, common.legend = T)


scoreplot(read_in[[4]][[660]], comps = 4, labels="names")
library(pls)
plot(read_in[[4]][[660]], "biplot")


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


pred_sites <- read_in[[2]]
pred_sites_plsr <- read_in[[3]]
# sanity check!
summary(lm(pred_sites_plsr$Median ~ pred_sites_plsr$plsr_pred_observed_4comp))
plot(pred_sites_plsr$Median ~ pred_sites_plsr$plsr_pred_observed_4comp); abline(0,1)
ggplot(pred_sites_plsr, aes(x = Median, y = plsr_pred_observed_4comp)) +
	geom_point() +
	facet_grid(~model_name)  +
	stat_smooth(method = "lm") +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))	)
ggplot(pred_sites, aes(x = Median, y = pred)) +
	geom_point() +
	facet_grid(~model_name)  +
	stat_smooth(method = "lm") +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))	)


summary(lm(pred_sites$Median ~ pred_sites$pred))
plot(pred_sites$Median ~ pred_sites$pred); abline(0,1)

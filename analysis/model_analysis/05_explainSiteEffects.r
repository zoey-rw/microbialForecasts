# Determine important predictors of site effects (random effects)
# Creates the "site_effects_dredged.rds" files
# Has to be run before creating any taxonomic forecasts (03a_hindcast_single_taxon)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn)

# Read in site effect estimates from model
site_effects_all <- readRDS(here("data/summary/site_effects.rds"))
site_effects_refit <- site_effects_all %>% filter(time_period == "2015-11_2020-01" | time_period == "2016-01_2018-01"  &
																										model_name == "all_covariates")

site_effects_calibration <- site_effects_all %>% filter(time_period == "2015-11_2018-01" &
																										model_name == "all_covariates")
site_effects_calibration <- site_effects_calibration %>% mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))

df_predictors <- readRDS(here("data/summary/site_effect_predictors.rds"))

# Site effects for explaining
df_merged <- site_effects_calibration %>%
	select(siteID,Mean,rank,rank_only,group,taxon) %>%
	merge(df_predictors, all=T)

# Unobserved sites for predicting
df_unobserved <- df_predictors %>% filter(!siteID %in% site_effects_calibration$siteID)

by_site <- df_merged %>%
	pivot_wider(id_cols = c("rank","taxon"), names_from = "siteID", values_from = "Mean") %>%
	select(-c(1:2))


# Loop through each taxon and see which factors best explain its site effects
options(na.action = "na.fail")

dredge.out <- list()
weird <- list()

# TODO: ADD functional groups back in once estimates are avail
#tax_only <- site_effects_refit %>% filter(fcast_type=="Taxonomic")

tax_list <- unique(site_effects_calibration$taxon)

#tax_list=tax_list[223:302]
#tax_list <- tax_list[!tax_list %in% c("other",NA)]#,"proteobacteria")]
s = "cyanobacteriia"
s = "other_order_bac_lachnospirales"
#for (s in tax_list){
pacman::p_load(doParallel)
cl <- makeCluster(27, type="FORK", outfile="")
registerDoParallel(cl)

clusterExport(cl,varlist=list("site_effects_calibration", "df_unobserved", "df_merged"),envir = environment())


#Run for multiple groups, in parallel (via PSOCK)
# Why does this work with "do" but not "dopar"!!!!!
dredge.out = foreach(s=tax_list, .errorhandling = 'remove') %dopar% {
	pacman::p_load(dplyr, MuMIn)
	options(na.action = "na.fail")

		taxon_eff <- site_effects_calibration %>% filter(taxon == !!s)

	taxon_rank <- taxon_eff %>% select(rank) %>% unique() %>% unlist()
	#species_dat <- df_merged %>% filter(taxon==s) %>% select(-c(1,3:6))
	#	species_dat <- df_merged %>% filter(taxon==s)  %>% select_if(is.numeric)
	#fm <-lm(Mean ~ . -ecoregion -nlcd, data = species_dat)

	pred_list = c("MAT", "MAP",
								"caNh4d", "kNh4d", "mgNh4d",
								"naNh4d",
								"cecdNh4", "alKcl",
								"feOxalate", "mnOxalate", "pOxalate", "siOxalate",
								#"estimatedOC",
								"no2Satx", "no3Satx", "so4Satx",
								"totalP")
	species_dat <- as.data.frame(df_merged) %>% filter(taxon==!!s | is.na(taxon))  %>% select(c("Mean","siteID", !!pred_list))
	species_dat <- species_dat[complete.cases(species_dat),]
	siteID_vec <- species_dat$siteID #unique(df_predictors$siteID)


	species_dat$siteID <- NULL
	#corrplot::corrplot(cor(species_dat), type="upper")
	fm <-lm(Mean ~ ., data = species_dat)
	#print(summary(fm))

	# Up to 5 predictors allowed for explaining each taxon's site effect
	#models <- lapply(MuMIn:::.dredge.par(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)
	#models <- lapply(dredge(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)

	temp <- MuMIn:::.dredge.par(fm, m.lim = c(1,6))

	# Average top models and get importance
	models <- get.models(temp,  subset = delta <= 3)
	ma <- model.avg(models)
	imp1 <- cbind.data.frame(predictor = names(ma$sw),
													 values = ma$sw)
	attr(imp1$values, "n.models") <- NULL
	imp <- imp1 %>% mutate(taxon_rank = taxon_rank,
												 taxon = !!s)
	rownames(imp) <- NULL

	# Subset prediction dataset to the important predictors
	#df_unobserved_taxon <- df_unobserved %>% select(siteID, !!names(ma$sw)) %>% na.omit
	df_unobserved_taxon <- df_unobserved %>% select(siteID, !!pred_list) %>% na.omit
	Weights(ma) <- cos2Weights(models, data = df_unobserved_taxon)
	# Predict using covariates for sites not in calibration
	df_unobserved_taxon <- df_unobserved_taxon %>% mutate(taxon_rank = taxon_rank, taxon = !!s,
																												pred = predict(ma,
																																			 newdata = df_unobserved_taxon))

	# Predict using covariates for sites in calibration
	modeled = species_dat %>% mutate(taxon_rank = taxon_rank,
																	 taxon = !!s, siteID = siteID_vec,
																	 pred = predict(ma,
																	 							 newdata = species_dat))


	# Merge with fitted values once sites are observed
	# newsites_observed <- site_effects_refit %>% filter(taxon == !!s) %>% select(siteID, Mean, SD, model_name)
	# predicted <- merge(df_unobserved_taxon, newsites_observed, all.x=T, by=c("siteID"))

	# Just testing whether the refit site effects estimates are similar to the calibration (they are)
	# newsites_observed$refit_Mean <- newsites_observed$Mean
	# modeled <- merge(modeled, newsites_observed[,c("siteID", "refit_Mean")])
	#plot(predicted$Mean ~ predicted$pred)
	#plot(modeled$Mean ~ modeled$pred)

	# dredged <- as.data.frame(avg_results$coefmat.subset) %>%
	# 	rownames_to_column("predictor") %>% mutate(taxon = !!s,
	# 																						 taxon_rank = taxon_rank)

	out <- list(imp = imp,
							unobserved = df_unobserved_taxon,
							modeled = modeled)
	return(out)
}

dredge.out
dredged_sites <- map_df(dredge.out, 1) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))
unobs_sites <- map_df(dredge.out, 2) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))
pred_sites <- map_df(dredge.out, 3) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))


site_effects_refit <- site_effects_refit %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))
#plot(unobs_sites$Mean ~ unobs_sites$pred)
#plot(pred_sites$Mean ~ pred_sites$pred)

all.out <- left_join(dredged_sites, unique(site_effects_calibration[,c("taxon","rank_only","pretty_group")]))
unobs_sites <- left_join(unobs_sites, unique(site_effects_calibration[,c("taxon","rank_only","pretty_group")]))

pred_sites <- left_join(pred_sites, unique(site_effects_calibration[,c("taxon","rank_only","pretty_group")]))
#all.out$`effect size` <- abs(all.out$Estimate)

# all.out_sorted <- all.out %>% #group_by(pretty_group) %>%
# 	arrange(`effect size`, predictor)

saveRDS(list(all.out, pred_sites), here("data/summary/site_effects_dredged.rds"))
saveRDS(unobs_sites, here("data/summary/site_effects_unobserved.rds"))






site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))

all.out <- site_eff_dredged_in[[1]]
pred_sites <- site_eff_dredged_in[[2]]

### Overall predictability at unobserved sites ###
ggscatter(unobs_sites, x = "Mean", y = "pred",color = "pretty_group",
					add = "reg.line",                                 # Add regression line
					conf.int = TRUE,                                  # Add confidence interval
					add.params = list(color = "pretty_group"))+
	stat_cor(method = "pearson", label.x = -2,  p.digits = 2) +
	geom_abline(slope = 1, intercept = 0, linetype=2) +
	xlab("Observed") + ylab("Predicted") + theme(
		text = element_text(size = 16)) +
	ggtitle("Bacterial site effects are predictable from soil chemistry and climate") + facet_grid(~pretty_group)




all.out <- all.out %>%
	group_by(predictor) %>%
	mutate(importance = round(mean(values), 2)) %>%
	ungroup() %>%
	#group_by(pretty_group, only_rank) %>%
	mutate(predictor = factor(predictor),
				 predictor = fct_reorder(predictor, importance))

all.out$predictor_category_nutrient = recode(all.out$predictor, "siOxalate"  = "micronutrient",
																								 "no2Satx" = "macronutrient",
																								 "totalP" = "macronutrient",
																								 "no3Satx" = "macronutrient",
																								 "cecdNh4" = "cations",
																								 "pOxalate"= "macronutrient",
																								 "so4Satx"= "micronutrient",
																								 "naNh4d" = "micronutrient",
																								 "mnOxalate" = "micronutrient",
																								 "MAT"	     = "climate",
																								 "feOxalate" = "micronutrient",
																								 "MAP"      = "climate",
																								 "kNh4d"  = "micronutrient",
																								 "caNh4d" = "micronutrient")


all.out$predictor_category_other = recode(all.out$predictor,
																						 "cecdNh4" = "cations",
																						 "MAT"	     = "climate",
																						 "MAP"      = "climate")


all.out$predictor_category_metal = recode(all.out$predictor, "siOxalate"  = "micronutrient",
																						 "no2Satx" = "macronutrient",
																						 "totalP" = "macronutrient",
																						 "no3Satx" = "macronutrient",
																						 "cecdNh4" = "cations",
																						 "pOxalate"= "macronutrient",
																						 "so4Satx"= "micronutrient",
																						 "naNh4d" = "salt",
																						 "mnOxalate" = "salt",
																						 "MAT"	     = "climate",
																						 "feOxalate" = "metal",
																					"alKcl" = "metal",
																						 "MAP"      = "climate",
																						 "kNh4d"  = "salt",
																						 "caNh4d" = "salt")


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

# FACET BY RANK
library(tidytext)
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


library(ggExtra)
ggMarginal(p, type = "density")



site_effects_all <- readRDS(here("data/summary/site_effects.rds"))
site_effects_refit_only <- site_effects_all %>% filter(time_period == "2015-11_2020-01" &
																										model_name == "all_covariates")
site_effects_refit_only <- site_effects_refit_only %>%
	filter(!site_effects_refit_only$siteID %in% site_effects_refit$siteID) %>%
	rename(obs_site_effect = Mean)

predicting_site_effects <- merge(unobs_sites, site_effects_refit_only)


ggscatter(pred_sites, x = "Mean", y = "pred",color = "only_rank",
					add = "reg.line",                                 # Add regression line
					conf.int = TRUE,                                  # Add confidence interval
					add.params = list(color = "only_rank"))+
	stat_cor(method = "pearson", label.x = -2, label.y = 1, p.digits = 2) +
	geom_abline(slope = 1, intercept = 0, linetype=2) +
	xlab("Observed") + ylab("Predicted") + theme(
		text = element_text(size = 16)) + ggtitle("Site effects are predictable from soil chemistry and climate")







all.out[all.out$taxon_rank=="functional_group",]

ggplot(all.out %>% filter(predictor != "(Intercept)")) + geom_point(aes(x = only_rank,
																 y = `effect size`,
																 fill = pretty_group),
														 size=3, alpha = .5,
														 shape=21,
														 position=position_jitterdodge(dodge.width = 1, jitter.width = .3, jitter.height = 0)) +
	facet_wrap(~predictor, scale="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")) + xlab(NULL) + ggtitle("Variables best explaining site random effects")


# Visualize absolute size of site effects
ggplot(data=site_effects_refit,
			 aes(x = only_rank,y = abs(Mean)))+
	geom_jitter(aes(color = as.factor(siteID)), size = 4, width=.2) +
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



stopCluster(cl)

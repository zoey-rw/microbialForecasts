# Determine important predictors of site effects (random effects)
# Creates the "site_effect_predictors.rds" and "site_effects_dredged.rds" files

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn)

# Read in site effect estimates from model
site_effects_all <- readRDS(here("data/summary/site_effects.rds"))
site_effects_refit <- site_effects_all %>% filter(time_period == "2015-11_2020-01" | time_period == "2016-01_2018-01"  &
																										model_name == "all_covariates")

## Get site climate data from NEON
fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
fieldsites <- fieldsites_raw %>% filter(!grepl("Aquatic", field_site_type)) %>%
	select(siteID = field_site_id,
				 MAT = field_mean_annual_temperature_C,
				 MAP = field_mean_annual_precipitation_mm,
				 nlcd = field_dominant_nlcd_classes,
				 #ecoregion = field_domain_id)
				 ecoregion = field_domain_id)

fieldsites$evergreen   <- ifelse(grepl("Evergreen", fieldsites$nlcd), "Evergreen", "")
fieldsites$deciduous   <- ifelse(grepl("Deciduous", fieldsites$nlcd), "Deciduous", "")
fieldsites$wetlands 	 <- ifelse(grepl("Wetlands", fieldsites$nlcd), "Wetlands", "")
fieldsites$herbaceous <- ifelse(grepl("Herbaceous", fieldsites$nlcd), "Herbaceous", "")
fieldsites$shrub      <- ifelse(grepl("Grassland|Shrub", fieldsites$nlcd), "Shrub", "")
fieldsites$crops      <- ifelse(grepl("Crops", fieldsites$nlcd), "Crops", "")
fieldsites$nlcd_clean <- apply(fieldsites[,c(6:11) ], 1, paste0, collapse = "")

# fieldsites$evergreen   <- ifelse(grepl("Evergreen", fieldsites$nlcd), 1, 0)
# fieldsites$deciduous   <- ifelse(grepl("Deciduous", fieldsites$nlcd), 1, 0)
# fieldsites$wetlands 	 <- ifelse(grepl("Wetlands", fieldsites$nlcd), 1, 0)
# fieldsites$herbaceous <- ifelse(grepl("Crops", fieldsites$nlcd), 1, 0)
# fieldsites$shrub      <- ifelse(grepl("Grassland|Shrub", fieldsites$nlcd), 1, 0)
# fieldsites$crops      <- ifelse(grepl("Crops", fieldsites$nlcd), 1, 0)
#
# classes <- expand(fieldsites, nesting(evergreen,
# 																			deciduous,
# 																			wetlands,
# 																			herbaceous,
# 																			shrub,
# 																			crops))
# classes1 <- classes %>%
# 	mutate(across(1:6, function(x)
# 		factor(x, c(0, 1),  c("", cur_column()))))
# classes1 <- lapply(classes1[,c(1:6)], as.character) %>% as.data.frame()
# classes1$cat  <- apply(classes1[,c(1:6) ], 1, paste0, collapse = "")
#
#
# fieldsites$nlcd_clean <- ifelse(grepl("Crops", fieldsites$nlcd), "Crops", fieldsites$nlcd)
# fieldsites$nlcd_clean <- ifelse(grepl("Evergreen Forest", fieldsites$nlcd), "Evergreen Forest", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Deciduous|Mixed", fieldsites$nlcd), "Deciduous/Mixed Forest", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Grassland|Scrub", fieldsites$nlcd), "Grassland/Scrub", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Wetlands", fieldsites$nlcd), "Wetlands", fieldsites$nlcd_clean)


### Read in bulk density data
bd_raw <- readRDS(here("data/raw/bulkDensity_allsites.rds"))
bd <- bd_raw %>% filter(bulkDensTopDepth < 25) %>%
	group_by(siteID) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) #%>%
bd$bulkDens <- ifelse(is.na(bd$bulkDensThirdBar), bd$bulkDensFieldMoist, bd$bulkDensThirdBar)
bd <- bd %>% select(siteID, bulkDens)

### Read in micronutrient data
bgc_raw <- readRDS(here("data/raw/dp1.10047.rds"))
# remove deep horizons.
bgc <- bgc_raw %>% filter(biogeoTopDepth < 25)
most_relevant <- c("siteID", "plotID",
									 #"biogeoTopDepth", "biogeoBottomDepth",
									 #"gypsumConc",
									 #"caco3Conc",
									 "caNh4d",
									 #"caSatx",
									 "kNh4d",
									 #"kSatx",
									 "mgNh4d", #"mgSatx",
									 "naNh4d",
									 #"naSatx",
									 "cecdNh4",
									 #"alSatCecd33",
									 "alKcl",
									 "feKcl", "feOxalate",
									 #"mnKcl",
									 "mnOxalate",
									 "pOxalate", "siOxalate",
									 "estimatedOC", "sulfurTot",
									 #"pSatx",
									 #"ececCecd33",
									 "no2Satx", "no3Satx", "so4Satx",
									 "OlsenPExtractable", "MehlichIIITotP", "Bray1PExtractable")
bgc_relevant <- bgc %>% select(!!!most_relevant)
bgc_relevant$totalP <- ifelse(is.na(bgc_relevant$MehlichIIITotP),
															bgc_relevant$OlsenPExtractable,
															bgc_relevant$MehlichIIITotP)
bgc_relevant[,c("MehlichIIITotP","OlsenPExtractable","Bray1PExtractable")] <-NULL
bgc_relevant <- bgc_relevant %>% group_by(siteID) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# Merge predictor dfs together
df_predictors <- fieldsites %>%
	merge(bgc_relevant, all=T) %>%
	merge(bd, all=T) %>%
	# Center and scale predictor values
	mutate(across(where(is.numeric), function(x) as.vector(scale(x))))

df_predictors <- merge(df_predictors, fieldsites_raw %>% select(siteID = field_site_id, latitude = field_latitude))

saveRDS(df_predictors, here("data/summary/site_effect_predictors.rds"))

df_predictors <- readRDS(here("data/summary/site_effect_predictors.rds"))

# Site effects for explaining
df_merged <- site_effects_refit %>%
	select(siteID,Mean,rank,rank_only,group,taxon) %>%
	merge(df_predictors)

# Unobserved sites for predicting
df_unobserved <- df_predictors %>% filter(!siteID %in% site_effects_refit$siteID)

by_site <- df_merged %>%
	pivot_wider(id_cols = c("rank","taxon"), names_from = "siteID", values_from = "Mean") %>%
	select(-c(1:2))


# Loop through each taxon and see which factors best explain its site effects
options(na.action = "na.fail")

dredge.out <- list()
weird <- list()
tax_list <- unique(site_effects_refit$taxon)
tax_list <- tax_list[!tax_list %in% c("other",NA)]#,"proteobacteria")]
s = "cyanobacteriia"
#for (s in tax_list){
pacman::p_load(doParallel)
cl <- makeCluster(28, type="FORK", outfile="")
registerDoParallel(cl)

clusterExport(cl,varlist=list("site_effects_refit", "df_unobserved", "df_merged"),envir = environment())


#Run for multiple groups, in parallel (via PSOCK)
# Why does this work with "do" but not "dopar"!!!!!
dredge.out = foreach(s=tax_list, .errorhandling = 'remove') %dopar% {
	pacman::p_load(dplyr, MuMIn)
	options(na.action = "na.fail")

		taxon_eff <- site_effects_refit %>% filter(taxon == !!s)
	df_unobserved_taxon <- df_unobserved[complete.cases(df_unobserved),]

	taxon_rank <- taxon_eff %>% select(rank_name) %>% unique() %>% unlist()
	#species_dat <- df_merged %>% filter(taxon==s) %>% select(-c(1,3:6))
	#	species_dat <- df_merged %>% filter(taxon==s)  %>% select_if(is.numeric)
	#fm <-lm(Mean ~ . -ecoregion -nlcd, data = species_dat)

	species_dat <- as.data.frame(df_merged) %>% filter(taxon==!!s)  %>% select(c("Mean", "MAT", "MAP",
																															"caNh4d", "kNh4d", #"mgNh4d",
																															"naNh4d",
																															"cecdNh4", #"alKcl",
																															"feOxalate", "mnOxalate", "pOxalate", "siOxalate",
																															#"estimatedOC",
																															"no2Satx", "no3Satx", "so4Satx",
																															"totalP","siteID"))

	species_dat <- species_dat[complete.cases(species_dat),]
	siteID_vec <- species_dat$siteID
	species_dat$siteID <- NULL
	#corrplot::corrplot(cor(species_dat), type="upper")
	fm <-lm(Mean ~ ., data = species_dat)
	#print(summary(fm))

	# Up to 5 predictors allowed for explaining each taxon's site effect
	#models <- lapply(MuMIn:::.dredge.par(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)
	#models <- lapply(dredge(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)

	temp <- MuMIn:::.dredge.par(fm, m.lim = c(1,5))
	models <- get.models(temp,  subset = delta <= 2)

	#if (length(models = 1))
	ma <- model.avg(models)
	imp1 <- cbind.data.frame(predictor = names(ma$sw), values = ma$sw)
	attr(imp1$values, "n.models") <- NULL
	imp <- imp1 %>% mutate(taxon_rank = taxon_rank, taxon = !!s)
	rownames(imp) <- NULL
	# imp <- ma$sw %>%
	# 	rownames_to_column("predictor") %>% as_tibble()   %>% as.data.frame()
	# attr(MyData$VAR, "ATT") <- NULL

	# imp <- cbind.data.frame(importance = sw(ma)) %>% as.data.frame() %>%
	# 	rownames_to_column("predictor") %>% mutate(taxon_rank = taxon_rank,
	# 																						 taxon = !!s)
	#avg_results <- summary(ma)

	Weights(ma) <- cos2Weights(models, data = df_unobserved_taxon)
	#df_unobserved_taxon$pred <- predict(ma, newdata = df_unobserved_taxon)
	df_unobserved_taxon <- df_unobserved_taxon %>% mutate(taxon_rank = taxon_rank, taxon = !!s,
																												pred = predict(ma,
																																			 newdata = df_unobserved_taxon))
	modeled = species_dat %>% mutate(taxon_rank = taxon_rank,
																	 taxon = !!s, siteID = siteID_vec,
											 pred = predict(ma,
											 							 newdata = species_dat))

	# dredged <- as.data.frame(avg_results$coefmat.subset) %>%
	# 	rownames_to_column("predictor") %>% mutate(taxon = !!s,
	# 																						 taxon_rank = taxon_rank)

	out <- list(imp = imp,
							unobserved = df_unobserved_taxon,
							modeled = modeled)
	return(out)
}

dredge.out
dredged_sites <- map_df(dredge.out, 1)
unobs_sites <- map_df(dredge.out, 2)
pred_sites <- map_df(dredge.out, 3)


all.out <- merge(dredged_sites, unique(site_effects_refit[,c("taxon","only_rank","pretty_group")]))
unobs_sites <- merge(unobs_sites, unique(site_effects_refit[,c("taxon","only_rank","pretty_group")]))

pred_sites <- merge(pred_sites, unique(site_effects_refit[,c("taxon","only_rank","pretty_group")]))
#all.out$`effect size` <- abs(all.out$Estimate)

# all.out_sorted <- all.out %>% #group_by(pretty_group) %>%
# 	arrange(`effect size`, predictor)

saveRDS(list(all.out, pred_sites), here("data/summary/site_effects_dredged.rds"))
saveRDS(unobs_sites, here("data/summary/site_effects_unobserved.rds"))






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

all.out$predictor_category = recode(all.out$predictor, "siOxalate"  = "micronutrient",
																								 "no2Satx" = "macronutrient",
																								 "totalP" = "macronutrient",
																								 "no3Satx" = "macronutrient",
																								 "cecdNh4" = "cations",
																								 "pOxalate"= "macronutrient",
																								 "so4Satx"= "macronutrient",
																								 "naNh4d" = "micronutrient",
																								 "mnOxalate" = "micronutrient",
																								 "MAT"	     = "climate",
																								 "feOxalate" = "micronutrient",
																								 "MAP"      = "climate",
																								 "kNh4d"  = "macronutrient",
																								 "caNh4d" = "macronutrient")

# FACET BY PREDICTOR
ggplot(all.out) + geom_point(aes(x = only_rank,
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
ggplot(all.out) + geom_point(aes(x = reorder_within(predictor, values, only_rank),
																 y = values,
																 color = pretty_group),
														 size=3, alpha = .2,
														 position=position_jitterdodge(dodge.width = 1, jitter.width = .1, jitter.height = .1)) +
	facet_wrap(~only_rank, scale="free") +
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




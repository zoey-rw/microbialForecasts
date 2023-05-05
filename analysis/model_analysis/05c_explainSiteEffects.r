# Determine important predictors of site effects (random effects) for env_cycl forecasts (instead of adding another loop to this script... just duplicating it. Lazy but it works)
# Creates the "site_effects_dredged_env_cycl.rds" files
# Has to be run before creating any taxonomic forecasts (06_createHindcasts)

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn)


# Created by running "./data_construction/covariate_prep/combine_site_descriptors.r"
# Requires extra NEON data downloads; use existing file if possible!
df_predictors <- readRDS(here("data/clean/site_effect_predictors.rds"))


# Read in site effect estimates from model; split into calibration timeperiod and full dataset
site_effects_all <- readRDS(here("data/summary/site_effects.rds"))
site_effects_refit <- site_effects_all %>% filter(time_period == "2015-11_2020-01" &
																										model_name == "env_cycl")
site_effects_calibration <- site_effects_all %>% filter(time_period == "2015-11_2018-01" &
																													model_name == "env_cycl")
site_effects_calibration <- site_effects_calibration %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))


# Observed site effects for explaining
df_observed <- site_effects_calibration %>%
	select(siteID,Mean,rank,rank_only,group,taxon) %>%
	merge(df_predictors, all.x=T)
by_site <- df_observed %>%
	pivot_wider(id_cols = c("rank","taxon"),
							names_from = "siteID", values_from = "Mean") %>%
	select(-c(1:2))

# Unobserved sites for predicting
df_unobserved <- df_predictors %>%
	filter(!siteID %in% site_effects_calibration$siteID)



# Loop through each taxon and see which factors best explain its site effects
options(na.action = "na.fail")

dredge.out <- list()
weird <- list()

# TODO: ADD functional groups back in once estimates are avail
#tax_only <- site_effects_refit %>% filter(fcast_type=="Taxonomic")

tax_list <- unique(site_effects_calibration$taxon)

s = "acidobacteriota"
s = "cyanobacteriia"
s = "other_order_bac_lachnospirales"
#for (s in tax_list){
pacman::p_load(doParallel)
cl <- makeCluster(27, type="FORK", outfile="")
registerDoParallel(cl)

clusterExport(cl,varlist=list("df_unobserved", "df_observed"),envir = environment())


#Run for multiple groups, in parallel
dredge.out = foreach(s=tax_list, .errorhandling = 'remove') %dopar% {
	pacman::p_load(dplyr, MuMIn)
	options(na.action = "na.fail")

	pred_list = c("MAT", "MAP",
								"caNh4d", "kNh4d", "mgNh4d",
								"naNh4d",
								"cecdNh4", "alKcl",#"latitude",
								"feOxalate", "mnOxalate", "pOxalate", "siOxalate",
								#"estimatedOC",
								"no2Satx", "no3Satx", "so4Satx",
								"totalP")
	species_dat <- as.data.frame(df_observed) %>%
		filter(taxon==!!s | is.na(taxon))  %>%
		select(c("Mean","siteID", !!pred_list))
	species_dat <- species_dat[complete.cases(species_dat),]
	siteID_vec <- species_dat$siteID

	taxon_rank <- df_observed %>% filter(taxon == !!s) %>% select(rank) %>% unique() %>% unlist()


	species_dat$siteID <- NULL
	#corrplot::corrplot(cor(species_dat), type="upper")
	fm <-lm(Mean ~ ., data = species_dat)
	#print(summary(fm))

	# Up to 6 predictors allowed for explaining each taxon's site effect
	#models <- lapply(MuMIn:::.dredge.par(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)
	#models <- lapply(dredge(fm, evaluate = FALSE, m.lim = c(1,5), cluster = cl), eval)

	temp <- MuMIn:::.dredge.par(fm, m.lim = c(1,6), extra="adjR^2", beta="none")
	
	# Average top models and get importance
	models <- get.models(temp,  subset = delta <= 3)
	ma <- model.avg(models)
	predictor_importance <- cbind.data.frame(predictor = names(ma$sw),
																					 values = ma$sw) %>% 
		mutate(taxon_rank = taxon_rank,
					 taxon = !!s)
	rownames(predictor_importance) <- NULL
	attr(predictor_importance$values, "n.models") <- NULL
	
	
	# Subset unobserved prediction dataset to the important predictors
	df_unobserved_taxon <- df_unobserved %>% select(siteID, !!pred_list) %>% 
		na.omit %>% 
		mutate(taxon_rank = taxon_rank, taxon = !!s)
	Weights(ma) <- cos2Weights(models, data = df_unobserved_taxon)
	
	# Predict using covariates for unobserved sites (not in calibration)
	pred = predict(ma, newdata = df_unobserved_taxon, se.fit=T)
	new_site_effects = cbind.data.frame(df_unobserved_taxon, pred)
	sample_size = nrow(df_unobserved_taxon)
	new_site_effects$sd.fit = new_site_effects$se.fit * sqrt(sample_size)
	new_site_effects$ci_lo = new_site_effects$fit - (new_site_effects$se.fit * 1.96)
	new_site_effects$ci_hi = new_site_effects$fit + (new_site_effects$se.fit * 1.96)
	
	# Predict using covariates for sites in calibration (checking in-sample model accuracy)
	modeled = species_dat %>% mutate(taxon_rank = taxon_rank,
																	 taxon = !!s, siteID = siteID_vec,
																	 pred = predict(ma,
																	 							 newdata = species_dat))
	summary(lm(modeled$Mean ~ modeled$pred))
	
	out <- list(predictor_importance = predictor_importance,
							new_site_effects = new_site_effects,
							modeled_values = modeled)
	return(out)
}

# Fix some taxon names and recombine dfs
dredged_predictor_importance <- map_df(dredge.out, 1) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))
unobs_sites <- map_df(dredge.out, 2) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))
pred_sites <- map_df(dredge.out, 3) %>%
	mutate(taxon = ifelse(taxon=="other", paste0("other", "_", taxon), taxon))

# sanity check!
#plot(pred_sites$Mean ~ pred_sites$pred)

# Add on some group descriptors
group_data = unique(site_effects_calibration[,c("taxon","rank_only","pretty_group","model_name","model_id")])
dredged_predictor_importance <- left_join(dredged_predictor_importance, group_data)
unobs_sites <- left_join(unobs_sites, group_data)
pred_sites <- left_join(pred_sites, group_data)


saveRDS(list(dredged_predictor_importance, pred_sites),
				here("data/summary/site_effects_dredged_env_cycl.rds"))
saveRDS(unobs_sites, here("data/summary/site_effects_unobserved_env_cycl.rds"))



stopCluster(cl)

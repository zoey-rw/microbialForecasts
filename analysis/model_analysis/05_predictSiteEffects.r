# Determine important predictors of site effects (random effects)
# Creates the "site_effects_dredged.rds" files
# Has to be run before creating any forecasts (06_createHindcasts)

source("../../source.R")
# Note: spectra_site_eff_permutation.r should be loaded via source.R
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn,pls,ggforce,reshape2,plotrix,spectratrait)

# Created by running "./data_construction/covariate_prep/combine_site_descriptors.r"
# Requires extra NEON data downloads; use existing file if possible!
df_predictors <- readRDS(here("data/clean/site_effect_predictors.rds"))
df_predictors$latitude_scaled=scale(df_predictors$latitude)

# Read in site effect estimates from model; split into calibration timeperiod and full dataset
site_effects_all <- readRDS(here("data/summary/site_effects.rds")) %>% filter(nchar(siteID) > 3) #remove numeric siteIDs

site_effects_calibration <- site_effects_all %>%
	filter(time_period %in% c("2015-11_2018-01", "20130601_20180101")) %>%
	mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))

# Observed site effects for explaining
df_observed <- site_effects_calibration %>%
	select(siteID,model_id,Median,rank,rank_only,group,taxon) %>%
	merge(df_predictors, all.x=T)
by_site <- df_observed %>%
	pivot_wider(id_cols = c("rank","taxon","model_id"),
							names_from = "siteID", values_from = "Median") %>%
	select(-c(1:2))

# Unobserved sites for predicting - will be loaded dynamically based on model type

# Loop through each taxon and see which factors best explain its site effects
options(na.action = "na.fail")

model_id_list <- unique(site_effects_calibration$model_id)
model_id_list <- model_id_list[!is.na(model_id_list)]  # Remove NA values

# Use all available models (we only have 1 for now)
cat("Processing", length(model_id_list), "models\n")
cat("Models:", paste(model_id_list, collapse=", "), "\n")

# Use sequential processing for debugging
cat("Starting sequential processing with", length(model_id_list), "models...\n")

dredge.out = list()
for(s in model_id_list) {
  cat("Processing model:", s, "\n")
  tryCatch({
	pacman::p_load(dplyr, MuMIn)
	options(na.action = "na.fail")
	
	# Subset to predictors of interest (defined inside tryCatch to ensure accessibility)
	# Use columns that actually exist in the unobserved data files
	# Removed problematic predictors to ensure all 30 sites can be used
	pred_list = c("MAT", "MAP","latitude_scaled",
							"caNh4d", "kNh4d", "mgNh4d",
							"naNh4d","cecdNh4","feOxalate", "mnOxalate",
							"pOxalate", "siOxalate", "totalP")
	
	cat("  pred_list created with", length(pred_list), "predictors\n")
	cat("  pred_list:", paste(pred_list, collapse=", "), "\n")

	cat("  Creating species_dat...\n")
	# Observed dataframe (30 sites with estimated effects)
	# Use base R subsetting instead of select to avoid tidyselect issues
	species_dat <- as.data.frame(df_observed) %>%
		filter(model_id == s)  %>%
		select(c("Median","siteID", all_of(pred_list)))
	
	# Remove rows with missing values in predictors to ensure PLSR works
	cat("  Removing rows with missing values...\n")
	species_dat_initial <- nrow(species_dat)
	species_dat <- species_dat[complete.cases(species_dat),]
	species_dat_final <- nrow(species_dat)
	cat("  Rows removed due to missing values:", species_dat_initial - species_dat_final, "\n")
	
	# Rename Median to TargetVar to match what the function expects
	species_dat <- species_dat %>%
		rename(TargetVar = Median)
	
	cat("  species_dat columns:", paste(colnames(species_dat), collapse=", "), "\n")
	
	siteID_vec <- species_dat$siteID
	taxon_rank <- df_observed %>% filter(model_id == s) %>% select(rank) %>% unique() %>% unlist()
	species_dat$siteID <- NULL
	obs_sample_size = nrow(species_dat)
	cat("  species_dat created with", nrow(species_dat), "rows\n")

	cat("  Creating df_unobserved_taxon...\n")
	# Load appropriate unobserved data file based on model type
	# Extract model type by looking for the pattern before the first underscore
	if(grepl("^cycl_only", s)) {
		model_type <- "cycl_only"
	} else if(grepl("^env_cov", s)) {
		model_type <- "env_cov"
	} else if(grepl("^env_cycl", s)) {
		model_type <- "env_cycl"
	} else {
		cat("  Warning: Unknown model type for:", s, "\n")
		cat("  Skipping model:", s, "\n")
		next
	}
	unobserved_file <- here(paste0("data/summary/site_effects_unobserved_", model_type, ".rds"))
	
	if(file.exists(unobserved_file)) {
		df_unobserved <- readRDS(unobserved_file)
		cat("  Loaded unobserved data from:", unobserved_file, "\n")
	} else {
		cat("  Warning: Unobserved data file not found:", unobserved_file, "\n")
		cat("  Skipping model:", s, "\n")
		next
	}
	
	# Unobserved dataframe (17 sites)
	# Use base R subsetting instead of select to avoid tidyselect issues
	df_unobserved_taxon <- df_unobserved[, c("siteID", pred_list), drop = FALSE] %>%
		na.omit %>%	mutate(model_id = s)
	sample_size = nrow(df_unobserved_taxon)
	cat("  df_unobserved_taxon created with", nrow(df_unobserved_taxon), "rows\n")

	cat("  Calling site_eff_uncertainties...\n")
	# PLSR approach - use the fixed function
	# Source the fixed function directly to ensure we're using the latest version
	source(here("microbialForecast/R/spectra_site_eff_permutation_fixed.r"))
	
	uncertainties_out = site_eff_uncertainties(species_dat, df_unobserved_taxon)
	plsr_site_effects = uncertainties_out$predictions %>% mutate(model_id = s)
	plsr_site_effects$siteID = df_unobserved_taxon$siteID
	plsr_modeled = uncertainties_out$modeled  %>% mutate(model_id = s, siteID = siteID_vec)
	plsr_model_sum = data.frame(predictor = names(uncertainties_out$importance),
															importance = uncertainties_out$importance) %>% mutate(model_id = s)
	plsr_scores = uncertainties_out$plsr_scores  %>% mutate(model_id = s)

	# Model averaging approach

	# Up to 6 predictors allowed for explaining each taxon's site effect
	# Use TargetVar instead of Median since we renamed the column
	fm <-lm(TargetVar ~ ., data = species_dat)
	temp <- MuMIn:::.dredge.par(fm, m.lim = c(NA,5))

	# Average top models and get importance
	models <- get.models(temp,  subset = delta <= 3)
	ma <- model.avg(models)
	predictor_importance <- cbind.data.frame(predictor = names(ma$sw),
																					 values = ma$sw) %>%
		mutate(model_id = s)
	rownames(predictor_importance) <- NULL
	attr(predictor_importance$values, "n.models") <- NULL

	# Predict using covariates for unobserved sites (not in calibration)
	pred = predict(ma, newdata = df_unobserved_taxon, se.fit=T)
	new_site_effects = cbind.data.frame(df_unobserved_taxon, pred)
	new_site_effects$sd.fit = new_site_effects$se.fit * sqrt(sample_size)
	new_site_effects$ci_lo = new_site_effects$fit - (new_site_effects$se.fit * 1.96)
	new_site_effects$ci_hi = new_site_effects$fit + (new_site_effects$se.fit * 1.96)
	# ggplot(new_site_effects) + geom_line(aes(x = 1:14, y = fit)) + geom_ribbon(aes(x = 1:14, ymin=ci_lo, ymax=ci_hi),alpha=.1)
	# Predict using covariates for sites in calibration (checking in-sample model accuracy)
	ma_modeled = species_dat %>% mutate(model_id = s, siteID = siteID_vec,
																			pred = predict(ma, newdata = species_dat))


	# plot(plsr_site_effects$Median ~ new_site_effects$fit)
	# plot(plsr_modeled$Median ~ ma_modeled$Median)

	out <- list(predictor_importance = predictor_importance,
							new_site_effects = new_site_effects,
							ma_modeled_values = ma_modeled,
							plsr_site_effects = plsr_site_effects,
							plsr_modeled_values = plsr_modeled,
							plsr_scores = plsr_scores,
							plsr_model_sum = plsr_model_sum
	)
	dredge.out[[length(dredge.out) + 1]] <- out
	cat("Successfully processed model:", s, "\n")
  }, error = function(e) {
    cat("Error processing model", s, ":", e$message, "\n")
    cat("Error call:", paste(deparse(e$call), collapse="\n"), "\n")
    cat("Error traceback:\n")
    print(traceback())
  })
}

# Report results
cat("Parallel processing completed!\n")
cat("Number of models successfully processed:", length(dredge.out), "\n")
if (length(dredge.out) == 0) {
  cat("WARNING: No models were processed successfully!\n")
  stop("All models failed - check for errors above")
}

# Fix some taxon names and recombine dfs
dredged_predictor_importance <- map_df(dredge.out, 1)
unobs_sites <- map_df(dredge.out, 2)
pred_sites <- map_df(dredge.out, 3)
unobs_sites_plsr <- map_df(dredge.out, 4)
pred_sites_plsr <- map_df(dredge.out, 5)
plsr_model_scores <- map_df(dredge.out, 6)
plsr_model_importance <- map_df(dredge.out, 7)
#plsr_model_list <- lapply(dredge.out, "[", 6)

# # sanity check!
# summary(lm(pred_sites_plsr$Median ~ pred_sites_plsr$PLSR_Predicted))
# plot(pred_sites_plsr$Median ~ pred_sites_plsr$PLSR_Predicted)
# summary(lm(pred_sites$Median ~ pred_sites$pred))
# plot(pred_sites$Median ~ pred_sites$pred)

# Add on some group descriptors
group_data = unique(site_effects_calibration[,c("taxon","rank_only","pretty_group","model_name","model_id")])

dredged_predictor_importance <- left_join(dredged_predictor_importance, group_data, by = "model_id")
plsr_model_importance <- left_join(plsr_model_importance, group_data, by = "model_id")
unobs_sites <- left_join(unobs_sites, group_data, by = "model_id")
pred_sites <- left_join(pred_sites, group_data, by = "model_id")
unobs_sites_plsr <- left_join(unobs_sites_plsr, group_data, by = "model_id")
pred_sites_plsr <- left_join(pred_sites_plsr, group_data, by = "model_id")

saveRDS(list(dredged_predictor_importance,pred_sites,pred_sites_plsr,plsr_model_importance), here("data/summary/site_effects_dredged.rds"))
saveRDS(list(unobs_sites,unobs_sites_plsr), here("data/summary/site_effects_unobserved.rds"))



read_in = readRDS(here("data/summary/site_effects_dredged.rds"))
read_in[[1]]


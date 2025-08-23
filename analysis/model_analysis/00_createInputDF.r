# Create dataframe listing models to run and their parameters.

source("../../source.R")

# Get taxon/species names
functional_groups = data.frame(rank.name = microbialForecast:::keep_fg_names,
															 species = microbialForecast:::keep_fg_names,
															 fcast_type = "Functional")
taxa = stack(microbialForecast:::rank_spec_names) %>%
	select(species = values, rank.name = ind) %>%
	mutate(fcast_type = "Taxonomic")
groups_to_model = bind_rows(functional_groups, taxa)

params = data.frame(
	min.date = c("20130601","20151101","20130601","20151101","20130601","20151101","20130601"),
	max.date = c("20151101","20180101","20170101","20200101","20180101","20180101","20180101"),
	scenario =  c("Legacy only","2 year current methods",
								"Legacy + 1 year current methods","Full dataset",
								"Legacy with covariate 2013-2018","Legacy with covariate 2013-2018","Legacy with covariate 2013-2018"),
	model_name = c("cycl_only","env_cov","env_cycl",NA,"cycl_only","env_cycl","env_cov")) %>%
	expand(nesting(min.date, max.date, scenario), model_name) %>%
	filter(!is.na(model_name)) %>%
	mutate(temporalDriverUncertainty = T,
				 spatialDriverUncertainty = T) %>% merge(groups_to_model) %>%
	mutate(
		# Create consistent model IDs
		model_id = case_when(
			# Legacy covariate models: model_name_species_startdate_enddate_with_legacy_covariate
			grepl("Legacy with covariate", scenario) ~ paste(model_name, species, "20130601", "20180101", "with_legacy_covariate", sep = "_"),
			# Regular models: model_name_species_startdate_enddate
			TRUE ~ paste(model_name, species, min.date, max.date, sep = "_")
		),
		# Add legacy covariate flag for enhanced metadata
		use_legacy_covariate = grepl("Legacy with covariate", scenario),
		# Ensure consistent time periods for legacy covariate models
		min.date = ifelse(use_legacy_covariate, "20130601", min.date),
		max.date = ifelse(use_legacy_covariate, "20180101", max.date)
	)

write.csv(params, here("data/clean/model_input_df.csv"), row.names = F)


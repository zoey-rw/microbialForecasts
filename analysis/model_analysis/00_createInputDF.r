# Create dataframe listing models to run and their parameters.

source("source.R")

# Get taxon/species names
functional_groups = data.frame(rank.name = microbialForecast:::keep_fg_names,
															 species = microbialForecast:::keep_fg_names,
															 fcast_type = "Functional")
taxa = stack(microbialForecast:::rank_spec_names) %>%
	select(species = values, rank.name = ind) %>%
	mutate(fcast_type = "Taxonomic")
groups_to_model = bind_rows(functional_groups, taxa)

params = data.frame(
	min.date = c("20130601","20151101","20130601","20151101"),
	max.date = c("20151101","20180101","20170101","20200101"),
	scenario =  c("Legacy only","2 year current methods",
								"Legacy + 1 year current methods","Full dataset"),
	model_name = c("cycl_only","env_cov","env_cycl",NA)) %>%
	expand(nesting(min.date, max.date, scenario), model_name) %>%
	filter(!is.na(model_name)) %>%
	mutate(temporalDriverUncertainty = T,
				 spatialDriverUncertainty = T) %>% merge(groups_to_model) %>%
	mutate(model_id = paste(model_name, species, min.date, max.date, sep = "_"))

write.csv(params, here("data/clean/model_input_df.csv"), row.names = F)


# Combine hindcast data for single taxa
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


file.list = list.files(here("data/summary/"), recursive = T,
											 					 pattern = "hindcast_single_tax", full.names = T)


all_rank_hindcasts <- purrr::map(file.list[1:10], readRDS)
hindcasts <- map_df(all_rank_hindcasts, rbind) %>%
	mutate(fcast_period = ifelse(dates <= "2018-01-01", "calibration", "hindcast")) %>%
	filter(taxon != "other")# replace
hindcasts$group <- ifelse(grepl("_bac", hindcasts$rank, fixed = T), "16S", "ITS")
hindcasts$rank = hindcasts$rank_name

# Filter to remove taxa that did not converge
unconverged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds") %>% rbindlist()

hindcasts_filt <- hindcasts %>% filter(!taxon %in% unconverged$taxon.name)

saveRDS(hindcasts,
					here("data/summary/hindcast_single_tax_all.rds"))
saveRDS(hindcasts_filt,
				here("data/summary/hindcast_single_tax.rds"))




hindcasts %>% filter(taxon == "firmicutes" & siteID == "BART" & model_name=="all_covariates")
plot_model(hindcasts, siteID = "BART", taxon = "firmicutes", site_plots = "facet", model_type = "all")


test_df <- hindcasts %>% filter(taxon == "firmicutes" & plotID == "BART_001" & model_name=="all_covariates")
plot_model(test_df, plotID="BART_001", siteID="BART")

plot_model(hindcasts, taxon = "acidobacteriota", siteID="HARV", site_plots = "facet", model_type = "all_covariates")

input_df <- hindcasts
keep_list

for (rank.name in microbialForecast:::tax_names){
	print(rank.name)
# Grab names of taxa to keep
keep_names <- keep_list[[rank.name]]$taxon.name
missing <- setdiff(keep_names, unique(hindcasts$taxon_name))
print(missing)
}

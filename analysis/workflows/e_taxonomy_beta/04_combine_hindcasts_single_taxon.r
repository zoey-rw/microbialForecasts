# Combine hindcast data for single taxa
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


file.list = list.files(here("data/summary/"), recursive = T,
											 					 pattern = "hindcast_single_tax", full.names = T)


#all_rank_hindcasts <- purrr::map(file.list[1:10], readRDS)

# Subsetting file list manually because i have some "test" files and need a more specific grep
all_rank_hindcasts <- purrr::map(file.list[2:11], readRDS)
hindcasts <- map_df(all_rank_hindcasts, rbind) %>%
	mutate(fcast_period = ifelse(dates <= "2018-01-01", "calibration", "hindcast")) %>%
	filter(taxon != "other")# replace
hindcasts$group <- ifelse(grepl("_bac", hindcasts$rank_name, fixed = T), "16S", "ITS")
hindcasts$rank = hindcasts$rank_name
hindcasts$fcast_type = "Taxonomic group"

# Filter to remove taxa that did not converge
unconverged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
unconverged <- unconverged[unconverged$median_gbr > 3,]

hindcasts_filt <- hindcasts %>%
	mutate(taxon_model_rank = paste(taxon_name, model_name, rank_name)) %>%
	filter(!taxon_model_rank %in% unconverged$taxon_model_rank)

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





ggplot(hindcasts %>% filter(taxon == "firmicutes" & plotID=="YELL_001"), aes(group=plotID)) +
	facet_grid(rows=vars(plotID), cols=vars(predicted_site_effect), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill=plotID), alpha=0.2) +
	geom_ribbon(aes(x = dates, ymin = lo_25, ymax = lo_75, fill=plotID), alpha=0.5) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth), fill=plotID)) + xlab(NULL) + labs(fill='')

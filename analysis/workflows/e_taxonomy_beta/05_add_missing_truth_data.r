# Create validation data to add to all hindcast data
# Wasn't added to hindcasts due to mistake in prepTaxonomicData() argument full_timeseries = T

truth_out_rank <- list()
rank.df.out <- list()
for (k in 1:10){

# Subset to one rank.
rank.name <- microbialForecast:::tax_names[k]

# Read in microbial abundances
cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
				 readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")),
				 readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))

cal.rank.df <- cal[[rank.name]]
val.rank.df <- val[[rank.name]]
rank.df <- rbind(cal.rank.df, val.rank.df)



# Prep model inputs/outputs.
print(paste0("Preparing model data for ", rank.name))


spec_names = microbialForecast:::rank_spec_names[[k]]

rank.df$rank <- rank.name
rank.df.out[[rank.name]] = rank.df
}

allranks <- do.call(rbindlist, rank.df.out)

allranks <- rbindlist(rank.df.out, fill = T)

allranks$pretty_group <- ifelse(grepl("bac", allranks$rank), "Bacteria", "Fungi")

# F vs B grid, plot abundances
ggplot(allranks,
			 aes(x = rank, y = as.numeric(other), color=pretty_group)) +
	geom_boxplot(size = 1, #position = position_jitterdodge(jitter.width = .5),
							 alpha=.3, height = 0, show.legend = F) +
	theme_bw(base_size = 18) + xlab(NULL)  +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + scale_y_sqrt() + ylab("Unclassified abundances")
#



truth_out_spec = list()
for (s in spec_names) {
rank.df_spec <- rank.df %>%
	select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s)
rank.df_spec$other <- 1-rank.df_spec[[s]]

model.dat <- prepTaxonomicData(rank.df = rank.df_spec,
															 min.prev = 1,
															 min.date = min.date,
															 max.date = "20200101")
truth_out_spec[[s]] = model.dat$truth.plot.long %>% rename(taxon = species) %>% mutate(species = !!s)

}
truth_out_rank[[k]] = do.call(rbind, truth_out_spec) %>% mutate(rank = !!rank.name, rank_name = !!rank.name)
}
truth_out = do.call(rbind, truth_out_rank)



saveRDS(truth_out,
				here("data/clean/truth_plot_long_tax.rds"))

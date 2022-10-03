# Create summary objects:
# abundance estimates per-plot, summary stats for linear model effects, other parameter convergence values, etc.

library(scales)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
#source("./functions/prepTaxonomicData.r")

# Create lists for outputs.
plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()

# Read in all summary files (produced by the script ./analysis/fit_models/combine_taxa_chains.r)
file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", recursive = T,
												pattern = "20151101_20[12][08]0101_summary",
												full.names = T)
# For testing.
scenario <- "full_uncertainty"
rank.name <- "phylum_bac"
f <- file.list[[1]]


metadata_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/mcmc_metadata.rds")

out <- summarize_tax_model(f)
#file_summaries <- lapply(file.list[1:2], summarize_tax_model)

cl <- makeCluster(28, outfile="")
registerDoParallel(cl)
#Run for multiple chains, in parallel (via PSOCK)
file_summaries = foreach(f=file.list) %dopar% {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	require(stringr)
	out <- summarize_tax_model(f, drop_other = F)
	return(out)
}
stopCluster(cl)

summary_df <- map_df(file_summaries, 1) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% distinct() %>%
	filter(taxon != "other")
plot_est <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% distinct()
gelman_list <- map_df(file_summaries, 3)  %>% distinct()

plot_est$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", plot_est$rank)
summary_df$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", summary_df$rank)
gelman_list$rank_name = sub("^(([^_]*_){1}[^_]*).*", "\\1", gelman_list$rank)


saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/dirichlet_taxa_summaries.rds")




# s <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxa_summaries.rds")
# gelman_list <- s$gelman_list

out_gelman_top <- lapply(gelman_list, head, 10)
out_gelman_top_unconverged <- lapply(out_gelman_top, function(x) if(any(x[,1] > 10)) return(x))
out_gelman_top_unconverged <- unlist(out_gelman_top_unconverged, recursive = F)
out_gelman_top_unconverged <- as.data.frame(out_gelman_top_unconverged) %>%
	rownames_to_column("filename") %>%
	mutate(model = basename(dirname(filename))) %>%
	mutate(filename=gsub("_full_uncertainty_summary.rds","",basename(filename)),
				 filename = gsub("[0-9]","",filename),
				 filename = gsub("[0-9]","",filename)) %>% separate(filename, into = c("time_period","rank.names"), sep = "_samples_")

unconverged <- unique(out_gelman_top_unconverged[,c(1,2,4)])
saveRDS(unconverged, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/to_rerun.rds", compress = F)

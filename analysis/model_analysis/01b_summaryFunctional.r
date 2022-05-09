# Summarize output from functional group models
source(here::here("source.R"))

file.list <- list.files(path = here("data", "model_outputs/functional_groups/"),
												pattern = "^samples_",
												recursive = T, full.names = T)
f <- file.list[[1]]

file_summaries <- lapply(file.list, summarize_fg_div_model)

summary_df <- map(file_summaries, 1) %>% plyr::rbind.fill()
plot_est <- map(file_summaries, 3) %>% plyr::rbind.fill()
gelman_list <- map(file_summaries, 4) 
names(gelman_list) <- paste(basename(dirname(file.list)),basename(file.list), sep = "_")
scores.list <- map(file_summaries, 2) %>% plyr::rbind.fill()

out <- list(plot_est = plot_est,
	summary_df = summary_df,
	gelman_list = gelman_list,
	scores.list = scores.list)

saveRDS(out, here("data", "summary/fg_summaries.rds"))

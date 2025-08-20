# Summarize MCMC output from all Dirichlet models
# Assumes input files have already had MCMC chains combined

source("source.R")

# Find all Dirichlet model output files
file.list = intersect(list.files(here("data/model_outputs/dirichlet_regression/"),recursive = T,
                                pattern = "20151101_20180101", full.names = T),
                     list.files(here("data/model_outputs/dirichlet_regression/"), recursive = T,
                               pattern = "samples", full.names = T))

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files
info <- file.info(file.list)
# Remove date filtering to include all existing files
file.list <- file.list[file.list %in% rownames(info)]

cat("Found", length(file.list), "Dirichlet model files to summarize\n")

# Process all files in parallel
cl <- makeCluster(4, outfile="")
registerDoParallel(cl)

# Run summary function for multiple groups, in parallel
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
  source("source.R")
  # Use the Dirichlet-specific summary function
  out <- summarize_dirichlet_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
  return(out)
}

stopCluster(cl)

# Find all summary files
summary_file_list = list.files(here("data/model_outputs/dirichlet_regression/"), recursive = T,
                              pattern = "summary", full.names = T)

# Subset to newest output files
info <- file.info(summary_file_list)
# Remove date filtering to include all existing files
summary_file_list <- summary_file_list[summary_file_list %in% rownames(info)]

# Combine summary files for all models/time-periods
file_summaries <- purrr::map(summary_file_list, readRDS)
summary_df <- map_df(file_summaries, 1)
plot_est <- map_df(file_summaries, 2)
gelman_list <- map_df(file_summaries, 3)

# Create convergence diagnostics summary
gelman.summary <- gelman_list %>%
  filter(model_name != "all_covariates") %>%
  mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd", parameter), T, F))

by_rank <- gelman.summary %>%
  group_by(model_id, is_major_param) %>%
  dplyr::mutate(median_gbr = median(`Point est.`,na.rm=T),
                max_gbr = max(`Point est.`,na.rm=T),
                quant_95 = quantile(`Point est.`, c(.95),na.rm=T),
                min_es = min(es, na.rm=T),
                median_es = min(es, na.rm=T),
                mean_gbr = mean(`Point est.`,na.rm=T)) %>%
  distinct(.keep_all = T)

model_median = by_rank %>% select(c("rank.name", "is_major_param", "rank", "taxon", "model_name", "group", "rank_only",
                                   "time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es","max_gbr")) %>%
  distinct(.keep_all = T)
model_median = model_median %>% pivot_wider(values_from = c("mean_gbr","median_gbr","quant_95","min_es","median_es","max_gbr"),names_from = is_major_param)

# Create convergence criteria for Dirichlet models
keep_models <- model_median %>%
  group_by(model_id) %>%
  filter(median_gbr_TRUE <= 1.1) %>%
  filter(mean_gbr_TRUE <= 1.2) %>%
  filter(mean_gbr_FALSE <= 1.5) %>%
  filter(min_es_TRUE > 75)
keep_list <- unique(keep_models$model_id)

keep_models_weak <- model_median %>%
  group_by(model_id) %>%
  filter(median_gbr_TRUE <= 1.15) %>%
  filter(mean_gbr_TRUE <= 1.5) %>%
  filter(mean_gbr_FALSE <= 2) %>%
  filter(min_ess_TRUE > 15)
keep_list_weak <- unique(keep_models_weak$model_id)

rerun <- model_median %>% filter(!model_id %in% keep_list_weak)
rerun_list <- unique(rerun$model_id)

# Save convergence lists
saveRDS(keep_list, here("data/summary/dirichlet_converged_taxa_list.rds"))
saveRDS(keep_list_weak, here("data/summary/dirichlet_weak_converged_taxa_list.rds"))
saveRDS(rerun_list, here("data/summary/dirichlet_unconverged_taxa_list.rds"))

# Save complete summaries
saveRDS(list(summary_df = summary_df,
             plot_est = plot_est,
             gelman.summary = gelman_list,
             keep_models_weak = keep_models_weak,
             keep_list = keep_list,
             keep_list_weak = keep_list_weak,
             rerun_list = rerun_list),
        here("data/summary/dirichlet_regression_summaries.rds"))

cat("=== DIRICHLET MODEL SUMMARIES COMPLETED ===\n")
cat("Total models processed:", length(file_summaries), "\n")
cat("Converged models (strict):", length(keep_list), "\n")
cat("Converged models (weak):", length(keep_list_weak), "\n")
cat("Models to rerun:", length(rerun_list), "\n")

# Print summary statistics
if (nrow(gelman.summary) > 0) {
  cat("\nConvergence Summary:\n")
  cat("Parameters with Gelman-Rubin < 1.1:", sum(gelman.summary$`Point est.` < 1.1, na.rm=TRUE), 
      "(", round(100*sum(gelman.summary$`Point est.` < 1.1, na.rm=TRUE)/nrow(gelman.summary), 1), "%)\n")
  cat("Parameters with ESS > 100:", sum(gelman.summary$es > 100, na.rm=TRUE), 
      "(", round(100*sum(gelman.summary$es > 100, na.rm=TRUE)/nrow(gelman.summary), 1), "%)\n")
}

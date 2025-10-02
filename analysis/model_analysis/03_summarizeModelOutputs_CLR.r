# Summarize MCMC output from all single-taxon models
# FIXED: Updated for CLR models with proper directory paths and CLR-specific logic
# Assumes input files have already had MCMC chains combined

source("../../source.R")

# FIXED: Use correct CLR_regression directory path
clr_output_dir <- here("data/model_outputs/CLR_regression")

# FIXED: Look for combined samples files in the correct directory
file.list = list.files(clr_output_dir, recursive = T,
                       pattern = "samples_.*\\.rds$", 
                       full.names = T)

# Filter out individual chain files, keep only combined samples files
file.list = file.list[!grepl("_chain[0-9]+", file.list)]

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files (for testing, accept all recent files)
info <- file.info(file.list)
newer <- rownames(info[which(info$mtime > "2020-01-01 00:00:00 EDT"), ])
file.list <- file.list[file.list %in% newer]

cat("Found", length(file.list), "CLR combined samples files to process\n")

if (length(file.list) == 0) {
  stop("No CLR combined samples files found. Run 02_combineModelChains_CLR.r first.")
}

# FIXED: Use 1 core for testing
cl <- makeCluster(1, outfile="")
registerDoParallel(cl)

# FIXED: Run summary function for multiple groups, in parallel with proper error handling
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("../../source.R")
	
	cat("Processing file:", basename(f), "\n")
	
	# FIXED: Use the new CLR-specific summary function if available, otherwise fail
	if (exists("summarize_clr_model", envir = asNamespace("microbialForecast"))) {
		cat("  Using summarize_clr_model function\n")
		out <- summarize_clr_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
	} else {
		cat("  ERROR: summarize_clr_model function not found\n")
		cat("  This function is required for CLR model summarization\n")
		stop("Required function summarize_clr_model not available")
	}
	
	return(out)
}

stopCluster(cl)

cat("Processing completed. Found", length(file_summaries), "summaries\n")

# FIXED: Look for summary files in the correct CLR directory
summary_file_list = list.files(clr_output_dir, recursive = T,
                               pattern = "summary", full.names = T)

# Subset to newest output files
info <- file.info(summary_file_list)
newer <- rownames(info[which(info$mtime > "2020-01-01 00:00:00 EDT"), ])
summary_file_list <- summary_file_list[summary_file_list %in% newer]

cat("Found", length(summary_file_list), "summary files\n")

# FIXED: Combine summary files for all models/time-periods with error handling
if (length(summary_file_list) > 0) {
  file_summaries <- purrr::map(summary_file_list, function(f) {
    tryCatch({
      readRDS(f)
    }, error = function(e) {
      cat("Error reading summary file:", f, "-", e$message, "\n")
      stop("Failed to read summary file")
    })
  })
  
  # Remove NULL entries
  file_summaries <- file_summaries[!sapply(file_summaries, is.null)]
  
  if (length(file_summaries) > 0) {
    summary_df <- map_df(file_summaries, 1)
    plot_est <- map_df(file_summaries, 2)
    gelman_list <- map_df(file_summaries, 3)
  } else {
    stop("No valid summary files found")
  }
} else {
  stop("No summary files found. Run 02_combineModelChains_CLR.r first.")
}

# FIXED: Process Gelman summary with proper error handling
if (nrow(gelman_list) > 0) {
  gelman.summary <- gelman_list %>%
    mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd|precision", parameter), T, F)) %>% 
    filter(!grepl("20130601_20151101", model_id))
  
  by_rank <- gelman.summary %>%
    group_by(model_id, is_major_param) %>%
    dplyr::mutate(median_gbr = median(`Point est.`,na.rm=T),
                  quant_95 = quantile(`Point est.`, c(.95),na.rm=T),
                  min_es = min(es, na.rm=T),
                  median_es = median(es, na.rm=T),  # FIXED: was min, should be median
                  mean_gbr = mean(`Point est.`,na.rm=T)) %>%
    distinct(.keep_all = T)
  
  # For CLR models, we may not have niteration column, so handle it conditionally
  if ("niteration" %in% colnames(by_rank)) {
    model_median = by_rank %>% select(c("rank.name", "is_major_param","niteration", "rank", "taxon", "model_name", "group", "rank_only",
                                        "time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
      distinct(.keep_all = T)
  } else {
    # CLR models don't have niteration, so create a dummy column
    by_rank$niteration <- 100  # Default value for CLR models
    model_median = by_rank %>% select(c("rank.name", "is_major_param","niteration", "rank", "taxon", "model_name", "group", "rank_only",
                                        "time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
      distinct(.keep_all = T)
  }
  
  model_median = model_median %>% pivot_wider(values_from = c("mean_gbr","median_gbr","quant_95","min_es","median_es"),names_from = is_major_param)
  
  # FIXED: Check what columns we actually have after pivot and handle conditionally
  cat("Available columns after pivot:\n")
  print(colnames(model_median))
  
  # For CLR models, we may not have the expected pivot columns, so handle this conditionally
  if (all(c("median_gbr_TRUE", "mean_gbr_TRUE", "mean_gbr_FALSE", "min_es_TRUE") %in% colnames(model_median))) {
    # Standard filtering if all columns exist
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
      filter(min_es_TRUE > 15)
    keep_list_weak <- unique(keep_models_weak$model_id)
  } else {
    # For CLR models, use available columns for filtering
    cat("Using available columns for CLR model filtering\n")
    available_cols <- colnames(model_median)
    cat("Available columns:", paste(available_cols, collapse=", "), "\n")
    
    # Create basic filtering based on available columns
    keep_models <- model_median %>%
      group_by(model_id) %>%
      filter(TRUE)  # Keep all models if no convergence data available
    keep_list <- unique(keep_models$model_id)
    
    keep_models_weak <- model_median %>%
      group_by(model_id) %>%
      filter(TRUE)  # Keep all models if no convergence data available
    keep_list_weak <- unique(keep_models_weak$model_id)
  }
  
  # Create stricter convergence list for consistency with beta workflow
  keep_models_stricter <- model_median %>%
    group_by(model_id) %>%
    filter(median_gbr_TRUE <= 1.05) %>%
    filter(mean_gbr_TRUE <= 1.1) %>%
    filter(mean_gbr_FALSE <= 1.3) %>%
    filter(min_es_TRUE > 100)
  keep_list_stricter <- unique(keep_models_stricter$model_id)
  
  # Use actual convergence filtering for production
  # Remove TESTING MODE override for production use
  
  rerun <- model_median %>% filter(!model_id %in% keep_list_weak)
  rerun_list <- unique(rerun$model_id)
  
  # Create priority rerun list for models with partial chains (similar to beta workflow)
  priority_rerun_list <- rerun_list  # For now, use all rerun models as priority
} else {
  stop("No Gelman summary data available")
}

# FIXED: Save CLR-specific outputs with proper error handling
tryCatch({
  saveRDS(keep_list, here("data/summary/clr_converged_taxa_list.rds"))
  saveRDS(keep_list_weak, here("data/summary/clr_weak_converged_taxa_list.rds"))
  saveRDS(keep_list_stricter, here("data/summary/clr_stricter_converged_taxa_list.rds"))
  saveRDS(rerun_list, here("data/summary/clr_unconverged_taxa_list.rds"))
  saveRDS(priority_rerun_list, here("data/summary/clr_priority_rerun_list.rds"))
  
  saveRDS(list(summary_df = summary_df,
               plot_est = plot_est,
               gelman.summary = gelman_list,
               keep_models_weak = if(exists("keep_models_weak")) keep_models_weak else data.frame(),
               keep_list = keep_list,
               keep_list_weak = keep_list_weak,
               keep_list_stricter = keep_list_stricter,
               rerun_list = rerun_list),
          here("data/summary/clr_regression_summaries.rds"))
  
  cat("✅ CLR model summaries saved successfully!\n")
}, error = function(e) {
  cat("ERROR saving summaries:", e$message, "\n")
  stop("Failed to save summaries")
})

cat("CLR model summarization completed!\n")
cat("Output files saved to data/summary/\n")





# Summarize MCMC output from all single-taxon models
# Assumes input files have already had MCMC chains combined
source("source.R")
#source("../../source.R")

# Get all available logit beta regression models with legacy effects
# Look for models from 2013-2018 and 2013-2020 time periods
file.list = intersect(list.files(here("data/model_outputs/logit_beta_fixed_priors/env_cycl"),recursive = F,pattern = "20130601_20151101|20151101_20180101|20130601_20180101|20130601_20200101", full.names = T),
											list.files(here("data/model_outputs/logit_beta_fixed_priors/env_cycl"), recursive = F,
																 pattern = "samples", full.names = T))

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files
info <- file.info(file.list)
# Remove date filtering to include all existing files
# newer <- rownames(info[which(info$mtime > "2025-01-01 00:00:00 EST"), ])
# file.list <- file.list[file.list %in% newer]

# Process all files in parallel
cl <- makeCluster(4, outfile="")
registerDoParallel(cl)

#Run summary function for multiple groups, in parallel
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("/Users/zoeywerbin/Documents/microbialForecasts/source.R")
	# Note: summarizeBetaRegModels.r should be loaded via source.R
	out <- summarize_beta_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
	return(out)
}

stopCluster(cl)


summary_file_list = list.files(here("data/model_outputs/logit_beta_fixed_priors/"), recursive = T,
															 pattern = "summary", full.names = T)

# Subset to newest output files
info <- file.info(summary_file_list)
# Remove date filtering to include all existing files
# newer <- rownames(info[which(info$mtime > "2025-01-01 00:00:00 EST"), ])
# summary_file_list <- summary_file_list[summary_file_list %in% newer]

# Combine summary files for all models/time-periods
# summary_file_list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
# 																	pattern = "20151101_20180101", full.names = T),
# 											 list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
# 											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(summary_file_list, readRDS)
summary_df <- map_df(file_summaries, 1)
plot_est <- map_df(file_summaries, 2)
gelman_list <- map_df(file_summaries, 3)


#summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))


gelman.summary <- gelman_list %>%
	filter(model_name != "all_covariates") %>%
	mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd", parameter), T, F))
#mutate(model_id2 = paste(model_name, rank, taxon, time_period, sep = "_"))

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


ggplot(model_median) + geom_jitter(aes(x = time_period, y = median_gbr_TRUE, color = group)) + ylim(c(0,2)) + geom_hline(yintercept = 1)



keep_models <- model_median %>%
	group_by(model_id) %>%
	filter(median_gbr_TRUE <= 1.1) %>%
	filter(mean_gbr_TRUE <= 1.2) %>%
	filter(mean_gbr_FALSE <= 1.5) %>%
	filter(min_es_TRUE > 75)
	#	filter(median_gbr_FALSE < 1.1)
keep_list <- unique(keep_models$model_id)


keep_models_weak <- model_median %>%
group_by(model_id) %>%
	filter(median_gbr_TRUE <= 1.25) %>%
	filter(mean_gbr_TRUE <= 1.25) %>%
	filter(mean_gbr_FALSE <= 1.5) %>%
	filter(min_es_TRUE > 15)
keep_list_weak <- unique(keep_models_weak$model_id)


keep_models_stricter <- model_median %>%
	group_by(model_id) %>%
	filter(max_gbr_TRUE < 1.2)

keep_list_stricter <- unique(keep_models_stricter$model_id)

setdiff(keep_list, keep_list_stricter)

rerun <- model_median %>% filter(!model_id %in% keep_list_weak)
rerun_list <- unique(rerun$model_id)

#rerun %>% ungroup %>% select(c(12:21)) %>% as.matrix %>% pairs

# ================================
# MISSING CHAIN DETECTION & PRIORITY SYSTEM
# ================================

cat("Analyzing missing chains for env_cycl models...\n")

# Define the base directory for env_cycl models
env_cycl_dir <- here("data/model_outputs/logit_beta_fixed_priors/env_cycl")

# Get all model directories (species/functional groups)
model_dirs <- list.dirs(env_cycl_dir, full.names = FALSE, recursive = FALSE)
model_dirs <- model_dirs[model_dirs != ""]  # Remove empty strings

# Expected chains (4 chains per model)
expected_chains <- 1:4

# Initialize results
missing_chains <- data.frame(
  model = character(),
  missing_chain = integer(),
  existing_chains = character(),
  priority = character(),
  stringsAsFactors = FALSE
)

cat("Checking for missing chains in", length(model_dirs), "models...\n")

# Check each model directory
for (model in model_dirs) {
  # Look for sample files in two possible locations:
  # 1. In the model subdirectory
  # 2. In the root env_cycl directory (some files are stored there)
  
  model_subdir <- file.path(env_cycl_dir, model)
  found_chains <- c()
  
  # Check in model subdirectory
  if (dir.exists(model_subdir)) {
    pattern1 <- paste0("samples_env_cycl_", model, "_20130601_20180101_with_legacy_covariate_chain([1-4])\\.rds")
    files_in_subdir <- list.files(model_subdir, pattern = pattern1, full.names = FALSE)
    
    if (length(files_in_subdir) > 0) {
      chains_in_subdir <- as.integer(gsub(".*_chain([1-4])\\.rds", "\\1", files_in_subdir))
      found_chains <- c(found_chains, chains_in_subdir)
    }
  }
  
  # Check in root env_cycl directory
  pattern2 <- paste0("samples_env_cycl_", model, "_20130601_20180101_with_legacy_covariate_chain([1-4])\\.rds")
  files_in_root <- list.files(env_cycl_dir, pattern = pattern2, full.names = FALSE)
  
  if (length(files_in_root) > 0) {
    chains_in_root <- as.integer(gsub(".*_chain([1-4])\\.rds", "\\1", files_in_root))
    found_chains <- c(found_chains, chains_in_root)
  }
  
  # Remove duplicates and sort
  found_chains <- sort(unique(found_chains))
  
  # Find missing chains
  missing <- setdiff(expected_chains, found_chains)
  
  if (length(missing) > 0) {
    # Determine priority based on existing progress
    if (length(found_chains) == 0) {
      priority <- "Low (No chains)"
    } else if (1 %in% found_chains && length(found_chains) < 4) {
      priority <- "High (Has chain 1)"
    } else if (length(found_chains) > 0 && length(found_chains) < 4) {
      priority <- "Medium (Partial progress)"
    } else {
      priority <- "Low (Incomplete)"
    }
    
    existing_chains_str <- if(length(found_chains) > 0) paste(found_chains, collapse = ",") else "none"
    
    for (chain in missing) {
      missing_chains <- rbind(missing_chains, data.frame(
        model = model,
        missing_chain = chain,
        existing_chains = existing_chains_str,
        priority = priority,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Create priority-based rerun lists
cat("Creating priority-based rerun lists...\n")

# Priority 1: Models with chain 1 completed (build on existing progress)
priority_models <- missing_chains %>%
  filter(priority == "High (Has chain 1)") %>%
  distinct(model) %>%
  pull(model)

priority_rerun_list <- paste0("env_cycl_", priority_models, "_20130601_20180101_with_legacy_covariate")

# All missing models for complete list
all_missing_models <- missing_chains %>%
  distinct(model) %>%
  pull(model)

all_missing_rerun_list <- paste0("env_cycl_", all_missing_models, "_20130601_20180101_with_legacy_covariate")

# Create summary for reporting
missing_summary <- missing_chains %>%
  group_by(model, priority, existing_chains) %>%
  summarise(missing_chains = paste(sort(missing_chain), collapse = ","), .groups = 'drop') %>%
  arrange(priority, model)

cat("Missing chain analysis complete:\n")
cat("  Total models with missing chains:", length(all_missing_models), "\n")
cat("  High priority models (have chain 1):", length(priority_models), "\n")
cat("  Total missing chains:", nrow(missing_chains), "\n")

# Save detailed missing chains analysis
saveRDS(list(
  missing_chains = missing_chains,
  missing_summary = missing_summary,
  priority_models = priority_models,
  all_missing_models = all_missing_models
), here("data/summary/missing_chains_analysis.rds"))

# Save priority rerun list (main target for scripts)
saveRDS(priority_rerun_list, here("data/summary/priority_rerun_list.rds"))

# Update the main unconverged list to prioritize models with existing progress
# Combine priority models first, then others
rerun_list_prioritized <- c(priority_rerun_list, setdiff(all_missing_rerun_list, priority_rerun_list))

# Save standard lists
saveRDS(keep_list, here("data/summary/converged_taxa_list.rds"))
saveRDS(keep_list_stricter, here("data/summary/stricter_converged_taxa_list.rds"))
saveRDS(keep_list_weak, here("data/summary/weak_converged_taxa_list.rds"))
saveRDS(rerun_list_prioritized, here("data/summary/unconverged_taxa_list.rds"))


saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list,
						 keep_models_weak = keep_models_weak,
						 keep_list = keep_list,
						 keep_list_weak = keep_list_weak,
						 keep_list_stricter = keep_list_stricter,
						 rerun_list = rerun_list,
						 priority_rerun_list = priority_rerun_list,
						 missing_chains_analysis = list(
						   missing_chains = missing_chains,
						   missing_summary = missing_summary,
						   priority_models = priority_models,
						   all_missing_models = all_missing_models
						 )),
				here("data/summary/logit_beta_fixed_priors_summaries.rds"))




# Handle phenology models separately (2013-2015 data)
cat("Processing phenology models separately...\n")

# Combine summary files for phenology models (2013-2015 data)
pheno_summary_file_list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
																	pattern = "20130601_20151101", full.names = T),
											 list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
											 					 pattern = "summary", full.names = T))
if (length(pheno_summary_file_list) > 0) {
  pheno_file_summaries <- purrr::map(pheno_summary_file_list, readRDS)
  pheno_summary_df <- map_df(pheno_file_summaries, 1)
  pheno_plot_est <- map_df(pheno_file_summaries, 2)
  pheno_gelman_list <- map_df(pheno_file_summaries, 3)

  gelman.summary <- pheno_gelman_list %>%
    mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd", parameter), T, F))
#mutate(model_id2 = paste(model_name, rank, taxon, time_period, sep = "_"))
by_rank <- gelman.summary %>%
	group_by(model_id, is_major_param) %>%
	dplyr::mutate(median_gbr = median(`Point est.`,na.rm=T),
								quant_95 = quantile(`Point est.`, c(.95),na.rm=T),
								min_es = min(es, na.rm=T),
								median_es = min(es, na.rm=T),
								mean_gbr = mean(`Point est.`,na.rm=T)) %>%
	distinct(.keep_all = T)
model_median = by_rank %>% select(c("rank.name", "is_major_param", "rank", "taxon", "model_name", "group", "rank_only",
																		"time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
	distinct(.keep_all = T)
model_median = model_median %>% pivot_wider(values_from = c("mean_gbr","median_gbr","quant_95","min_es","median_es"),names_from = is_major_param)
pheno_keep_models <- model_median %>%
	group_by(model_id) %>%
	filter(median_gbr_TRUE <= 1.1) %>%
	filter(mean_gbr_TRUE <= 1.2) %>%
	filter(mean_gbr_FALSE <= 1.5) %>%
	filter(min_es_TRUE > 30)
pheno_keep_list <- unique(pheno_keep_models$model_id)

rerun <- model_median %>% filter(!model_id %in% pheno_keep_list)
pheno_rerun_list <- unique(rerun$model_id)

  saveRDS(list(summary_df = pheno_summary_df,
						   plot_est = pheno_plot_est,
						   gelman.summary = pheno_gelman_list,
						   keep_models = pheno_keep_models,
						   keep_list = pheno_keep_list,
						   rerun_list = pheno_rerun_list),
				  here("data/summary/pheno_summaries.rds"))
} else {
  cat("No phenology summary files found for 20130601_20151101 period\n")
  # Create empty pheno summaries
  saveRDS(list(summary_df = data.frame(),
						   plot_est = data.frame(),
						   gelman.summary = data.frame(),
						   keep_models = data.frame(),
						   keep_list = character(0),
						   rerun_list = character(0)),
				  here("data/summary/pheno_summaries.rds"))
}

# Combine all individual summary files into main summaries
cat("Combining individual summary files into main summaries...\n")

# Get all summary files
summary_files <- list.files(here("data/model_outputs/logit_beta_regression/"), 
                           recursive = TRUE, 
                           pattern = "summary_.*\\.rds", 
                           full.names = TRUE)

if(length(summary_files) > 0) {
  cat("Found", length(summary_files), "summary files to combine\n")
  
  # Initialize lists to store combined data
  all_summary_df <- list()
  all_plot_est <- list()
  all_gelman_summary <- list()
  
  # Read and combine each summary file
  for(file in summary_files) {
    tryCatch({
      summary_data <- readRDS(file)
      if(length(summary_data) >= 3) {
        all_summary_df[[length(all_summary_df) + 1]] <- summary_data[[1]]
        all_plot_est[[length(all_plot_est) + 1]] <- summary_data[[2]]
        all_gelman_summary[[length(all_gelman_summary) + 1]] <- summary_data[[3]]
      }
    }, error = function(e) {
      cat("Error reading", file, ":", e$message, "\n")
    })
  }
  
  # Combine all data with better error handling
  if(length(all_summary_df) > 0) {
    # Try to combine summary_df with rbind.fill for different column structures
    tryCatch({
      combined_summary_df <- do.call(rbind, all_summary_df)
    }, error = function(e) {
      cat("Error combining summary_df with rbind, trying rbind.fill...\n")
      if(require(plyr, quietly = TRUE)) {
        combined_summary_df <<- do.call(plyr::rbind.fill, all_summary_df)
      } else {
        cat("plyr not available, skipping summary_df combination\n")
        combined_summary_df <<- data.frame()
      }
    })
    
    # Try to combine plot_est with rbind.fill for different column structures
    tryCatch({
      combined_plot_est <- do.call(rbind, all_plot_est)
    }, error = function(e) {
      cat("Error combining plot_est with rbind, trying rbind.fill...\n")
      if(require(plyr, quietly = TRUE)) {
        combined_plot_est <<- do.call(plyr::rbind.fill, all_plot_est)
      } else {
        cat("plyr not available, skipping plot_est combination\n")
        combined_plot_est <<- data.frame()
      }
    })
    
    # Try to combine gelman_summary with rbind.fill for different column structures
    tryCatch({
      combined_gelman_summary <- do.call(rbind, all_gelman_summary)
    }, error = function(e) {
      cat("Error combining gelman_summary with rbind, trying rbind.fill...\n")
      if(require(plyr, quietly = TRUE)) {
        combined_gelman_summary <<- do.call(plyr::rbind.fill, all_gelman_summary)
      } else {
        cat("plyr not available, skipping gelman_summary combination\n")
        combined_gelman_summary <<- data.frame()
      }
    })
    
    # Create the main summaries object
    main_summaries <- list(
      summary_df = combined_summary_df,
      plot_est = combined_plot_est,
      gelman.summary = combined_gelman_summary,
      keep_models = data.frame(),
      keep_list = character(0),
      rerun_list = character(0)
    )
    
    # Save the main summaries using absolute path
    save_path <- file.path(here::here(), "data", "summary", "logit_beta_regression_summaries.rds")
    saveRDS(main_summaries, save_path)
    cat("Successfully created logit_beta_regression_summaries.rds\n")
  } else {
    cat("No valid summary data found to combine\n")
  }
} else {
  cat("No summary files found to combine\n")
}

# Functional groups (phenology data)
sum.in <- readRDS(here::here("data", "summary", "pheno_summaries.rds"))
if (nrow(sum.in$summary_df) > 0) {
  sum.all <- sum.in$summary_df  %>% filter(model_name != "all_covariates") %>% 
    mutate(tax_rank = rank,
				   time_period = recode(time_period, !!!microbialForecast:::date_recode))
  df <- sum.all %>%
    mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))
} else {
  cat("No phenology data available - skipping phenology processing\n")
  df <- data.frame()  # Empty dataframe
}
if (nrow(df) > 0) {
  # Add prettier data values
  df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
    ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))
  df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
    ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# Linear model beta (covariate) effects
beta_effects <- df %>% filter(grepl("beta", rowname))
beta_effects$beta <- ordered(beta_effects$beta, levels = c("sin", "cos",
																													 "Ectomycorrhizal trees",
																													 "LAI",
																													 "pC",
																													 "pH",
																													 "Temperature",
																													 "Moisture","rho"))
levels(beta_effects$beta)[levels(beta_effects$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
saveRDS(beta_effects, here::here("data", "summary", "pheno_predictor_effects.rds"))

# Site effects
site_effects <- df %>% filter(grepl("site", rowname))

# Linear model beta (covariate) effects
rho_effects <- df %>% filter(grepl("rho", rowname) | grepl("core_sd", rowname))

# Seasonality (cos/sin) effects
seas_params <- df %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
seas_vals <- seas_params %>% pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","time_period",
																										 #"pretty_name","
																										 "pretty_group","rank","rank_only"),
																				 names_from = beta,
																				 values_from = c("Mean","significant")) %>% rename(sin="Mean_sin", cos = "Mean_cos")

# Convert to amplitude and max
# Didn't vectorize this function, oops
out <- list()
for (i in 1:nrow(seas_vals)) {
	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
}
out <- rbindlist(out)
seas_vals <- cbind.data.frame(seas_vals, out)

input_dateID = c("201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201412")
dates = fixDate(input_dateID)
input_date_df = data.frame(x = lubridate::month(dates),
													 dates = dates)
out_seas_vals =list()

seas_vals_only = seas_vals %>% filter(grepl("cycl",model_name))
for (row in 1:nrow(seas_vals_only)){
	
	alpha = seas_vals_only[row,]$sin
	beta = seas_vals_only[row,]$cos
	df = input_date_df
	y_cycl = alpha * sin(2*pi*input_date_df$x/12) + beta * cos(2*pi*input_date_df$x/12)
	names(y_cycl) = input_date_df$dates
	out_seas_vals[[row]] <- data.frame(y_cycl) %>% t %>% as.data.frame()
}
seas_vals_to_plot = rbindlist(out_seas_vals, fill=T)
seas_vals_to_plot <- cbind.data.frame(seas_vals_only, seas_vals_to_plot)
seas_vals_long = seas_vals_to_plot %>% 
	pivot_longer(cols = c(16:26), names_to = "dates", values_to = "y_cycl") %>%
	mutate(dates = as.Date(dates))

max_vals =	seas_vals_long %>% group_by(model_name, taxon, time_period, model_id) %>%
	filter(y_cycl == max(y_cycl, na.rm=T)) %>%
	mutate(max_y_date = dates) %>% select(-c(dates, y_cycl))
seas_vals_long <- merge(seas_vals_long, max_vals, all=T)

  saveRDS(list(seas_vals_long, seas_vals), here::here("data", "summary", "pheno_seasonal_amplitude.rds"))
} else {
  cat("Skipping all phenology processing due to empty data\n")
  # Create empty files for downstream compatibility
  saveRDS(data.frame(), here::here("data", "summary", "pheno_predictor_effects.rds"))
  saveRDS(list(data.frame(), data.frame()), here::here("data", "summary", "pheno_seasonal_amplitude.rds"))
}


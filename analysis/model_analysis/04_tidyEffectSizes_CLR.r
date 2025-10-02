# Combine & separately save model parameter and effect size estimates (beta covariates) from all models
# FIXED: Updated for CLR models with proper file paths and CLR-specific logic

source("../../source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects
# FIXED: Use correct CLR summary file path

# Check if CLR summaries exist, otherwise fail
clr_summary_file <- here("data", "summary", "clr_regression_summaries.rds")
if (!file.exists(clr_summary_file)) {
  stop("CLR summary file not found: ", clr_summary_file, "\n",
       "This script should be run after 03_summarizeModelOutputs_CLR.r")
}

cat("Reading CLR summaries from:", clr_summary_file, "\n")
sum.in <- readRDS(clr_summary_file)

# FIXED: Process summary data with proper error handling
if (nrow(sum.in$summary_df) == 0) {
  stop("No summary data available in CLR summaries")
}

sum.all <- sum.in$summary_df %>% 
  filter(model_name != "all_covariates") %>% 
  mutate(tax_rank = rank,
         time_period = recode(time_period, !!!microbialForecast:::date_recode))

df <- sum.all %>%
  mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
  ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
  ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

cat("Processed", nrow(df), "CLR model summaries\n")

# FIXED: For saving: filter by effect type with proper error handling
if (nrow(df) > 0) {
  # Linear model beta (covariate) effects
  beta_effects <- df %>% filter(grepl("beta", rowname))
  
  if (nrow(beta_effects) > 0) {
    # FIXED: Handle CLR-specific beta parameter names
    beta_effects$beta <- ordered(beta_effects$beta, levels = c("sin", "cos",
                                                               "Ectomycorrhizal trees",
                                                               "LAI",
                                                               "pC",
                                                               "pH",
                                                               "Temperature",
                                                               "Moisture","rho"))
    levels(beta_effects$beta)[levels(beta_effects$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
    
    # Save beta effects
    saveRDS(beta_effects, here("data", "summary/clr_predictor_effects.rds"))
    cat("✅ CLR predictor effects saved:", nrow(beta_effects), "effects\n")
  } else {
    cat("WARNING: No beta effects found in CLR summaries\n")
  }
  
  # FIXED: Create CLR-specific site effects summary if available
  if ("site_effect" %in% colnames(df)) {
    site_effects <- df %>% filter(grepl("site", rowname))
    if (nrow(site_effects) > 0) {
      saveRDS(site_effects, here("data", "summary/clr_site_effects.rds"))
      cat("✅ CLR site effects saved:", nrow(site_effects), "effects\n")
    }
  }
  
  # FIXED: Create CLR-specific precision/rho effects summary if available
  if ("precision" %in% colnames(df) || any(grepl("rho", colnames(df)))) {
    precision_effects <- df %>% filter(grepl("precision|rho", rowname))
    if (nrow(precision_effects) > 0) {
      saveRDS(precision_effects, here("data", "summary/clr_precision_effects.rds"))
      cat("✅ CLR precision effects saved:", nrow(precision_effects), "effects\n")
    }
  }
  
  # FIXED: Create CLR-specific seasonality effects summary
  seas_params <- df %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
  if (nrow(seas_params) > 0) {
    seas_vals <- seas_params %>% 
      pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","time_period",
                              "pretty_name","pretty_group","rank","rank_only"),
                  names_from = beta,
                  values_from = c("Mean","significant")) %>% 
      rename(sin="Mean_sin", cos = "Mean_cos")
    
    # Convert to amplitude and max
    # FIXED: Check if sin_cos_to_seasonality function exists
    if (exists("sin_cos_to_seasonality", envir = asNamespace("microbialForecast"))) {
      out <- list()
      for (i in 1:nrow(seas_vals)) {
        out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
      }
      out <- rbindlist(out)
      seas_vals <- cbind.data.frame(seas_vals, out)
      
      # Save seasonality values
      saveRDS(seas_vals, here("data", "summary/clr_seasonality_effects.rds"))
      cat("✅ CLR seasonality effects saved:", nrow(seas_vals), "effects\n")
    } else {
      cat("WARNING: sin_cos_to_seasonality function not found, saving basic seasonality data\n")
      saveRDS(seas_vals, here("data", "summary/clr_seasonality_effects.rds"))
    }
  } else {
    cat("WARNING: No seasonality parameters found in CLR summaries\n")
  }
  
} else {
  stop("No data to process")
}

cat("✅ CLR effect size tidying completed!\n")
cat("Output files saved to data/summary/\n")

# FIXED: Remove commented-out code and provide clean summary
cat("\n=== CLR Effect Size Summary ===\n")
cat("Total models processed:", nrow(df), "\n")
cat("Model types:", paste(unique(df$model_name), collapse=", "), "\n")
cat("Taxonomic ranks:", paste(unique(df$rank_only), collapse=", "), "\n")
cat("Groups:", paste(unique(df$pretty_group), collapse=", "), "\n")


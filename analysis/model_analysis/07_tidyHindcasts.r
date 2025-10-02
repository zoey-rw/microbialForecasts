# Combine hindcasts from all workflows
library(here)
source(here("source.R"))
library(data.table)
library(dplyr)

cat("=== COMBINING HINDCAST DATA ===\n")

# Load observed sites - use most recent files
observed_files <- c(
  here("data/summary/combined_observed_hindcasts_corrected_taxa.rds"),
  here("data/summary/combined_observed_hindcasts_clean.rds"),
  here("data/summary/combined_observed_hindcasts.rds")
)

observed_hindcasts <- data.frame()
for (file in observed_files) {
  if (file.exists(file)) {
    observed_hindcasts <- readRDS(file)
    cat(sprintf("Loaded %d rows from: %s\n", nrow(observed_hindcasts), basename(file)))
    break
  }
}

# Load unobserved sites - use combined files only
unobserved_files <- c(
  here("data/summary/combined_unobserved_hindcasts_optimized.rds"),
  here("data/summary/combined_unobserved_hindcasts.rds")
)

unobserved_hindcasts <- data.frame()
for (file in unobserved_files) {
  if (file.exists(file)) {
    unobserved_hindcasts <- readRDS(file)
    cat(sprintf("Loaded %d rows from: %s\n", nrow(unobserved_hindcasts), basename(file)))
    break
  }
}

# Combine datasets
if (nrow(observed_hindcasts) > 0 && nrow(unobserved_hindcasts) > 0) {
  common_cols <- intersect(colnames(observed_hindcasts), colnames(unobserved_hindcasts))
  all_hindcasts <- rbindlist(list(
    as.data.table(observed_hindcasts)[, ..common_cols],
    as.data.table(unobserved_hindcasts)[, ..common_cols]
  ), fill = TRUE)
} else if (nrow(observed_hindcasts) > 0) {
  all_hindcasts <- as.data.table(observed_hindcasts)
} else if (nrow(unobserved_hindcasts) > 0) {
  all_hindcasts <- as.data.table(unobserved_hindcasts)
} else {
  stop("No hindcast files found!")
}

rm(observed_hindcasts, unobserved_hindcasts)
gc()

cat(sprintf("Combined: %d rows\n", nrow(all_hindcasts)))

# Add dates if missing
if (!"dates" %in% colnames(all_hindcasts)) {
  all_hindcasts[, dates := fixDate(dateID)]
}

# Separate taxonomic and functional groups
tax_mask <- grepl("_bac|_fun", all_hindcasts$rank_name)
tax <- all_hindcasts[tax_mask]
fg <- all_hindcasts[!tax_mask]

# Process taxonomic data
if (nrow(tax) > 0) {
  tax[, `:=`(
    group = ifelse(grepl("16S|bac", rank_name), "16S", "ITS"),
    category = "taxonomic_rank",
    fcast_type = "Taxonomic",
    rank_only = sapply(strsplit(rank_name, "_", fixed = TRUE), `[`, 1)
  )]
}

# Process functional groups
if (nrow(fg) > 0) {
  fg[, `:=`(
    category = assign_fg_categories(species),
    fcast_type = "Functional",
    rank_only = "functional"
  )]
  fg[, group := assign_fg_kingdoms(category)]
}

# Combine
hindcast_data <- rbindlist(list(tax, fg), fill = TRUE)
rm(tax, fg, all_hindcasts)
gc()

# Add derived columns
hindcast_data[, `:=`(
  pretty_group = ifelse(grepl("16S|bac", group), "Bacteria", "Fungi"),
  taxon_name = species,
  taxon = species
)]

# Set rank order
hindcast_data[, rank_only := ordered(rank_only, 
  levels = c("genus", "family", "order", "class", "phylum", "functional"))]

# Fix newsite column
if (!"newsite" %in% colnames(hindcast_data)) {
  hindcast_data[, newsite := ifelse(
    siteID %in% c("ABBY", "BARR", "BONA", "DEJU", "HEAL", "KONA", 
                  "LAJA", "LENO", "MLBS", "RMNP", "SOAP", "TOOL", 
                  "WREF", "YELL"),
    "New site", "Observed site"
  )]
}

hindcast_data[, pretty_name := rank_only]

# Remove "other" taxa
hindcast_data <- hindcast_data[!grepl("other", taxon)]

# Add site_prediction if missing
if (!"site_prediction" %in% colnames(hindcast_data)) {
  hindcast_data[, site_prediction := case_when(
    predicted_site_effect == TRUE & newsite == "New site" ~ "New time x site (modeled effect)",
    predicted_site_effect == FALSE & newsite == "New site" ~ "New time x site (random effect)",
    TRUE ~ "New time (observed site)"
  )]
}

# Add timepoint if missing
if (!"timepoint" %in% colnames(hindcast_data)) {
  hindcast_data[, timepoint := date_num]
}

# CRITICAL: Load plot_start and site_start data for proper calibration filtering
cat("\n=== LOADING MODEL START DATES ===\n")
# Try to get start dates from the model configuration used in Script 6
model_inputs_exist <- FALSE

# Check if we can load from a sample model to get start dates
sample_model_paths <- list.files(
  here("data/model_outputs/reconstructed_from_checkpoints"),
  pattern = "summary_.*_beta_regression.rds",
  recursive = TRUE,
  full.names = TRUE
)[1:2]  # Just check first two

plot_starts_bac <- list()
site_starts_bac <- list()
plot_starts_fun <- list()
site_starts_fun <- list()

for (model_file in sample_model_paths) {
  if (file.exists(model_file)) {
    tryCatch({
      summary_data <- readRDS(model_file)
      if (is.list(summary_data) && "model_data" %in% names(summary_data)) {
        model_data <- summary_data$model_data
        
        # Check if this is bacteria or fungi based on file path
        is_bac <- grepl("16S|bac", model_file, ignore.case = TRUE)
        
        if ("plot_start" %in% names(model_data)) {
          if (is_bac) {
            plot_starts_bac <- model_data$plot_start
          } else {
            plot_starts_fun <- model_data$plot_start
          }
        }
        
        if ("site_start" %in% names(model_data)) {
          if (is_bac) {
            site_starts_bac <- model_data$site_start
          } else {
            site_starts_fun <- model_data$site_start
          }
        }
        
        if (length(plot_starts_bac) > 0 && length(plot_starts_fun) > 0) {
          model_inputs_exist <- TRUE
          break
        }
      }
    }, error = function(e) {
      cat("Could not load start dates from", basename(model_file), "\n")
    })
  }
}

# If we found start dates, merge them
if (model_inputs_exist && (length(plot_starts_bac) > 0 || length(plot_starts_fun) > 0)) {
  cat("Found start dates from model files\n")
  
  # Create data frames from start dates
  if (length(site_starts_bac) > 0) {
    bac_site_start <- data.frame(
      siteID = names(site_starts_bac),
      site_start_date = unlist(site_starts_bac),
      pretty_group = "Bacteria"
    )
  } else {
    bac_site_start <- data.frame()
  }
  
  if (length(plot_starts_bac) > 0) {
    bac_plot_start <- data.frame(
      plotID = names(plot_starts_bac),
      plot_start_date = unlist(plot_starts_bac),
      pretty_group = "Bacteria"
    )
  } else {
    bac_plot_start <- data.frame()
  }
  
  if (length(site_starts_fun) > 0) {
    fun_site_start <- data.frame(
      siteID = names(site_starts_fun),
      site_start_date = unlist(site_starts_fun),
      pretty_group = "Fungi"
    )
  } else {
    fun_site_start <- data.frame()
  }
  
  if (length(plot_starts_fun) > 0) {
    fun_plot_start <- data.frame(
      plotID = names(plot_starts_fun),
      plot_start_date = unlist(plot_starts_fun),
      pretty_group = "Fungi"
    )
  } else {
    fun_plot_start <- data.frame()
  }
  
  # Merge start dates into hindcast_data
  if (nrow(bac_site_start) > 0) {
    hindcast_data <- merge(hindcast_data, bac_site_start, 
                          by = c("siteID", "pretty_group"), all.x = TRUE)
  }
  if (nrow(bac_plot_start) > 0) {
    hindcast_data <- merge(hindcast_data, bac_plot_start, 
                          by = c("plotID", "pretty_group"), all.x = TRUE)
  }
  if (nrow(fun_site_start) > 0) {
    hindcast_data <- merge(hindcast_data, fun_site_start, 
                          by = c("siteID", "pretty_group"), all.x = TRUE)
  }
  if (nrow(fun_plot_start) > 0) {
    hindcast_data <- merge(hindcast_data, fun_plot_start, 
                          by = c("plotID", "pretty_group"), all.x = TRUE)
  }
  
  cat("Merged start dates into hindcast data\n")
} else {
  cat("Warning: Could not load start dates from model files\n")
  cat("Calibration filtering may not exclude first observation dates\n")
}

# Add required columns for script 8 if they don't exist
# Since the combined files were generated with old code, we need to add these columns
if (!"is_site_start_date" %in% colnames(hindcast_data)) {
  hindcast_data[, is_site_start_date := FALSE]
}
if (!"is_plot_start_date" %in% colnames(hindcast_data)) {
  hindcast_data[, is_plot_start_date := FALSE]
}
if (!"is_any_start_date" %in% colnames(hindcast_data)) {
  hindcast_data[, is_any_start_date := FALSE]
}
if (!"site_start_date" %in% colnames(hindcast_data)) {
  hindcast_data[, site_start_date := NA_real_]
}
if (!"plot_start_date" %in% colnames(hindcast_data)) {
  hindcast_data[, plot_start_date := NA_real_]
}

# Add the new columns that were added to hindcast functions
if (!"plot_start" %in% colnames(hindcast_data)) {
  hindcast_data[, plot_start := NA_real_]
}
if (!"site_start" %in% colnames(hindcast_data)) {
  hindcast_data[, site_start := NA_real_]
}

# Sort
setorder(hindcast_data, model_name, fcast_type, pretty_group, taxon, plotID, date_num)

# Save
saveRDS(hindcast_data, here("data/summary/all_hindcasts_plsr2.rds"))
cat("\n✓ Saved: all_hindcasts_plsr2.rds\n")
cat(sprintf("Final output: %d rows, %d columns\n", nrow(hindcast_data), ncol(hindcast_data)))

# Save parquet if arrow available
if (require(arrow, quietly = TRUE)) {
  dir.create(here("data/summary/parquet"), showWarnings = FALSE, recursive = TRUE)
  write_parquet(hindcast_data, here("data/summary/parquet/all_hindcasts_plsr2.parquet"))
  cat("✓ Saved: all_hindcasts_plsr2.parquet\n")
}

cat("\n=== HINDCAST TIDYING COMPLETE ===\n")
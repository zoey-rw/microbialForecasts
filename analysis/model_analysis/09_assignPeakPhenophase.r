# read in seasonal values
library(lubridate)
library(here)
library(dplyr)
source(here("source.R"))
# Phenology functions (assign_pheno_*) loaded via microbialForecast package

# Load data from 03_summarizeModelOutputs.r outputs
# Step 3 creates summary_df, plot_est, and convergence lists
seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))

# Note: Using newest summary file for consistency (though not directly used here)
# The script loads plot estimates directly from parquet chunks, which are the newest data
# sum.in_pheno <- readRDS(here("data", "summary/pheno_summaries.rds"))  # OLD - not used
# seas_in_pheno = readRDS(here("data/summary/pheno_seasonal_amplitude.rds"))  # OLD - not used

# Load weakly converged models list from step 3 output first
keep_models_weak <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))

# Normalize converged model IDs to handle _beta_regression suffix mismatch
keep_models_weak_normalized <- gsub("_beta_regression$", "", keep_models_weak)

# Create a data frame with model_id column for filtering
keep_models_weak_df <- data.frame(model_id = keep_models_weak, stringsAsFactors = FALSE)

seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% keep_models_weak_df$model_id)
seas_vals_wide <- seas_in[[2]] %>% filter(model_id %in% keep_models_weak_df$model_id)

# Load plot_estimates from the outputs of 03c_extract_plot_summaries.r and 03d_combine_plot_chunks.r
# Use the package function if available, otherwise load from combined file
cat("Loading plot_estimates from summary files...\n")

# Prefer combined file (from 03d_combine_plot_chunks.r) as it has properly joined Mean and metadata
# Chunk files have Mean and metadata in separate rows, which is incorrect
plot_estimates_file <- here("data/summary/plot_estimates.rds")
if(!file.exists(plot_estimates_file)) {
  plot_estimates_file <- here("data/summary/plot_estimates.parquet")
}

use_chunk_loading <- FALSE
if(file.exists(plot_estimates_file)) {
  cat("Found combined plot_estimates file - using it (has properly joined Mean and metadata)\n")
  use_chunk_loading <- FALSE
} else {
  # Fallback to chunks only if combined file doesn't exist
  chunk_dirs <- c("plot_estimates_env_cycl", "plot_estimates_cycl_only", "plot_estimates_env_cov")
  chunk_dirs_exist <- sapply(chunk_dirs, function(d) dir.exists(here("data/summary", d)))
  
  if(any(chunk_dirs_exist)) {
    cat("⚠️ Combined file not found - using chunks (note: chunks may have Mean/metadata separation issue)\n")
    cat("Found chunk directories, using chunk-based loading...\n")
    use_chunk_loading <- TRUE
  
  # load_plot_estimates() is provided by the microbialForecast package
  if(!exists("load_plot_estimates", mode="function")) {
    stop("load_plot_estimates() not found. Ensure microbialForecast package is loaded via source.R")
  }
  
  # Load from all available model types
  all_plot_data <- list()
  for(model_type in c("env_cycl", "cycl_only", "env_cov")) {
    chunk_dir <- here("data/summary", paste0("plot_estimates_", model_type))
    if(dir.exists(chunk_dir)) {
      cat("Loading", model_type, "from chunks...\n")
      tryCatch({
        # Load ALL columns (not subset) to ensure we get all available data including Mean
        # The subset_cols approach was filtering out rows that don't have those columns
        if(exists("load_plot_estimates", mode="function")) {
          cat("  Loading all columns for", model_type, "(not subsetting)...\n")
          model_data <- load_plot_estimates(model_type = model_type, subset_cols = NULL)
          cat("  Loaded", nrow(model_data), "rows with", ncol(model_data), "columns\n")
        } else {
          # Manual loading from chunks
          chunk_files <- list.files(chunk_dir, pattern = "\\.parquet$", full.names = TRUE)
          if(length(chunk_files) > 0) {
            cat("  Loading", length(chunk_files), "chunks manually...\n")
            chunk_list <- list()
            for(chunk_file in chunk_files) {
              if(requireNamespace("arrow", quietly = TRUE)) {
                chunk_data <- arrow::read_parquet(chunk_file)
              } else if(requireNamespace("nanoparquet", quietly = TRUE)) {
                chunk_data <- nanoparquet::read_parquet(chunk_file)
              } else {
                stop("Neither arrow nor nanoparquet available")
              }
              chunk_list[[length(chunk_list) + 1]] <- data.table::as.data.table(chunk_data)
            }
            model_data <- rbindlist(chunk_list, fill = TRUE, use.names = TRUE)
          } else {
            model_data <- data.table()
          }
        }
        
        if(nrow(model_data) > 0) {
          all_plot_data[[model_type]] <- model_data
          cat("  ✅ Loaded", nrow(model_data), "rows for", model_type, "\n")
        }
      }, error = function(e) {
        cat("  ❌ Error loading", model_type, ":", e$message, "\n")
      })
    }
  }
  
    if(length(all_plot_data) > 0) {
      plot_estimates <- rbindlist(all_plot_data, fill = TRUE, use.names = TRUE)
      cat("✅ Combined", length(all_plot_data), "model types into", nrow(plot_estimates), "rows\n")
      cat("   Mean non-NA:", sum(!is.na(plot_estimates$Mean)), "out of", nrow(plot_estimates), "\n")
      cat("   siteID non-NA:", sum(!is.na(plot_estimates$siteID)), "out of", nrow(plot_estimates), "\n")
      cat("   dates non-NA:", sum(!is.na(plot_estimates$dates)), "out of", nrow(plot_estimates), "\n")
      cat("   model_id non-NA:", sum(!is.na(plot_estimates$model_id)), "out of", nrow(plot_estimates), "\n")
    } else {
      cat("⚠️ No data loaded from chunks, falling back to combined file...\n")
      use_chunk_loading <- FALSE
    }
  }
}

# Load from combined file (preferred - has properly joined data)
if(!use_chunk_loading) {
  if(!file.exists(plot_estimates_file)) {
    stop("Plot estimates file not found. Please run 03c_extract_plot_summaries.r and 03d_combine_plot_chunks.r first.\n",
         "Expected files: data/summary/plot_estimates.parquet or data/summary/plot_estimates.rds")
  }
  
  # Check file size - if too small, it's probably incomplete
  file_size_mb <- file.info(plot_estimates_file)$size / 1048576
  cat("Loading combined plot_estimates file (", round(file_size_mb, 2), "MB)...\n", sep="")
  if(file_size_mb < 10) {
    cat("Warning: plot_estimates file is very small. It may be incomplete.\n")
    cat("Consider running 03d_combine_plot_chunks.r to regenerate the combined file.\n")
  }
  
  # Load plot_estimates
  if(grepl("\\.parquet$", plot_estimates_file)) {
    if(requireNamespace("arrow", quietly = TRUE)) {
      plot_estimates <- arrow::read_parquet(plot_estimates_file)
      plot_estimates <- data.table::as.data.table(plot_estimates)
    } else if(requireNamespace("nanoparquet", quietly = TRUE)) {
      plot_estimates <- nanoparquet::read_parquet(plot_estimates_file)
      plot_estimates <- data.table::as.data.table(plot_estimates)
    } else {
      stop("Neither arrow nor nanoparquet available to read parquet files")
    }
  } else {
    plot_estimates <- readRDS(plot_estimates_file)
    plot_estimates <- data.table::as.data.table(plot_estimates)
  }
  
  cat("✅ Loaded", nrow(plot_estimates), "plot estimate rows from combined file\n")
  cat("   Rows with Mean and model_id:", sum(!is.na(plot_estimates$Mean) & !is.na(plot_estimates$model_id)), "\n")
}

# Normalize model_ids for matching
if(data.table::is.data.table(plot_estimates)) {
  plot_estimates[, model_id_normalized := gsub("_beta_regression$|_combined$", "", model_id)]
} else {
  plot_estimates$model_id_normalized <- gsub("_beta_regression$|_combined$", "", plot_estimates$model_id)
}

cat("Before filtering to converged models:\n")
cat("  Total rows:", nrow(plot_estimates), "\n")
cat("  Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("  Unique normalized model_ids:", length(unique(plot_estimates$model_id_normalized)), "\n")

# Filter to converged models
if(data.table::is.data.table(plot_estimates)) {
  plot_estimates <- plot_estimates[model_id_normalized %in% keep_models_weak_normalized]
} else {
  plot_estimates <- plot_estimates[plot_estimates$model_id_normalized %in% keep_models_weak_normalized, ]
}

cat("After filtering to converged models:\n")
cat("  Total rows:", nrow(plot_estimates), "rows\n")
cat("  Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("  Unique normalized model_ids:", length(unique(plot_estimates$model_id_normalized)), "\n")

# Ensure Mean column exists (use 50% or med if Mean is missing or all NA)
if(!"Mean" %in% names(plot_estimates)) {
  # Mean column doesn't exist - try to create it from quantiles
  if("50%" %in% names(plot_estimates)) {
    if(data.table::is.data.table(plot_estimates)) {
      plot_estimates[, Mean := `50%`]
    } else {
      plot_estimates$Mean <- plot_estimates$`50%`
    }
    cat("Created Mean from '50%' column\n")
  } else if("med" %in% names(plot_estimates)) {
    if(data.table::is.data.table(plot_estimates)) {
      plot_estimates[, Mean := med]
    } else {
      plot_estimates$Mean <- plot_estimates$med
    }
    cat("Created Mean from 'med' column\n")
  } else {
    # Create Mean column filled with NA if no source available
    if(data.table::is.data.table(plot_estimates)) {
      plot_estimates[, Mean := NA_real_]
    } else {
      plot_estimates$Mean <- NA_real_
    }
    cat("Created Mean column filled with NA (no quantile source available)\n")
  }
} else {
  # Mean column exists - check if it's all NA and try to fill from quantiles
  mean_non_na_count <- sum(!is.na(plot_estimates$Mean), na.rm=TRUE)
  if(mean_non_na_count == 0) {
    cat("Mean column exists but all values are NA. Attempting to fill from quantiles...\n")
    if("50%" %in% names(plot_estimates)) {
      non_na_50 <- sum(!is.na(plot_estimates$`50%`), na.rm=TRUE)
      if(non_na_50 > 0) {
        if(data.table::is.data.table(plot_estimates)) {
          plot_estimates[is.na(Mean), Mean := `50%`]
        } else {
          plot_estimates$Mean[is.na(plot_estimates$Mean)] <- plot_estimates$`50%`[is.na(plot_estimates$Mean)]
        }
        cat("Filled", sum(!is.na(plot_estimates$Mean)), "Mean values from '50%' column\n")
      }
    } else if("med" %in% names(plot_estimates)) {
      non_na_med <- sum(!is.na(plot_estimates$med), na.rm=TRUE)
      if(non_na_med > 0) {
        if(data.table::is.data.table(plot_estimates)) {
          plot_estimates[is.na(Mean), Mean := med]
        } else {
          plot_estimates$Mean[is.na(plot_estimates$Mean)] <- plot_estimates$med[is.na(plot_estimates$Mean)]
        }
        cat("Filled", sum(!is.na(plot_estimates$Mean)), "Mean values from 'med' column\n")
      }
    }
  }
}

# Check metadata completeness
required_metadata <- c("siteID", "plotID", "dates", "model_id")
missing_metadata <- sapply(required_metadata, function(col) {
  if(col %in% names(plot_estimates)) {
    sum(is.na(plot_estimates[[col]]))
  } else {
    nrow(plot_estimates)  # Column doesn't exist
  }
})

cat("Metadata completeness:\n")
cat("  Total rows:", nrow(plot_estimates), "\n")
for(i in seq_along(required_metadata)) {
  col <- required_metadata[i]
  if(col %in% names(plot_estimates)) {
    cat("  ", col, ": ", sum(!is.na(plot_estimates[[col]])), " non-NA out of ", nrow(plot_estimates), "\n", sep="")
  } else {
    cat("  ", col, ": MISSING COLUMN\n", sep="")
  }
}

# Filter out rows with missing critical metadata (but keep rows with NA Mean if metadata is present)
# For phenology analysis, we need at least siteID, dates, and model_id
if(data.table::is.data.table(plot_estimates)) {
  plot_estimates <- plot_estimates[!is.na(siteID) & !is.na(dates) & !is.na(model_id)]
} else {
  plot_estimates <- plot_estimates[!is.na(plot_estimates$siteID) & !is.na(plot_estimates$dates) & !is.na(plot_estimates$model_id), ]
}

cat("After filtering to rows with complete metadata (siteID, dates, model_id):\n")
cat("  Total rows:", nrow(plot_estimates), "\n")
cat("  Rows with non-NA Mean:", sum(!is.na(plot_estimates$Mean)), "\n")
cat("  Rows with NA Mean (but metadata present):", sum(is.na(plot_estimates$Mean)), "\n")

# If all Mean values are NA, we can't do peak detection
# This suggests the plot_estimates from summary files don't have Mean populated
# We need Mean values for peak detection, so check if we have any
if(sum(!is.na(plot_estimates$Mean)) == 0) {
  cat("\n⚠️ WARNING: All Mean values are NA in plot_estimates!\n")
  cat("   This means we cannot perform peak detection.\n")
  cat("   The plot_estimates from summary files may not have Mean values populated.\n")
  cat("   Consider using hindcast data or regenerating plot_estimates with Mean values.\n")
  cat("   For now, we will proceed but peak detection will find 0 models.\n")
}

# Ensure dates column exists
if(!"dates" %in% names(plot_estimates)) {
  if("dateID" %in% names(plot_estimates)) {
    if(exists("fixDate", mode="function")) {
      plot_estimates[, dates := fixDate(dateID)]
    } else {
      plot_estimates[, dates := as.Date(paste0(as.character(dateID), "01"), format = "%Y%m%d")]
    }
  } else if("date_num" %in% names(plot_estimates)) {
    plot_estimates[, dates := as.Date(date_num, origin = "1970-01-01")]
  } else {
    stop("No date column found (dates, dateID, or date_num)")
  }
}

cat("✅ Loaded plot_estimates:", nrow(plot_estimates), "rows\n")
cat("Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("Non-NA Mean:", sum(!is.na(plot_estimates$Mean)), "rows\n")

# Plot estimates are already filtered to converged models and have Mean from hindcast data
# Add month and year columns for phenology analysis
cat("\nPreparing plot estimates for phenology analysis...\n")
cat("Total rows:", nrow(plot_estimates), "\n")
cat("Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("Non-NA Mean:", sum(!is.na(plot_estimates$Mean)), "rows\n")

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>%
	select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)
cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")

# Convert to data.frame for dplyr operations if needed
if(data.table::is.data.table(plot_estimates)) {
  plot_estimates <- as.data.frame(plot_estimates)
}

plot_estimates = plot_estimates %>% filter(model_name != "all_covariates")

# Extract month and year from dates column
plot_estimates <- plot_estimates %>%
  mutate(
    month = lubridate::month(dates),
    year = lubridate::year(dates)
  )

cat("After filtering and adding month/year:\n")
cat("  Total rows:", nrow(plot_estimates), "\n")
cat("  Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("  Years available:", paste(sort(unique(plot_estimates$year)), collapse=", "), "\n")
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cov_est = plot_estimates %>% filter(grepl("env_cov",model_name))


pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
# Use all siteIDs from plot_estimates, not just cycl_only_est
pheno_categories_long = pheno_categories_in[[2]] %>% filter(ID %in% unique(plot_estimates$siteID))
cat("Phenology categories filtered to", length(unique(pheno_categories_long$ID)), "sites from plot_estimates\n")


max_dates = seas_vals_short %>% mutate(dates = max_y_date)


# Process plot estimates in smaller chunks to avoid memory issues
cat("Processing plot estimates for peak abundance analysis in chunks...\n")

# Get unique years available in the data and process them
# Use ALL years that have Mean data (not just hardcoded years)
available_years <- sort(unique(lubridate::year(plot_estimates$dates[!is.na(plot_estimates$Mean)])))
available_years <- as.character(available_years[!is.na(available_years)])
cat("Available years with Mean data:", paste(available_years, collapse=", "), "\n")
# Process ALL available years to get data for all converged models
target_years <- if(length(available_years) > 0) {
  available_years
} else {
  c("2014", "2016", "2017")  # Fallback to default if no years found
}
cat("Processing years:", paste(target_years, collapse=", "), "\n")
cat("Total years to process:", length(target_years), "\n")
# Assign phenophase to every monthly row, then derive per-site-year peak rows by
# slicing. Running the per-site loop inside the year loop bounds memory: only
# one year of monthly data is labeled at a time.
assign_pheno_to_df <- function(df, pheno_categories_long) {
  df <- df %>%
    mutate(dates = ymd(paste0(year, "-", sprintf("%02d", month), "-15"))) %>%
    filter(siteID %in% pheno_categories_long$ID)

  cat("Processing", nrow(df), "rows with vectorized phenology assignment...\n")

  df$site_cat <- NA_character_
  unique_sites <- unique(df$siteID)
  cat("Processing", length(unique_sites), "unique sites...\n")

  df$dates_ymd <- ymd(df$dates)

  for(i in seq_along(unique_sites)) {
    site <- unique_sites[i]
    site_rows <- which(df$siteID == site)

    if(length(site_rows) > 0) {
      site_dates_ymd <- df$dates_ymd[site_rows]
      site_pheno <- pheno_categories_long %>% filter(ID == site)

      if(nrow(site_pheno) > 0) {
        valid_dates <- !is.na(site_dates_ymd)
        if(sum(valid_dates) > 0) {
          site_cats <- character(length(site_dates_ymd))
          for(j in seq_along(site_dates_ymd)) {
            if(valid_dates[j]) {
              keep_index <- site_dates_ymd[j] %within% site_pheno$value
              out_cats <- site_pheno[keep_index, ]$name
              if(length(out_cats) > 0) {
                out_cats <- gsub("_interval", "", out_cats)
                out_cats <- gsub("dormancy_interval[12]", "dormancy", out_cats)
                unique_cats <- unique(out_cats)

                # If a date falls in overlapping intervals, prefer the
                # narrower phenophase. Peak > greenup > greendown > dormancy
                # avoids labeling a summer date "dormancy" when both apply.
                if(length(unique_cats) > 1) {
                  priority_order <- c("peak", "greenup", "greendown", "dormancy")
                  selected_cat <- NULL
                  for(priority_cat in priority_order) {
                    if(priority_cat %in% unique_cats) {
                      selected_cat <- priority_cat
                      break
                    }
                  }
                  site_cats[j] <- if(!is.null(selected_cat)) selected_cat else unique_cats[1]
                } else {
                  site_cats[j] <- unique_cats[1]
                }
              } else {
                site_cats[j] <- NA_character_
              }
            } else {
              site_cats[j] <- NA_character_
            }
          }
          df$site_cat[site_rows] <- site_cats
        }
      }
    }

    if(i %% 10 == 0) {
      cat("  Processed site", i, "of", length(unique_sites), "\n")
    }
  }

  df$dates_ymd <- NULL
  df %>% filter(!is.na(site_cat))
}

# Process one year at a time: summarize to monthly, assign phenophase, store.
# This caps memory at ~1/N_years of the full monthly dataset.
monthly_chunks <- list()

for(year in target_years) {
  cat("\nProcessing year:", year, "\n")

  median_col <- if("50%" %in% names(plot_estimates)) "50%" else "Mean"
  year_num <- as.numeric(year)

  year_rows_before <- plot_estimates %>% filter(year == year_num)
  cat("  Rows before NA filter:", nrow(year_rows_before), "\n")

  year_data <- plot_estimates %>%
    filter(year == year_num) %>%
    filter(!is.na(.data[[median_col]]))

  cat("  Rows after NA filter:", nrow(year_data), "\n")
  cat("  Unique model_ids:", length(unique(year_data$model_id)), "\n")
  cat("  Unique siteIDs:", length(unique(year_data$siteID)), "\n")

  monthly_year <- year_data %>%
    group_by(model_id, siteID, year, month, dates) %>%
    summarize(mean_modeled_abun = mean(.data[[median_col]], na.rm = TRUE),
              .groups = "drop") %>%
    filter(!is.na(mean_modeled_abun))

  cat("  Rows after summarization (all months):", nrow(monthly_year), "\n")

  if(nrow(monthly_year) > 0) {
    monthly_year <- assign_pheno_to_df(monthly_year, pheno_categories_long)
    cat("  Rows after phenology assignment:", nrow(monthly_year), "\n")
    monthly_chunks[[year]] <- monthly_year
  }

  rm(year_data, monthly_year)
  gc()
}

if(length(monthly_chunks) > 0) {
  all_monthly_abun <- bind_rows(monthly_chunks)
  rm(monthly_chunks); gc()
  cat("\nCombined monthly abundance data (labeled):",
      nrow(all_monthly_abun), "rows\n")

  # Derive per-site-year peak rows by slicing the labeled monthly data.
  max_abun <- all_monthly_abun %>%
    group_by(model_id, siteID, year) %>%
    filter(mean_modeled_abun == max(mean_modeled_abun, na.rm = TRUE)) %>%
    ungroup()
  cat("Derived peak rows from labeled monthly data:", nrow(max_abun), "rows\n")
} else {
  all_monthly_abun <- data.frame()
  max_abun <- data.frame()
}

# About 6 sites are missing all MODIS phenophase data; these "NA" values will be excluded
# About 7000 additional dates have some missing data
pheno_levels <- c("dormancy", "greenup", "peak", "greendown")
if(nrow(max_abun) > 0) {
  max_abun$sampling_season <- factor(max_abun$site_cat, ordered = TRUE, levels = pheno_levels)
}
if(nrow(all_monthly_abun) > 0) {
  all_monthly_abun$sampling_season <- factor(all_monthly_abun$site_cat, ordered = TRUE, levels = pheno_levels)
}



# Build per-model metadata for downstream figures (group columns, amplitude, and
# seasonality-significance flags). The peak-of-peaks per-guild "peak phenophase"
# that used to be computed here -- an abundance-weighted argmax over per-site-year
# peak months (seasonality_mode2), plus a >50%-frequency variant (seasonality_mode3)
# -- has been REMOVED. It was not robust (it disagreed with the mean-abundance peak
# definition and with itself across sites) and no downstream script used it;
# consumers read only metadata from element [[1]].
if(nrow(max_abun) > 0) {
  cat("Building per-model metadata table...\n")

  # Amplitude + significance flags keyed by normalized model_id
  max_dates_for_merge <- max_dates %>%
    mutate(model_id_normalized = gsub("_beta_regression$|_combined$", "", model_id)) %>%
    select(model_id_normalized, fcast_type, pretty_group, rank_only, model_name, taxon, amplitude, significant_sin, significant_cos) %>%
    distinct(model_id_normalized, .keep_all = TRUE)

  cat("Extracting metadata from plot_estimates...\n")
  if(!data.table::is.data.table(plot_estimates)) {
    plot_estimates <- data.table::as.data.table(plot_estimates)
  }
  plot_estimates_unique <- unique(plot_estimates[, .(model_id, model_id_normalized, model_name, taxon, fcast_type, pretty_group, rank_only, species)])
  plot_estimates_df <- as.data.frame(plot_estimates_unique)

  fg_names <- microbialForecast:::keep_fg_names
  all_bacteria <- unlist(microbialForecast:::rank_spec_names[grepl("_bac$", names(microbialForecast:::rank_spec_names))])
  all_fungi <- unlist(microbialForecast:::rank_spec_names[grepl("_fun$", names(microbialForecast:::rank_spec_names))])

  plot_estimates_df <- plot_estimates_df %>%
    mutate(
      fcast_type = ifelse(is.na(fcast_type) & !is.na(taxon),
                         ifelse(taxon %in% fg_names, "Functional", "Taxonomic"),
                         fcast_type),
      pretty_group = ifelse(is.na(pretty_group) & !is.na(taxon),
                           ifelse(taxon %in% fg_names,
                                 ifelse(taxon %in% all_fungi, "Fungi", "Bacteria"),
                                 ifelse(taxon %in% all_fungi, "Fungi",
                                       ifelse(taxon %in% all_bacteria, "Bacteria", NA))),
                           pretty_group)
    )

  plot_estimates_metadata <- plot_estimates_df %>%
    mutate(model_id_normalized = gsub("_beta_regression$|_combined$", "", model_id)) %>%
    filter(!is.na(fcast_type) | !is.na(pretty_group) | !is.na(model_name)) %>%
    select(model_id_normalized, fcast_type, pretty_group, rank_only, model_name, taxon) %>%
    arrange(model_id_normalized, !is.na(fcast_type), !is.na(pretty_group), !is.na(model_name)) %>%
    distinct(model_id_normalized, .keep_all = TRUE)

  all_metadata <- plot_estimates_metadata %>%
    full_join(max_dates_for_merge %>% select(model_id_normalized, amplitude, significant_sin, significant_cos),
              by = "model_id_normalized") %>%
    mutate(
      fcast_type = ifelse(is.na(fcast_type) & model_id_normalized %in% max_dates_for_merge$model_id_normalized,
                         max_dates_for_merge$fcast_type[match(model_id_normalized, max_dates_for_merge$model_id_normalized)],
                         fcast_type),
      pretty_group = ifelse(is.na(pretty_group) & model_id_normalized %in% max_dates_for_merge$model_id_normalized,
                           max_dates_for_merge$pretty_group[match(model_id_normalized, max_dates_for_merge$model_id_normalized)],
                           pretty_group),
      rank_only = ifelse(is.na(rank_only) & model_id_normalized %in% max_dates_for_merge$model_id_normalized,
                        max_dates_for_merge$rank_only[match(model_id_normalized, max_dates_for_merge$model_id_normalized)],
                        rank_only),
      model_name = ifelse(is.na(model_name) & model_id_normalized %in% max_dates_for_merge$model_id_normalized,
                         max_dates_for_merge$model_name[match(model_id_normalized, max_dates_for_merge$model_id_normalized)],
                         model_name),
      taxon = ifelse(is.na(taxon) & model_id_normalized %in% max_dates_for_merge$model_id_normalized,
                    max_dates_for_merge$taxon[match(model_id_normalized, max_dates_for_merge$model_id_normalized)],
                    taxon)
    ) %>%
    distinct(model_id_normalized, .keep_all = TRUE)

  # Element [[1]]: per-model metadata (model_id + group cols + amplitude + sig flags)
  model_metadata <- all_metadata %>% rename(model_id = model_id_normalized)

  # Element [[4]]: per-site-year peak month rows + metadata
  max_abun_normalized <- max_abun %>%
    mutate(model_id_normalized = gsub("_beta_regression$|_combined$", "", model_id))
  max_abun_to_plot <- merge(max_abun_normalized, all_metadata,
                            by.x = "model_id_normalized", by.y = "model_id_normalized", all.x = TRUE) %>%
    select(-model_id_normalized)

  # Element [[6]]: every model x site x month estimate + phenophase + metadata
  if(nrow(all_monthly_abun) > 0) {
    all_monthly_abun_normalized <- all_monthly_abun %>%
      mutate(model_id_normalized = gsub("_beta_regression$|_combined$", "", model_id))
    all_monthly_abun_to_plot <- merge(all_monthly_abun_normalized, all_metadata,
                                      by.x = "model_id_normalized", by.y = "model_id_normalized", all.x = TRUE) %>%
      select(-model_id_normalized)
  } else {
    all_monthly_abun_to_plot <- data.frame()
  }
} else {
  cat("No max_abun data available\n")
  model_metadata <- data.frame()
  max_abun_to_plot <- data.frame()
  all_monthly_abun_to_plot <- data.frame()
}

cat("\n=== SAVING RESULTS ===\n")
# Element layout. The former peak-of-peaks summaries were removed; [[3]] and [[5]]
# are kept as empty placeholders so positional indexing of [[4]] and [[6]] in
# downstream scripts is unchanged.
#   [[1]] model_metadata            -- one row per model: group cols, amplitude, sig flags
#   [[2]] max_abun                  -- per site-year peak month with phenophase
#   [[3]] (removed: seasonality_mode3)
#   [[4]] max_abun_to_plot          -- max_abun + metadata (peak-month view)
#   [[5]] (removed: seasonality_mode_to_plot)
#   [[6]] all_monthly_abun_to_plot  -- every model x site x month estimate +
#                                      phenophase + metadata (full seasonal profile)
saveRDS(list(model_metadata, max_abun, data.frame(),
             max_abun_to_plot, data.frame(),
             all_monthly_abun_to_plot),
        here("data/clean/pheno_group_peak_phenophases.rds"))

cat("\n=== SUMMARY ===\n")
cat("✅ Used hindcast data for calibration period plot estimates\n")
cat("   Total rows:", nrow(plot_estimates), "\n")
cat("   Non-NA Mean:", sum(!is.na(plot_estimates$Mean)), "rows\n")
cat("   Unique model_ids:", length(unique(plot_estimates$model_id)), "\n")
cat("\nScript completed successfully!\n")

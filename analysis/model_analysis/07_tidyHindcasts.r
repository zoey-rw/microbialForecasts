#!/usr/bin/env Rscript

# ==============================================================================
# HINDCAST PROCESSING PIPELINE (OPTIMIZED)
# ==============================================================================
# Drop-in replacement for previous processing script.
# Uses Arrow Datasets to handle memory-efficient, resumable processing.
# ==============================================================================

# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
tryCatch({
  mem.maxVSize(Inf)
  cat("Memory limit increased to unlimited\n")
}, error = function(e) {
  cat("Note: Could not increase memory limit (OS dependent)\n")
})

# Ensure nanoparquet is available (Lightweight, memory-efficient)
if (!requireNamespace("nanoparquet", quietly = TRUE)) {
  cat("nanoparquet package not available. Attempting installation...\n")
  tryCatch({
    install.packages("nanoparquet", repos = "https://cran.rstudio.com/", dependencies = FALSE, quiet = TRUE)
    if (!requireNamespace("nanoparquet", quietly = TRUE)) {
      stop("nanoparquet installation failed")
    }
  }, error = function(e) {
    cat("⚠️  nanoparquet installation failed.\n")
    cat("   Please install manually: install.packages('nanoparquet')\n")
    stop("nanoparquet package required but not available")
  })
}

library(here)
library(data.table)
library(dplyr)
library(nanoparquet)

# Source external helper if it exists (for fixDate), otherwise define safe fallback
if (file.exists(here("source.R"))) source(here("source.R"))

# --- CONSTANTS & PATHS ---
BATCH_SIZE   <- 50L  # Number of files to process before flushing to disk
SUMMARY_DIR  <- here("data/summary")
PARQUET_DIR  <- here("data/summary/parquet")
DATASET_DIR  <- here("data/summary/parquet/_dataset_staging") # Intermediate parquet storage
LOG_FILE     <- here("data/summary/parquet/_processed_log.rds") # Resume capability

# Final Output Targets (Matches original script)
OUT_RDS_COMBINED   <- here("data/summary/all_hindcasts_plsr2.rds")
OUT_RDS_OBSERVED   <- here("data/summary/all_hindcasts_observed_plsr2.rds")
OUT_RDS_UNOBSERVED <- here("data/summary/all_hindcasts_unobserved_plsr2.rds")
OUT_PARQUET_FULL   <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")

# Ensure directories exist
dir.create(DATASET_DIR, recursive = TRUE, showWarnings = FALSE)

# 2. HELPER FUNCTIONS (CLEANING & RESUME LOGIC)
# ------------------------------------------------------------------------------

# --- Resume Logic ---
get_processed_log <- function() {
  if (!file.exists(LOG_FILE)) {
    return(data.table(filename = character(0), mtime = as.POSIXct(character(0))))
  }
  log <- readRDS(LOG_FILE)
  # Backwards compatibility: convert old character vector format to data.table
  if (is.character(log)) {
    cat("Converting old log format (character vector) to new format (filename + mtime).\n")
    log <- data.table(filename = log, mtime = as.POSIXct("2000-01-01"))
  }
  if (!is.data.table(log)) log <- as.data.table(log)
  log
}

update_processed_log <- function(filenames, mtimes) {
  current <- get_processed_log()
  new_entries <- data.table(filename = filenames, mtime = mtimes)
  # Remove any existing entries for these filenames, then add new ones
  updated <- rbindlist(list(current[!filename %in% filenames], new_entries))
  saveRDS(updated, LOG_FILE)
}

# --- Data Cleaning ---
# Logic ported from original script to ensure consistent schema
clean_hindcast_dt <- function(dt, filename) {
  # 1. Normalize Site Prediction (vectorized so each row uses its own predicted_site_effect and newsite)
  if (!"site_prediction" %in% names(dt) || any(is.na(dt$site_prediction))) {
    ns  <- if ("newsite" %in% names(dt)) dt$newsite else NA
    pse <- if ("predicted_site_effect" %in% names(dt)) dt$predicted_site_effect else NA
    n   <- nrow(dt)
    if (length(ns) == 1L && is.na(ns)) ns <- rep(NA, n)
    if (length(pse) == 1L && is.na(pse)) pse <- rep(NA, n)
    ns_norm <- data.table::fcase(
      is.character(ns) & tolower(ns) %in% c("new site","new_site","newsite","true"), "New site",
      is.character(ns), "Observed site",
      !is.na(ns) & (ns == TRUE | ns == 1L), "New site",
      !is.na(ns) & (ns == FALSE | ns == 0L), "Observed site",
      default = NA_character_
    )
    if (!"site_prediction" %in% names(dt)) dt[, site_prediction := NA_character_]
    idx <- is.na(dt$site_prediction)
    if (any(idx)) {
      dt[idx, site_prediction := data.table::fcase(
        pse[idx] == TRUE  & ns_norm[idx] == "New site", "New time x site (modeled effect)",
        pse[idx] == FALSE & ns_norm[idx] == "New site", "New time x site (random effect)",
        default = "New time (observed site)"
      )]
    }
  }

  # 2. Fix Dates
  if (!"dates" %in% names(dt) && "dateID" %in% names(dt)) {
    if (exists("fixDate", mode="function")) {
      dt[, dates := fixDate(dateID)]
    } else {
      # Fallback: convert dateID (YYYYMM) to Date (first of month)
      # dateID is numeric like 201801 for Jan 2018
      dt[, dates := as.Date(paste0(as.character(dateID), "01"), format = "%Y%m%d")]
    }
  }
  if (!"timepoint" %in% names(dt) && "date_num" %in% names(dt)) dt[, timepoint := date_num]

  # 3. Model Name Extraction
  if (!"model_name" %in% names(dt) && "model_id" %in% names(dt)) {
    # Extract model_name from model_id with proper handling of multi-part names
    dt[, model_name := {
      parts_list <- strsplit(model_id, "_", fixed = TRUE)
      sapply(parts_list, function(parts) {
        if (length(parts) >= 2) {
          if (parts[1] == "cycl" && parts[2] == "only") return("cycl_only")
          if (parts[1] == "env" && parts[2] == "cov") return("env_cov")
          if (parts[1] == "env" && parts[2] == "cycl") return("env_cycl")
        }
        if (length(parts) >= 1) return(parts[1])
        return(NA_character_)
      })
    }]
  }

  # 4. Taxonomy / Pretty Groups — use the canonical package function.
  dt <- fill_pretty_group(dt)
  
  if (!"rank_only" %in% names(dt) && "rank_name" %in% names(dt)) {
    dt[, rank_only := tstrsplit(rank_name, "_", keep = 1)]
    rank_levels <- c("genus", "family", "order", "class", "phylum", "functional", "diversity")
    # Only set factor if valid levels exist
    valid_ranks <- intersect(unique(dt$rank_only), rank_levels)
    if (length(valid_ranks) > 0) {
      dt[, rank_only := factor(rank_only, levels = rank_levels, ordered = TRUE)]
    }
  }

  # 5. Metadata
  dt$.__source_file__ <- filename
  # Determine mode from filename for partitioning
  dt$.__mode__ <- ifelse(grepl("unobserved", filename, ignore.case = TRUE), "unobserved", "observed")

  # 6. Remove list columns (Arrow doesn't like mixed list columns)
  list_cols <- sapply(dt, is.list)
  if (any(list_cols)) {
    dt[, (names(dt)[list_cols]) := NULL]
  }

  return(dt)
}

# 3. MAIN PROCESSING LOOP (STREAMING)
# ------------------------------------------------------------------------------
process_streaming <- function() {
  cat("=== STARTING PROCESSING ===\n")
  
  # A. Discovery
  root_dir <- here("data/hindcasts/driver_uncertainty")
  cat("Scanning", root_dir, "...\n")
  
  # Find all relevant RDS files
  all_files <- list.files(root_dir, pattern = "_(observed|unobserved)\\.rds$", full.names = TRUE, recursive = TRUE)
  
  if (length(all_files) == 0) {
    cat("No files found. Exiting.\n")
    return(invisible(NULL))
  }
  
  # B. Resume Logic (with mtime tracking)
  processed_log <- get_processed_log()

  # Get current modification times for all files
  file_info <- data.table(
    full_path = all_files,
    filename = basename(all_files),
    current_mtime = file.mtime(all_files)
  )

  # A file needs processing if it's new OR its mtime is newer than what we logged
  file_info <- merge(file_info, processed_log, by = "filename", all.x = TRUE)
  files_to_do_idx <- is.na(file_info$mtime) | file_info$current_mtime > file_info$mtime
  files_to_do <- file_info$full_path[files_to_do_idx]
  n_new <- sum(is.na(file_info$mtime[files_to_do_idx]))
  n_updated <- sum(!is.na(file_info$mtime[files_to_do_idx]))

  cat(sprintf("Total files: %d | Already processed: %d | To process: %d (new: %d, updated: %d)\n",
              length(all_files), nrow(processed_log), length(files_to_do), n_new, n_updated))

  # If any previously-processed files need re-processing, we must do a full rebuild
  # because parquet partitions are batched and can't be selectively updated
  if (n_updated > 0 && nrow(processed_log) > 0) {
    cat("Detected updated files. Clearing staging parquet files and log for full rebuild.\n")
    existing_parquets <- list.files(DATASET_DIR, pattern = "\\.parquet$", full.names = TRUE)
    if (length(existing_parquets) > 0) {
      file.remove(existing_parquets)
      cat(sprintf("  Removed %d stale parquet partition files.\n", length(existing_parquets)))
    }
    # Clear the log so all files get re-processed
    saveRDS(data.table(filename = character(0), mtime = as.POSIXct(character(0))), LOG_FILE)
    processed_log <- get_processed_log()
    # Now all files need processing
    files_to_do <- all_files
    file_info <- data.table(
      full_path = all_files,
      filename = basename(all_files),
      current_mtime = file.mtime(all_files)
    )
    cat(sprintf("  Full rebuild: %d files to process.\n", length(files_to_do)))
  }
  
  if (length(files_to_do) == 0) {
    cat("All files up to date. Skipping processing loop.\n")
    return(invisible(NULL))
  }

  # C. Batch Processing
  # Split files into chunks of BATCH_SIZE
  batches <- split(files_to_do, ceiling(seq_along(files_to_do) / BATCH_SIZE))
  
  for (i in seq_along(batches)) {
    batch_files <- batches[[i]]
    cat(sprintf("  Processing batch %d/%d (%d files)...\n", i, length(batches), length(batch_files)))
    
    batch_data <- vector("list", length(batch_files))
    valid_indices <- integer(0)
    
    # Read and clean files in this batch
    for (j in seq_along(batch_files)) {
      f <- batch_files[j]
      tryCatch({
        raw <- readRDS(f)
        if (!is.null(raw) && nrow(raw) > 0) {
          dt <- as.data.table(raw)
          dt <- clean_hindcast_dt(dt, basename(f))
          batch_data[[j]] <- dt
          valid_indices <- c(valid_indices, j)
        }
      }, error = function(e) {
        cat(sprintf("    Error reading %s: %s\n", basename(f), e$message))
      })
    }
    
    # If we have data, write to Parquet Dataset
    if (length(valid_indices) > 0) {
      combined_batch <- rbindlist(batch_data[valid_indices], use.names = TRUE, fill = TRUE)
      
      # Write partition
      # We use a timestamp to avoid overwriting if run multiple times quickly
      part_name <- sprintf("part-%s-%04d.parquet", format(Sys.time(), "%Y%m%d%H%M%S"), i)
      nanoparquet::write_parquet(combined_batch, file.path(DATASET_DIR, part_name))
      
      # Update Log (Save success with mtimes)
      batch_basenames <- basename(batch_files)
      batch_mtimes <- file_info$current_mtime[match(batch_basenames, file_info$filename)]
      update_processed_log(batch_basenames, batch_mtimes)

      rm(combined_batch)
    } else {
      # Even if empty/error, mark as processed so we don't loop forever on bad files
      batch_basenames <- basename(batch_files)
      batch_mtimes <- file_info$current_mtime[match(batch_basenames, file_info$filename)]
      update_processed_log(batch_basenames, batch_mtimes)
    }
    
    # Garbage Collection
    rm(batch_data)
    gc(verbose = FALSE)
  }
}

# Run the streaming processor
process_streaming()


# 4. FINAL OUTPUT GENERATION (COMBINE & SAVE)
# ------------------------------------------------------------------------------
# Local aliases removed — use fill_pretty_group() from microbialForecast package.
# fill_pretty_group() derives kingdom from rank_name suffix (_bac/_fun) and
# assign_fg_kingdoms() for functional groups. It is the single canonical source.

cat("\n=== GENERATING FINAL OUTPUTS ===\n")

# Read all parquet partitions ONCE, then split by mode
parquet_files <- list.files(DATASET_DIR, pattern = "\\.parquet$", full.names = TRUE)

if (length(parquet_files) == 0) {
  cat("No parquet files found in dataset directory. Nothing to combine.\n")
} else {
  cat(sprintf("Reading %d parquet partition files (single pass)...\n", length(parquet_files)))

  # Read all partitions in batches and combine into one data.table
  BATCH_COMBINE <- 50L
  all_data <- vector("list", ceiling(length(parquet_files) / BATCH_COMBINE))
  batch_idx <- 0L

  for (i in seq(1, length(parquet_files), by = BATCH_COMBINE)) {
    j <- min(i + BATCH_COMBINE - 1L, length(parquet_files))
    batch_chunks <- lapply(parquet_files[i:j], function(f) {
      as.data.table(nanoparquet::read_parquet(f))
    })
    batch_idx <- batch_idx + 1L
    all_data[[batch_idx]] <- rbindlist(batch_chunks, use.names = TRUE, fill = TRUE)
    rm(batch_chunks); gc(verbose = FALSE)
    if (batch_idx %% 4 == 0) cat(sprintf("  Read %d/%d files...\n", j, length(parquet_files)))
  }

  df_full <- rbindlist(all_data[1:batch_idx], use.names = TRUE, fill = TRUE)
  rm(all_data); gc(verbose = FALSE)
  df_full <- fill_pretty_group(df_full)
  cat(sprintf("  Combined: %d rows\n", nrow(df_full)))

  # Save combined parquet
  cat("Saving combined parquet...\n")
  tryCatch({
    nanoparquet::write_parquet(df_full, OUT_PARQUET_FULL)
    cat(sprintf("  ✓ Saved %s (%d rows)\n", basename(OUT_PARQUET_FULL), nrow(df_full)))
  }, error = function(e) cat("  Failed:", e$message, "\n"))

  # Split by mode and save observed/unobserved RDS
  mode_col <- ".__mode__"
  if (mode_col %in% names(df_full)) {
    cat("Saving Observed RDS...\n")
    tryCatch({
      df_obs <- df_full[get(mode_col) == "observed"]
      if (nrow(df_obs) > 0) {
        saveRDS(df_obs, OUT_RDS_OBSERVED)
        cat(sprintf("  ✓ Saved %s (%d rows)\n", basename(OUT_RDS_OBSERVED), nrow(df_obs)))
      } else {
        cat("  Warning: No observed data found.\n")
      }
      rm(df_obs); gc(verbose = FALSE)
    }, error = function(e) cat("  Failed:", e$message, "\n"))

    cat("Saving Unobserved RDS...\n")
    tryCatch({
      df_unobs <- df_full[get(mode_col) == "unobserved"]
      if (nrow(df_unobs) > 0) {
        saveRDS(df_unobs, OUT_RDS_UNOBSERVED)
        cat(sprintf("  ✓ Saved %s (%d rows)\n", basename(OUT_RDS_UNOBSERVED), nrow(df_unobs)))
      } else {
        cat("  Warning: No unobserved data found.\n")
      }
      rm(df_unobs); gc(verbose = FALSE)
    }, error = function(e) cat("  Failed:", e$message, "\n"))
  }

  # Save combined RDS
  cat("Saving Combined RDS...\n")
  tryCatch({
    cat(sprintf("  Total Rows: %d\n", nrow(df_full)))
    if ("truth" %in% names(df_full)) {
      cat(sprintf("  Truth Coverage: %.1f%%\n", 100 * sum(!is.na(df_full$truth)) / nrow(df_full)))
    }
    saveRDS(df_full, OUT_RDS_COMBINED)
    cat(sprintf("  ✓ Saved %s\n", basename(OUT_RDS_COMBINED)))
    rm(df_full); gc(verbose = FALSE)
  }, error = function(e) {
    cat("\n⚠️  Could not save combined RDS (memory). Parquet and split files ARE saved.\n")
  })
}

cat("\n=== PROCESSING COMPLETE ===\n")

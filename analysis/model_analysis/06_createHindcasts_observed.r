#!/usr/bin/env Rscript

cat("=== HINDCASTS FOR OBSERVED SITES ===\n")

## --- Setup ---

options(bitmapType = "cairo")

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
Sys.setenv("OMP_THREAD_LIMIT" = "1")  # avoids OpenMP device contention in parallel workers

source("../../source.R")
source(here::here("microbialForecast/R/prepBetaRegData.r"))
source(here::here("microbialForecast/R/run_hindcast.r"))

# ---- HARDENED PLOTTER: override generate_hindcast_diagnostics ----------------
# Never pivots truth; plots ribbons/medians from predictions and points/lines from truth separately
generate_hindcast_diagnostics <- function(df, model_id, taxon, out_dir, per_site = TRUE) {
  stopifnot(dir.exists(out_dir) || dir.create(out_dir, recursive = TRUE, showWarnings = FALSE))
  
  # 1) Make dates sane and numeric columns truly numeric
  if ("dateID" %in% names(df) && "dates" %in% names(df) && !inherits(df$dates, "Date")) {
    df$dates <- as.Date(paste0(df$dateID, "01"), "%Y%m%d")
  }
  num_cols <- intersect(c("lo","lo_25","med","hi_75","hi","mean","sd","truth","truth_obs"), names(df))
  for (nm in num_cols) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  
  # 2) Create a single truth column we'll actually plot
  # CRITICAL FIX: Check if truth columns exist before using them
  if (!"truth_obs" %in% names(df)) df$truth_obs <- NA_real_
  if (!"truth" %in% names(df)) df$truth <- NA_real_
  df$truth_plot <- dplyr::coalesce(df$truth_obs, df$truth)
  
  # 3) Basic invariants before we touch anything
  pre_truth_n <- sum(!is.na(df$truth_plot))
  if (pre_truth_n == 0) {
    warning("No non-NA truth to plot for ", taxon, " (model ", model_id, ")")
  }
  
  # 4) Split prediction vs truth frames (truth never gets pivoted)
  pred_df <- df %>%
    dplyr::filter(!is.na(med) | !is.na(lo) | !is.na(hi)) %>%
    dplyr::select(dplyr::any_of(c(
      "siteID","plotID","dates","dateID","fcast_period","model_id","species","taxon","rank_name",
      "site_prediction","newsite","predicted_site_effect","lo","lo_25","med","hi_75","hi"
    )))
  
  pred_df <- dedup_key(pred_df, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
  
  truth_df <- df %>%
    dplyr::filter(!is.na(truth_plot)) %>%
    dplyr::select(dplyr::any_of(c(
      "siteID","plotID","dates","dateID","truth_plot","fcast_period","model_id","species","taxon"
    )))
  
  truth_df <- dedup_key(truth_df, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
  
  # 5) Plotting helper (one PNG per site or one big)
  period_colours <- c("calibration" = "grey50", "hindcast" = "steelblue")
  period_fills   <- c("calibration" = "grey70", "hindcast" = "steelblue")

  make_plot <- function(p_pred, p_truth, title) {
    g <- ggplot2::ggplot()

    # Ensure fcast_period exists and is a factor with consistent levels
    if (nrow(p_pred)) {
      if (!"fcast_period" %in% names(p_pred)) p_pred$fcast_period <- "hindcast"
      p_pred$fcast_period <- factor(p_pred$fcast_period, levels = c("calibration", "hindcast"))

      # Deduplicate by dates (x-axis) to prevent visual duplicates
      if ("dates" %in% names(p_pred) && nrow(p_pred) > 0) {
        p_pred <- p_pred %>%
          dplyr::arrange(dates) %>%
          dplyr::distinct(plotID, dates, .keep_all = TRUE)
        if ("dateID" %in% names(p_pred)) {
          p_pred <- p_pred %>%
            dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
            dplyr::arrange(dates)
        }
      }

      g <- g +
        ggplot2::geom_ribbon(
          data = p_pred,
          ggplot2::aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period, group = plotID),
          alpha = 0.20
        ) +
        ggplot2::geom_line(
          data = p_pred,
          ggplot2::aes(x = dates, y = med, colour = fcast_period, group = plotID),
          linewidth = 0.35
        ) +
        ggplot2::scale_colour_manual(values = period_colours, drop = FALSE) +
        ggplot2::scale_fill_manual(values = period_fills, drop = FALSE)

      # Add vertical line at calibration/hindcast boundary
      if ("fcast_period" %in% names(p_pred) && all(c("calibration", "hindcast") %in% p_pred$fcast_period)) {
        cal_dates <- p_pred$dates[p_pred$fcast_period == "calibration"]
        if (length(cal_dates) > 0) {
          boundary_date <- max(cal_dates, na.rm = TRUE)
          g <- g + ggplot2::geom_vline(xintercept = boundary_date, linetype = "dashed",
                                       colour = "grey30", linewidth = 0.4)
        }
      }
    }

    if (nrow(p_truth)) {
      p_truth <- p_truth %>% dplyr::arrange(dates)
      if ("dates" %in% names(p_truth) && nrow(p_truth) > 0) {
        p_truth <- dedup_key(p_truth, c("plotID", "dates"), arrange_by = "dates")
        if ("dateID" %in% names(p_truth) && anyDuplicated(p_truth[, c("plotID", "dates")])) {
          p_truth <- dedup_key(p_truth, c("plotID", "dateID"), arrange_by = "dates")
        }
      }

      g <- g +
        ggplot2::geom_point(
          data = p_truth,
          ggplot2::aes(x = dates, y = truth_plot, group = plotID),
          size = 0.9, stroke = 0
        )
    }

    # Facet by plot within a site; or all plots if site is missing
    facet_var <- if ("plotID" %in% names(p_pred)) "plotID" else if ("plotID" %in% names(p_truth)) "plotID" else NULL
    if (!is.null(facet_var)) g <- g + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_x")

    g +
      ggplot2::labs(title = title, x = NULL, y = "Relative abundance (0-1)") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     legend.position = "bottom", legend.title = ggplot2::element_blank())
  }
  
  # 6) Write figures - per plotID to match original format
  if ("plotID" %in% names(df)) {
    for (pid in sort(unique(df$plotID))) {
      plot_pred  <- pred_df  %>% dplyr::filter(plotID == pid)
      plot_truth <- truth_df %>% dplyr::filter(plotID == pid)
      
      # CRITICAL FIX: Deduplicate by dateID to prevent visual duplicates in plots
      if ("dateID" %in% names(plot_pred) && nrow(plot_pred) > 0) {
        plot_pred <- plot_pred %>%
          dplyr::arrange(dateID) %>%
          dplyr::distinct(plotID, dateID, .keep_all = TRUE)
      }
      if ("dateID" %in% names(plot_truth) && nrow(plot_truth) > 0) {
        plot_truth <- plot_truth %>%
          dplyr::arrange(dateID) %>%
          dplyr::distinct(plotID, dateID, .keep_all = TRUE)
      }
      
      # Get siteID for this plot
      site_id <- if ("siteID" %in% names(df)) unique(df$siteID[df$plotID == pid])[1] else "unknown"
      
      title <- paste0("Observed sites – ", taxon, " (", model_id, ") | site=", site_id)
      p <- make_plot(plot_pred, plot_truth, title)
      
      # Use original filename format: hindcast_<plotID>_<taxon>.png
      out_png <- file.path(out_dir, paste0("hindcast_", pid, "_", gsub("[^A-Za-z0-9]+","_", taxon), ".png"))
      tryCatch({
        ggplot2::ggsave(out_png, p, width = 10, height = 6, dpi = 150, limitsize = FALSE)
      }, error = function(e) {
        warning("Failed to save figure ", out_png, ": ", conditionMessage(e))
      })
    }
  } else if (per_site && "siteID" %in% names(df)) {
    for (sid in sort(unique(df$siteID))) {
      sid_pred  <- pred_df  %>% dplyr::filter(siteID == sid)
      sid_truth <- truth_df %>% dplyr::filter(siteID == sid)
      
      sid_pred <- dedup_key(sid_pred, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
      sid_truth <- dedup_key(sid_truth, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
      
      title <- paste0("Observed sites – ", taxon, " (", model_id, ") | site=", sid)
      p <- make_plot(sid_pred, sid_truth, title)
      
      out_png <- file.path(out_dir, paste0("hindcast_", gsub("[^A-Za-z0-9]+","_", taxon),
                                           "_", gsub("[^A-Za-z0-9]+","_", model_id),
                                           "_site_", sid, ".png"))
      tryCatch({
        ggplot2::ggsave(out_png, p, width = 10, height = 6, dpi = 150, limitsize = FALSE)
      }, error = function(e) {
        warning("Failed to save figure ", out_png, ": ", conditionMessage(e))
      })
    }
  } else {
    pred_df <- dedup_key(pred_df, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
    truth_df <- dedup_key(truth_df, c("plotID", "dateID"), arrange_by = c("plotID", "dateID"))
    
    title <- paste0("Observed sites – ", taxon, " (", model_id, ")")
    p <- make_plot(pred_df, truth_df, title)
    out_png <- file.path(out_dir, paste0("hindcast_", gsub("[^A-Za-z0-9]+","_", taxon),
                                         "_", gsub("[^A-Za-z0-9]+","_", model_id), ".png"))
    tryCatch({
      ggplot2::ggsave(out_png, p, width = 12, height = 7, dpi = 150, limitsize = FALSE)
    }, error = function(e) {
      warning("Failed to save figure ", out_png, ": ", conditionMessage(e))
    })
  }
  
  # 7) Post-condition: Truth must not be lost by anything we did
  post_truth_n <- sum(!is.na(df$truth_plot))
  if (post_truth_n < pre_truth_n) {
    stop("Plotting prep lost truth rows: ", pre_truth_n, " → ", post_truth_n)
  }
  invisible(TRUE)
}
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(stringr); library(lubridate); library(tidyr)
})

# Force all common threaders to single-thread (prevent runaway threading)
# 1) data.table
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)

# 2) RcppParallel (TBB users)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  try(RcppParallel::setThreadOptions(numThreads = 1), silent=TRUE)
}

# 3) RhpcBLASctl (if available) – belt & suspenders for BLAS/OpenMP
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try({
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }, silent=TRUE)
}

project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
setwd(project_root)
cat("PROJECT ROOT:", project_root, "\n")

## --- CLI Argument Parsing ---

args <- commandArgs(trailingOnly=TRUE)
get_flag <- function(key, default=NULL) {
  hit <- grep(paste0("^--", key, "="), args, value=TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

pattern_model  <- get_flag("model",  "")
pattern_taxon  <- get_flag("taxon",  "")
pattern_mname  <- get_flag("mname",  "")  # env_cycl / env_cov / cycl_only
models_file    <- get_flag("models", "")   # optional file with one model_id per line (e.g. env_cycl_missing_observed_32.txt)
LOCAL_TEST     <- identical(tolower(get_flag("local", Sys.getenv("LOCAL_TEST", "false"))), "true")
make_figs      <- identical(tolower(get_flag("figs", "false")), "true")
force_rerun    <- identical(tolower(get_flag("force", "false")), "true")  # Force re-run even if files exist
sites_arg      <- get_flag("sites", "")  # optional comma-separated list, e.g. CPER,HARV,WOOD,BART
figs_only      <- identical(tolower(get_flag("figs-only", "false")), "true")  # generate figures from existing hindcast RDS only (main process)
sequential     <- identical(tolower(get_flag("sequential", "false")), "true")  # force sequential processing for debugging

## --- Constants ---

# CRITICAL: Only use cloglog_beta_driver_uncertainty models (exclude logit, fixed_priors, fixed_drivers, and regression)
# Override with MODEL_DIRS env var to read from a different output directory (e.g., rerun)
env_model_dirs <- Sys.getenv("MODEL_DIRS", "")
if (nzchar(env_model_dirs)) {
  MODEL_DIRS <- strsplit(env_model_dirs, ",")[[1]]
  cat("MODEL_DIRS override:", paste(MODEL_DIRS, collapse=", "), "\n")
} else {
  MODEL_DIRS <- c("cloglog_beta_driver_uncertainty")
}

TRAIN_MIN <- "20130601"
ACCEPT_MAX <- "20180101"  # Only 2018 is the accepted max

# observed sites from rank tables (exclude unobserved)
UNOBSERVED <- c("ABBY","BARR","BONA","DCFS","DEJU","HEAL","KONA","LAJA","LENO","MLBS","RMNP","SOAP","TOOL","WREF","YELL")

## --- Reproducibility ---

set.seed(42L)

## --- Helper Functions ---

as_char_cols <- function(df, cols) {
  for (nm in cols) if (nm %in% names(df)) df[[nm]] <- as.character(df[[nm]])
  df
}

# Deduplicate by keys with optional prioritization
dedup_key <- function(df, keys, prefer = NULL, arrange_by = NULL) {
  keys <- keys[keys %in% names(df)]
  if (!length(keys) || !nrow(df)) return(df)
  
  if (!is.null(prefer) && prefer %in% names(df)) {
    df <- df %>%
      dplyr::mutate(.prefer = !is.na(.data[[prefer]]) & is.finite(.data[[prefer]])) %>%
      dplyr::arrange(dplyr::desc(.prefer))
  }
  if (!is.null(arrange_by) && all(arrange_by %in% names(df))) {
    df <- df %>% dplyr::arrange(dplyr::across(dplyr::all_of(arrange_by)))
  }
  df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(keys)), .keep_all = TRUE) %>%
    dplyr::select(-dplyr::any_of(".prefer"))
}

# Ensure dateID and dates columns are properly formatted
ensure_dates <- function(df) {
  if ("dateID" %in% names(df)) {
    df$dateID <- suppressWarnings(as.numeric(as.character(df$dateID)))
  }
  if ("dates" %in% names(df) && !inherits(df$dates, "Date")) {
    df$dates <- tryCatch(
      as.Date(df$dates),
      error = function(e) {
        if ("dateID" %in% names(df)) {
          as.Date(paste0(df$dateID, "01"), "%Y%m%d")
        } else {
          as.Date(NA)
        }
      }
    )
  }
  if (!"dates" %in% names(df) && "dateID" %in% names(df)) {
    df$dates <- as.Date(paste0(df$dateID, "01"), "%Y%m%d")
  }
  df
}

# Ensure columns exist with default values
ensure_cols <- function(df, cols, default = NA_real_) {
  for (nm in cols) {
    if (!nm %in% names(df)) df[[nm]] <- default
  }
  df
}

# Debug printing helper
dbg <- function(..., .on = FALSE) {
  if (.on) {
    cat(...)
    cat("\n")
    try(flush.console(), silent = TRUE)
  }
}

# Safe readRDS with retry logic
safe_readRDS <- function(path, retry = 1, sleep = 0.5) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      if (retry > 0 && grepl("connection|cannot open", e$message, ignore.case = TRUE)) {
        Sys.sleep(sleep)
        safe_readRDS(path, retry = retry - 1, sleep = sleep)
      } else {
        stop(e)
      }
    }
  )
}

first_match_col <- function(df, patterns) {
  for (p in patterns) {
    hits <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hits)) return(hits[1])
  }
  NA_character_
}

# detect driver flag from model_id or dir
is_driver_model <- function(samples_path, model_id) {
  grepl("driver_uncertainty", samples_path) || grepl("driver_uncertainty", model_id)
}

# checks for existing per-site outputs for a model (tolerant of filename variants)
has_all_site_outputs <- function(model_id, samples_path, required_sites, project_root) {
  dir <- file.path(project_root, "data", "hindcasts",
                   if (is_driver_model(samples_path, model_id)) "driver_uncertainty" else "standard")
  if (!dir.exists(dir)) return(FALSE)
  # accept any file that matches model_id + site + observed.rds, regardless of extra middle bits
  # e.g., hindcasts_<model_id>(_anything_optional)_<SITE>_observed.rds
  # escape special regex characters in model_id (fixed, valid, and simple)
  model_id_escaped <- gsub("([.^$*+?()|[\\\\]])", "\\\\\\1", model_id)
  model_id_escaped <- gsub("(\\{|\\})", "\\\\\\1", model_id_escaped)
  site_done <- vapply(required_sites, function(sid) {
    pat <- paste0("^hindcasts_", model_id_escaped, "(_.*)?_", sid, "_observed\\.rds$")
    any(grepl(pat, list.files(dir)))
  }, logical(1))
  all(site_done)
}

# lightweight run log
runlog_path <- file.path(project_root, "data", "hindcasts", "run_log.csv")
append_log <- function(model_id, taxon, status, msg="", nrows=NA_integer_) {
  line <- data.frame(ts=as.character(Sys.time()), model_id, taxon, status, nrows, msg, stringsAsFactors=FALSE)
  if (!file.exists(runlog_path)) {
    write.table(line, runlog_path, sep=",", row.names=FALSE, col.names=TRUE, append=FALSE)
  } else {
    write.table(line, runlog_path, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}

## --- Parse Model Files ---

all_summary <- character(0)
for (sub in MODEL_DIRS) {
  base <- file.path(project_root, "data/model_outputs", sub)
  if (!dir.exists(base)) next
  hits <- list.files(base, pattern="^summary_.*_beta_regression\\.rds$", recursive=TRUE, full.names=TRUE)
  all_summary <- c(all_summary, hits)
}
if (length(all_summary) == 0) stop("No summary_*.rds files found under data/model_outputs/*")

parse_model <- function(p) {
  id <- basename(p)
  id <- sub("^summary_", "", id)
  id <- sub("_beta_regression\\.rds$", "", id)
  
  # CRITICAL FIX: Handle double-prefix models (e.g., cycl_only_cycl_only_...)
  # First, identify the actual model_name by checking for known patterns
  if (grepl("^cycl_only_cycl_only_", id)) {
    model_name <- "cycl_only"
    taxon <- sub("^cycl_only_cycl_only_", "", id)
  } else if (grepl("^env_cycl_env_cycl_", id)) {
    model_name <- "env_cycl"
    taxon <- sub("^env_cycl_env_cycl_", "", id)
  } else if (grepl("^env_cov_env_cov_", id)) {
    model_name <- "env_cov"
    taxon <- sub("^env_cov_env_cov_", "", id)
  } else {
    # Standard parsing for normal models
    model_name <- sub("^([^_]+_[^_]+)_.*$", "\\1", id)
    taxon <- sub("^[^_]+_[^_]+_", "", id)
  }
  
  taxon <- sub("_[0-9]{8}_.*$", "", taxon)
  min_date <- sub("^.*_([0-9]{8})_[0-9]{8}_.*$", "\\1", id)
  max_date <- sub("^.*_[0-9]{8}_([0-9]{8})_.*$", "\\1", id)
  
  # Get file modification time for deduplication
  file_info <- file.info(p)
  mtime <- file_info$mtime
  
  data.frame(path=p, model_id=id, model_name, taxon, min_date, max_date, mtime=mtime, stringsAsFactors=FALSE)
}

mi <- do.call(rbind, lapply(all_summary, parse_model))
mi <- mi[mi$min_date %in% TRAIN_MIN & mi$max_date %in% ACCEPT_MAX & mi$model_name %in% c("env_cycl","env_cov","cycl_only"), , drop=FALSE]
if (nrow(mi) == 0) stop("No candidate models (after filter).")

# CRITICAL: Exclude reconstructed_from_checkpoints models
n_before <- nrow(mi)
mi <- mi[!grepl("reconstructed_from_checkpoints", mi$path, ignore.case=TRUE) & 
         !grepl("reconstructed_from_checkpoints", mi$model_id, ignore.case=TRUE), , drop=FALSE]
if (n_before > nrow(mi)) {
  cat("Excluded", n_before - nrow(mi), "reconstructed_from_checkpoints models\n")
}

# CRITICAL: Ensure all models are from driver_uncertainty directories only
n_before <- nrow(mi)
mi <- mi[grepl("driver_uncertainty", mi$path, ignore.case=TRUE), , drop=FALSE]
if (n_before > nrow(mi)) {
  cat("Excluded", n_before - nrow(mi), "non-driver_uncertainty models\n")
}

if (nrow(mi) == 0) stop("No driver_uncertainty models found (after filters).")

# CRITICAL: Deduplicate by model_id (keep newest file based on modification time)
# There may be multiple summary files for the same model_id (e.g., from different chains or runs)
# Newer files are typically in parent directories, older ones in taxon-specific subdirectories
n_before_dedup <- nrow(mi)
if (n_before_dedup > 0) {
  # Sort by model_id and mtime (newest first), then keep first occurrence of each model_id
  mi <- mi[order(mi$model_id, -as.numeric(mi$mtime)), ]
  mi <- mi[!duplicated(mi$model_id), , drop=FALSE]
  # Remove mtime column as it's no longer needed
  mi$mtime <- NULL
}
if (n_before_dedup > nrow(mi)) {
  cat("Deduplicated", n_before_dedup - nrow(mi), "duplicate model_ids (keeping newest file by modification time)\n")
  cat("Final unique models:", nrow(mi), "\n")
}

# CLI filters
if (nzchar(pattern_model)) {
  mi <- mi[grepl(pattern_model, mi$model_id, ignore.case=TRUE), , drop=FALSE]
  cat("Filtered to models matching pattern:", pattern_model, "\n")
}
if (nzchar(pattern_taxon)) {
  mi <- mi[grepl(pattern_taxon, mi$taxon, ignore.case=TRUE), , drop=FALSE]
  cat("Filtered to taxa matching pattern:", pattern_taxon, "\n")
}
if (nzchar(pattern_mname)) {
  mi <- mi[mi$model_name == pattern_mname, , drop=FALSE]
  cat("Filtered to model name:", pattern_mname, "\n")
}
if (nzchar(models_file)) {
  mpath <- if (file.exists(models_file)) models_file else file.path(project_root, "data/hindcasts", models_file)
  if (!file.exists(mpath)) stop("Models file not found: ", models_file, " or ", mpath)
  want <- trimws(readLines(mpath))
  want <- want[nzchar(want) & !grepl("^#", want)]
  n_before <- nrow(mi)
  mi <- mi[mi$model_id %in% want, , drop=FALSE]
  cat("Filtered to", nrow(mi), "model_ids from", mpath, "(requested", length(want), ")\n")
  if (nrow(mi) == 0) stop("No models in mi match the model_ids in ", mpath)
}

if (nrow(mi) == 0) stop("No models remaining after filters.")

# Limit models in LOCAL_TEST mode (final limit applied later at line ~2675)
if (LOCAL_TEST && nrow(mi) > 5) {
  cat("LOCAL_TEST mode: Limiting to first 5 models (out of", nrow(mi), "total)\n")
  mi <- mi[1:5, , drop=FALSE]
}

# Get required sites (pre-load rank data once)
r16 <- readRDS(file.path(project_root, "data/clean/groupAbundances_16S_2023.rds"))
rits <- readRDS(file.path(project_root, "data/clean/groupAbundances_ITS_2023.rds"))
all_sites <- unique(unlist(lapply(c(r16, rits), function(x) unique(x$siteID))))
required_sites <- sort(setdiff(all_sites, c(UNOBSERVED, NA)))
if (nzchar(sites_arg)) {
  sites_filter <- trimws(strsplit(sites_arg, ",")[[1]])
  sites_filter <- sites_filter[nzchar(sites_filter)]
  if (length(sites_filter) > 0) {
    required_sites <- intersect(required_sites, sites_filter)
    cat("Restricted to sites:", paste(required_sites, collapse = ", "), "\n")
  }
}
cat("Required sites:", length(required_sites), "\n")

# Load predictor data once (avoid re-reading per model)
env_data <- readRDS(file.path(project_root, "data/clean/all_predictor_data.rds"))
cat("Loaded predictor data\n")

## --- Main Processing Function ---

run_one_model <- function(row, project_root, required_sites, env_data=NULL, r16=NULL, rits=NULL, LOCAL_TEST=FALSE, force_rerun=FALSE, make_figs=FALSE) {
  # In workers, r16/rits/env_data should already be loaded via clusterEvalQ
  # In sequential mode, we use the passed-in versions
  # Fallback: load if not provided (shouldn't happen in parallel)
  if (is.null(env_data)) {
    if (exists("env_data", envir=.GlobalEnv)) {
      env_data <- get("env_data", envir=.GlobalEnv)
    } else {
      env_data <- readRDS(file.path(project_root, "data/clean/all_predictor_data.rds"))
    }
  }
  if (is.null(r16)) {
    if (exists("r16", envir=.GlobalEnv)) {
      r16 <- get("r16", envir=.GlobalEnv)
    } else {
      r16 <- readRDS(file.path(project_root, "data/clean/groupAbundances_16S_2023.rds"))
    }
  }
  if (is.null(rits)) {
    if (exists("rits", envir=.GlobalEnv)) {
      rits <- get("rits", envir=.GlobalEnv)
    } else {
      rits <- readRDS(file.path(project_root, "data/clean/groupAbundances_ITS_2023.rds"))
    }
  }
  # Ensure single values to avoid vector condition errors
  model_id <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
  taxon <- if (length(row$taxon) == 1) row$taxon else row$taxon[1]
  
  cat("\n=== Processing:", model_id, "taxon:", taxon, "model:", row$model_name, "===\n")
  
    tryCatch({
      # Locate paired samples_ file in same tree
      model_root <- dirname(dirname(row$path))
      # Escape regex metacharacters in model_id for safe pattern matching
      model_id_escaped <- gsub("([.^$*+?()|[\\\\]])", "\\\\\\1", model_id)
      model_id_escaped <- gsub("(\\{|\\})", "\\\\\\1", model_id_escaped)
      samples_guess <- list.files(model_root, pattern=paste0("^samples_", model_id_escaped, ".*\\.rds$"), recursive=TRUE, full.names=TRUE)
      samples_guess <- samples_guess[!grepl("_chain[0-9]", samples_guess)]
    if (length(samples_guess) < 1) {
      msg <- "No samples file found"
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    samples_path <- samples_guess[[1]]
    
    # Detect link function from directory name for metadata
    link_function <- if (grepl("cloglog", basename(model_root), ignore.case=TRUE)) {
      "cloglog"
    } else if (grepl("logit", basename(model_root), ignore.case=TRUE)) {
      "logit"
    } else {
      "cloglog"  # default fallback
    }
    
    # Check if already complete (unless force rerun)
    model_id_for_check <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
    if (!force_rerun && has_all_site_outputs(model_id_for_check, samples_path, required_sites, project_root)) {
      cat("⏭️  SKIP (complete):", model_id, "\n")
      append_log(model_id, taxon, "SKIP", "already complete", NA_integer_)
      return(list(success=TRUE, model_id=model_id, taxon=taxon, skipped=TRUE))
    }
    if (force_rerun) {
      cat("⚠️  FORCE RERUN:", model_id, "(ignoring existing files)\n")
    }
    
    # Load model + data
    # CRITICAL FIX: Add retry logic for connection errors
    model_summary <- tryCatch({
      readRDS(row$path)
    }, error = function(e) {
      if (grepl("connection|cannot open", e$message, ignore.case=TRUE)) {
        # Retry once after a short delay
        Sys.sleep(0.5)
        tryCatch({
          readRDS(row$path)
        }, error = function(e2) {
          stop(paste("Failed to read model summary after retry:", e2$message))
        })
      } else {
        stop(e)
      }
    })
    # CRITICAL FIX: model_summary structure is list(summary_df, pred.means, pred.quantiles, gd)
    # We need pred.quantiles (element 3) which has columns like "2.5%", "50%", "97.5%"
    # NOT pred.means (element 2) which has Mean, SD, etc.
    stopifnot(is.list(model_summary), length(model_summary) >= 3L, is.data.frame(model_summary[[3]]))
    cal_quants <- model_summary[[3]]  # pred.quantiles with quantile columns
    
    # CRITICAL FIX: Check if cal_quants is empty (some models have no prediction data)
    if (nrow(cal_quants) == 0) {
      msg <- "model_summary[[3]] (pred.quantiles) is empty - model has no prediction data"
      append_log(model_id, taxon, "ERROR", msg, 0L)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    model_dat <- safe_readRDS(samples_path)
    param_samples <- model_dat$samples
    
    # Find rank table containing this taxon
    all_ranks <- c(r16, rits)
    rank_name <- NULL
    taxon_for_search <- if (length(row$taxon) == 1) row$taxon else row$taxon[1]
    for (nm in names(all_ranks)) {
      condition_result <- taxon_for_search %in% colnames(all_ranks[[nm]])
      if (length(condition_result) == 1 && condition_result) { 
        rank_name <- nm
        break 
      }
    }
    if (is.null(rank_name)) {
      msg <- paste("Could not find rank table containing taxon:", row$taxon)
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    rank.df <- all_ranks[[rank_name]]
    rank.df <- as_char_cols(rank.df, c("siteID","plotID"))
    taxon_for_select <- if (length(row$taxon) == 1) row$taxon else row$taxon[1]
    rank.df_spec <- rank.df %>% select(siteID, plotID, dateID, sampleID, dates, plot_date, all_of(taxon_for_select))
    rank.df_spec$other <- 1 - rank.df_spec[[taxon_for_select]]
    
    # Apply same horizon filtering as during model fitting to ensure truth data matches training data
    dom_soil_horizons <- readRDS(file.path(project_root, "data/clean/dominantHorizonsSite.rds"))
    keep_hor <- paste0(dom_soil_horizons$siteID, dom_soil_horizons$horizon)
    rank.df_spec <- rank.df_spec %>%
      mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O"),
             site_hor = paste0(siteID, horizon)) %>%
      dplyr::filter(site_hor %in% keep_hor) %>%
      select(-c(horizon, site_hor))
    
    # Time mapping from cal_quants
    required_cols <- c("timepoint")
    missing_cols <- setdiff(required_cols, names(cal_quants))
    if (length(missing_cols) > 0) {
      msg <- paste("model_summary[[3]] (pred.quantiles) is missing required columns:", paste(missing_cols, collapse=", "))
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    # Reconstruct plotID if missing (from model data)
    if (!"plotID" %in% names(cal_quants)) {
      # Try to reconstruct plotID from plot_num using model data
      if ("plot_num" %in% names(cal_quants)) {
        # Get plotID mapping from model data (truth.plot.long should always have both plot_num and plotID)
        tryCatch({
          if (!is.null(model_dat$metadata) && !is.null(model_dat$metadata$model_data)) {
            truth_data <- if ("truth.plot.long" %in% names(model_dat$metadata$model_data)) {
              model_dat$metadata$model_data$truth.plot.long
            } else {
              model_dat$metadata$model_data
            }
            if (is.data.frame(truth_data) && "plot_num" %in% names(truth_data) && "plotID" %in% names(truth_data)) {
              plot_mapping <- truth_data %>% 
                select(plot_num, plotID) %>% 
                distinct() %>%
                filter(!is.na(plot_num) & !is.na(plotID))
              if (nrow(plot_mapping) > 0) {
                plot_key <- plot_mapping$plotID
                names(plot_key) <- as.character(plot_mapping$plot_num)
                cal_quants$plotID <- plot_key[as.character(cal_quants$plot_num)]
                cat("  ✓ Reconstructed plotID from model data\n")
              } else {
                msg <- "model_summary[[3]] (pred.quantiles) is missing plotID column - no plot mapping found in model data"
                append_log(model_id, taxon, "ERROR", msg, NA_integer_)
                return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
              }
            } else {
              msg <- "model_summary[[3]] (pred.quantiles) is missing plotID column - model data structure invalid"
              append_log(model_id, taxon, "ERROR", msg, NA_integer_)
              return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
            }
          } else {
            msg <- "model_summary[[3]] (pred.quantiles) is missing plotID column - no model metadata available"
            append_log(model_id, taxon, "ERROR", msg, NA_integer_)
            return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
          }
        }, error = function(e) {
          msg <- paste("model_summary[[3]] (pred.quantiles) is missing plotID column - reconstruction failed:", e$message)
          append_log(model_id, taxon, "ERROR", msg, NA_integer_)
          return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
        })
      } else {
        msg <- "model_summary[[3]] (pred.quantiles) is missing plotID column and plot_num"
        append_log(model_id, taxon, "ERROR", msg, NA_integer_)
        return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
      }
    }
    
    # Reconstruct siteID if missing (from model data or plotID)
    if (!"siteID" %in% names(cal_quants)) {
      # First try to reconstruct from site_num using model data
      if ("site_num" %in% names(cal_quants)) {
        # Get siteID mapping from model data (truth.plot.long should always have both site_num and siteID)
        tryCatch({
          if (!is.null(model_dat$metadata) && !is.null(model_dat$metadata$model_data)) {
            truth_data <- if ("truth.plot.long" %in% names(model_dat$metadata$model_data)) {
              model_dat$metadata$model_data$truth.plot.long
            } else {
              model_dat$metadata$model_data
            }
            if (is.data.frame(truth_data) && "site_num" %in% names(truth_data) && "siteID" %in% names(truth_data)) {
              site_mapping <- truth_data %>% 
                select(site_num, siteID) %>% 
                distinct() %>%
                filter(!is.na(site_num) & !is.na(siteID))
              if (nrow(site_mapping) > 0) {
                site_key <- site_mapping$siteID
                names(site_key) <- as.character(site_mapping$site_num)
                cal_quants$siteID <- site_key[as.character(cal_quants$site_num)]
                cat("  ✓ Reconstructed siteID from model data\n")
              } else {
                # Fallback: extract from plotID if available
                if ("plotID" %in% names(cal_quants)) {
                  cal_quants$siteID <- substr(cal_quants$plotID, 1, 4)
                  cat("  ✓ Reconstructed siteID from plotID\n")
                } else {
                  msg <- "model_summary[[3]] (pred.quantiles) is missing siteID column - no site mapping found"
                  append_log(model_id, taxon, "ERROR", msg, NA_integer_)
                  return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
                }
              }
            } else {
              # Fallback: extract from plotID if available
              if ("plotID" %in% names(cal_quants)) {
                cal_quants$siteID <- substr(cal_quants$plotID, 1, 4)
                cat("  ✓ Reconstructed siteID from plotID\n")
              } else {
                msg <- "model_summary[[3]] (pred.quantiles) is missing siteID column - model data structure invalid"
                append_log(model_id, taxon, "ERROR", msg, NA_integer_)
                return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
              }
            }
          } else {
            # Fallback: extract from plotID if available
            if ("plotID" %in% names(cal_quants)) {
              cal_quants$siteID <- substr(cal_quants$plotID, 1, 4)
              cat("  ✓ Reconstructed siteID from plotID\n")
            } else {
              msg <- "model_summary[[3]] (pred.quantiles) is missing siteID column - no model metadata available"
              append_log(model_id, taxon, "ERROR", msg, NA_integer_)
              return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
            }
          }
        }, error = function(e) {
          # Fallback: extract from plotID if available
          if ("plotID" %in% names(cal_quants)) {
            cal_quants$siteID <<- substr(cal_quants$plotID, 1, 4)
            cat("  ✓ Reconstructed siteID from plotID (fallback)\n")
          } else {
            msg <- paste("model_summary[[3]] (pred.quantiles) is missing siteID column - reconstruction failed:", e$message)
            append_log(model_id, taxon, "ERROR", msg, NA_integer_)
            return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
          }
        })
      } else {
        # Fallback: extract from plotID if available
        if ("plotID" %in% names(cal_quants)) {
          cal_quants$siteID <- substr(cal_quants$plotID, 1, 4)
          cat("  ✓ Reconstructed siteID from plotID\n")
        } else {
          msg <- "model_summary[[3]] (pred.quantiles) is missing siteID column and site_num"
          append_log(model_id, taxon, "ERROR", msg, NA_integer_)
          return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
        }
      }
    }
    
    cal_quants <- as_char_cols(cal_quants, c("siteID","plotID"))
    
    # Create plot mapping (filter out NA values)
    required_plot_cols <- c("plotID", "siteID", "plot_num")
    missing_plot_cols <- setdiff(required_plot_cols, names(cal_quants))
    if (length(missing_plot_cols) > 0) {
      msg <- paste("model_summary[[3]] (pred.quantiles) is missing required columns for plot mapping:", paste(missing_plot_cols, collapse=", "))
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    cal_quants_for_plot_map <- cal_quants %>%
      filter(!is.na(plotID) & !is.na(siteID) & !is.na(plot_num))
    
    if (nrow(cal_quants_for_plot_map) == 0) {
      msg <- "model_summary[[3]] has no valid plot mapping data (all plotID, siteID, or plot_num are NA)"
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    # Create unique plot mapping (one row per plotID)
    plot_map <- cal_quants_for_plot_map %>% 
      select(plotID, siteID, plot_num, any_of("site_num")) %>% 
      arrange(plotID) %>%
      distinct(plotID, .keep_all = TRUE)
    
    if (nrow(plot_map) == 0) {
      msg <- "model_summary[[3]] has no plot mapping data after filtering and distinct()"
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    # Reconstruct site_num if missing
    if (!"site_num" %in% names(plot_map)) {
      tryCatch({
        if (!is.null(model_dat$metadata) && !is.null(model_dat$metadata$model_data)) {
          truth_data <- if ("truth.plot.long" %in% names(model_dat$metadata$model_data)) {
            model_dat$metadata$model_data$truth.plot.long
          } else {
            model_dat$metadata$model_data
          }
          if (is.data.frame(truth_data) && "plotID" %in% names(truth_data) && "site_num" %in% names(truth_data)) {
            site_num_mapping <- truth_data %>% 
              select(plotID, any_of("site_num")) %>% 
              distinct() %>%
              filter(!is.na(plotID))
            if (nrow(site_num_mapping) > 0) {
              plot_map <- plot_map %>% 
                left_join(site_num_mapping, by = "plotID")
              # Don't filter out NA site_num - keep all rows and handle per-plot
              cat("  ✓ Reconstructed site_num from model data\n")
            }
          }
        }
        # If still missing, try to get from full_in later (will be available after model inputs are loaded)
      }, error = function(e) {
        cat("  ⚠️  Warning: Could not reconstruct site_num:", e$message, "\n")
      })
    }
    
    # Validate timepoint and dateID exist
    if (!"timepoint" %in% names(cal_quants)) {
      msg <- "model_summary[[3]] (pred.quantiles) is missing required columns: timepoint"
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    if (!"dateID" %in% names(cal_quants)) {
      # Try to reconstruct dateID from dates if available
      if ("dates" %in% names(cal_quants)) {
        cal_quants$dateID <- as.numeric(format(as.Date(cal_quants$dates), "%Y%m"))
        cat("  ✓ Reconstructed dateID from dates\n")
      } else {
        msg <- "model_summary[[3]] (pred.quantiles) is missing required columns: dateID (and dates)"
        append_log(model_id, taxon, "ERROR", msg, NA_integer_)
        return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
      }
    }
    
    # Build time mapping: timepoint is global (consistent across plots) from prepBetaRegData
    # Each timepoint should map to one dateID; handle duplicates if present
    time_train_raw <- cal_quants %>%
      select(timepoint, dateID) %>%
      filter(!is.na(timepoint) & !is.na(dateID))
    
    # Check for duplicate timepoint->dateID mappings (shouldn't happen if timepoint is global)
    dup_check <- time_train_raw %>%
      group_by(timepoint) %>%
      summarise(n_dateids = n_distinct(dateID), .groups="drop") %>%
      filter(n_dateids > 1)
    
    if (nrow(dup_check) > 0) {
      cat("  ⚠️  DIAGNOSTIC: Found timepoints with multiple dateID values (data integrity issue):\n")
      print(head(dup_check, 10))
      cat("  ⚠️  This suggests the merge in summarizeBetaRegModels created incorrect mappings.\n")
      cat("  ⚠️  Using most common dateID for each timepoint.\n")
      fflush()
      
      # For each timepoint with multiple dateIDs, use the most common one
      time_train <- time_train_raw %>%
        group_by(timepoint, dateID) %>%
        summarise(n = n(), .groups="drop") %>%
        group_by(timepoint) %>%
        slice_max(n, n=1, with_ties=FALSE) %>%
        ungroup() %>%
        select(timepoint, dateID) %>%
        arrange(timepoint)
    } else {
      # No duplicates - use distinct
      time_train <- time_train_raw %>%
        distinct(timepoint, dateID) %>%
        arrange(timepoint)
    }
    
    time_train <- time_train %>%
      mutate(
        dateID = as.numeric(dateID),
        trained_date_num = timepoint,
        dates = {
          first_dateID <- dateID[1]
          if (length(first_dateID) == 1 && !is.na(first_dateID)) {
            as.Date(paste0(dateID, "01"), "%Y%m%d")
          } else {
            base_date <- as.Date(row$min_date, "%Y%m%d")
            # CRITICAL FIX: Use lubridate::months() explicitly or calculate manually
            # months() function may not be exported in newer lubridate versions
            tryCatch({
              base_date %m+% lubridate::months(0:(n()-1))
            }, error = function(e) {
              # Fallback: calculate months manually
              seq(from = base_date, by = "1 month", length.out = n())
            })
          }
        }
      )
    
    cal_end_dateID <- suppressWarnings(max(time_train$dateID, na.rm=TRUE))
    if (!is.finite(cal_end_dateID)) {
      base <- as.Date(row$min_date, "%Y%m%d")
      # CRITICAL FIX: Use lubridate::months() explicitly or calculate manually
      tryCatch({
        seq_dates <- base %m+% lubridate::months(0:(nrow(time_train)-1))
      }, error = function(e) {
        # Fallback: calculate months manually
        seq_dates <- seq(from = base, by = "1 month", length.out = nrow(time_train))
      })
      time_train$dateID <- as.numeric(format(seq_dates, "%Y%m"))
      time_train$dates  <- as.Date(paste0(time_train$dateID, "01"), "%Y%m%d")
      cal_end_dateID <- max(time_train$dateID, na.rm=TRUE)
    }
    max_trained_num <- max(time_train$trained_date_num, na.rm=TRUE)
    
    # Build hindcast horizon (month-by-month, robust across years)
    hindcast_end <- "20200101"
    # starting month is the month *after* cal_end_dateID
    start_mon <- as.Date(paste0(cal_end_dateID, "01"), "%Y%m%d")
    # CRITICAL FIX: Use lubridate::months() explicitly
    tryCatch({
      start_mon <- lubridate::floor_date(start_mon, unit = "month") %m+% lubridate::months(1)
    }, error = function(e) {
      # Fallback: add 1 month manually
      start_mon <- lubridate::floor_date(start_mon, unit = "month")
      start_mon <- seq(from = start_mon, by = "1 month", length.out = 2)[2]
    })
    end_mon <- as.Date(hindcast_end, "%Y%m%d")
    # monthly sequence as Date, then back to YYYYMM numeric
    mon_seq <- seq(from = start_mon, to = end_mon, by = "1 month")
    hind_dates <- data.frame(
      dateID = as.numeric(format(mon_seq, "%Y%m")),
      trained_date_num = max_trained_num + seq_along(mon_seq),
      dates = as.Date(paste0(format(mon_seq, "%Y%m"), "01"), "%Y%m%d")
    )
    
    # DIAGNOSTIC: Check for duplicate dateIDs BEFORE creating trained_time_map
    # This will help us understand the root cause
    if (anyDuplicated(time_train$dateID)) {
      dups_time_train <- time_train$dateID[duplicated(time_train$dateID)]
      cat("  ⚠️  DIAGNOSTIC: time_train has duplicate dateID values:", paste(unique(dups_time_train), collapse=", "), "\n")
      cat("  ⚠️  This means multiple timepoints map to the same dateID\n")
      # Show which timepoints map to duplicate dateIDs
      dup_mappings <- time_train %>%
        filter(dateID %in% unique(dups_time_train)) %>%
        arrange(dateID, timepoint)
      cat("  ⚠️  Duplicate mappings:\n")
      print(head(dup_mappings, 10))
      fflush()
    }
    
    # Check for overlap between time_train and hind_dates
    overlap <- intersect(time_train$dateID, hind_dates$dateID)
    if (length(overlap) > 0) {
      cat("  ⚠️  DIAGNOSTIC: hind_dates overlaps with time_train dateIDs:", paste(overlap, collapse=", "), "\n")
      cat("  ⚠️  This suggests the calibration period extends into the hindcast period\n")
      cat("  ⚠️  Filtering out overlapping dateIDs from hind_dates\n")
      fflush()
      # CRITICAL FIX: Remove overlapping dateIDs from hind_dates to prevent duplicates
      hind_dates <- hind_dates %>% filter(!dateID %in% overlap)
      if (nrow(hind_dates) == 0) {
        warning("All hindcast dateIDs were filtered out due to overlap with calibration period")
      }
    }
    
    trained_time_map <- bind_rows(
      time_train %>% select(dateID, trained_date_num, dates),
      hind_dates %>% select(dateID, trained_date_num, dates)
    ) %>% arrange(dateID) %>% mutate(trained_date_num = seq_len(nrow(.)))
    
    # Enhanced error message with diagnostics
    if (anyDuplicated(trained_time_map$dateID)) {
      dups <- trained_time_map$dateID[duplicated(trained_time_map$dateID)]
      dup_details <- trained_time_map %>%
        filter(dateID %in% unique(dups)) %>%
        arrange(dateID)
      cat("  ❌ CRITICAL: trained_time_map has duplicate dateID values:\n")
      print(dup_details)
      fflush()
      stop("CRITICAL DATA INTEGRITY ERROR: trained_time_map has duplicate dateID values: ", 
           paste(unique(dups), collapse=", "), 
           ". This indicates multiple timepoints map to the same dateID, or hind_dates overlaps with time_train.")
    }
    
    # Full inputs
    keep_vec <- c("siteID","plotID","dateID","sampleID","dates","plot_date", row$taxon)
    full_in <- prepBetaRegData(
      rank.df = rank.df_spec,
      min.prev = 0,
      min.date = TRAIN_MIN,
      max.date = as.character(hindcast_end),
      full_timeseries = TRUE,
      keep_vec = keep_vec,
      predictor_data = env_data
    )
    
    full_in$truth.plot.long <- as_char_cols(full_in$truth.plot.long, c("siteID","plotID"))
    # CRITICAL FIX: site_num may not exist, so use any_of() to make it optional
    plot_site_key <- full_in$truth.plot.long %>% select(siteID, plotID, dateID, date_num, plot_num, any_of("site_num")) %>% distinct()
    
    # Reconstruct site_num from plot_site_key if still missing
    if (!"site_num" %in% names(plot_map) && "site_num" %in% names(plot_site_key)) {
      plot_map <- plot_map %>% 
        left_join(plot_site_key %>% select(plotID, any_of("site_num")) %>% distinct(), by = "plotID")
      # Only filter out NA site_num if the column exists and we have other rows
      if ("site_num" %in% names(plot_map)) {
        # Don't filter out rows - keep them even if site_num is NA (we'll handle it per-plot)
        cat("  ✓ Reconstructed site_num from plot_site_key\n")
      }
    }
    
    # Overlap check (GUARDRAIL: stop if overlap == 0)
    trained_plots <- unique(plot_map$plotID)
    trained_plots <- as.character(trained_plots)
    in_inputs <- unique(plot_site_key$plotID)
    ovl <- intersect(trained_plots, in_inputs)
    cat("DEBUG overlap | trained:", length(trained_plots), " in_inputs:", length(in_inputs), " overlap:", length(ovl), "\n")
    fflush()
    if (length(ovl) == 0) {
      msg <- "No overlap between trained plots and input plots"
      cat("DEBUG sample trained:", paste(head(trained_plots,5), collapse=", "), "\n")
      cat("DEBUG sample inputs :", paste(head(in_inputs,5), collapse=", "), "\n")
      fflush()
      append_log(model_id, taxon, "ERROR", msg, NA_integer_)
      stop(msg)
    }
    
    # Restrict to required_sites (e.g. when --sites=CPER,HARV,... is set)
    if (length(required_sites) > 0 && "siteID" %in% names(plot_site_key)) {
      plot_sites <- plot_site_key %>% distinct(plotID, siteID)
      ovl_in_sites <- (plot_sites %>% filter(siteID %in% required_sites))$plotID
      ovl <- intersect(ovl, ovl_in_sites)
      if (length(ovl) == 0) {
        msg <- paste("No overlap plots in required sites:", paste(required_sites, collapse = ", "))
        append_log(model_id, taxon, "ERROR", msg, NA_integer_)
        return(list(success = FALSE, model_id = model_id, taxon = taxon, error = msg))
      }
    }

    # TEST MODE: Limit to first 2 plots for faster testing
    if (LOCAL_TEST && length(ovl) > 2) {
      ovl <- sort(ovl)[1:2]
      cat("TEST MODE: Limiting to", length(ovl), "plots for testing:", paste(ovl, collapse=", "), "\n")
      fflush()
    }
    
    # Extract taxon_name once at start (ensure single value to avoid vector condition errors)
    taxon_name <- as.character(row$taxon)[1]
    
    # Process each plot
    out_list <- list()
    for (plotID in ovl) {
      plot_info <- plot_map %>% filter(plotID == !!plotID) %>% distinct()
      # Ensure single row per plotID
      if (nrow(plot_info) > 1) {
        cat("  WARNING: plot_info has", nrow(plot_info), "rows for plotID", plotID, "- taking first\n")
        plot_info <- plot_info[1, , drop=FALSE]
      }
      if (nrow(plot_info) < 1) {
        cat("  WARNING: No plot_info found for", plotID, "- skipping\n")
        next
      }
      pnum <- plot_info$plot_num[1]
      snum <- if ("site_num" %in% names(plot_info)) plot_info$site_num[1] else NA_integer_
      sid <- as.character(plot_info$siteID[1])
      
      # Get site_num from plot_site_key if still missing
      if (is.na(snum) && "site_num" %in% names(plot_site_key)) {
        site_num_row <- plot_site_key %>% filter(plotID == !!plotID) %>% head(1)
        if (nrow(site_num_row) > 0) {
          site_num_val <- site_num_row$site_num[1]
          if (length(site_num_val) == 1 && !is.na(site_num_val)) {
            snum <- as.integer(site_num_val)
          }
        }
      }
      
      # Calibration truth for this plot (filter to calibration period only)
      condition_check <- (length(taxon_name) == 1L && taxon_name %in% names(rank.df_spec))
      if (!condition_check) {
        if (LOCAL_TEST) {
          cat("  WARN: Taxon", taxon_name, "not found in rank.df_spec after horizon filtering - creating empty cal_truth\n")
          fflush()
        }
        cal_truth <- data.frame(siteID = character(0), plotID = character(0), dateID = numeric(0), 
                                dates = as.Date(character(0)), truth = numeric(0), plot_num = integer(0), 
                                date_num = numeric(0), site_num = integer(0), timepoint = numeric(0),
                                stringsAsFactors = FALSE)
      } else {
        cal_truth <- rank.df_spec %>%
          filter(plotID == !!plotID, dateID <= cal_end_dateID) %>%
          mutate(dateID = as.numeric(as.character(dateID))) %>%
          select(siteID, plotID, dateID, all_of(taxon_name)) %>%
          rename(truth = !!taxon_name) %>%
          left_join(trained_time_map, by="dateID") %>%
          mutate(plot_num = pnum, site_num = snum, timepoint = NA_real_, date_num = trained_date_num) %>%
          select(siteID, plotID, dateID, dates, truth, plot_num, date_num, site_num, timepoint)
      }
      
      # Hindcast truth (may be empty)
      condition_check <- (length(taxon_name) == 1L && taxon_name %in% names(rank.df_spec))
      if (!condition_check) {
        hind_truth <- data.frame(siteID = character(0), plotID = character(0), dateID = numeric(0), 
                                  dates = as.Date(character(0)), truth = numeric(0), plot_num = integer(0), 
                                  date_num = numeric(0), site_num = integer(0), timepoint = numeric(0),
                                  stringsAsFactors = FALSE)
      } else {
        hind_truth <- rank.df_spec %>%
          filter(plotID == !!plotID, dateID > cal_end_dateID) %>%
          ensure_dates() %>%
          select(siteID, plotID, dateID, all_of(taxon_name)) %>%
          rename(truth = !!taxon_name) %>%
          mutate(species = taxon_name) %>%
          left_join(trained_time_map, by="dateID") %>%
          mutate(plot_num = pnum, site_num = snum, timepoint = NA_real_, date_num = trained_date_num) %>%
          select(siteID, plotID, dateID, dates, truth, plot_num, date_num, site_num, timepoint)
      }
      
      # Truth sanity (bounds)
      bound01 <- function(x) { x[!is.finite(x)] <- NA_real_; x }
      if (nrow(cal_truth) > 0) {
        cal_truth$truth <- bound01(as.numeric(cal_truth$truth))
        if (any(cal_truth$truth < 0 | cal_truth$truth > 1, na.rm=TRUE)) {
          stop("Calibration truth outside [0,1] for ", plotID)
        }
      }
      if (nrow(hind_truth) > 0) {
        hind_truth$truth <- bound01(as.numeric(hind_truth$truth))
        if (any(hind_truth$truth < 0 | hind_truth$truth > 1, na.rm=TRUE)) {
          stop("Hindcast truth outside [0,1] for ", plotID)
        }
      }
      
      # Filter truth data to dateIDs in trained_time_map (prevents NA date_num from left_join)
      valid_dateIDs <- trained_time_map$dateID
      
      if (nrow(cal_truth) > 0 && "truth" %in% names(cal_truth)) {
        cal_truth_valid <- cal_truth %>% filter(dateID %in% valid_dateIDs)
      } else {
        cal_truth_valid <- cal_truth  # Keep empty structure
      }
      if (nrow(hind_truth) > 0 && "truth" %in% names(hind_truth)) {
        hind_truth_valid <- hind_truth %>% filter(dateID %in% valid_dateIDs)
      } else {
        hind_truth_valid <- hind_truth  # Keep empty structure
      }
      
      # Warn if we're dropping truth data
      dropped_cal <- nrow(cal_truth) - nrow(cal_truth_valid)
      dropped_hind <- nrow(hind_truth) - nrow(hind_truth_valid)
      if (dropped_cal > 0 || dropped_hind > 0) {
        if (LOCAL_TEST) {
          cat("  WARNING: Dropping", dropped_cal, "calibration and", dropped_hind, 
              "hindcast truth rows with dateIDs not in trained_time_map\n")
          fflush()
        }
      }
      
      # Combine calibration and hindcast truth, fill missing dates
      comb <- bind_rows(cal_truth_valid, hind_truth_valid)
      comb <- ensure_cols(comb, "truth", default = NA_real_)
      missing_ids <- setdiff(trained_time_map$dateID, unique(comb$dateID))
      if (length(missing_ids)) {
        comb <- bind_rows(
          comb,
          trained_time_map %>% filter(dateID %in% missing_ids) %>%
            transmute(siteID = sid, plotID = plotID, dateID, dates, truth = NA_real_,
                      plot_num = pnum, date_num = trained_date_num, site_num = snum, timepoint = NA_real_)
        )
      }
      
      # Validate date_num (re-join if needed)
      na_date_num_count <- sum(is.na(comb$date_num))
      if (na_date_num_count > 0) {
        # Try to fix by re-joining with trained_time_map
        comb <- comb %>%
          select(-any_of(c("date_num", "trained_date_num", "dates"))) %>%
          left_join(trained_time_map %>% select(dateID, trained_date_num, dates), by="dateID") %>%
          mutate(date_num = trained_date_num)
        
        # If still NA, filter them out with warning
        na_date_num_count <- sum(is.na(comb$date_num))
        if (na_date_num_count > 0) {
          if (LOCAL_TEST) {
            cat("  WARNING: Removing", na_date_num_count, "rows with NA date_num after re-join\n")
            fflush()
          }
          comb <- comb %>% filter(!is.na(date_num))
        }
      }
      
      # Final check - if we have no valid rows, skip this plot
      if (nrow(comb) == 0) {
        cat("  WARNING: No valid truth data for plot", plotID, "after filtering - skipping\n")
        fflush()
        next
      }
      
      # Guard against nested parallel (if we're already in parallel mode)
      if (!LOCAL_TEST && requireNamespace("foreach", quietly=TRUE)) {
        workers <- try(foreach::getDoParWorkers(), silent=TRUE)
        if (!inherits(workers, "try-error") && workers > 1) {
          # If inner code tries to use foreach, make it sequential to avoid nested storms
          foreach::registerDoSEQ()
          on.exit({
            try(foreach::registerDoParallel(cl), silent=TRUE)
          }, add=TRUE)
        }
      }
      
      # Reduce Nmc to prevent memory limit errors (stay under 16GB)
      model_name_val <- if (length(row$model_name) == 1) row$model_name else row$model_name[1]
      if (model_name_val == "env_cov") {
        Nmc_to_use <- if (LOCAL_TEST) 150 else 250
      } else if (model_name_val == "env_cycl") {
        Nmc_to_use <- if (LOCAL_TEST) 150 else 300
      } else {
        Nmc_to_use <- if (LOCAL_TEST) 200 else 350
      }
      
      gc(verbose = FALSE, full = TRUE)
      
      # forecast
      hp <- fcast_logit_beta(
        plotID, full_in, param_samples, comb,
        plot_summary = cal_quants,
        Nmc = Nmc_to_use,
        rank.name = rank_name,
        model_id = row$model_id,
        metadata = list(species = row$taxon,
                        cal_end_dateID = cal_end_dateID,
                        trained_time_map = trained_time_map,
                        link_function = link_function)
      )
      
      # CRITICAL: Immediately clear forecast output and force garbage collection
      # This prevents memory accumulation across multiple plots/models
      if (!is.null(hp) && nrow(hp) > 0) {
        # Save hp to a temporary variable, then clear workspace
        hp_result <- hp
        rm(hp)
        gc(verbose = FALSE, full = TRUE)
        hp <- hp_result
      }
      
      # Aggressive garbage collection after forecast to free memory
      gc(verbose = FALSE, full = TRUE)
      
      if (is.null(hp) || nrow(hp) == 0) {
        cat("WARN: empty forecast for plot", plotID, "\n"); next
      }
      
      # Remove duplicates from fcast_logit_beta output
      if ("dateID" %in% names(hp) && "plotID" %in% names(hp)) {
        n_before <- nrow(hp)
        hp <- hp %>%
          mutate(dateID = as.numeric(as.character(dateID))) %>%
          arrange(dateID) %>%
          distinct(plotID, dateID, .keep_all = TRUE)
        if (LOCAL_TEST && nrow(hp) < n_before) {
          cat("  ⚠️  Removed", n_before - nrow(hp), "duplicate rows from fcast_logit_beta output\n")
          fflush()
        }
      }
      
      # ensure columns and annotate
      for (k in c("lo","lo_25","med","hi_75","hi")) if (!k %in% names(hp)) hp[[k]] <- NA_real_
      if (!"dates" %in% names(hp) && "dateID" %in% names(hp)) hp$dates <- as.Date(paste0(hp$dateID, "01"), "%Y%m%d")
      
      # Merge truth data from comb (average per dateID to get plot-level means)
      if ("truth" %in% names(comb) && nrow(comb) > 0) {
        truth_from_comb_early <- comb %>%
          filter(!is.na(truth), is.finite(truth)) %>%
          mutate(dateID = as.numeric(as.character(dateID))) %>%
          group_by(plotID, dateID) %>%
          summarise(
            truth = mean(truth, na.rm = TRUE),
            dates = dplyr::first(dates),  # Keep first date for each dateID
            .groups = "drop"
          )
      } else {
        # Create empty structure if no truth data
        truth_from_comb_early <- data.frame(plotID = character(0), dateID = numeric(0), 
                                            truth = numeric(0), dates = as.Date(character(0)),
                                            stringsAsFactors = FALSE)
      }
      
      # Merge truth to hp
      if (nrow(truth_from_comb_early) > 0 && "dateID" %in% names(hp)) {
        if ("truth" %in% names(hp)) {
          hp <- hp %>%
            mutate(dateID = as.numeric(as.character(dateID))) %>%
            left_join(truth_from_comb_early %>% select(plotID, dateID, dates, truth),
                      by = c("plotID", "dateID"), suffix = c("", "_from_comb")) %>%
            mutate(
              truth = dplyr::coalesce(truth_from_comb, truth),
              dates = dplyr::coalesce(dates_from_comb, dates)
            ) %>%
            select(-ends_with("_from_comb"))
        } else {
          # If hp doesn't have truth column, add it
          hp <- hp %>%
            mutate(dateID = as.numeric(as.character(dateID))) %>%
            left_join(truth_from_comb_early %>% select(plotID, dateID, truth), by=c("plotID", "dateID"))
        }
        
        # Safety check: ensure truth values are valid
        safety_condition <- is.na(hp$truth) | hp$dateID %in% comb$dateID[!is.na(comb$truth)]
        if (length(safety_condition) > 0 && !all(safety_condition)) {
          stop("Safety check failed: some truth values don't match comb$dateID")
        }
        
        # Add truth rows that don't have matching forecasts
        missing_truth <- truth_from_comb_early %>%
          filter(!dateID %in% hp$dateID)
        if (nrow(missing_truth) > 0) {
          # Ensure dates column is Date type (same as hp)
          missing_truth$dates <- as.Date(missing_truth$dates)
          
          # Create rows for truth without forecasts
          missing_rows <- missing_truth %>%
            mutate(
              plotID = plotID,
              siteID = sid,
              med = NA_real_, lo = NA_real_, hi = NA_real_,
              lo_25 = NA_real_, hi_75 = NA_real_,
              mean = NA_real_, sd = NA_real_,
              fcast_period = ifelse(!is.na(dateID) & dateID > cal_end_dateID, "hindcast", "calibration"),
              species = row$taxon,
              taxon = row$taxon,
              model_name = row$model_name,
              time_period = paste0(row$min_date, "_", row$max_date),
              rank_name = rank_name,
              predicted_site_effect = FALSE,
              newsite = "Observed site",
              model_id = row$model_id,
              site_prediction = "New time (observed site)"
            )
          
          # Add missing columns to match hp structure
          missing_cols <- setdiff(names(hp), names(missing_rows))
          if (length(missing_cols) > 0) {
            for (col in missing_cols) {
              # Use same type as in hp
              if (col %in% names(hp)) {
                missing_rows[[col]] <- switch(class(hp[[col]][1]),
                  "character" = NA_character_,
                  "numeric" = NA_real_,
                  "integer" = NA_integer_,
                  "Date" = as.Date(NA),
                  "logical" = NA,
                  NA
                )
              }
            }
          }
          
          missing_rows <- missing_rows[, names(hp), drop=FALSE]
          # Ensure dates are same type before binding
          if ("dates" %in% names(hp) && "dates" %in% names(missing_rows)) {
            # Get the type of dates in hp
            if (inherits(hp$dates, "Date")) {
              missing_rows$dates <- as.Date(missing_rows$dates)
            } else if (is.numeric(hp$dates)) {
              # If hp has numeric dates, convert missing_rows dates to numeric
              missing_rows$dates <- as.numeric(missing_rows$dates)
            }
          }
          
          hp <- data.table::rbindlist(list(hp, missing_rows), fill=TRUE, use.names=TRUE)
          
          # Ensure dates remain as Date objects after rbindlist
          if ("dateID" %in% names(hp) && "dates" %in% names(hp)) {
            if (!inherits(hp$dates, "Date")) {
              dateID_val <- hp$dateID[1]
              if (length(dateID_val) == 1 && !is.na(dateID_val)) {
                hp$dates <- as.Date(paste0(hp$dateID, "01"), "%Y%m%d")
              }
            }
          }
        }
      }
      
      # Ensure dates are Date objects
      hp <- ensure_dates(hp)
      
      hp <- hp %>%
        mutate(model_name = row$model_name,
               time_period = paste0(row$min_date, "_", row$max_date),
               species = row$taxon,
               taxon = row$taxon,
               rank_name = rank_name,
               predicted_site_effect = FALSE,
               newsite = "Observed site",
               model_id = row$model_id,
               site_prediction = "New time (observed site)",
               siteID = sid)
      
      # Add calibration ribbons
      if (!"plotID" %in% names(cal_quants)) {
        cat("  ERROR: cal_quants missing plotID column - cannot filter for plot", plotID, "\n")
        fflush()
        next
      }
      
      pc_raw <- cal_quants %>% filter(plotID == !!plotID)
      
      if (nrow(pc_raw) == 0) {
        if (LOCAL_TEST) {
          cat("  WARNING: No calibration quantile data for plot", plotID, "- skipping calibration ribbons\n")
          fflush()
        }
        pc <- data.frame()
      } else {
      if (LOCAL_TEST && nrow(pc_raw) > 0) {
        cat("  cal_quants filtered by plotID:", nrow(pc_raw), "rows\n")
        if ("timepoint" %in% names(pc_raw)) {
          pc_timepoints <- sort(unique(pc_raw$timepoint))
          cat("  cal_quants timepoint range: [", min(pc_timepoints, na.rm=TRUE), ", ", 
              max(pc_timepoints, na.rm=TRUE), "] (", length(pc_timepoints), " unique timepoints)\n", sep="")
          time_train_timepoints <- sort(unique(time_train$timepoint))
          cat("  time_train timepoint range: [", min(time_train_timepoints, na.rm=TRUE), ", ", 
              max(time_train_timepoints, na.rm=TRUE), "] (", length(time_train_timepoints), " unique timepoints)\n", sep="")
          missing_in_train <- setdiff(pc_timepoints, time_train_timepoints)
          if (length(missing_in_train) > 0) {
            cat("  WARNING:", length(missing_in_train), "timepoints in cal_quants not in time_train:",
                paste(head(missing_in_train, 5), collapse=", "), ifelse(length(missing_in_train) > 5, "...", ""), "\n")
          }
        }
        fflush()
      }
      
      # Handle duplicate dateID/dates columns from join (coalesce to use merged result)
      if (!"dates" %in% names(time_train)) {
        # Reconstruct dates from dateID if missing
        time_train$dates <- as.Date(paste0(time_train$dateID, "01"), "%Y%m%d")
      }
      
      pc <- pc_raw %>% left_join(
        time_train %>% select(timepoint, dateID, dates), by="timepoint", suffix = c("_from_quants", "_from_train")
      )
      
      # Ensure plotID and siteID are preserved after join
      if (!"plotID" %in% names(pc) && "plotID" %in% names(pc_raw)) {
        pc$plotID <- pc_raw$plotID
      }
      if (!"siteID" %in% names(pc) && "siteID" %in% names(pc_raw)) {
        pc$siteID <- pc_raw$siteID
      }
      
      # Coalesce dateID and dates columns after join (prefer time_train's version)
      if ("dateID_from_train" %in% names(pc)) {
        pc <- pc %>% mutate(dateID_from_train = suppressWarnings(as.numeric(as.character(dateID_from_train))))
        if ("dateID_from_quants" %in% names(pc)) {
          pc <- pc %>% mutate(
            dateID_from_quants = suppressWarnings(as.numeric(as.character(dateID_from_quants)))
          )
          if ("dateID" %in% names(pc)) {
            pc <- pc %>% mutate(dateID = dplyr::coalesce(dateID_from_train, dateID_from_quants, dateID))
          } else {
            pc <- pc %>% mutate(dateID = dplyr::coalesce(dateID_from_train, dateID_from_quants))
          }
        } else {
          if ("dateID" %in% names(pc)) {
            pc <- pc %>% mutate(dateID = dplyr::coalesce(dateID_from_train, dateID))
          } else {
            pc <- pc %>% mutate(dateID = dateID_from_train)
          }
        }
        pc <- pc %>% select(-any_of(c("dateID_from_quants", "dateID_from_train")))
      } else if ("dateID.y" %in% names(pc)) {
        pc <- pc %>% mutate(
          dateID.y = suppressWarnings(as.numeric(as.character(dateID.y)))
        )
        if ("dateID.x" %in% names(pc)) {
          pc <- pc %>% mutate(
              dateID.x = suppressWarnings(as.numeric(as.character(dateID.x))),
              dateID = dplyr::coalesce(dateID.y, dateID.x)
          )
        } else {
          pc <- pc %>% mutate(dateID = dateID.y)
        }
        pc <- pc %>% select(-any_of(c("dateID.x", "dateID.y")))
        } else if ("dateID" %in% names(pc)) {
          # Ensure existing dateID is numeric
          pc <- pc %>% mutate(dateID = suppressWarnings(as.numeric(as.character(dateID))))
      }
      
      # Coalesce dates columns (prefer time_train's version)
      # CRITICAL FIX: Ensure all dates columns are Date objects before coalescing to avoid type mismatch
      # IMPORTANT: We must preserve the original dates values, not reconstruct from dateID (which could be different)
      if ("dates_from_train" %in% names(pc)) {
        # Convert to Date if not already - preserve original values, don't reconstruct
        if (!inherits(pc$dates_from_train, "Date")) {
          # Check what type it actually is and convert appropriately
          if (is.character(pc$dates_from_train)) {
            pc <- pc %>% mutate(
              dates_from_train = as.Date(dates_from_train)
            )
          } else if (is.numeric(pc$dates_from_train)) {
            # If numeric, it might be dateID format (YYYYMM) or POSIX timestamp
            # Try dateID format first (most common)
            pc <- pc %>% mutate(
              dates_from_train = tryCatch(
                as.Date(paste0(as.character(dates_from_train), "01"), "%Y%m%d"),
                error = function(e) {
                  # If that fails, try as POSIX timestamp
                  tryCatch(
                    as.Date(as.POSIXct(dates_from_train, origin = "1970-01-01")),
                    error = function(e2) {
                      # Last resort: reconstruct from dateID if available, but WARN
                      if ("dateID" %in% names(pc)) {
                        warning("CRITICAL: Could not convert dates_from_train (numeric) to Date. ",
                                "Reconstructing from dateID - this may not match original value! ",
                                "Original value: ", paste(head(unique(dates_from_train), 3), collapse=", "))
                        as.Date(paste0(dateID, "01"), "%Y%m%d")
                      } else {
                        stop("CRITICAL: Cannot convert dates_from_train to Date and no dateID available. ",
                             "Original value type: ", class(dates_from_train), 
                             " Sample values: ", paste(head(unique(dates_from_train), 3), collapse=", "))
                      }
                    }
                  )
                }
              )
            )
          } else {
            # Unknown type - try direct conversion, but error if it fails
            pc <- pc %>% mutate(
              dates_from_train = tryCatch(
                as.Date(dates_from_train),
                error = function(e) {
                  stop("CRITICAL: Cannot convert dates_from_train to Date. ",
                       "Type: ", class(dates_from_train),
                       " Error: ", e$message,
                       " Sample values: ", paste(head(unique(dates_from_train), 3), collapse=", "))
                }
              )
            )
          }
        }
        if ("dates_from_quants" %in% names(pc)) {
          if (!inherits(pc$dates_from_quants, "Date")) {
            # Same conversion logic as above
            if (is.character(pc$dates_from_quants)) {
              pc <- pc %>% mutate(
                dates_from_quants = as.Date(dates_from_quants)
              )
            } else if (is.numeric(pc$dates_from_quants)) {
              pc <- pc %>% mutate(
                dates_from_quants = tryCatch(
                  as.Date(paste0(as.character(dates_from_quants), "01"), "%Y%m%d"),
                  error = function(e) {
                    tryCatch(
                      as.Date(as.POSIXct(dates_from_quants, origin = "1970-01-01")),
                      error = function(e2) {
                        if ("dateID" %in% names(pc)) {
                          warning("CRITICAL: Could not convert dates_from_quants (numeric) to Date. ",
                                  "Reconstructing from dateID - this may not match original value!")
                          as.Date(paste0(dateID, "01"), "%Y%m%d")
                        } else {
                          stop("CRITICAL: Cannot convert dates_from_quants to Date and no dateID available.")
                        }
                      }
                    )
                  }
                )
              )
            } else {
              pc <- pc %>% mutate(
                dates_from_quants = tryCatch(
                  as.Date(dates_from_quants),
                  error = function(e) {
                    stop("CRITICAL: Cannot convert dates_from_quants to Date. Type: ", class(dates_from_quants))
                  }
                )
              )
            }
          }
          pc <- pc %>% mutate(dates = dplyr::coalesce(dates_from_train, dates_from_quants))
        } else {
          pc <- pc %>% mutate(dates = dates_from_train)
        }
        pc <- pc %>% select(-any_of(c("dates_from_quants", "dates_from_train")))
      } else if ("dates.y" %in% names(pc)) {
        # Same conversion logic
        if (!inherits(pc$dates.y, "Date")) {
          if (is.character(pc$dates.y)) {
            pc <- pc %>% mutate(dates.y = as.Date(dates.y))
          } else if (is.numeric(pc$dates.y)) {
            pc <- pc %>% mutate(
              dates.y = tryCatch(
                as.Date(paste0(as.character(dates.y), "01"), "%Y%m%d"),
                error = function(e) {
                  tryCatch(
                    as.Date(as.POSIXct(dates.y, origin = "1970-01-01")),
                    error = function(e2) {
                      if ("dateID" %in% names(pc)) {
                        warning("CRITICAL: Could not convert dates.y (numeric) to Date. Reconstructing from dateID!")
                        as.Date(paste0(dateID, "01"), "%Y%m%d")
                      } else {
                        stop("CRITICAL: Cannot convert dates.y to Date and no dateID available.")
                      }
                    }
                  )
                }
              )
            )
          } else {
            pc <- pc %>% mutate(
              dates.y = tryCatch(
                as.Date(dates.y),
                error = function(e) {
                  stop("CRITICAL: Cannot convert dates.y to Date. Type: ", class(dates.y))
                }
              )
            )
          }
        }
        if ("dates.x" %in% names(pc)) {
          if (!inherits(pc$dates.x, "Date")) {
            if (is.character(pc$dates.x)) {
              pc <- pc %>% mutate(dates.x = as.Date(dates.x))
            } else if (is.numeric(pc$dates.x)) {
              pc <- pc %>% mutate(
                dates.x = tryCatch(
                  as.Date(paste0(as.character(dates.x), "01"), "%Y%m%d"),
                  error = function(e) {
                    tryCatch(
                      as.Date(as.POSIXct(dates.x, origin = "1970-01-01")),
                      error = function(e2) {
                        if ("dateID" %in% names(pc)) {
                          warning("CRITICAL: Could not convert dates.x (numeric) to Date. Reconstructing from dateID!")
                          as.Date(paste0(dateID, "01"), "%Y%m%d")
                        } else {
                          stop("CRITICAL: Cannot convert dates.x to Date and no dateID available.")
                        }
                      }
                    )
                  }
                )
              )
            } else {
              pc <- pc %>% mutate(
                dates.x = tryCatch(
                  as.Date(dates.x),
                  error = function(e) {
                    stop("CRITICAL: Cannot convert dates.x to Date. Type: ", class(dates.x))
                  }
                )
              )
            }
          }
          pc <- pc %>% mutate(dates = dplyr::coalesce(dates.y, dates.x))
        } else {
          pc <- pc %>% mutate(dates = dates.y)
        }
        pc <- pc %>% select(-any_of(c("dates.x", "dates.y")))
      } else if ("dates" %in% names(pc)) {
        # Ensure existing dates column is Date type
        if (!inherits(pc$dates, "Date")) {
          if (is.character(pc$dates)) {
            pc <- pc %>% mutate(dates = as.Date(dates))
          } else if (is.numeric(pc$dates)) {
            pc <- pc %>% mutate(
              dates = tryCatch(
                as.Date(paste0(as.character(dates), "01"), "%Y%m%d"),
                error = function(e) {
                  tryCatch(
                    as.Date(as.POSIXct(dates, origin = "1970-01-01")),
                    error = function(e2) {
                      if ("dateID" %in% names(pc)) {
                        warning("CRITICAL: Could not convert dates (numeric) to Date. Reconstructing from dateID!")
                        as.Date(paste0(dateID, "01"), "%Y%m%d")
                      } else {
                        stop("CRITICAL: Cannot convert dates to Date and no dateID available.")
                      }
                    }
                  )
                }
              )
            )
          } else {
            pc <- pc %>% mutate(
              dates = tryCatch(
                as.Date(dates),
                error = function(e) {
                  stop("CRITICAL: Cannot convert dates to Date. Type: ", class(dates))
                }
              )
            )
          }
        }
      }
      
      # Ensure dateID and dates are correct types
      if ("dateID" %in% names(pc)) {
        pc <- pc %>% mutate(dateID = as.numeric(as.character(dateID)))
      }
      if ("dates" %in% names(pc) && !all(is.na(pc$dates))) {
        pc <- pc %>% mutate(dates = as.Date(dates))
      }
      
      # CRITICAL DEBUG: Check join results
      if (LOCAL_TEST && nrow(pc) > 0) {
        pc_with_dateID <- sum(!is.na(pc$dateID))
        pc_without_dateID <- sum(is.na(pc$dateID))
        cat("  After join: ", pc_with_dateID, " rows with dateID, ", pc_without_dateID, " rows without dateID\n", sep="")
        if (pc_without_dateID > 0 && "timepoint" %in% names(pc)) {
          missing_dateID_timepoints <- unique(pc$timepoint[is.na(pc$dateID)])
          cat("  Timepoints missing dateID:", paste(head(missing_dateID_timepoints, 5), collapse=", "),
              ifelse(length(missing_dateID_timepoints) > 5, "...", ""), "\n")
        }
        fflush()
      }
      
      # CRITICAL FIX: Remove duplicates - ensure one row per dateID
      # The join might create duplicates if there are multiple timepoints per dateID
      # CRITICAL: Prioritize rows with quantile data when deduplicating
      if (nrow(pc) > 0 && "dateID" %in% names(pc)) {
        # Check if quantile columns exist before deduplication
        quant_cols <- grep("^2\\.5%|^25%|^50%|^75%|^97\\.5%", names(pc), value=TRUE)
        if (length(quant_cols) > 0) {
          # Create a helper column to prioritize rows with non-NA quantiles
          # Use the median (50%) column if available, otherwise use the first quantile column
          priority_col <- if("50%" %in% quant_cols) "50%" else quant_cols[1]
          pc <- pc %>%
            mutate(
              dateID = as.numeric(as.character(dateID)),
              has_quantile = !is.na(.data[[priority_col]]) & is.finite(.data[[priority_col]])
            ) %>%
            # Sort so rows with quantiles come first, then by timepoint to ensure consistent ordering
            arrange(desc(has_quantile), timepoint) %>%
            distinct(dateID, .keep_all = TRUE) %>%
            select(-has_quantile) %>%
            arrange(dateID)  # Sort by dateID for clean plotting
        } else {
          # No quantile columns found - use standard deduplication
        pc <- pc %>%
          mutate(dateID = as.numeric(as.character(dateID))) %>%
          distinct(dateID, .keep_all = TRUE) %>%
          arrange(dateID)  # Sort by dateID for clean plotting
        }
      }
      
      # CRITICAL DEBUG: Check if early timepoints have dateID after deduplication
      if (LOCAL_TEST && nrow(pc) > 0 && "dateID" %in% names(pc)) {
        pc_by_dateID <- pc %>% arrange(dateID) %>% filter(!is.na(dateID))
        if (nrow(pc_by_dateID) > 0) {
          early_dateIDs <- head(sort(pc_by_dateID$dateID), 5)
          late_dateIDs <- tail(sort(pc_by_dateID$dateID), 5)
          cat("  After dedup: ", nrow(pc_by_dateID), " rows with dateID\n", sep="")
          cat("  Early dateIDs:", paste(early_dateIDs, collapse=", "), "\n")
          cat("  Late dateIDs:", paste(late_dateIDs, collapse=", "), "\n")
        }
        fflush()
      }
        
        if (nrow(pc) > 0) {
        # Debug: show available columns in cal_quants before filtering
        if (LOCAL_TEST) {
          cat("cal_quants column names:\n")
          print(names(cal_quants))
          cat("\n")
        }
        
        # Debug: show available columns in pc
        cat("  Available columns in pc:", paste(colnames(pc), collapse=", "), "\n")
        fflush()
        
        # CRITICAL FIX: Use robust pick_q() function that recognizes all common quantile name variants
        # This handles R's make.names() style (e.g., X2.5., X25., X50., X75., X97.5.) as well as other formats
        pick_q <- function(df, want = c("2.5","25","50","75","97.5")) {
          nms <- names(df)
          
          # patterns for each target quantile (case-insensitive, handles multiple naming conventions)
          pats <- list(
            "2.5"   = c("^x?2\\.5%?$", "^x?2\\.5\\.?$", "^q0?\\.025$","^lo$","^lo_?0?25$","^x?0?2_?5%?$", "^2\\.5%?$", "^2\\.5\\.?$"),
            "25"    = c("^x?25%?$",    "^x?25\\.?$",    "^q0?25$","^lo_?25$", "^25%?$", "^25\\.?$"),
            "50"    = c("^x?50%?$",    "^x?50\\.?$",    "^median$","^q0?5$","^med$", "^50%?$", "^50\\.?$"),
            "75"    = c("^x?75%?$",    "^x?75\\.?$",    "^q0?75$","^hi_?75$", "^75%?$", "^75\\.?$"),
            "97.5"  = c("^x?97\\.5%?$","^x?97\\.5\\.?$","^q0?975$","^hi$","^hi_?97_?5$", "^97\\.5%?$", "^97\\.5\\.?$")
          )
          
          finder <- function(keys) {
            hit <- NA_character_
            for (p in keys) {
              g <- grep(p, nms, ignore.case = TRUE, value = TRUE)
              if (length(g)) { hit <- g[1]; break }
            }
            hit
          }
          
          # Extract quantiles with safe fallback to NA if column not found
          get_col <- function(col_name) {
            found_col <- finder(pats[[col_name]])
            if (is.na(found_col)) {
              if (LOCAL_TEST) {
                cat("  WARNING: Could not find column matching", col_name, "quantile pattern\n")
                fflush()
              }
              rep(NA_real_, nrow(df))
            } else {
              if (LOCAL_TEST) {
                cat("  Found", col_name, "quantile in column:", found_col, "\n")
                fflush()
              }
                # CRITICAL FIX: Handle columns with special characters properly
                # Access column using double brackets for columns with special chars like "2.5%"
                raw_values <- df[[found_col]]
                
                # Debug: Check what type we got
                if (LOCAL_TEST && nrow(df) > 0) {
                  cat("    Column type:", class(raw_values), "\n")
                  cat("    First 3 values:", paste(head(raw_values, 3), collapse=", "), "\n")
                  cat("    NA count:", sum(is.na(raw_values)), "out of", length(raw_values), "\n")
                  fflush()
                }
                
                # Convert to numeric - handle character/character-like columns
                if (is.character(raw_values) || is.factor(raw_values)) {
                  # Try to convert character to numeric
                  numeric_values <- suppressWarnings(as.numeric(raw_values))
                  if (LOCAL_TEST && any(is.na(numeric_values) & !is.na(raw_values))) {
                    cat("    WARNING: Some values failed to convert from", class(raw_values), "to numeric\n")
                    cat("    First failing value:", head(raw_values[is.na(numeric_values) & !is.na(raw_values)], 1), "\n")
                    fflush()
                  }
                  numeric_values
                } else if (is.numeric(raw_values)) {
                  raw_values
                } else if (is.list(raw_values)) {
                  # Handle list columns (shouldn't happen but be safe)
                  if (LOCAL_TEST) {
                    cat("    WARNING: Column is a list, extracting first element\n")
                    fflush()
                  }
                  suppressWarnings(as.numeric(sapply(raw_values, function(x) if (is.null(x)) NA else x[1])))
                } else {
                  suppressWarnings(as.numeric(raw_values))
                }
            }
          }
          
          list(
            lo    = get_col("2.5"),
            lo_25 = get_col("25"),
            med   = get_col("50"),
            hi_75 = get_col("75"),
            hi    = get_col("97.5")
          )
        }
        
        qcols <- pick_q(pc)
        pc$lo    <- qcols$lo
        pc$lo_25 <- qcols$lo_25
        pc$med   <- qcols$med
        pc$hi_75 <- qcols$hi_75
        pc$hi    <- qcols$hi
        
        # CRITICAL DEBUG: Check if quantiles are NA in the source data
        # NOTE: If cal_quants (model_summary[[3]]) has NA quantiles for some timepoints,
        # ribbons will be missing for those timepoints. This is expected when the model
        # summary doesn't have quantile data for those timepoints (e.g., no data during
        # training period, or convergence issues).
        if (LOCAL_TEST && sum(is.na(pc$med)) > 0) {
          na_rows <- which(is.na(pc$med))
          cat("  NOTE: ", length(na_rows), " rows have NA medians (timepoints:", 
              paste(unique(pc$timepoint[na_rows]), collapse=", "), 
              ") - ribbons will be missing for these timepoints\n", sep="")
          if (length(na_rows) <= 10) {
            cat("  NA quantile rows (first", min(5, length(na_rows)), "):\n")
            for (i in head(na_rows, 5)) {
              cat("    timepoint", pc$timepoint[i], "dateID", pc$dateID[i], 
                  "- 2.5%:", pc[["2.5%"]][i], "50%:", pc[["50%"]][i], "\n")
            }
          }
          fflush()
        }
        
          # CRITICAL GUARD: Check if we have any finite medians for calibration ribbons
          # NOTE: Even if ALL calibration ribbons are NA, the hindcast can still run because
          # fcast_logit_beta uses samples2_matrix/plot_summary (not cal_quants) for initial conditions
          # and finds the last non-NA value. So we only warn, but don't skip the hindcast.
          # CRITICAL FIX: Check if med column exists before checking finite values
          if ("med" %in% names(pc)) {
            finite_med_count <- sum(is.finite(pc$med), na.rm=TRUE)
            if (finite_med_count == 0) {
              # No calibration ribbons available, but hindcast can still proceed
              cat("  WARNING: No finite calibration quantiles found for plotting ribbons!\n")
              cat("    Available columns:", paste(names(pc), collapse=", "), "\n")
              cat("    NOTE: Hindcast will still proceed - it uses last non-NA value from samples2_matrix/plot_summary\n")
              fflush()
              pc <- data.frame()  # Set to empty so we skip calibration ribbon plotting only
            } else {
              # Some medians are finite - we can plot calibration ribbons
              na_med_count <- sum(is.na(pc$med), na.rm=TRUE)
              if (LOCAL_TEST && na_med_count > 0) {
                cat("  NOTE: Plot", plotID, "has", finite_med_count, "finite and", na_med_count, 
                    "NA medians - will plot ribbons for finite timepoints only\n")
                fflush()
              }
            }
          } else {
            cat("  WARNING: 'med' column missing from pc - calibration ribbons will not be available\n")
            fflush()
            pc <- data.frame()  # Set to empty so we skip calibration ribbon plotting only
          }
        }
      }
        
      # Only process calibration quantiles if we have data
      if (nrow(pc) > 0) {
        # CRITICAL DEBUG: Check if early timepoints have NA quantiles
        if (LOCAL_TEST && "dateID" %in% names(pc)) {
          pc_by_dateID <- pc %>% arrange(dateID) %>% filter(!is.na(dateID))
          if (nrow(pc_by_dateID) > 0) {
            # Check quantiles for early vs late timepoints
            early_rows <- head(pc_by_dateID, 5)
            late_rows <- tail(pc_by_dateID, 5)
            
            early_med_na <- sum(is.na(early_rows$med))
            late_med_na <- sum(is.na(late_rows$med))
            
            cat("  Quantile check:\n")
            cat("    Early rows (first 5):", early_med_na, "NA medians out of", nrow(early_rows), "\n")
            cat("    Late rows (last 5):", late_med_na, "NA medians out of", nrow(late_rows), "\n")
            if (early_med_na > 0) {
              cat("    Early rows with NA med - dateIDs:", 
                  paste(early_rows$dateID[is.na(early_rows$med)], collapse=", "), "\n")
              cat("    Early rows with NA med - timepoints:", 
                  paste(early_rows$timepoint[is.na(early_rows$med)], collapse=", "), "\n")
            }
            fflush()
          }
        }
        
        # Debug: show quantile values after extraction
        if (LOCAL_TEST) {
          cat("  Summary of pc quantile columns:\n")
          print(summary(pc[, c("lo","med","hi")]))
          cat("  Quantile range: [", sprintf("%.6f", min(pc$med, na.rm=TRUE)), ", ", 
              sprintf("%.6f", max(pc$med, na.rm=TRUE)), "]\n", sep="")
          fflush()
        }
        
        # cal_quants from model_summary[[3]] are already on [0,1] scale (from plot_mu)
        # No transformation needed - use Mean and SD from pred.quantiles
        pc <- pc %>% mutate(
          mean  = if("Mean" %in% names(.)) as.numeric(Mean) else NA_real_,
          sd    = if("SD" %in% names(.)) as.numeric(SD) else NA_real_,
          fcast_period = "calibration",
          model_name = if (length(row$model_name) == 1) row$model_name else row$model_name[1],
          time_period = paste0(if (length(row$min_date) == 1) row$min_date else row$min_date[1], "_", if (length(row$max_date) == 1) row$max_date else row$max_date[1]),
          species = taxon_name, taxon = taxon_name, rank_name = rank_name,
          predicted_site_effect = FALSE, newsite = "Observed site", model_id = if (length(row$model_id) == 1) row$model_id else row$model_id[1],
          site_prediction = "New time (observed site)", siteID = sid
        ) %>%
        select(-any_of(c("Mean", "SD")))  # drop uppercase to avoid duplicates in rbindlist
        
        # Ensure all timepoints have dateID and dates from time_train
        if ("timepoint" %in% names(pc)) {
          pc <- pc %>%
            select(-any_of(c("dateID", "dates"))) %>%
            left_join(
              time_train %>% select(timepoint, dateID, dates) %>% mutate(dateID = as.numeric(as.character(dateID))),
              by="timepoint"
            )
        }
        
        # Add calibration truth data to pc (average per dateID to get plot-level means)
        if (nrow(cal_truth) > 0 && "dateID" %in% names(pc) && "dateID" %in% names(cal_truth)) {
          pc <- pc %>% mutate(dateID = as.numeric(as.character(dateID)))
          pc_has_dateID <- sum(!is.na(pc$dateID))
          
          if (!"truth" %in% names(cal_truth)) {
            if (LOCAL_TEST) {
              cat("WARN: cal_truth missing 'truth' column - skipping truth join\n")
              fflush()
            }
            cal_truth_for_join <- data.frame(dateID = numeric(0), truth = numeric(0))
          } else {
            cal_truth_for_join <- cal_truth %>%
              mutate(dateID = as.numeric(as.character(dateID))) %>%
              filter(!is.na(dateID), !is.na(truth), is.finite(truth)) %>%
              group_by(dateID) %>%
              summarise(truth = mean(truth, na.rm = TRUE), .groups = "drop")
          }
          
          if (nrow(cal_truth_for_join) > 0 && pc_has_dateID > 0) {
            pc <- pc %>%
              left_join(cal_truth_for_join, by="dateID")
            
            # CRITICAL FIX: Check for dates with truth but zero predictions
            # This happens when plot_start_date is set after the first observation
            # Extract from samples2 (plot_mu) for those timepoints if available
            # CRITICAL FIX: Check if truth column exists before filtering
            if ("truth" %in% names(pc)) {
              dates_with_truth_but_zero_pred <- pc %>%
                filter(!is.na(truth), !is.na(dateID), 
                       (is.na(med) | med == 0) & (is.na(lo) | lo == 0) & (is.na(hi) | hi == 0))
            } else {
              dates_with_truth_but_zero_pred <- pc[0, ]
            }
            
            if (nrow(dates_with_truth_but_zero_pred) > 0 && exists("model_dat") && 
                !is.null(model_dat$samples2) && exists("plot_map")) {
              plot_info <- plot_map %>% filter(plotID == !!plotID) %>% distinct()
              if (nrow(plot_info) > 1) {
                plot_info <- plot_info[1, , drop=FALSE]
              }
              if (nrow(plot_info) > 0) {
                plot_idx <- plot_info$plot_num[1]
                samples2_matrix <- NULL
                if (inherits(model_dat$samples2, "mcmc.list")) {
                  samples2_matrix <- as.matrix(do.call(rbind, model_dat$samples2))
                } else if (is.matrix(model_dat$samples2)) {
                  samples2_matrix <- model_dat$samples2
                }
                
                if (!is.null(samples2_matrix) && !is.null(colnames(samples2_matrix))) {
                  # Extract plot_mu columns for this plot
                  cal_mu_cols <- grep(paste0("^plot_mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
                  if (length(cal_mu_cols) == 0) {
                    cal_mu_cols <- grep(paste0("^mu\\[", plot_idx, ","), colnames(samples2_matrix), value = TRUE)
                  }
                  
                  if (length(cal_mu_cols) > 0) {
                    # Extract timepoint numbers from column names
                    timepoint_nums <- as.numeric(gsub(".*,\\s*([^]]+)\\]", "\\1", cal_mu_cols))
                    
                    # Map dateIDs to timepoints using trained_time_map
                    for (i in 1:nrow(dates_with_truth_but_zero_pred)) {
                      dateID_val <- dates_with_truth_but_zero_pred$dateID[i]
                      tp_match <- which(trained_time_map$dateID == dateID_val)
                      if (length(tp_match) > 0) {
                        tp <- trained_time_map$trained_date_num[tp_match[1]]
                        col_match <- which(timepoint_nums == tp)
                        if (length(col_match) > 0) {
                          # Extract samples for this timepoint
                          col_name <- cal_mu_cols[col_match[1]]
                          if (col_name %in% colnames(samples2_matrix)) {
                            tp_samples <- samples2_matrix[, col_name, drop = FALSE]
                            if (ncol(tp_samples) > 0 && sum(!is.na(tp_samples)) > 0) {
                              # Calculate quantiles
                              tp_quantiles <- quantile(tp_samples[, 1], c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
                              tp_mean <- mean(tp_samples[, 1], na.rm = TRUE)
                              tp_sd <- sd(tp_samples[, 1], na.rm = TRUE)
                              
                              # Update pc for this dateID
                              pc_idx <- which(pc$dateID == dateID_val)
                              if (length(pc_idx) > 0) {
                                pc$lo[pc_idx[1]] <- tp_quantiles[1]
                                pc$lo_25[pc_idx[1]] <- tp_quantiles[2]
                                pc$med[pc_idx[1]] <- tp_quantiles[3]
                                pc$hi_75[pc_idx[1]] <- tp_quantiles[4]
                                pc$hi[pc_idx[1]] <- tp_quantiles[5]
                                pc$mean[pc_idx[1]] <- tp_mean
                                pc$sd[pc_idx[1]] <- tp_sd
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            
          }
        }
        
        # Preserve truth data before removing calibration rows from hp
        truth_data <- NULL
        if ("truth" %in% names(hp) && sum(!is.na(hp$truth)) > 0) {
          truth_data <- hp %>%
            filter(!is.na(truth), is.finite(truth)) %>%
            select(plotID, dateID, truth, dates) %>%
            group_by(plotID, dateID) %>%
            summarise(
              truth = mean(truth, na.rm = TRUE),
              dates = dplyr::first(dates),
              .groups = "drop"
            )
          if (LOCAL_TEST && nrow(truth_data) > 0) {
            cat("  Preserving", nrow(truth_data), "truth rows from hp before deduplication\n")
            fflush()
          }
        }
        
        # Remove calibration period rows from hp (use pc for calibration quantiles instead)
        if ("fcast_period" %in% names(hp)) {
          hp_cal <- hp %>% filter(fcast_period == "calibration")
          hp_hind <- hp %>% filter(fcast_period == "hindcast")
          if (LOCAL_TEST && nrow(hp_cal) > 0) {
            cat("  Removing", nrow(hp_cal), "calibration rows from hp (will use pc instead)\n")
            fflush()
          }
          hp <- hp_hind  # Use only hindcast period from hp
        } else {
          # If fcast_period not set, filter by dateID to remove calibration rows
          if ("dateID" %in% names(hp)) {
            hp <- hp %>% mutate(dateID = as.numeric(as.character(dateID)))
            if (LOCAL_TEST) {
              removed <- sum(hp$dateID <= cal_end_dateID, na.rm=TRUE)
              if (removed > 0) {
                cat("  Removed", removed, "calibration rows from hp (filtered by dateID >", cal_end_dateID, ")\n")
                fflush()
              }
            }
            hp <- hp %>% filter(dateID > cal_end_dateID)  # Only keep hindcast period (after calibration end)
          }
          hp$fcast_period <- "hindcast"
        }
        
        # Filter hp by dateID to remove calibration period rows (definitive filter)
        if ("dateID" %in% names(hp)) {
          hp <- hp %>% mutate(dateID = as.numeric(as.character(dateID)))
          n_before <- nrow(hp)
          # Remove ALL rows with dateID <= cal_end_dateID (calibration period)
          hp <- hp %>% filter(dateID > cal_end_dateID)
          if (LOCAL_TEST && nrow(hp) < n_before) {
            removed <- n_before - nrow(hp)
            cat("  ⚠️  Removed", removed, "rows from hp with dateID <=", cal_end_dateID, "(calibration period)\n")
            fflush()
          }
        }
        
        # Deduplicate hp by dateID (prefer rows with truth if available)
        if ("dateID" %in% names(hp) && nrow(hp) > 0) {
          n_before <- nrow(hp)
          hp <- dedup_key(hp, c("plotID", "dateID"), prefer = "truth", arrange_by = "dateID")
          dbg("  ⚠️  Deduplicated", n_before - nrow(hp), "duplicate rows from hp (multiple rows per dateID)", .on = LOCAL_TEST)
        }
        
        # Preserve truth from pc
        if ("truth" %in% names(pc) && sum(!is.na(pc$truth)) > 0) {
          pc_truth_data <- pc %>%
            filter(!is.na(truth), is.finite(truth)) %>%
            select(plotID, dateID, truth, dates) %>%
            group_by(plotID, dateID) %>%
            summarise(
              truth = mean(truth, na.rm = TRUE),
              dates = dplyr::first(dates),
              .groups = "drop"
            )
          if (!is.null(truth_data)) {
            truth_data <- bind_rows(truth_data, pc_truth_data) %>%
              group_by(plotID, dateID) %>%
              summarise(
                truth = mean(truth, na.rm = TRUE),
                dates = dplyr::first(dates),
                .groups = "drop"
              )  # Average any duplicates from combining hp and pc truth
          } else {
            truth_data <- pc_truth_data
          }
          if (LOCAL_TEST && nrow(pc_truth_data) > 0) {
            cat("  Preserving", nrow(pc_truth_data), "truth rows from pc\n")
            fflush()
          }
        }
        
        # Combine pc (calibration) with hp (hindcast only) and deduplicate
        # Ensure both are data.frames before combining
        if (nrow(pc) == 0 && nrow(hp) == 0) {
          hp <- data.frame()
        } else {
          combine_list <- list()
          if (nrow(pc) > 0) {
            if (!"fcast_period" %in% names(pc)) pc$fcast_period <- "calibration"
            combine_list <- c(combine_list, list(pc))
          }
          if (nrow(hp) > 0) combine_list <- c(combine_list, list(hp))
          
          if (length(combine_list) > 0) {
            hp <- data.table::rbindlist(combine_list, fill = TRUE, use.names = TRUE)
            hp <- as.data.frame(hp)  # Convert to data.frame for consistent dplyr behavior
            hp <- hp %>%
              dplyr::mutate(dateID = as.numeric(as.character(dateID))) %>%
              dplyr::arrange(dateID)
            
            # Deduplicate: prioritize calibration rows for calibration period, hindcast rows for hindcast period, rows with truth
            if ("truth" %in% names(hp)) {
              hp <- hp %>%
                dplyr::mutate(
                  is_calibration_period = dateID <= cal_end_dateID,
                  has_truth = !is.na(truth) & is.finite(truth),
                  priority = dplyr::case_when(
                    is_calibration_period & fcast_period == "calibration" ~ 1,
                    !is_calibration_period & fcast_period == "hindcast" ~ 2,
                    has_truth ~ 3,
                    TRUE ~ 4
                  )
                ) %>%
                dplyr::arrange(priority, dplyr::desc(has_truth), dateID) %>%
                dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                dplyr::select(-priority, -is_calibration_period, -has_truth) %>%
                dplyr::arrange(dateID)
            } else {
              hp <- hp %>%
                dplyr::mutate(
                  is_calibration_period = dateID <= cal_end_dateID,
                  priority = dplyr::case_when(
                    is_calibration_period & fcast_period == "calibration" ~ 1,
                    !is_calibration_period & fcast_period == "hindcast" ~ 2,
                    TRUE ~ 4
                  )
                ) %>%
                dplyr::arrange(priority, dateID) %>%
                dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                dplyr::select(-priority, -is_calibration_period) %>%
                dplyr::arrange(dateID)
            }
            # Force continuity at calibration-hindcast boundary: first hindcast point = last calibration
            # (avoids disconnect from IC sampling or one-step-ahead vs IC display)
            if (nrow(hp) > 0 && "fcast_period" %in% names(hp) && "dateID" %in% names(hp)) {
              cal_rows <- hp %>% dplyr::filter(fcast_period == "calibration") %>% dplyr::arrange(dateID)
              hind_rows <- hp %>% dplyr::filter(fcast_period == "hindcast") %>% dplyr::arrange(dateID)
              if (nrow(cal_rows) > 0 && nrow(hind_rows) > 0) {
                last_cal <- cal_rows %>% dplyr::slice_tail(n = 1)
                first_hind_dateID <- min(hind_rows$dateID, na.rm = TRUE)
                qcols <- c("med", "mean", "lo", "lo_25", "hi_75", "hi")
                qcols <- qcols[qcols %in% names(hp)]
                for (qc in qcols) {
                  if (qc %in% names(last_cal) && !all(is.na(last_cal[[qc]]))) {
                    hp[hp$fcast_period == "hindcast" & hp$dateID == first_hind_dateID, qc] <- last_cal[[qc]][1]
                  }
                }
              }
            }
          } else {
            hp <- data.frame()
          }
        }
        
        # CRITICAL FIX: Restore truth data after deduplication
        # Merge truth back in case it was lost during distinct()
        if (!is.null(truth_data) && nrow(truth_data) > 0 && "dateID" %in% names(hp)) {
          hp <- hp %>%
            select(-any_of("truth")) %>%  # Remove existing truth (might be incomplete)
            left_join(
              truth_data %>% 
                mutate(dateID = as.numeric(as.character(dateID))) %>%
                select(plotID, dateID, truth),
              by = c("plotID", "dateID")
            )
          if (LOCAL_TEST) {
            truth_count <- sum(!is.na(hp$truth))
            cat("  Restored truth data: ", truth_count, " rows with truth\n")
            fflush()
          }
        }
        
        # Merge truth from comb (ensures all truth data is included)
        if (!is.null(truth_from_comb_early) && nrow(truth_from_comb_early) > 0 && "dateID" %in% names(hp)) {
          truth_from_comb_final <- truth_from_comb_early %>%
            select(plotID, dateID, truth) %>%
            group_by(plotID, dateID) %>%
            summarise(truth = mean(truth, na.rm = TRUE), .groups = "drop")
          
          hp <- hp %>%
            left_join(truth_from_comb_final, by = c("plotID", "dateID"), suffix = c("", "_from_comb")) %>%
            mutate(truth = dplyr::coalesce(truth_from_comb, truth)) %>%
            select(-ends_with("_from_comb"))
          dbg("  After comb merge: ", sum(!is.na(hp$truth)), " rows with truth", .on = LOCAL_TEST)
        }
        
        # Ensure dates remain as Date objects after rbindlist
        if ("dateID" %in% names(hp) && "dates" %in% names(hp)) {
          # If dates are numeric but dateID exists, reconstruct dates from dateID
          if (is.numeric(hp$dates) && any(hp$dates > 10000, na.rm=TRUE)) {
            # Looks like dateID values, reconstruct from dateID
            hp$dates[!is.na(hp$dateID)] <- as.Date(paste0(hp$dateID[!is.na(hp$dateID)], "01"), "%Y%m%d")
          } else if (!inherits(hp$dates, "Date") && "dateID" %in% names(hp)) {
            # Convert dates to Date if they're not already
            hp$dates <- as.Date(paste0(hp$dateID, "01"), "%Y%m%d")
          } else if (is.numeric(hp$dates)) {
            # If numeric dates, try converting from dateID
            hp$dates <- as.Date(paste0(hp$dateID, "01"), "%Y%m%d")
          }
        }
      } else {
        # Assign fcast_period based on dateID relative to calibration end, not blindly as "hindcast"
        if ("dateID" %in% names(hp) && exists("cal_end_dateID")) {
          hp$fcast_period <- ifelse(is.na(hp$fcast_period),
            ifelse(!is.na(hp$dateID) & as.numeric(as.character(hp$dateID)) <= cal_end_dateID,
                   "calibration", "hindcast"),
            hp$fcast_period)
        } else {
          hp$fcast_period <- ifelse(is.na(hp$fcast_period), "hindcast", hp$fcast_period)
        }
      }
      
      # Final deduplication before adding to out_list
      if ("plotID" %in% names(hp) && "dateID" %in% names(hp) && nrow(hp) > 0) {
        n_before <- nrow(hp)
        hp <- ensure_dates(hp)
        if ("truth" %in% names(hp)) {
          # Prioritize: rows with medians, then rows with truth, then by fcast_period
          hp <- hp %>%
            dplyr::mutate(
              has_med = !is.na(med) & is.finite(med),
              has_truth = !is.na(truth) & is.finite(truth)
            ) %>%
            dplyr::arrange(desc(has_med), desc(has_truth), fcast_period, dateID)
          hp <- dedup_key(hp, c("plotID", "dateID"), arrange_by = "dateID") %>%
            dplyr::select(-has_med, -has_truth)
        } else {
          # Prioritize: rows with medians, then by fcast_period
          hp <- hp %>%
            dplyr::mutate(has_med = !is.na(med) & is.finite(med)) %>%
            dplyr::arrange(desc(has_med), fcast_period, dateID)
          hp <- dedup_key(hp, c("plotID", "dateID"), arrange_by = "dateID") %>%
            dplyr::select(-has_med)
        }
        dbg("  ⚠️  Final deduplication removed", n_before - nrow(hp), "duplicate rows before adding to out_list", .on = LOCAL_TEST)
      }
      
      out_list[[plotID]] <- hp
      dbg("OK plot", plotID, "rows:", nrow(hp), .on = TRUE)
      
      # Clear plot-specific data to prevent memory accumulation
      rm(list = c("hp", "comb", "cal_truth", "hind_truth"))
      # Check if full_in_plot exists before removing
      if (exists("full_in_plot")) rm(full_in_plot)
      gc(verbose = FALSE, full = TRUE)
    }
    
    # Clear model-specific data after processing all plots
    rm(list = c("full_in", "param_samples", "model_dat", "cal_quants", "time_train", "hind_dates", "trained_time_map"))
    gc(verbose = FALSE, full = TRUE)
    
    # GUARDRAIL: if no forecasts, write marker and fail
    if (!length(out_list)) {
      driver_flag <- is_driver_model(samples_path, row$model_id)
      hind_dir <- file.path(project_root, "data", "hindcasts", if (driver_flag) "driver_uncertainty" else "standard")
      dir.create(hind_dir, showWarnings=FALSE, recursive=TRUE)  # ensure exists
      empty_marker <- file.path(hind_dir, paste0("hindcasts_", row$model_id, "_EMPTY.txt"))
      writeLines(paste("No plots produced forecasts for", row$model_id), empty_marker)
      msg <- "No plots produced forecasts"
      append_log(model_id, taxon, "ERROR", msg, 0L)
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
    tax_output <- data.table::rbindlist(out_list, fill = TRUE)
    # Convert to data.frame for consistent dplyr behavior
    tax_output <- as.data.frame(tax_output)
    
    # Numeric guarantees
    must_num <- c("lo","lo_25","med","hi_75","hi")
    for (nm in must_num) if (nm %in% names(tax_output)) tax_output[[nm]] <- suppressWarnings(as.numeric(tax_output[[nm]]))

    # Label each row as Bacteria or Fungi
    tax_output <- fill_pretty_group(data.table::as.data.table(tax_output))
    tax_output <- as.data.frame(tax_output)

    # Write outputs (defensive I/O)
    driver_flag <- is_driver_model(samples_path, row$model_id)
    hind_dir <- file.path(project_root, "data", "hindcasts", if (driver_flag) "driver_uncertainty" else "standard")
    sum_dir  <- file.path(project_root, "data", "summary",   if (driver_flag) "driver_uncertainty" else "standard")
    dir.create(hind_dir, showWarnings=FALSE, recursive=TRUE)
    dir.create(sum_dir,  showWarnings=FALSE, recursive=TRUE)
    
    # Sanity check: ensure no duplicates by (model_id, plotID, dateID) before writing
    if ("model_id" %in% names(tax_output) && "plotID" %in% names(tax_output) && "dateID" %in% names(tax_output)) {
      dup_check <- tax_output %>%
        dplyr::count(.data[["model_id"]], .data[["plotID"]], .data[["dateID"]]) %>%
        dplyr::filter(n > 1)
      if (nrow(dup_check) > 0) {
        warning("Duplicates remain in tax_output: ", nrow(dup_check), " (model_id, plotID, dateID) keys")
      }
    }
    
    sites_written <- character(0)
    if ("siteID" %in% names(tax_output)) {
      for (sid in unique(tax_output$siteID)) {
        site_df <- tax_output[tax_output$siteID == sid, ]
        fp <- file.path(hind_dir, paste0("hindcasts_", row$model_id, "_", sid, "_observed.rds"))
        try(saveRDS(site_df, fp), silent=FALSE)
        if (file.exists(fp)) {
          sites_written <- c(sites_written, sid)
          cat("WROTE site file:", fp, "rows:", nrow(site_df), "\n")
          fflush()
        } else {
          warning("Failed to write site file: ", fp)
        }
      }
    }
    
    sp_out <- file.path(sum_dir, paste0("hindcasts_", row$taxon, "_observed.rds"))
    
    # Combine with existing taxon summary if it exists (multiple models can share same taxon)
    if (file.exists(sp_out)) {
      existing_data <- tryCatch({
        readRDS(sp_out)
      }, error = function(e) {
        cat("⚠️  Warning: Could not read existing taxon summary:", conditionMessage(e), "\n")
        fflush()
        NULL
      })
      
      if (!is.null(existing_data) && nrow(existing_data) > 0) {
        combined <- data.table::rbindlist(list(existing_data, tax_output), fill=TRUE, use.names=TRUE)
        
        # Remove duplicates by model_id + plotID + dateID
        if ("model_id" %in% names(combined) && "plotID" %in% names(combined) && "dateID" %in% names(combined)) {
          combined <- combined %>%
            arrange(model_id, plotID, dateID) %>%
            distinct(model_id, plotID, dateID, .keep_all = TRUE)
        }
        
        tax_output <- combined
        
        if (LOCAL_TEST) {
          cat("  Merged with existing taxon summary:", 
              nrow(existing_data), "existing rows +", 
              nrow(tax_output) - nrow(existing_data), "new rows =", 
              nrow(tax_output), "total rows\n")
          fflush()
        }
      }
    }
    
    try(saveRDS(tax_output, sp_out), silent=FALSE)
    if (file.exists(sp_out)) {
      # Result checksum (quick integrity check)
      ck <- tryCatch({ sum(!is.na(tax_output$med)); }, error=function(e) NA_integer_)
      unique_models <- if ("model_id" %in% names(tax_output)) length(unique(tax_output$model_id)) else 1
      cat("WROTE taxon summary:", sp_out, "rows:", nrow(tax_output), "models:", unique_models, "nonNA_med:", ck, "\n")
      fflush()
      
      # Generate diagnostic figures (throttled: only in LOCAL_TEST or when --figs=true)
      if (make_figs || LOCAL_TEST) {
        if (exists("generate_hindcast_diagnostics")) {
          tryCatch({
            fig_root <- if (driver_flag) {
              file.path(project_root, "figures", "hindcast_diagnostics", "observed_sites", "driver_uncertainty", model_id)
            } else {
              file.path(project_root, "figures", "hindcast_diagnostics", "observed_sites", "standard", model_id)
            }
            dir.create(fig_root, recursive=TRUE, showWarnings=FALSE)
            
            # Check for required columns
            need <- c("plotID","siteID","dates","med","lo","hi","fcast_period")
            missing_cols <- setdiff(need, names(tax_output))
            if (length(missing_cols) == 0) {
              # Filter tax_output to only current model_id (may contain multiple models from taxon summary)
              if ("model_id" %in% names(tax_output)) {
                tax_output_filtered <- tax_output %>% dplyr::filter(model_id == !!model_id)
                if (nrow(tax_output_filtered) == 0) {
                  cat("⚠️  No data found for model_id", model_id, "in tax_output (contains", length(unique(tax_output$model_id)), "models)\n")
                  fflush()
                } else {
                  tax_output <- tax_output_filtered
                }
              }
              
              # Create stable copy of truth for plotting
              plot_df <- tax_output %>%
                dplyr::mutate(truth_obs = truth) %>%
                dplyr::relocate(truth_obs, .after = truth)
              
              # Deduplicate plot_df before plotting
              if ("plotID" %in% names(plot_df) && "dateID" %in% names(plot_df) && nrow(plot_df) > 0) {
                n_before <- nrow(plot_df)
                plot_df <- plot_df %>%
                  dplyr::mutate(
                    dateID = as.numeric(as.character(dateID)),
                    has_med = !is.na(med) & is.finite(med),
                    has_truth = !is.na(truth_obs) & is.finite(truth_obs)
                  ) %>%
                  dplyr::arrange(desc(has_med), desc(has_truth), fcast_period, plotID, dateID) %>%
                  dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                  dplyr::select(-has_med, -has_truth) %>%
                  dplyr::arrange(plotID, dateID)
                
                if (LOCAL_TEST && nrow(plot_df) < n_before) {
                  cat("  ⚠️  Deduplicated plot_df: removed", n_before - nrow(plot_df), "duplicate rows before plotting\n")
                  fflush()
                }
              }
              
              # Check truth preservation at (plotID, dateID) level (accounts for deduplication)
              truth_tax <- tax_output %>%
                dplyr::filter(!is.na(truth)) %>%
                dplyr::distinct(plotID, dateID)
              truth_plot <- plot_df %>%
                dplyr::filter(!is.na(truth_obs)) %>%
                dplyr::distinct(plotID, dateID)
              if (nrow(truth_plot) != nrow(truth_tax)) {
                warning("Truth (plotID, dateID) count differs after dedup: tax_output ", nrow(truth_tax), " vs plot_df ", nrow(truth_plot), " — continuing with figure generation")
              }
              
              # Debug: log truth/pred counts
              msg <- plot_df %>%
                dplyr::summarise(n_truth = sum(!is.na(truth_obs)),
                                n_pred  = sum(!is.na(med)),
                                n_rows = dplyr::n()) %>% as.list()
              cat("DEBUG truth/pred counts:", paste(names(msg), unlist(msg), collapse="  "), "\n")
              fflush()
              
              # Guard: warn if truth is missing before plotting
              pre_truth <- sum(!is.na(plot_df$truth_obs))
              if (!"truth_obs" %in% names(plot_df) || pre_truth == 0) {
                dbg("⚠️  WARNING: Truth data missing before plotting (rows with truth:", pre_truth, ")", .on = TRUE)
              } else {
                dbg("Preparing", pre_truth, "truth observations for plotting", .on = TRUE)
              }
              
              generate_hindcast_diagnostics(plot_df, model_id, taxon, out_dir=fig_root)
              
              # Verify truth survived plotting prep
              # (Note: this check is in the diagnostic function itself)
              
              pngs <- list.files(fig_root, pattern="\\.png$", full.names=TRUE)
              cat("WROTE", length(pngs), "diagnostic figure(s) to:", fig_root, "\n")
              fflush()
            } else {
              cat("⚠️  Skipping diagnostic figures - missing columns:", paste(missing_cols, collapse=", "), "\n")
              fflush()
            }
          }, error=function(e) {
            cat("⚠️  Error generating diagnostic figures:", conditionMessage(e), "\n")
            fflush()
          })
        }
      }
      
      append_log(model_id, taxon, "SUCCESS", paste0("nonNA_med=", ck), nrow(tax_output))
      return(list(success=TRUE, model_id=model_id, taxon=taxon, rows=nrow(tax_output), sites=length(sites_written)))
    } else {
      msg <- "Failed to write taxon summary file"
      append_log(model_id, taxon, "ERROR", msg, nrow(tax_output))
      return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
    }
    
  }, error = function(e) {
    msg <- conditionMessage(e)
    append_log(model_id, taxon, "ERROR", msg, NA_integer_)
    return(list(success=FALSE, model_id=model_id, taxon=taxon, error=msg))
  })
}

## --- Helper for log flushing ---

fflush <- function() try({ flush(stdout()); flush(stderr()) }, silent=TRUE)

## --- Parallel Setup ---

# Dependencies to export to workers (exclude big objects - load them in workers instead)
deps_export <- c(
  "run_one_model", "has_all_site_outputs", "is_driver_model", "append_log",
  "as_char_cols", "first_match_col",
  "TRAIN_MIN", "UNOBSERVED", "project_root", "required_sites", "force_rerun", "make_figs",
  "generate_hindcast_diagnostics"  # CRITICAL: Export the diagnostic function to workers
)

if (!LOCAL_TEST && !sequential) {
  if (!requireNamespace("doParallel", quietly=TRUE)) stop("doParallel package required for parallel execution")
  if (!requireNamespace("foreach", quietly=TRUE)) stop("foreach package required for parallel execution")
  library(doParallel)
  library(foreach)
  
  # Get cores from CLI or use safe default
  cores_arg <- as.integer(get_flag("cores", NA))
  if (!is.na(cores_arg) && cores_arg >= 1) {
    n_cores <- cores_arg
    cat("Using", n_cores, "cores from --cores flag\n")
  } else {
    # CRITICAL: Limit to max 2 cores to reduce RAM usage
    n_cores <- min(2, max(1, parallel::detectCores(logical = FALSE) - 1))
    cat("Using", n_cores, "cores (limited to max 2 to reduce RAM usage)\n")
  }
  
  cat("Creating parallel cluster with", n_cores, "workers...\n")
  fflush()
  cl <- parallel::makeCluster(n_cores, type="PSOCK")
  cat("Cluster created, registering...\n")
  fflush()
  doParallel::registerDoParallel(cl)
  cat("Exporting functions to workers...\n")
  fflush()
  parallel::clusterExport(cl, deps_export, envir = environment())
  parallel::clusterExport(cl, c("fflush", "project_root"), envir = environment())
  cat("Initializing workers (loading data - this may take a minute)...\n")
  fflush()
  parallel::clusterEvalQ(cl, {
    setwd(project_root)  # CRITICAL: ensure workers start in project root
    
    # CRITICAL: Force single-threading in workers too to prevent nested parallelization
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
    Sys.setenv("OMP_THREAD_LIMIT" = "1")
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
    
    # Load libraries first
    suppressPackageStartupMessages({
      library(dplyr)
      library(data.table)
      library(lubridate)
      library(tidyr)
      library(here)
      # CRITICAL: Load ggplot2 for generate_hindcast_diagnostics
      library(ggplot2)
      # CRITICAL: Load truncnorm for run_hindcast.r (needed for initial conditions)
      if (requireNamespace("truncnorm", quietly = TRUE)) {
        library(truncnorm)
      }
    })
    
    # Load big data once per worker (avoid exporting/duplicating)
    # CRITICAL: These files are large and loading can take time
    cat("Worker: Loading r16...\n")
    r16 <<- readRDS(file.path(project_root, "data/clean/groupAbundances_16S_2023.rds"))
    cat("Worker: Loading rits...\n")
    rits <<- readRDS(file.path(project_root, "data/clean/groupAbundances_ITS_2023.rds"))
    cat("Worker: Loading env_data...\n")
    env_data <<- readRDS(file.path(project_root, "data/clean/all_predictor_data.rds"))
    
    # Source required files (use absolute paths to avoid issues)
    cat("Worker: Sourcing source.R...\n")
    try({
      source(file.path(project_root, "source.R"), local=TRUE)
    }, silent=TRUE)
    cat("Worker: Sourcing prepBetaRegData.r...\n")
    try({
      source(file.path(project_root, "microbialForecast/R/prepBetaRegData.r"), local=TRUE)
    }, silent=TRUE)
    cat("Worker: Sourcing run_hindcast.r...\n")
    try({
      source(file.path(project_root, "microbialForecast/R/run_hindcast.r"), local=TRUE)
    }, silent=TRUE)
    
    # CRITICAL: Load ggplot2 library needed for generate_hindcast_diagnostics
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      library(ggplot2)
    }
    cat("Worker: Initialization complete\n")
  })
  cat("Workers initialized successfully\n")
  fflush()
  
  # CRITICAL: Re-export generate_hindcast_diagnostics after workers are set up
  # This ensures the local override (defined at top of script) is available in workers
  # even if run_hindcast.r defines its own version
  parallel::clusterExport(cl, "generate_hindcast_diagnostics", envir = environment())
  
  cat("Using", n_cores, "cores for parallel processing\n")
  on.exit(try(parallel::stopCluster(cl), silent=TRUE), add=TRUE)
} else {
  cat("LOCAL_TEST mode: sequential processing\n")
}

## --- Iterate Models ---

mi <- mi[order(mi$model_name, mi$taxon, mi$max_date), , drop=FALSE]

# NOTE: Process all models - don't filter to Fusarium onwards anymore
# (This was skipping many missing models that need processing)
cat("Processing all", nrow(mi), "models\n")

if (LOCAL_TEST) {
  mi <- mi[1:min(3, nrow(mi)), , drop=FALSE]   # limit to 3 for quick test
  cat("LOCAL_TEST: Processing", nrow(mi), "models\n")
} else {
  cat("Processing", nrow(mi), "models\n")
}

if (figs_only) {
  cat("FIGURES-ONLY: Generating diagnostic figures from existing hindcast RDS (main process)\n")
  for (k in seq_len(nrow(mi))) {
    row <- mi[k, , drop = FALSE]
    model_id <- if (length(row$model_id) == 1) row$model_id else row$model_id[[1]]
    taxon <- if (length(row$taxon) == 1) row$taxon else row$taxon[[1]]
    model_root <- dirname(dirname(row$path))
    model_id_escaped <- gsub("([.^$*+?()|[\\\\]])", "\\\\\\1", model_id)
    model_id_escaped <- gsub("(\\{|\\})", "\\\\\\1", model_id_escaped)
    samples_guess <- list.files(model_root, pattern = paste0("^samples_", model_id_escaped, ".*\\.rds$"), recursive = TRUE, full.names = TRUE)
    samples_guess <- samples_guess[!grepl("_chain[0-9]", samples_guess)]
    if (!length(samples_guess)) {
      cat("SKIP no samples for", model_id, "\n")
      next
    }
    samples_path <- samples_guess[[1]]
    driver_flag <- is_driver_model(samples_path, model_id)
    hind_dir <- file.path(project_root, "data", "hindcasts", if (driver_flag) "driver_uncertainty" else "standard")
    pat <- paste0("^hindcasts_", model_id_escaped, "(_[^_]*)?_([^_]+)_observed\\.rds$")
    all_site_files <- list.files(hind_dir, pattern = pat, full.names = TRUE)
    site_ids <- sub(pat, "\\2", basename(all_site_files))
    site_files <- all_site_files[site_ids %in% required_sites]
    if (length(site_files) == 0) {
      cat("SKIP no site files for", model_id, "in required_sites (", paste(required_sites, collapse = ", "), ")\n")
      next
    }
    tax_output <- data.table::rbindlist(lapply(site_files, readRDS), fill = TRUE)
    if ("plotID" %in% names(tax_output) && "dateID" %in% names(tax_output) && nrow(tax_output) > 0) {
      if ("truth" %in% names(tax_output)) {
        tax_output <- tax_output %>%
          dplyr::mutate(
            dateID = as.numeric(as.character(dateID)),
            has_med = !is.na(med) & is.finite(med),
            has_truth = !is.na(truth) & is.finite(truth)
          ) %>%
          dplyr::arrange(dplyr::desc(has_med), dplyr::desc(has_truth), fcast_period, plotID, dateID) %>%
          dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
          dplyr::select(-has_med, -has_truth) %>%
          dplyr::arrange(plotID, dateID)
      } else {
        tax_output <- tax_output %>%
          dplyr::mutate(
            dateID = as.numeric(as.character(dateID)),
            has_med = !is.na(med) & is.finite(med)
          ) %>%
          dplyr::arrange(dplyr::desc(has_med), fcast_period, plotID, dateID) %>%
          dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
          dplyr::select(-has_med) %>%
          dplyr::arrange(plotID, dateID)
      }
    }
    fig_root <- file.path(project_root, "figures", "hindcast_diagnostics", "observed_sites",
                          if (driver_flag) "driver_uncertainty" else "standard", model_id)
    dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)
    need <- c("plotID", "siteID", "dates", "med", "lo", "hi", "fcast_period")
    missing_cols <- setdiff(need, names(tax_output))
    if (length(missing_cols) > 0) {
      cat("SKIP", model_id, "- missing columns:", paste(missing_cols, collapse = ", "), "\n")
      next
    }
    if ("model_id" %in% names(tax_output)) {
      tax_output <- tax_output %>% dplyr::filter(model_id == !!model_id)
    }
    plot_df <- tax_output %>%
      dplyr::mutate(truth_obs = .data[["truth"]])
    if ("truth" %in% names(plot_df)) {
      plot_df <- plot_df %>% dplyr::relocate(truth_obs, .after = truth)
    }
    if ("plotID" %in% names(plot_df) && "dateID" %in% names(plot_df) && nrow(plot_df) > 0) {
      plot_df <- plot_df %>%
        dplyr::mutate(
          dateID = as.numeric(as.character(dateID)),
          has_med = !is.na(med) & is.finite(med),
          has_truth = !is.na(truth_obs) & is.finite(truth_obs)
        ) %>%
        dplyr::arrange(dplyr::desc(has_med), dplyr::desc(has_truth), fcast_period, plotID, dateID) %>%
        dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
        dplyr::select(-has_med, -has_truth) %>%
        dplyr::arrange(plotID, dateID)
    }
    tryCatch({
      generate_hindcast_diagnostics(plot_df, model_id, taxon, out_dir = fig_root)
      pngs <- list.files(fig_root, pattern = "\\.png$", full.names = TRUE)
      cat("WROTE", length(pngs), "diagnostic figure(s) for", model_id, "to", fig_root, "\n")
    }, error = function(e) {
      cat("Error generating figures for", model_id, ":", conditionMessage(e), "\n")
    })
  }
  cat("FIGURES-ONLY done\n")
  quit(save = "no", status = 0)
}

# Iteration function (shared by sequential and parallel)
iter_fun <- function(row) {
  tryCatch({
    # Re-locate samples
    model_root <- dirname(dirname(row$path))
    # Escape special regex characters in model_id
    model_id_escaped <- gsub("([.^$*+?()|[\\\\]])", "\\\\\\1", row$model_id)
  model_id_escaped <- gsub("(\\{|\\})", "\\\\\\1", model_id_escaped)
  samples_guess <- list.files(model_root, pattern=paste0("^samples_", model_id_escaped, ".*\\.rds$"),
                              recursive=TRUE, full.names=TRUE)
  samples_guess <- samples_guess[!grepl("_chain[0-9]", samples_guess)]
  if (!length(samples_guess)) {
    cat("SKIP no samples for", row$model_id, "\n")
    fflush()
    return(list(success=FALSE, model_id=row$model_id, error="no samples"))
  }
  samples_path <- samples_guess[[1]]
  
  # Check if already complete (unless force rerun)
  is_complete <- has_all_site_outputs(row$model_id, samples_path, required_sites, project_root)
  if (!force_rerun && is_complete) {
    # If figures requested, still generate them even if model is complete
    if (make_figs) {
      dbg("⏭️  Model complete, generating figures only:", row$model_id, .on = TRUE)
      tryCatch({
        # Load existing output from site files and generate figures
        driver_flag <- is_driver_model(samples_path, row$model_id)
        hind_dir <- file.path(project_root, "data", "hindcasts", 
                             if (driver_flag) "driver_uncertainty" else "standard")
        
        # Load all site files for this model
        model_id_escaped <- gsub("([.^$*+?()|[\\\\]])", "\\\\\\1", row$model_id)
        model_id_escaped <- gsub("(\\{|\\})", "\\\\\\1", model_id_escaped)
        pat <- paste0("^hindcasts_", model_id_escaped, "(_[^_]*)?_([^_]+)_observed\\.rds$")
        site_files <- list.files(hind_dir, pattern=pat, full.names=TRUE)
        
        if (length(site_files) > 0) {
          # Load and combine all site files
          tax_output <- data.table::rbindlist(lapply(site_files, readRDS), fill=TRUE)
          
          # Deduplicate after combining site files
          if ("plotID" %in% names(tax_output) && "dateID" %in% names(tax_output) && nrow(tax_output) > 0) {
            n_before <- nrow(tax_output)
            if ("truth" %in% names(tax_output)) {
              tax_output <- tax_output %>%
                dplyr::mutate(
                  dateID = as.numeric(as.character(dateID)),
                  has_med = !is.na(med) & is.finite(med),
                  has_truth = !is.na(truth) & is.finite(truth)
                ) %>%
                dplyr::arrange(desc(has_med), desc(has_truth), fcast_period, plotID, dateID) %>%
                dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                dplyr::select(-has_med, -has_truth) %>%
                dplyr::arrange(plotID, dateID)
            } else {
              tax_output <- tax_output %>%
                dplyr::mutate(
                  dateID = as.numeric(as.character(dateID)),
                  has_med = !is.na(med) & is.finite(med)
                ) %>%
                dplyr::arrange(desc(has_med), fcast_period, plotID, dateID) %>%
                dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                dplyr::select(-has_med) %>%
                dplyr::arrange(plotID, dateID)
            }
            
            if (nrow(tax_output) < n_before) {
              cat("  ⚠️  Deduplicated combined site files: removed", n_before - nrow(tax_output), "duplicate rows\n")
              fflush()
            }
          }
          
          fig_root <- if (driver_flag) {
            file.path(project_root, "figures", "hindcast_diagnostics", "observed_sites", 
                     "driver_uncertainty", row$model_id)
          } else {
            file.path(project_root, "figures", "hindcast_diagnostics", "observed_sites", 
                     "standard", row$model_id)
          }
          dir.create(fig_root, recursive=TRUE, showWarnings=FALSE)
          
          need <- c("plotID","siteID","dates","med","lo","hi","fcast_period")
          missing_cols <- setdiff(need, names(tax_output))
          if (length(missing_cols) == 0) {
            # Filter tax_output to only current model_id (may contain multiple models)
            if ("model_id" %in% names(tax_output)) {
              model_id_for_filter <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
              tax_output_filtered <- tax_output %>% dplyr::filter(model_id == !!model_id_for_filter)
              if (nrow(tax_output_filtered) == 0) {
                cat("⚠️  No data found for model_id", model_id_for_filter, "in tax_output (contains", length(unique(tax_output$model_id)), "models)\n")
                fflush()
              } else {
                tax_output <- tax_output_filtered
              }
            }
            
            plot_df <- tax_output %>%
              dplyr::mutate(truth_obs = truth) %>%
              dplyr::relocate(truth_obs, .after = truth)
            
            # Deduplicate plot_df before plotting
            if ("plotID" %in% names(plot_df) && "dateID" %in% names(plot_df) && nrow(plot_df) > 0) {
              n_before <- nrow(plot_df)
              plot_df <- plot_df %>%
                dplyr::mutate(
                  dateID = as.numeric(as.character(dateID)),
                  has_med = !is.na(med) & is.finite(med),
                  has_truth = !is.na(truth_obs) & is.finite(truth_obs)
                ) %>%
                dplyr::arrange(desc(has_med), desc(has_truth), fcast_period, plotID, dateID) %>%
                dplyr::distinct(plotID, dateID, .keep_all = TRUE) %>%
                dplyr::select(-has_med, -has_truth) %>%
                dplyr::arrange(plotID, dateID)
              
              if (nrow(plot_df) < n_before) {
                cat("  ⚠️  Deduplicated plot_df: removed", n_before - nrow(plot_df), "duplicate rows before plotting\n")
                fflush()
              }
            }
            
            model_id_for_diag <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
            taxon_for_diag <- if (length(row$taxon) == 1) row$taxon else row$taxon[1]
            generate_hindcast_diagnostics(plot_df, model_id_for_diag, taxon_for_diag, out_dir=fig_root)
            pngs <- list.files(fig_root, pattern="\\.png$", full.names=TRUE)
            dbg("WROTE", length(pngs), "diagnostic figure(s) to:", fig_root, .on = TRUE)
            model_id_log <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
            taxon_log <- if (length(row$taxon) == 1) row$taxon else row$taxon[1]
            append_log(model_id_log, taxon_log, "FIGURES", paste0("generated ", length(pngs), " figures"), nrow(plot_df))
          } else {
            dbg("⚠️  Missing columns for figures:", paste(missing_cols, collapse=", "), .on = TRUE)
          }
        } else {
          dbg("⚠️  No site files found for model", row$model_id, "in", hind_dir, .on = TRUE)
        }
      }, error = function(e) {
        dbg("⚠️  Error generating figures for complete model:", conditionMessage(e), .on = TRUE)
      })
    } else {
      dbg("⏭️  SKIP (complete):", row$model_id, .on = TRUE)
    }
    model_id_val <- if (length(row$model_id) == 1) row$model_id else row$model_id[1]
    return(list(success=TRUE, model_id=model_id_val, skipped=TRUE))
  }
  if (force_rerun) {
    dbg("⚠️  FORCE RERUN:", row$model_id, "(ignoring existing files)", .on = TRUE)
  }
  
  # Process model
  result <- tryCatch(
    {
      withCallingHandlers(
        {
          dbg("🔍 DEBUG: About to call run_one_model", .on = LOCAL_TEST)
          dbg("  row$model_id length:", length(row$model_id), .on = LOCAL_TEST)
          dbg("  row$taxon length:", length(row$taxon), .on = LOCAL_TEST)
          run_one_model(row, project_root, required_sites, env_data, r16, rits, LOCAL_TEST=LOCAL_TEST, force_rerun=force_rerun, make_figs=make_figs)
        },
        error = function(e) {
          dbg("❌ ERROR HANDLER in withCallingHandlers: ", conditionMessage(e), .on = TRUE)
          if (grepl("condition has length", conditionMessage(e), ignore.case=TRUE)) {
            dbg("❌ CONDITION LENGTH ERROR DETECTED - Full traceback:", .on = TRUE)
            tryCatch({
              traceback(20)
            }, error = function(e2) {
              dbg("  (traceback failed:", conditionMessage(e2), ")", .on = TRUE)
            })
            dbg("❌ Current row structure:", .on = TRUE)
            dbg("  row$model_id:", paste(row$model_id, collapse=", "), "length:", length(row$model_id), .on = TRUE)
            dbg("  row$taxon:", paste(row$taxon, collapse=", "), "length:", length(row$taxon), .on = TRUE)
            dbg("  row$model_name:", paste(row$model_name, collapse=", "), "length:", length(row$model_name), .on = TRUE)
          }
        }
      )
    },
    error = function(e) {
      dbg("❌ ERROR in outer tryCatch for model", if(length(row$model_id)==1) row$model_id else row$model_id[1], ":", conditionMessage(e), .on = TRUE)
      if (grepl("condition has length", conditionMessage(e), ignore.case=TRUE)) {
        dbg("❌ CONDITION LENGTH ERROR - Attempting traceback:", .on = TRUE)
        tryCatch({
          traceback(20)
        }, error = function(e2) {
          dbg("  (traceback failed:", conditionMessage(e2), ")", .on = TRUE)
        })
      }
      list(success=FALSE, model_id=if(length(row$model_id)==1) row$model_id else row$model_id[1], error=conditionMessage(e))
    }
  )
  return(result)
  }, error = function(e) {
    # Outer error handler for iter_fun itself
    dbg("❌ CRITICAL ERROR in iter_fun:", conditionMessage(e), .on = TRUE)
    if (grepl("condition has length", conditionMessage(e), ignore.case=TRUE)) {
      dbg("❌ CONDITION LENGTH ERROR in iter_fun - Full traceback:", .on = TRUE)
      traceback(10)
      if (exists("row")) {
        dbg("❌ Current row structure:", .on = TRUE)
        dbg("  row$model_id:", paste(row$model_id, collapse=", "), "length:", length(row$model_id), .on = TRUE)
        dbg("  row$taxon:", paste(row$taxon, collapse=", "), "length:", length(row$taxon), .on = TRUE)
        dbg("  row$model_name:", paste(row$model_name, collapse=", "), "length:", length(row$model_name), .on = TRUE)
      }
    }
    list(success=FALSE, model_id=if(exists("row") && length(row$model_id)==1) row$model_id else "unknown", error=conditionMessage(e))
  })
}

# Process models sequentially or in parallel
if (LOCAL_TEST || sequential) {
  cat("Using SEQUENTIAL processing (", if (sequential) "--sequential flag" else "LOCAL_TEST", ")\n")
  fflush()
  results <- lapply(seq_len(nrow(mi)), function(k) {
    row <- mi[k, , drop=FALSE]
    cat("\n=== RUN", k, "of", nrow(mi), ":", row$model_id, "taxon:", row$taxon, "model:", row$model_name, "===\n")
    fflush()
    iter_fun(row)
  })
} else {
  # Add interrupt handler to ensure cluster shuts down on ^C
  dbg("Starting parallel processing of", nrow(mi), "models...", .on = TRUE)
  
  # CRITICAL: Test worker access to variables before starting foreach
  cat("Testing worker access to variables...\n")
  test_result <- tryCatch({
    parallel::clusterEvalQ(cl, {
      list(
        has_mi = exists("mi", envir = .GlobalEnv),
        has_iter_fun = exists("iter_fun", envir = .GlobalEnv),
        has_project_root = exists("project_root", envir = .GlobalEnv),
        worker_pid = Sys.getpid()
      )
    })
  }, error = function(e) {
    dbg("ERROR testing worker access:", e$message, .on = TRUE)
    return(NULL)
  })
  
  if (!is.null(test_result)) {
    dbg("Worker test results:", .on = TRUE)
    for (i in seq_along(test_result)) {
      dbg("  Worker", i, ":", paste(names(test_result[[i]]), "=", unlist(test_result[[i]]), collapse=", "), .on = TRUE)
    }
  }
  
  results <- tryCatch({
    withCallingHandlers({
      foreach(k = seq_len(nrow(mi)),
              .packages = c("dplyr","data.table","lubridate","tidyr","ggplot2"),
              .export = c("mi", "iter_fun", "required_sites", "project_root", 
                         "force_rerun", "make_figs", "LOCAL_TEST"),
              .combine = function(...) unlist(list(...), recursive=FALSE),
              .errorhandling = "pass") %dopar% {
        # CRITICAL: Add file-based worker logging (cat() output may not appear in main log)
        worker_id <- Sys.getpid()
        log_file <- file.path(project_root, "data", "hindcasts", paste0("worker_", worker_id, ".log"))
        log_msg <- function(...) {
          msg <- paste(..., "\n", sep="")
          cat(msg, file=log_file, append=TRUE)
          cat(msg)  # Also try stdout
          fflush()
        }
        
        log_msg("Worker", worker_id, ": Starting iteration k =", k, "of", nrow(mi))
        
        # Test access to mi
        if (!exists("mi")) {
          log_msg("Worker", worker_id, ": ERROR - mi not found!")
          return(list(success=FALSE, model_id=paste0("unknown_", k), error="mi not found in worker"))
        }
        
        row <- mi[k, , drop=FALSE]
        log_msg("Worker", worker_id, ": === RUN", k, "of", nrow(mi), ":", row$model_id, "taxon:", row$taxon, "model:", row$model_name)
        
        # Test access to iter_fun
        if (!exists("iter_fun")) {
          log_msg("Worker", worker_id, ": ERROR - iter_fun not found!")
          return(list(success=FALSE, model_id=row$model_id, error="iter_fun not found in worker"))
        }
        
        log_msg("Worker", worker_id, ": Calling iter_fun for", row$model_id)
        result <- tryCatch({
          iter_fun(row)
        }, error = function(e) {
          log_msg("Worker", worker_id, ": ERROR in iter_fun:", e$message)
          list(success=FALSE, model_id=row$model_id, error=e$message)
        })
        
        log_msg("Worker", worker_id, ": Completed iteration", k)
        result
      }
    }, interrupt = function(e) {
      try(parallel::stopCluster(cl), silent = TRUE)
      stop("Interrupted", call. = FALSE)
    })
  }, error = function(e) {
    try(parallel::stopCluster(cl), silent = TRUE)
    stop(conditionMessage(e))
  })
}

## --- Summary ---

# Ensure results is a list
if (!is.list(results)) results <- as.list(results)

# Safely extract success/skipped status
get_status <- function(x) {
  if (is.list(x) && "success" %in% names(x)) {
    if (isTRUE(x$skipped)) return("skipped")
    if (isTRUE(x$success)) return("success")
    return("failed")
  }
  return("unknown")
}

statuses <- vapply(results, get_status, character(1))
ok <- sum(statuses == "success")
skipped <- sum(statuses == "skipped")
failed <- sum(statuses == "failed")

cat("\n=== SUMMARY ===\n")
cat("Total models:", length(results), "\n")
cat("Success:", ok, "\n")
cat("Skipped (already complete):", skipped, "\n")
cat("Failed:", failed, "\n")

if (failed > 0) {
  cat("\nFailed models:\n")
  for (i in seq_along(results)) {
    r <- results[[i]]
    if (is.list(r) && !isTRUE(r$success) && !isTRUE(r$skipped)) {
      model_id <- if ("model_id" %in% names(r)) r$model_id else paste0("model_", i)
      error_msg <- if ("error" %in% names(r)) r$error else "unknown error"
      cat("  ", model_id, ":", error_msg, "\n")
    }
  }
}

cat("\nRun log saved to:", runlog_path, "\n")
cat("=== DONE ===\n")

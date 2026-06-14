# Calculate scoring metrics for hindcasts — memory-safe refactor

tryCatch({
  mem.maxVSize(Inf)
  cat("Memory limit increased to unlimited\n")
}, error = function(e) {
  cat("Note: Could not increase memory limit:", e$message, "\n")
})

# Limit thread usage to prevent system overload
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1") 
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")

library(here)
source(here("source.R"))
# statsFunctions.r functions (add_scoring_metrics, etc.) loaded via microbialForecast package
# robust version now merged into microbialForecast::add_scoring_metrics
library(data.table)

cat("Loading hindcast data...\n")

# Load hindcasts via the package loader, which reads and unions the per-model
# parquet files (hindcasts_<model>.parquet) written by step 07.
hindcast_data <- as.data.table(load_hindcasts())
cat(sprintf("✓ Loaded hindcasts: %d rows, %d cols\n",
            nrow(hindcast_data), ncol(hindcast_data)))

# Assign pretty_group if missing or any NA (fill in missing values)
if (!"pretty_group" %in% names(hindcast_data) || any(is.na(hindcast_data$pretty_group))) {
  cat("Assigning pretty_group based on species taxonomy...\n")
  fg_names <- microbialForecast:::keep_fg_names
  rank_spec_names <- microbialForecast:::rank_spec_names
  all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
  all_fungi <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])
  
  # Initialize pretty_group column if missing
  if (!"pretty_group" %in% names(hindcast_data)) {
    hindcast_data[, pretty_group := NA_character_]
  }
  
  # Only assign for rows where pretty_group is NA
  na_rows <- is.na(hindcast_data$pretty_group)
  if (any(na_rows)) {
    cat("Filling", sum(na_rows), "rows with missing pretty_group\n")
    
    # Assign pretty_group based on species - use data.table syntax for efficiency
    # First, assign for functional groups
    fg_species <- hindcast_data$species[na_rows] %in% fg_names
    if (any(fg_species)) {
      fg_indices <- which(na_rows)[fg_species]
      fg_kingdoms <- microbialForecast:::assign_fg_kingdoms(microbialForecast:::assign_fg_categories(hindcast_data$species[fg_indices]))
      hindcast_data[fg_indices, pretty_group := ifelse(fg_kingdoms == "16S", "Bacteria", "Fungi")]
      na_rows[fg_indices] <- FALSE  # Update na_rows to exclude these
    }
    
    # Then assign for taxonomic bacteria
    bac_species <- hindcast_data$species[na_rows] %in% all_bacteria
    if (any(bac_species)) {
      hindcast_data[which(na_rows)[bac_species], pretty_group := "Bacteria"]
      na_rows[which(na_rows)[bac_species]] <- FALSE
    }
    
    # Then assign for taxonomic fungi
    fun_species <- hindcast_data$species[na_rows] %in% all_fungi
    if (any(fun_species)) {
      hindcast_data[which(na_rows)[fun_species], pretty_group := "Fungi"]
      na_rows[which(na_rows)[fun_species]] <- FALSE
    }
    
    # Fallback: use rank_name pattern if species-based assignment failed
    still_na <- na_rows
    if (any(still_na)) {
      hindcast_data[which(still_na), pretty_group := ifelse(grepl("_bac$|16S", rank_name), "Bacteria",
                                                             ifelse(grepl("_fun$|ITS", rank_name), "Fungi", NA_character_))]
    }
  }
  
  cat("Final pretty_group distribution:", table(hindcast_data$pretty_group, useNA = "ifany"), "\n")
}

# Add rank_only column if missing (first part of rank_name before underscore)
if (!"rank_only" %in% names(hindcast_data) && "rank_name" %in% names(hindcast_data)) {
  hindcast_data[, rank_only := sapply(strsplit(rank_name, "_", fixed = TRUE), function(x) {
    if (length(x) > 0 && x[1] != "") {
      # Handle special case: "functional" groups don't have _bac/_fun suffix
      if (x[1] %in% c("functional", "diversity")) {
        return(x[1])
      }
      # For taxonomic ranks, return first part (e.g., "phylum_bac" -> "phylum")
      return(x[1])
    } else {
      return(NA_character_)
    }
  })]
  # Set as ordered factor with standard levels
  rank_levels <- c("genus", "family", "order", "class", "phylum", "functional", "diversity")
  hindcast_data[, rank_only := factor(rank_only, levels = rank_levels, ordered = TRUE)]
  cat("Added rank_only column from rank_name\n")
}

# Handle nanoparquet column renaming: when the parquet has both uppercase (Mean, SD)
# and lowercase (mean, sd) columns, nanoparquet appends _1 to the lowercase duplicates.
# Standardize to lowercase (mean, sd) which is what run_hindcast.r produces.
col_renames <- c("mean_1" = "mean", "sd_1" = "sd")
for (old_name in names(col_renames)) {
  new_name <- col_renames[[old_name]]
  if (old_name %in% names(hindcast_data) && !new_name %in% names(hindcast_data)) {
    setnames(hindcast_data, old_name, new_name)
    cat(sprintf("Renamed '%s' -> '%s' (nanoparquet case-collision artifact)\n", old_name, new_name))
  }
}

# Keep only needed columns
needed_cols <- c("model_id","fcast_period","pretty_group","model_name",
                 "rank_name","rank_only","species","site_prediction","siteID","plotID",
                 "truth","mean","med","sd","date_num","plot_start_date")
drop_cols <- setdiff(names(hindcast_data), needed_cols)
if (length(drop_cols)) hindcast_data[, (drop_cols) := NULL]
cat(sprintf("Trimmed to %d cols\n", ncol(hindcast_data)))

# Convert truth column from character to numeric
hindcast_data[, truth := as.numeric(truth)]

# Set indexes for fast subsetting/grouping
setindex(hindcast_data, fcast_period)
setindex(hindcast_data, model_id)
setindex(hindcast_data, species)
setindex(hindcast_data, site_prediction)
setindex(hindcast_data, siteID)

hindcast_only <- hindcast_data[fcast_period == "hindcast" & !is.na(truth) & !is.na(med)]
calibration_only <- hindcast_data[fcast_period == "calibration" & !is.na(truth) & !is.na(med)]

# Filter calibration to exclude first dates if plot_start_date exists
if ("plot_start_date" %chin% names(calibration_only)) {
  calibration_only_not_first <- calibration_only[is.na(plot_start_date) | (date_num > plot_start_date)]
  cat("Filtered calibration by plot_start_date\n")
} else {
  cat("Warning: plot_start_date not found, using all calibration data\n")
  calibration_only_not_first <- copy(calibration_only)
}

cat(sprintf("Hindcast: %d, Calibration: %d, Calibration (not first): %d\n",
            nrow(hindcast_only), nrow(calibration_only), nrow(calibration_only_not_first)))

empty_metric_list <- function() list(
  RMSE=NA_real_, BIAS=NA_real_, MAE=NA_real_,
  CRPS=NA_real_, RSQ=NA_real_, RSQ.1=NA_real_,
  RMSE.iqr=NA_real_, RMSE.norm=NA_real_,
  CRPS_truncated=NA_real_, residual_variance=NA_real_,
  predictive_variance=NA_real_, total_PL=NA_real_
)

safe_metrics <- function(obs, mu, med, sdv) {
  if (length(obs) <= 1L) return(empty_metric_list())
  tryCatch(
    as.list(add_scoring_metrics(
      observed = obs, mean_predicted = mu, median_predicted = med, sd_predicted = sdv
    )),
    error = function(e) empty_metric_list()
  )
}

# Process in chunks by model_id
chunk_keys <- unique(hindcast_data$model_id)
cat("Processing", length(chunk_keys), "chunks by model_id...\n")

scoring_metrics_list <- list()
calibration_metrics_list <- list()
calibration_metrics_site_list <- list()
scoring_metrics_site_list <- list()

ctr <- 0L
for (m in chunk_keys) {
  ctr <- ctr + 1L
  if ((ctr %% 10L) == 0L) cat(".. chunk", ctr, "of", length(chunk_keys), "\n")

  hc <- hindcast_only[model_id == m]
  cb <- calibration_only_not_first[model_id == m]

  if (nrow(hc)) {
    by_cols <- c("model_id", "fcast_period", "pretty_group", "model_name", "rank_name", "species", "site_prediction")
    if ("rank_only" %in% names(hc)) by_cols <- c(by_cols, "rank_only")
    sm <- hc[!is.na(site_prediction),
      safe_metrics(truth, mean, med, sd),
      by = by_cols]
    scoring_metrics_list[[length(scoring_metrics_list)+1L]] <- sm

    model_name_val <- if (nrow(hc) > 0 && "model_name" %in% names(hc)) {
      unique(hc$model_name)[1]
    } else NA_character_
    
    sms <- hc[!is.na(site_prediction),
      {
        if (.N <= 1L) empty_metric_list()
        else safe_metrics(truth, mean, med, sd)
      },
      by = .(model_id, siteID, site_prediction)
    ]
    if (nrow(sms) > 0 && nrow(hc) > 0) {
      if (!is.na(model_name_val)) sms[, model_name := model_name_val]
      
      if ("pretty_group" %in% names(hc)) {
        pg_val <- unique(hc$pretty_group[!is.na(hc$pretty_group)])
        if (length(pg_val) > 0) sms[, pretty_group := pg_val[1]]
      }
      if ("rank_name" %in% names(hc)) {
        rn_val <- unique(hc$rank_name[!is.na(hc$rank_name)])
        if (length(rn_val) > 0) {
          sms[, rank_name := rn_val[1]]
          rn_first <- strsplit(rn_val[1], "_", fixed = TRUE)[[1]]
          sms[, rank_only := if (length(rn_first) > 0) rn_first[1] else NA_character_]
          sms[, pretty_name := rank_only]
        }
      }
      if ("species" %in% names(hc)) {
        sp_val <- unique(hc$species[!is.na(hc$species)])
        if (length(sp_val) > 0) sms[, species := sp_val[1]]
      }
      sms[, fcast_period := "hindcast"]
    }
    scoring_metrics_site_list[[length(scoring_metrics_site_list)+1L]] <- sms
  }

  if (nrow(cb)) {
    by_cal <- c("model_id", "fcast_period", "pretty_group", "model_name", "rank_name", "species")
    if ("rank_only" %in% names(cb)) by_cal <- c(by_cal, "rank_only")
    cm <- cb[, safe_metrics(truth, mean, med, sd), by = by_cal]
    calibration_metrics_list[[length(calibration_metrics_list)+1L]] <- cm

    by_cal_site <- c("model_id", "fcast_period", "pretty_group", "model_name", "rank_name", "species", "siteID")
    if ("rank_only" %in% names(cb)) by_cal_site <- c(by_cal_site, "rank_only")
    cms <- cb[, safe_metrics(truth, mean, med, sd), by = by_cal_site]
    calibration_metrics_site_list[[length(calibration_metrics_site_list)+1L]] <- cms
  }

  rm(hc, cb); gc(FALSE)
}

scoring_metrics <- rbindlist(scoring_metrics_list, use.names=TRUE, fill=TRUE)
scoring_metrics[, mean_crps_sample := CRPS]
if (!"pretty_name" %in% names(scoring_metrics)) {
  if ("rank_only" %in% names(scoring_metrics)) {
    scoring_metrics[, pretty_name := rank_only]
  } else if ("rank_name" %in% names(scoring_metrics)) {
    scoring_metrics[, pretty_name := sapply(strsplit(rank_name, "_", fixed = TRUE), function(x) if (length(x) > 0) x[1] else NA_character_)]
  }
}

calibration_metrics <- rbindlist(calibration_metrics_list, use.names=TRUE, fill=TRUE)
calibration_metrics_site <- rbindlist(calibration_metrics_site_list, use.names=TRUE, fill=TRUE)
scoring_metrics_site <- rbindlist(scoring_metrics_site_list, use.names=TRUE, fill=TRUE)
if (!"pretty_name" %in% names(calibration_metrics) && "rank_only" %in% names(calibration_metrics)) {
  calibration_metrics[, pretty_name := rank_only]
}
if (!"pretty_name" %in% names(calibration_metrics_site) && "rank_only" %in% names(calibration_metrics_site)) {
  calibration_metrics_site[, pretty_name := rank_only]
}
if ("rank_only" %in% names(scoring_metrics_site)) {
  scoring_metrics_site[, pretty_name := rank_only]
} else if ("rank_name" %in% names(scoring_metrics_site)) {
  scoring_metrics_site[, pretty_name := sapply(strsplit(rank_name, "_", fixed = TRUE), function(x) if (length(x) > 0) x[1] else NA_character_)]
}

# Pivot helpers in data.table
to_long <- function(dt, id_cols, measure_cols) {
  if (nrow(dt) == 0 || length(measure_cols) == 0) {
    result <- data.table(matrix(ncol = length(id_cols) + 2, nrow = 0))
    names(result) <- c(id_cols, "metric", "score")
    return(result)
  }
  id_cols <- intersect(id_cols, names(dt))
  melt(dt, id.vars = id_cols, measure.vars = measure_cols,
       variable.name = "metric", value.name = "score", variable.factor = FALSE)
}

metric_cols <- c("RMSE","BIAS","MAE","CRPS","RSQ","RSQ.1","RMSE.iqr",
                 "RMSE.norm","CRPS_truncated","residual_variance",
                 "predictive_variance","total_PL")

id_cols_long <- c("model_id", "fcast_period", "pretty_group", "model_name", "rank_name", "species", "site_prediction")
if ("rank_only" %in% names(scoring_metrics)) id_cols_long <- c(id_cols_long, "rank_only")
if ("pretty_name" %in% names(scoring_metrics)) id_cols_long <- c(id_cols_long, "pretty_name")
scoring_metrics_long <- to_long(
  scoring_metrics,
  id_cols = id_cols_long,
  measure_cols = intersect(metric_cols, names(scoring_metrics))
)

calibration_metrics_long <- to_long(
  calibration_metrics,
  id_cols = c("model_id","fcast_period","pretty_group","model_name",
              "rank_name","species"),
  measure_cols = intersect(metric_cols, names(calibration_metrics))
)

calibration_metrics_site_long <- to_long(
  calibration_metrics_site,
  id_cols = c("model_id","fcast_period","pretty_group","model_name",
              "rank_name","species","siteID"),
  measure_cols = intersect(metric_cols, names(calibration_metrics_site))
)[is.finite(score)]

id_site_long <- c("model_id", "fcast_period", "pretty_group", "model_name", "rank_name", "species", "siteID", "site_prediction")
if ("pretty_name" %in% names(scoring_metrics_site)) id_site_long <- c(id_site_long, "pretty_name")
scoring_metrics_site_long <- to_long(
  scoring_metrics_site,
  id_cols = id_site_long,
  measure_cols = intersect(metric_cols, names(scoring_metrics_site))
)[is.finite(score)]

# ---- CV calculations (env_cycl and driver uncertainty models) without dplyr ----
# Updated to handle new model types and naming patterns
truth_vals <- calibration_only[model_name == "env_cycl" | grepl("driver_uncertainty", model_id) | grepl("driver_uncertainty", model_name)]

cat("Starting CV calculations...\n")
calc_cv_dt <- function(x) {
  # same as your calc_cv(truth); if your calc_cv differs, call it here
  mu <- mean(x, na.rm=TRUE)
  s  <- sd(x, na.rm=TRUE)
  if (is.na(mu) || mu == 0) return(NA_real_)
  s / mu
}

cv_tax_per_plot <- truth_vals[, .(per_plot_cv = calc_cv_dt(truth)),
                              by = .(pretty_group, rank_name, species, siteID, plotID)]
cv_tax_per_plot_site <- cv_tax_per_plot[, .(mean_per_plot_site_cv = mean(per_plot_cv, na.rm=TRUE)),
                                        by = .(pretty_group, rank_name, species, siteID)]
cv_tax_per_plot_species <- cv_tax_per_plot[, .(mean_per_plot_cv = mean(per_plot_cv, na.rm=TRUE)),
                                           by = .(pretty_group, rank_name, species)]
cv_tax_per_site <- truth_vals[, .(per_site_cv = calc_cv_dt(truth)),
                              by = .(pretty_group, rank_name, species, siteID)]
cv_tax_per_site_species <- cv_tax_per_site[, .(mean_per_site_cv = mean(per_site_cv, na.rm=TRUE)),
                                           by = .(pretty_group, rank_name, species)]
cv_tax_overall <- truth_vals[, .(overall_cv = calc_cv_dt(truth)),
                             by = .(pretty_group, rank_name, species)]

# Remove duplicates before merging to prevent cartesian explosion
cv_tax_per_site_species <- unique(cv_tax_per_site_species, by = c("pretty_group","rank_name","species"))
cv_tax_per_site <- unique(cv_tax_per_site, by = c("pretty_group","rank_name","species","siteID"))

cv_tax <- merge(cv_tax_overall, cv_tax_per_site_species,
                by = c("pretty_group","rank_name","species"), all.x = TRUE)
# Don't merge cv_tax_per_site here as it has site-level data that would cause cartesian explosion
# cv_tax_per_site will be used separately for site-level analysis

# Merge CV with scoring metrics
scoring_metrics_cv <- merge(scoring_metrics_long, cv_tax, all.x = TRUE)
# Check if scoring_metrics_site_long has the required columns
if (all(c("pretty_group","rank_name","species","siteID") %in% names(scoring_metrics_site_long))) {
  scoring_metrics_cv_site <- merge(scoring_metrics_site_long, cv_tax_per_site,
                                   by = c("pretty_group","rank_name","species","siteID"), all.x = TRUE)
} else {
  cat("Warning: scoring_metrics_site_long missing required columns for CV merge\n")
  scoring_metrics_cv_site <- scoring_metrics_site_long
}

# Long CV (avoid pivot_longer)
# Only melt columns that actually exist
cv_cols <- c("overall_cv","per_site_cv","mean_per_site_cv")
existing_cv_cols <- intersect(cv_cols, names(scoring_metrics_cv))
cat("Available CV columns:", paste(existing_cv_cols, collapse=", "), "\n")

if (length(existing_cv_cols) > 0) {
  scoring_metrics_cv_long <- melt(
    scoring_metrics_cv,
    id.vars = setdiff(names(scoring_metrics_cv), existing_cv_cols),
    measure.vars = existing_cv_cols,
    variable.name = "cv_type", value.name = "cv"
  )
} else {
  cat("Warning: No CV columns found, creating empty long format\n")
  scoring_metrics_cv_long <- scoring_metrics_cv
  scoring_metrics_cv_long[, cv_type := NA_character_]
  scoring_metrics_cv_long[, cv := NA_real_]
}

# CV scaling per group without scale()
cat("Starting CV scaling...\n")
# Filter once to avoid Inf/NaN
scv <- scoring_metrics_cv_long[is.finite(cv) & is.finite(score)]

# Compute mean/sd by grouping, then join back to compute standardized values; avoids matrix allocs
grp_keys <- c("pretty_group","rank_name","cv_type","metric")
stats_cv <- scv[, .(cv_mean = mean(cv), cv_sd = sd(cv), sc_mean = mean(score), sc_sd = sd(score)),
               by = grp_keys][cv_sd > 0 & sc_sd > 0]
cv_metric_scaled <- stats_cv[scv, on = grp_keys][
  , `:=`(CV_scale = (cv - cv_mean)/cv_sd,
         metric_scale = (score - sc_mean)/sc_sd)][
  , c("cv_mean","cv_sd","sc_mean","sc_sd") := NULL][]

cat("CV scaling complete\n")
gc(FALSE)

# Skill scores (only need the two metrics; keep it small)
available_types <- unique(scoring_metrics_long$site_prediction)
cat("Site prediction types:", paste(available_types, collapse=", "), "\n")

skill_score_species <- data.table()
skill_score_species_RMSE <- data.table()
skill_score_rank <- data.table()

if (length(available_types) > 1L && "New time (observed site)" %chin% available_types) {
  # Helper to wide-ize a single metric with dcast
  make_skill <- function(metric_name) {
    wide <- scoring_metrics_long[metric == metric_name,
       dcast(.SD,
             model_id + fcast_period + pretty_group + model_name + rank_name + species ~ site_prediction,
             value.var = "score", fill = NA_real_)
    ]
    if ("New time (observed site)" %chin% names(wide)) {
      if ("New time x site (modeled effect)" %chin% names(wide))
        wide[, skill_score := 1 - `New time x site (modeled effect)`/`New time (observed site)`]
      else
        wide[, skill_score := NA_real_]  # Create column even if no modeled effect data
      
      if ("New time x site (random effect)" %chin% names(wide))
        wide[, skill_score_random := 1 - `New time x site (random effect)`/`New time (observed site)`]
      else
        wide[, skill_score_random := NA_real_]  # Create column even if no random effect data
    } else {
      wide[, skill_score := NA_real_]
      wide[, skill_score_random := NA_real_]
    }
    wide
  }

  skill_score_species <- make_skill("CRPS_truncated")
  skill_score_species_RMSE <- make_skill("RMSE.norm")

  if ("skill_score" %in% names(skill_score_species)) {
    skill_score_rank <- skill_score_species[, .(
      mean_skill_score = mean(skill_score, na.rm = TRUE)
    ), by = .(model_id, pretty_group, rank_name)]
  } else {
    cat("Warning: skill_score column not found, creating empty skill_score_rank\n")
    skill_score_rank <- data.table()
  }
} else {
  cat("Insufficient site_prediction types for skill score calculation\n")
}

# Convergence lists (unchanged)
unconverged <- if (file.exists(here("data/summary/unconverged_taxa_list.rds")))
                 readRDS(here("data/summary/unconverged_taxa_list.rds")) else NULL
converged <- if (file.exists(here("data/summary/weak_converged_taxa_list.rds")))
               readRDS(here("data/summary/weak_converged_taxa_list.rds")) else NULL
converged_strict <- if (file.exists(here("data/summary/converged_taxa_list.rds")))
                      readRDS(here("data/summary/converged_taxa_list.rds")) else NULL

# Save outputs (as before)
to_save <- list(
  cv_metric_scaled = cv_metric_scaled,
  scoring_metrics_cv = scoring_metrics_cv,
  scoring_metrics_cv_site = scoring_metrics_cv_site,
  scoring_metrics_cv_long = scoring_metrics_cv_long,
  calibration_truth_vals = truth_vals,

  scoring_metrics = scoring_metrics,
  scoring_metrics_long = scoring_metrics_long,
  scoring_metrics_site = scoring_metrics_site,
  scoring_metrics_site_long = scoring_metrics_site_long,

  calibration_metrics = calibration_metrics,
  calibration_metrics_long = calibration_metrics_long,
  calibration_metrics_site = calibration_metrics_site,
  calibration_metrics_site_long = calibration_metrics_site_long,

  skill_score_species = skill_score_species,
  skill_score_species_RMSE = skill_score_species_RMSE,
  skill_score_taxon = skill_score_species,
  skill_score_taxon_RMSE = skill_score_species_RMSE,
  skill_score_rank = skill_score_rank,
  unconverged_list = unconverged,
  converged_list = converged,
  converged_strict_list = converged_strict
)

saveRDS(to_save, here("data/summary/scoring_metrics_plsr2.rds"))
cat("Saved: scoring_metrics_plsr2.rds with", length(to_save), "components\n")
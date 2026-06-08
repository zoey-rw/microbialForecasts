#' Load hindcast data from parquet with pretty_group reconstruction
#'
#' Reads the hindcasts via duckdb (memory-efficient) and reconstructs
#' pretty_group from rank_name/species if it contains NAs. By default it reads
#' and unions the per-model-type files written by step 07
#' (hindcasts_env_cycl.parquet, hindcasts_cycl_only.parquet,
#' hindcasts_env_cov.parquet). When `model_names` is given, only the matching
#' per-model files are read. Falls back to the legacy combined
#' all_hindcasts_plsr2.parquet if the per-model files are absent.
#'
#' @param parquet_path Optional explicit path to a single parquet file. When
#'   supplied it overrides the per-model lookup (back-compatible).
#' @param fcast_period Filter to "hindcast", "calibration", or NULL for all.
#' @param observed_only If TRUE, filter to observed sites only (new_site == FALSE).
#' @param model_names Character vector of model_name values to keep, or NULL for all.
#' @param columns Character vector of columns to select, or NULL for all.
#' @param parquet_dir Directory holding the per-model parquet files. Defaults to
#'   here("data/summary/parquet").
#' @return data.table with pretty_group populated.
#' @export
load_hindcasts <- function(parquet_path = NULL,
                           fcast_period = NULL,
                           observed_only = FALSE,
                           model_names = NULL,
                           columns = NULL,
                           parquet_dir = NULL) {
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("duckdb package required for memory-efficient parquet reading")
  }

  all_models <- c("env_cycl", "cycl_only", "env_cov")
  if (is.null(parquet_dir)) {
    parquet_dir <- here::here("data/summary/parquet")
  }

  # Resolve which file(s) to read.
  if (!is.null(parquet_path)) {
    files <- parquet_path                       # explicit single file
  } else {
    wanted <- if (is.null(model_names)) all_models else intersect(model_names, all_models)
    per_model <- file.path(parquet_dir, sprintf("hindcasts_%s.parquet", wanted))
    per_model <- per_model[file.exists(per_model)]
    if (length(per_model) > 0L) {
      files <- per_model                        # per-model split (preferred)
    } else {
      combined <- file.path(parquet_dir, "all_hindcasts_plsr2.parquet")
      if (file.exists(combined)) {
        files <- combined                       # legacy combined fallback
      } else {
        stop("No hindcast parquet found in ", parquet_dir,
             "\nExpected hindcasts_<model>.parquet. Run ",
             "analysis/model_analysis/07_tidyHindcasts.r or download_data.R.")
      }
    }
  }

  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("Parquet file(s) not found: ", paste(missing, collapse = ", "))
  }

  # read_parquet over one or many files (union_by_name handles identical schemas).
  files_sql <- paste(sprintf("'%s'", files), collapse = ", ")
  read_expr <- if (length(files) > 1L) {
    sprintf("read_parquet([%s], union_by_name=true)", files_sql)
  } else {
    sprintf("read_parquet(%s)", files_sql)
  }

  # Build SQL query with optional filters
  col_clause <- if (is.null(columns)) "*" else paste(columns, collapse = ", ")
  where_clauses <- character(0)
  if (!is.null(fcast_period)) {
    where_clauses <- c(where_clauses,
                       sprintf("fcast_period = '%s'", fcast_period))
  }
  if (isTRUE(observed_only)) {
    where_clauses <- c(where_clauses,
                       "(new_site = FALSE OR new_site IS NULL)")
  }
  if (!is.null(model_names)) {
    quoted <- paste(sprintf("'%s'", model_names), collapse = ", ")
    where_clauses <- c(where_clauses,
                       sprintf("model_name IN (%s)", quoted))
  }

  where_sql <- if (length(where_clauses) > 0) {
    paste("WHERE", paste(where_clauses, collapse = " AND "))
  } else ""

  sql <- sprintf("SELECT %s FROM %s %s", col_clause, read_expr, where_sql)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  dt <- data.table::as.data.table(DBI::dbGetQuery(con, sql))

  # Reconstruct pretty_group if it has NAs
  dt <- fill_pretty_group(dt)

  dt
}


#' Fill NA values in pretty_group from rank_name and species
#'
#' Uses rank_name suffixes (_bac/_fun) for taxonomic groups and
#' assign_fg_kingdoms() for functional groups.
#'
#' @param dt data.table with rank_name, species, and optionally pretty_group columns
#' @return dt with pretty_group filled
#' @export
fill_pretty_group <- function(dt) {
  if (!"pretty_group" %in% names(dt)) {
    dt[, pretty_group := NA_character_]
  }

  n_na_before <- sum(is.na(dt$pretty_group))
  if (n_na_before == 0L) return(dt)

  # Step 1: rank_name suffix
  if ("rank_name" %in% names(dt)) {
    na_idx <- is.na(dt$pretty_group)
    rn <- dt$rank_name[na_idx]
    dt[na_idx & grepl("_bac$", rank_name), pretty_group := "Bacteria"]
    dt[na_idx & grepl("_fun$", rank_name), pretty_group := "Fungi"]
  }

  # Step 2: functional group species mapping
  if ("species" %in% names(dt)) {
    na_idx <- which(is.na(dt$pretty_group))
    if (length(na_idx) > 0L) {
      fg_names <- tryCatch(microbialForecast:::keep_fg_names, error = function(e) NULL)
      if (!is.null(fg_names)) {
        sp <- dt$species[na_idx]
        is_fg <- sp %in% fg_names
        if (any(is_fg)) {
          fg_kingdoms <- assign_fg_kingdoms(assign_fg_categories(sp[is_fg]))
          dt$pretty_group[na_idx[is_fg]] <- ifelse(
            fg_kingdoms == "16S", "Bacteria",
            ifelse(fg_kingdoms == "ITS", "Fungi", NA_character_)
          )
        }
      }
    }
  }

  n_na_after <- sum(is.na(dt$pretty_group))
  if (n_na_after > 0L && n_na_after < n_na_before) {
    message(sprintf("fill_pretty_group: filled %d/%d NAs (%d remain)",
                    n_na_before - n_na_after, n_na_before, n_na_after))
  }

  dt
}

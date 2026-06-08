# Custom NEON-store readers. These patch bugs in the neonstore package's file
# discovery / reading so that sensor products (e.g. soil moisture DP1.00094.001,
# soil temperature DP1.00041.001) with ragged column schemas can be read from a
# local NEON store. Migrated from functions/custom_neonstore.r for reproducibility.
#
# neon_read_custom() is the public entry point; the remaining functions are
# internal helpers it calls.


filter_duplicates_custom <- function (meta) {
  require(neonstore)
  meta$timestamp <- as.POSIXct(meta$timestamp, format = "%Y%m%dT%H%M%OS")
  meta_b <- meta[order(meta$timestamp, decreasing = TRUE),
                 ]
  meta_b$id <- paste(meta_b$product, meta_b$site, meta_b$table, meta_b$month,
                     sep = "-")
  out <- neonstore:::take_first_match(meta_b, "id")
  if (dim(out)[[1]] < dim(meta)[[1]])
    message(paste("Some raw files were detected with updated timestamps.\n",
                  "Using only most updated file to avoid duplicates."))
  out
}


meta_filter_custom <- function (meta, product = NA, table = NA, site = NA, start_date = NA,
                                end_date = NA, type = NA, timestamp = NA, ext = NA)
{
  if (!is.na(table)) {
    meta <- meta[grepl(table, meta$table), ]
  }
  if (!all(is.na(product))) {
    meta <- meta[meta$product %in% product, ]
  }
  if (!all(is.na(site))) {
    meta <- meta[meta$site %in% site, ]
  }
  if (!is.na(start_date)) {
    start_date <- as.Date(start_date)
    month <- as.Date(paste(meta$month,"-01",sep=""))
    keep <- month >= start_date
    keep[is.na(keep)] <- TRUE
    meta <- meta[keep, ]
  }
  if (!is.na(end_date)) {
    end_date <- as.Date(end_date)
    month <- as.Date(paste(meta$month,"-01",sep=""))
    keep <- month <= end_date
    keep[is.na(keep)] <- TRUE
    meta <- meta[keep, ]
  }
  if (!is.na(timestamp)) {
    meta <- meta[meta$timestamp < as.POSIXct(timestamp),
                 ]
  }
  if (!is.na(type)) {
    meta <- switch(type, basic = meta[meta$type != "expanded",
                                      ], expanded = meta[meta$type != "basic", ], meta)
  }
  if (!all(is.na(ext))) {
    meta <- meta[meta$ext %in% ext, ]
  }
  tibble::as_tibble(meta)
}


neon_index_custom <- function (product = NA, table = NA, site = NA, start_date = NA,
                        end_date = NA, type = NA, ext = NA, timestamp = NA, hash = NULL,
                        dir = neonstore::neon_dir())
{
  files <- list.files(dir)
  meta <- neonstore:::filename_parser(files)
  if (is.null(meta))
    return(NULL)
  meta$path <- file.path(dir, meta$path)
  meta$timestamp <- as.POSIXct(meta$timestamp, format = "%Y%m%dT%H%M%OS")
  meta <- meta_filter_custom(meta, product = product, table = table,
                      site = site, start_date = start_date, end_date = end_date,
                      type = type, timestamp = timestamp, ext = ext)
  meta$hash <- neonstore:::file_hash(meta$path, hash = hash)
  tibble::as_tibble(meta)
}

#' @title neon_read_custom
#' @description Read NEON data (including ragged-schema sensor products) from a local
#'   neonstore directory, de-duplicating files whose timestamps have been updated. A
#'   custom replacement for \code{neonstore::neon_read} that handles sensor metadata and
#'   inconsistent column schemas across months.
#' @param table,product,site,start_date,end_date,ext,timestamp Filters passed to the
#'   file index (see \code{neonstore::neon_index}).
#' @param dir Local neonstore directory. Defaults to \code{neonstore::neon_dir()}.
#' @param files Optional explicit vector of file paths to read instead of indexing \code{dir}.
#' @param sensor_metadata If TRUE, extract domain/site/position/date columns from sensor filenames.
#' @param altrep Passed to \code{vroom::vroom}.
#' @param ... Additional arguments passed to \code{vroom::vroom}.
#' @return A tibble of the combined NEON data, or NULL if no files are found.
#' @export
neon_read_custom <- function(table = NA,
                      product = NA,
                      site = NA,
                      start_date = NA,
                      end_date = NA,
                      ext = NA,
                      timestamp = NA,
                      dir = neonstore::neon_dir(),
                      files = NULL,
                      sensor_metadata = TRUE,
                      altrep = FALSE,
                      ...){

  if(is.null(files)){
    meta <- neon_index_custom(product = product,
                       table = table,
                       site = site,
                       start_date = start_date,
                       end_date = end_date,
                       ext = ext,
                       hash = NULL,
                       dir = dir)

    if(is.null(meta)) return(NULL)
    if(dim(meta)[[1]] == 0 )  return(NULL)

    ## If timestamp has changed but other metadata is the same, we only want the newer version
    meta <- filter_duplicates_custom(meta)
    files <- meta$path
  }

  if(length(files) == 0){
    if(is.null(table)) table <- "unspecified tables"
    warning(paste("no files found for", table, "in", dir, "\n",
                  "perhaps you need to download them first?"))
    return(NULL)
  }

  ## Handle the case of needing to add columns extracted from filenames
  if(neonstore:::is_sensor_data(files) && sensor_metadata){
    neon_read_sensor(meta, altrep = altrep, ...)
    ## Otherwise we can just read in:
  } else {
    read_csvs(files,  altrep = altrep, ...)
  }

}


neon_read_sensor <- function(meta, altrep = FALSE, ..., .id = "path") {
  suppressMessages({
    id <- unique(meta[[.id]])
    groups <-
      lapply(id,
             function(x){
               paths <- meta$path[meta[[.id]] == x]
               out <- read_csvs(paths, altrep = altrep, ...)
               out[.id] <- x
               out
             })
  })
  suppressWarnings({
    df <- ragged_bind(groups)
  })

  filename_meta <- neonstore:::neon_filename_parser(df$path)
  df$domainID <- filename_meta$DOM
  df$siteID <- filename_meta$SITE
  df$horizontalPosition <- filename_meta$HOR
  df$verticalPosition <- filename_meta$VER
  df$publicationDate <- as.POSIXct(filename_meta$GENTIME, format = "%Y%m%dT%H%M%OS")

  df
}

read_csvs <- function(files, altrep = FALSE, ...){
  ## vroom can read in a list of files, but only if columns are consistent
  tryCatch(vroom::vroom(files, altrep = altrep,  ...),
           error = function(e) vroom_ragged(files, altrep = altrep, ...),
           finally = NULL)
}


vroom_ragged <- function(files, altrep = FALSE, ...){

  ## We read the 1st line of every file to determine schema
  suppressMessages(
    schema <- lapply(files, vroom::vroom, n_max = 1, altrep = altrep, ...)
  )


  ## Now, we read in tables in groups of matching schema,
  ## filling in additional columns as in bind_rows.

  col_schemas <- lapply(schema, colnames)
  u_schemas <- unique(col_schemas)
  tbl_list <- vector("list", length=length(u_schemas))

  all_cols <- unique(unlist(u_schemas))

  i <- 1
  for(s in u_schemas){

    ## select tables that have matching schemas
    index <- vapply(col_schemas, identical, logical(1L), s)
    col_types <- vroom::spec(schema[index][[1]])

    ## Read in all those tables
    tbl <- vroom::vroom(files[index], col_types = col_types)

    ## append any columns missing from all_cols set
    missing <- all_cols[ !(all_cols %in% colnames(tbl)) ]
    tbl[ missing ] <- NA
    tbl_list[[i]] <- tbl
    i <- i+1

  }
  do.call(rbind, tbl_list)

}

## simpler case
ragged_bind <- function(x){

  col_schemas <- lapply(x, colnames)
  u_schemas <- unique(col_schemas)
  all_cols <- unique(unlist(u_schemas))
  for(i in seq_along(x)){
    ## append any columns missing from all_cols set
    missing <- all_cols[ !(all_cols %in% colnames(x[[i]])) ]
    x[[i]][ missing ] <- NA
  }
  do.call(rbind, x)

}

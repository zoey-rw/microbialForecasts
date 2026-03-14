#'  @title 			assign_pheno_date_vectorized_optimized
#'  @description Highly optimized vectorized version of assign_pheno_date for large datasets
#' @export
assign_pheno_date_vectorized_optimized <- function(dates, pheno_df = pheno_categories_long) {
  if (length(dates) == 0) return(character(0))
  
  # Convert dates to Date objects
  dates_ymd <- ymd(dates)
  
  # Remove NA dates
  valid_dates <- !is.na(dates_ymd)
  if (sum(valid_dates) == 0) return(rep(NA_character_, length(dates)))
  
  valid_dates_ymd <- dates_ymd[valid_dates]
  
  # Use data.table for maximum speed
  if (require(data.table, quietly = TRUE)) {
    # Convert pheno_df to data.table
    pheno_dt <- as.data.table(pheno_df)
    
    # Create result vector
    result <- character(length(valid_dates_ymd))
    
    # Process dates in chunks to balance memory and speed
    chunk_size <- 1000
    n_chunks <- ceiling(length(valid_dates_ymd) / chunk_size)
    
    for (chunk in 1:n_chunks) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, length(valid_dates_ymd))
      chunk_dates <- valid_dates_ymd[start_idx:end_idx]
      
      # Vectorized processing for this chunk
      for (i in seq_along(chunk_dates)) {
        date_i <- chunk_dates[i]
        keep_index <- date_i %within% pheno_dt$value
        out_categories <- pheno_dt[keep_index, name]
        
        if (length(out_categories) > 0) {
          # Clean up category names
          out_categories[out_categories == "dormancy_interval1"] <- "dormancy_interval"
          out_categories[out_categories == "dormancy_interval2"] <- "dormancy_interval"
          out_categories <- gsub("_interval", "", out_categories)
          
          # Handle multiple categories
          if (length(unique(out_categories)) > 1) {
            cat_table <- table(out_categories)
            cat_table <- cat_table / sum(cat_table)
            if (any(cat_table > 0.6)) {
              out_categories <- names(cat_table[cat_table > 0.6])
            } else {
              out_categories <- names(cat_table[cat_table > 0.3])
              out_categories <- paste(out_categories, collapse = "_")
            }
          } else {
            out_categories <- unique(out_categories)
          }
          
          result[start_idx + i - 1] <- out_categories
        } else {
          result[start_idx + i - 1] <- NA_character_
        }
      }
      
      # Progress indicator
      if (chunk %% 10 == 0) {
        cat("Processed phenology chunk", chunk, "of", n_chunks, "\n")
      }
    }
  } else {
    # Fallback to base R approach
    result <- character(length(valid_dates_ymd))
    
    for (i in seq_along(valid_dates_ymd)) {
      date_i <- valid_dates_ymd[i]
      keep_index <- date_i %within% pheno_df$value
      out_categories <- pheno_df[keep_index, ]$name
      
      if (length(out_categories) > 0) {
        # Clean up category names
        out_categories[out_categories == "dormancy_interval1"] <- "dormancy_interval"
        out_categories[out_categories == "dormancy_interval2"] <- "dormancy_interval"
        out_categories <- gsub("_interval", "", out_categories)
        
        # Handle multiple categories
        if (length(unique(out_categories)) > 1) {
          cat_table <- table(out_categories)
          cat_table <- cat_table / sum(cat_table)
          if (any(cat_table > 0.6)) {
            out_categories <- names(cat_table[cat_table > 0.6])
          } else {
            out_categories <- names(cat_table[cat_table > 0.3])
            out_categories <- paste(out_categories, collapse = "_")
          }
        } else {
          out_categories <- unique(out_categories)
        }
        
        result[i] <- out_categories
      } else {
        result[i] <- NA_character_
      }
    }
  }
  
  # Create final result vector with NAs for invalid dates
  final_result <- rep(NA_character_, length(dates))
  final_result[valid_dates] <- result
  
  return(final_result)
}

#'  @title 			assign_pheno_site_date_optimized
#'  @description Optimized version of assign_pheno_site_date for large datasets
#' @export
assign_pheno_site_date_optimized <- function(site = "HARV", date_ymd = "2015-01-01",
                                            pheno_df = pheno_categories_long) {
  
  # Handle vectorized inputs
  if (length(site) > 1 || length(date_ymd) > 1) {
    # Create data frame for vectorized processing
    site_date_df <- data.frame(siteID = site, dates = date_ymd, stringsAsFactors = FALSE)
    return(assign_pheno_category_optimized(site_date_df, pheno_df))
  }
  
  # Single site/date processing
  site_year_pheno <- pheno_df %>% filter(ID %in% site)
  
  keep_index <- date_ymd %within% site_year_pheno$value
  out_categories <- site_year_pheno[keep_index, ]$name
  
  if (length(out_categories) > 0) {
    out_categories[out_categories == "dormancy_interval1"] <- "dormancy_interval"
    out_categories[out_categories == "dormancy_interval2"] <- "dormancy_interval"
    out_categories <- gsub("_interval", "", out_categories)
    out_categories <- max(out_categories, na.rm = TRUE)
  } else {
    out_categories <- NA_character_
  }
  
  return(out_categories)
}

#'  @title 			assign_pheno_category_optimized
#'  @description Optimized version of assign_pheno_category for large datasets
#' @export
assign_pheno_category_optimized <- function(site_date_df,
                                           pheno_df = pheno_categories_long) {
  if ("siteID" %in% colnames(site_date_df)) {
    site <- site_date_df[["siteID"]]
    site_year_pheno <- pheno_df %>% filter(ID %in% site)
  } else {
    site_year_pheno <- pheno_df
  }
  
  date_ymd <- ymd(site_date_df[["dates"]])
  
  keep_index <- date_ymd %within% site_year_pheno$value
  out_categories <- site_year_pheno[keep_index, ]$name
  
  if (length(out_categories) > 0) {
    out_categories[out_categories == "dormancy_interval1"] <- "dormancy_interval"
    out_categories[out_categories == "dormancy_interval2"] <- "dormancy_interval"
    out_categories <- gsub("_interval", "", out_categories)
    
    # If date is at transition of two categories, use priority-based selection
    # Priority: peak > greenup > greendown > dormancy
    # This prevents dormancy from being incorrectly assigned to summer dates
    unique_cats <- unique(out_categories)
    if (length(unique_cats) > 1) {
      priority_order <- c("peak", "greenup", "greendown", "dormancy")
      # Find the highest priority category that matches
      selected_cat <- NULL
      for(priority_cat in priority_order) {
        if(priority_cat %in% unique_cats) {
          selected_cat <- priority_cat
          break
        }
      }
      # If no priority match, fall back to frequency-based selection
      if(is.null(selected_cat)) {
        cat_table <- table(out_categories)
        cat_table <- cat_table / sum(cat_table)
        if (any(cat_table > 0.6)) {
          selected_cat <- names(cat_table[cat_table > 0.7])[1]
        } else {
          selected_cat <- names(cat_table[cat_table > 0.3])[1]
        }
      }
      out_categories <- selected_cat
    } else {
      out_categories <- unique_cats[1]
    }
  } else {
    out_categories <- NA_character_
  }
  
  return(out_categories)
}

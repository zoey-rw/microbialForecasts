#'  @title 			assign_pheno_category
#'  @description Assign phenological category based on date and siteID dataframe
#' @export
#'
#'
# pheno_categories <- readRDS(here("data/clean/modis_greenup.rds"))
# # pheno_categories_long = pheno_categories %>% select(-c(3:10)) %>% pivot_longer(cols = 3:7) %>%# 	filter(ID %in% unique(cycl_only_est$siteID)

assign_pheno_category <- function(site_date_df,
																	pheno_df=pheno_categories_long) {
	if ("siteID" %in% colnames(site_date_df)) {
		site = site_date_df[["siteID"]]
		site_year_pheno = pheno_df %>% filter(ID %in% site)
	} else {
		site_year_pheno = pheno_df
	}
	date_ymd = ymd(site_date_df[["dates"]])

	keep_index <- date_ymd %within% site_year_pheno$value
	out_categories = site_year_pheno[keep_index,]$name
	#print(out_categories)

	out_categories[out_categories=="dormancy_interval1"] <- "dormancy_interval"
	out_categories[out_categories=="dormancy_interval2"] <- "dormancy_interval"
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
			cat_table = table(out_categories)
			cat_table = cat_table/sum(cat_table)
			if (any(cat_table > .6)) {
				selected_cat = names(cat_table[cat_table > .7])[1]
			} else {
				selected_cat = names(cat_table[cat_table > .3])[1]
			}
		}
		out_categories <- selected_cat
	} else {
		out_categories <- unique_cats[1]
	}
	return(out_categories)
}


#'  @title 			assign_pheno_site_date
#'  @description Assign phenological category based on date and siteID
#' @export
assign_pheno_site_date <- function(site = "HARV", date_ymd = "2015-01-01",
																	pheno_df=pheno_categories_long) {

	site_year_pheno = pheno_df %>% filter(ID %in% site)

	keep_index <- date_ymd %within% site_year_pheno$value
	out_categories = site_year_pheno[keep_index,]$name
	#print(out_categories)

	out_categories[out_categories=="dormancy_interval1"] <- "dormancy_interval"
	out_categories[out_categories=="dormancy_interval2"] <- "dormancy_interval"
	out_categories <- gsub("_interval", "", out_categories)
	# If date is at transition of two categories, keep the latter
	# Which is more likely to represent the sampling date (which can be anytime throughout the month,
	# we've just coded it as the first of the month)
		out_categories <- max(out_categories, na.rm=T)
	return(out_categories)
}

#'  @title 			assign_pheno_date
#'  @description Assign phenological category based only on date  dataframe
#' @export
assign_pheno_date  <- function(date,
															 pheno_df=pheno_categories_long) {
	# Handle single date (original behavior)
	if (length(date) == 1) {
		date_ymd = ymd(date)
		keep_index <- date_ymd %within% pheno_df$value
		out_categories = pheno_df[keep_index,]$name
		#print(out_categories)
		out_categories[out_categories=="dormancy_interval1"] <- "dormancy_interval"
		out_categories[out_categories=="dormancy_interval2"] <- "dormancy_interval"
		out_categories <- gsub("_interval", "", out_categories)
		# If date is at transition of two categories, keep the latter
		# Which is more likely to represent the sampling date (which can be anytime throughout the month,
		# we've just coded it as the first of the month)
		# If 90% of sites are in one interval, return that interval
		if (length(unique(out_categories)) > 1) {
			cat_table = table(out_categories)
			cat_table = cat_table/sum(cat_table)
			if (any(cat_table > .6)) {
				out_categories = names(cat_table[cat_table > .6])
			} else {
				out_categories = names(cat_table[cat_table > .3])
				out_categories <- paste(out_categories, collapse = "_")
			}
		} else {
			out_categories <- unique(out_categories)
		}
		return(out_categories)
	} else {
		# Handle multiple dates using optimized vectorized approach
		return(assign_pheno_date_vectorized(date, pheno_df))
	}
}

#'  @title 			assign_pheno_date_vectorized
#'  @description Optimized vectorized version of assign_pheno_date for multiple dates
#' @export
assign_pheno_date_vectorized <- function(dates, pheno_df = pheno_categories_long) {
  if (length(dates) == 0) return(character(0))
  
  # Convert dates to Date objects
  dates_ymd <- ymd(dates)
  
  # Remove NA dates
  valid_dates <- !is.na(dates_ymd)
  if (sum(valid_dates) == 0) return(rep(NA_character_, length(dates)))
  
  valid_dates_ymd <- dates_ymd[valid_dates]
  
  # Vectorized interval matching using data.table for speed
  if (require(data.table, quietly = TRUE)) {
    # Convert pheno_df to data.table for faster operations
    pheno_dt <- as.data.table(pheno_df)
    
    # Create a lookup table for each valid date
    result <- character(length(valid_dates_ymd))
    
    for (i in seq_along(valid_dates_ymd)) {
      date_i <- valid_dates_ymd[i]
      keep_index <- date_i %within% pheno_dt$value
      out_categories <- pheno_dt[keep_index, name]
      
      if (length(out_categories) > 0) {
        # Clean up category names
        out_categories[out_categories == "dormancy_interval1"] <- "dormancy_interval"
        out_categories[out_categories == "dormancy_interval2"] <- "dormancy_interval"
        out_categories <- gsub("_interval", "", out_categories)
        
        # Handle multiple categories - use priority-based selection
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
              selected_cat <- names(cat_table[cat_table > 0.6])[1]
            } else {
              selected_cat <- names(cat_table[cat_table > 0.3])[1]
            }
          }
          out_categories <- selected_cat
        } else {
          out_categories <- unique_cats[1]
        }
        
        result[i] <- out_categories
      } else {
        result[i] <- NA_character_
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
        
        # Handle multiple categories - use priority-based selection
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
              selected_cat <- names(cat_table[cat_table > 0.6])[1]
            } else {
              selected_cat <- names(cat_table[cat_table > 0.3])[1]
            }
          }
          out_categories <- selected_cat
        } else {
          out_categories <- unique_cats[1]
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

#'  @title 			assign_pheno_date_batch
#'  @description Batch processing version for very large datasets
#' @export
assign_pheno_date_batch <- function(dates, pheno_df = pheno_categories_long, batch_size = 1000) {
  if (length(dates) == 0) return(character(0))
  
  # Process in batches to avoid memory issues
  n_dates <- length(dates)
  n_batches <- ceiling(n_dates / batch_size)
  result <- character(n_dates)
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, n_dates)
    batch_dates <- dates[start_idx:end_idx]
    
    # Process this batch
    batch_result <- assign_pheno_date_vectorized(batch_dates, pheno_df)
    result[start_idx:end_idx] <- batch_result
    
    # Progress indicator
    if (batch %% 10 == 0) {
      cat("Processed batch", batch, "of", n_batches, "\n")
    }
  }
  
  return(result)
}


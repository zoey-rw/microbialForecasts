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
	# If date is at transition of two categories, keep the latter
	# Which is more likely to represent the sampling date (which can be anytime throughout the month,
	# we've just coded it as the first of the month)

	# If 90% of sites are in one interval, return that interval
	if (length(unique(out_categories)) > 1) {
		cat_table = table(out_categories)
		cat_table = cat_table/sum(cat_table)
		if (any(cat_table > .6)) {
			out_categories = names(cat_table[cat_table > .7])
		} else {
			out_categories = names(cat_table[cat_table > .3])
			out_categories <- paste(out_categories, collapse = "_")
		}
	} else {
		out_categories <- unique(out_categories)
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
}


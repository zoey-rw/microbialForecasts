# Combine hindcasts from all workflows
# This script must be run before other analysis scripts

source("../../source.R")

# Load required packages for memory efficiency
if (!require(arrow, quietly = TRUE)) {
  install.packages("arrow", repos = "https://cran.rstudio.com/")
  library(arrow)
}

# Check if Parquet file exists, otherwise use RDS
parquet_file <- here("data/summary/parquet/all_hindcasts_raw_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_raw_plsr2.rds")

if (file.exists(parquet_file)) {
  cat("Using Parquet file for memory efficiency...\n")
  all_rank_hindcasts <- read_parquet(parquet_file)
  cat(sprintf("Loaded %d rows from Parquet file\n", nrow(all_rank_hindcasts)))
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  cat("WARNING: This may cause memory issues!\n")
  all_rank_hindcasts <- readRDS(rds_file)
  cat(sprintf("Loaded %d rows from RDS file\n", nrow(all_rank_hindcasts)))
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}
hindcast_data <- all_rank_hindcasts %>%
	mutate(dates=fixDate(dateID),
				 fcast_period = ifelse(dates <= "2018-01-01" & (newsite == FALSE|is.na(newsite)),
				 											"calibration", "hindcast"))

# hindcasts$group <- ifelse(grepl("_bac", hindcasts$rank_name, fixed = T), "16S", "ITS")
# hindcasts$rank = hindcasts$rank_name
tax = hindcast_data %>% filter(grepl("_bac|_fun", rank_name)) %>%
	mutate(group = ifelse(grepl("16S|bac", rank_name), "16S", "ITS"),
				 category="taxonomic_rank",
				 fcast_type= "Taxonomic")

tax$rank_only= sapply(str_split(tax$rank_name, "_",  n = 2), `[`, 1)
fg = hindcast_data %>% filter(!grepl("_bac|_fun", rank_name)) %>%
	mutate(
	category = assign_fg_categories(species),
	group = assign_fg_kingdoms(category),
	fcast_type="Functional",
	rank_only="functional")
hindcast_data = rbind(tax, fg) %>%
		mutate(pretty_group = ifelse(grepl("16S|bac", group), "Bacteria", "Fungi"),
						taxon_name = species,
						taxon=species)



	hindcast_data$rank_only <- ordered(hindcast_data$rank_only, levels = c("genus",
																																			 "family",
																																			 "order",
																																			 "class",
																																			 "phylum", "functional"))



# Add prettier data values
hindcast_data$newsite <- ifelse(hindcast_data$new_site, "New site", "Observed site")
hindcast_data$pretty_name <- ordered(hindcast_data$rank_only, levels = c("genus",
																																					 "family",
																																					 "order",
																																					 "class",
																																					 "phylum", "functional"))

hindcast_data = hindcast_data %>%
	filter(!grepl("other", taxon)) %>%
	mutate(site_prediction = ifelse(predicted_site_effect==TRUE & new_site==TRUE,
																	"New time x site (modeled effect)",
																	ifelse(predicted_site_effect==FALSE & new_site==TRUE,
																				 "New time x site (random effect)",
																				 "New time (observed site)"))) %>% mutate(timepoint=date_num)

hindcast_only = hindcast_data %>% filter(fcast_period=="hindcast")
calibration_only = hindcast_data %>% filter(fcast_period=="calibration") %>%
	select(-c(predicted_site_effect, site_prediction)) %>% distinct(.keep_all = T)


# For the calibration, the first observed date per plot is always wonky due to model structure, so we leave it out of our scoring metrics
# To determine if an observation is the first date or not, we unfortunately have to pull in some data ("plot_start" and "site_start") used to fit the model
# pulled model dat from refit script with j=1 or j=6
model.dat_bac = quick_get_rank_df(k = 1)
model.dat_fun = quick_get_rank_df(k = 6)

# Add missing dates back in
dates_key <- model.dat_bac$truth.plot.long %>% select(date_num, dateID) %>% unique()
calibration_only$date_num <- dates_key[match(calibration_only$dateID, dates_key$dateID),]$date_num

bac_site_start_dates <- model.dat_bac$site_start %>% stack %>% mutate(pretty_group = "Bacteria")
bac_plot_start_dates <- model.dat_bac$plot_start %>% stack %>% mutate(pretty_group = "Bacteria")
colnames(bac_site_start_dates) <- c("site_start_date", "siteID", "pretty_group")
colnames(bac_plot_start_dates) <- c("plot_start_date", "plotID", "pretty_group")
bac_plot_start_dates$plotID <- as.character(bac_plot_start_dates$plotID)

fun_site_start_dates <- model.dat_fun$site_start %>% stack %>% mutate(pretty_group = "Fungi")
fun_plot_start_dates <- model.dat_fun$plot_start %>% stack %>% mutate(pretty_group = "Fungi")
colnames(fun_site_start_dates) <- c("site_start_date", "siteID", "pretty_group")
colnames(fun_plot_start_dates) <- c("plot_start_date", "plotID", "pretty_group")
fun_plot_start_dates$plotID <- as.character(fun_plot_start_dates$plotID)

# Merge this information back into the calibration dataset
calibration_bac <- merge(calibration_only %>% filter(pretty_group=="Bacteria"),
												 bac_site_start_dates, by=c("siteID","pretty_group"))
calibration_bac <- merge(calibration_bac,
												 bac_plot_start_dates, by=c("plotID","pretty_group"))
calibration_fun <- merge(calibration_only %>% filter(pretty_group=="Fungi"),
												 fun_site_start_dates, by=c("siteID","pretty_group"))
calibration_fun <- merge(calibration_fun,
												 fun_plot_start_dates, by=c("plotID","pretty_group"))

calibration_only = rbind.data.frame(calibration_bac, calibration_fun)

# Clean up memory after processing large datasets
cat("Cleaning up memory...\n")
rm(all_rank_hindcasts, bac_site_start_dates, bac_plot_start_dates, fun_site_start_dates, fun_plot_start_dates,
   calibration_bac, calibration_fun)
gc()

cat("✓ Hindcast processing complete!\n")
calibration_only$is_site_start_date = ifelse(calibration_only$date_num==calibration_only$site_start_date, T, F)
calibration_only$is_plot_start_date = ifelse(calibration_only$date_num==calibration_only$plot_start_date, T, F)
calibration_only$is_any_start_date = ifelse(calibration_only$is_plot_start_date|calibration_only$is_site_start_date, T, F)
calibration_only = calibration_only %>% filter(date_num >= plot_start_date)

hindcast_data_out <- rbindlist(list(hindcast_only, calibration_only), fill = T) %>%
	arrange(model_name, fcast_type, pretty_group, taxon, plotID, date_num)

hindcast_data_out <- hindcast_data_out %>% select(-c(is_any_start_date,is_site_start_date,site_start_date))

# Save as both RDS and Parquet for compatibility
saveRDS(hindcast_data_out, here("data/summary/all_hindcasts_plsr2.rds"))
cat("✓ Saved RDS file: all_hindcasts_plsr2.rds\n")

# Save as Parquet for memory efficiency
if (require(arrow, quietly = TRUE)) {
  arrow::write_parquet(hindcast_data_out, here("data/summary/parquet/all_hindcasts_plsr2.parquet"))
  cat("✓ Saved Parquet file: all_hindcasts_plsr2.parquet\n")
} else {
  cat("⚠️  Arrow package not available, skipping Parquet save\n")
}




# Filter to remove taxa that did not converge
# unconverged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
# unconverged <- unconverged[unconverged$median_gbr > 3,]

# hindcasts_filt <- hindcasts %>%
# 	mutate(taxon_model_rank = paste(taxon_name, model_name, rank_name)) %>%
# 	filter(!taxon_model_rank %in% unconverged$taxon_model_rank)
# saveRDS(hindcasts_filt,
# 				here("data/summary/hindcast_group_filt.rds"))

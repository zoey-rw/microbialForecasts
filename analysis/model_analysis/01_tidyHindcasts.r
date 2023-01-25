# Combine hindcasts from all workflows
# This script must be run before other analysis scripts

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Combine forecast types
div_hindcast_data <- readRDS(here("data/summary/hindcast_div.rds")) %>%
	mutate(fcast_type = "Diversity", rank = paste0("diversity_", group),
				 taxon = rank)

fg_hindcast_data <- readRDS(here("data/summary/beta_hindcast_fg_2015-11_2018-01.rds")) %>%
	mutate(fcast_type = "Functional group")

# tax_hindcast_data <- readRDS("./data/summary/hindcast_tax_test.rds") %>%
# 	mutate(fcast_type = "Taxonomic",
# 				 fcast_period = ifelse(dates <= "2017-01-01", "calibration", "hindcast")) %>%
# 	filter(species != "other")# replace

# Include unconverged models for hindcasts. Unfortunate.
tax_hindcast_data <- readRDS(here("data/summary/hindcast_single_tax_all.rds"))
# Nvm.
#tax_hindcast_data <- readRDS(here("data/summary/hindcast_single_tax.rds"))

hindcast_data <- rbindlist(list(div_hindcast_data,
																fg_hindcast_data,
																tax_hindcast_data), fill = T)

hindcast_data$rank_only <- sapply(str_split(hindcast_data$rank, "_",  n = 2), `[`, 1)
hindcast_data$rank_only <- ordered(hindcast_data$rank_only, levels = c("genus",
																																			 "family",
																																			 "order",
																																			 "class",
																																			 "phylum", "functional", "diversity"))
hindcast_data$rank <- ordered(hindcast_data$rank, levels = c("genus_bac","genus_fun",
																														 "family_bac","family_fun",
																														 "order_bac", "order_fun",
																														 "class_bac", "class_fun",
																														 "phylum_bac","phylum_fun",
																														 "functional_group", "diversity_16S", "diversity_ITS"))


# Add prettier data values
hindcast_data$newsite <- ifelse(hindcast_data$new_site, "New site", "Observed site")
hindcast_data$pretty_group <- ifelse(hindcast_data$group=="16S", "Bacteria", "Fungi")
hindcast_data$pretty_name <- recode(hindcast_data$rank, !!!microbialForecast:::pretty_rank_names)
hindcast_data$pretty_name <- ordered(hindcast_data$pretty_name, levels = c("Genus",
																																					 "Family",
																																					 "Order",
																																					 "Class",
																																					 "Phylum", "Functional group", "Diversity"))

hindcast_data <- hindcast_data %>% mutate(truth=as.numeric(truth),
					 mean=as.numeric(mean),
					 Mean=as.numeric(Mean))

hindcast_data = hindcast_data %>%
	filter(!grepl("other", taxon)) %>%
	mutate(site_prediction = ifelse(predicted_site_effect==TRUE & new_site==TRUE,
																	"New time x site (modeled effect)",
																	ifelse(predicted_site_effect==FALSE & new_site==TRUE,
																				 "New time x site (random effect)",
																				 "New time (observed site)")))


hindcast_data_simple <- hindcast_data %>% mutate(mean = ifelse(!is.na(Mean), Mean, mean),
																								 sd = ifelse(!is.na(SD), SD, sd),
																								 hi = ifelse(!is.na(`97.5%`), `97.5%`, hi),
																								 med = ifelse(!is.na(`50%`), `50%`, med),
																								 lo = ifelse(!is.na(`2.5%`), `2.5%`, lo),
																								 lo_25 = ifelse(!is.na(`25%`), `25%`, lo_25),
																								 hi_75 = ifelse(!is.na(`75%`), `75%`, lo_75)) %>%
	select(-c(lo_75, `2.5%`, `25%`, `50%`, `75%`, `97.5%`, Mean, SD, Shannon_orig, Shannon_scale_site, rank_name, species_num))

hindcast_only = hindcast_data_simple %>% filter(fcast_period=="hindcast")
calibration_only = hindcast_data_simple %>% filter(fcast_period=="calibration") %>%
	select(-c(predicted_site_effect, site_prediction)) %>% distinct(.keep_all = T)

# For the calibration, the first observed date per plot is always wonky due to model structure, so we leave it out of our scoring metrics
# To determine if an observation is the first date or not, we unfortunately have to pull in some data ("plot_start" and "site_start") used to fit the model
# pulled model dat from refit script with j=1 or j=6

model.dat_bac = quick_get_rank_df(k = 1)
model.dat_fun = quick_get_rank_df(k = 6)

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
calibration_only$is_site_start_date = ifelse(calibration_only$date_num==calibration_only$site_start_date, T, F)
calibration_only$is_plot_start_date = ifelse(calibration_only$date_num==calibration_only$plot_start_date, T, F)
calibration_only$is_any_start_date = ifelse(calibration_only$is_plot_start_date|calibration_only$is_site_start_date, T, F)
calibration_only = calibration_only %>% filter(timepoint >= site_start_date)

hindcast_data_out <- rbindlist(list(hindcast_only, calibration_only), fill = T) %>%
	arrange(model_name, fcast_type, pretty_group, taxon, plotID, timepoint)

saveRDS(hindcast_data_out, here("data/summary/all_hindcasts.rds"))



# truth_dat = readRDS(here("data/clean/truth_plot_long_tax.rds"))
# truth_dat %>% filter(siteID=="MLBS" & taxon == "actinobacteria" & !is.na(truth))


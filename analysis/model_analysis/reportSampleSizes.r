

fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")  %>% filter(!grepl("Aquatic", field_site_type)) %>%
	select(siteID = field_site_id,
				 MAT = field_mean_annual_temperature_C,
				 MAP = field_mean_annual_precipitation_mm,
				 latitude = field_latitude,
				 longitude = field_longitude,
				 nlcd = field_dominant_nlcd_classes,
				 ecoregion = field_domain_id)

hindcast_201511_201801 <- readRDS("./data/summary/beta_hindcast_fg_2015-11_2018-01.rds")

hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))

# Subset data from a bacterial group and fungal group since sample sizes differ slightly
hindcast_201511_201801 <- hindcast_in %>% filter(model_name=="env_cycl" & taxon %in% c("ectomycorrhizal", "copiotroph"))

# How many data points for forecast validation ? #424
n_newtime_plot_obs = hindcast_201511_201801 %>% filter(new_site==F & !is.na(truth)) %>%
	select(plotID, dateID) %>% distinct() %>% tally

# How many data points for out-of-sample (new-site) validation? #162
n_newsite_plot_obs = hindcast_201511_201801 %>% filter(new_site==T & !is.na(truth)) %>%
	select(plotID, dateID) %>% distinct() %>% tally

# How many samples in entire dataset? #6005 for fungi
abun_its = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/clean/groupAbundances_ITS_2023.rds")
unique(abun_its$phylum_fun$sampleID) %>% length

# How many samples in entire dataset? #6080 for bacteria
abun_16s = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/clean/groupAbundances_16S_2023.rds")
unique(abun_16s$phylum_bac$sampleID) %>% length

calibration = hindcast_201511_201801 %>% filter(fcast_period == "calibration" & !is.na(truth))


# How many samples are used for calibration? #3128 for bacteria
model.inputs_bac <- prepBetaRegData(rank.df = abun_16s$copiotroph,
														 min.prev = 3,
														 min.date = "20150101",
														 max.date = "20180101")
nrow(model.inputs_bac$y)

# How many unique plot/dates does this amount to? #1283
model.inputs_bac$truth.plot.long %>% filter(!is.na(truth)) %>%
	select(plotID, dateID) %>% distinct() %>% tally

# How many samples are used for calibration? #2705 for bacteria
model.inputs_fun <- prepBetaRegData(rank.df = abun_its$ectomycorrhizal,
																		min.prev = 3,
																		min.date = "20150101",
																		max.date = "20180101")
nrow(model.inputs_fun$y)

# How many unique plot/dates does this amount to? #1213
model.inputs_fun$truth.plot.long %>% filter(!is.na(truth)) %>%
	select(plotID, dateID) %>% distinct() %>% tally

# Now get sample sizes per site to create site table
fun_site_n = model.inputs_fun$sample_values %>% filter(!is.na(other)) %>% group_by(siteID) %>% distinct() %>% tally(name = "Calibration fungal samples")
bac_site_n = model.inputs_bac$sample_values %>% filter(!is.na(other)) %>% group_by(siteID) %>% distinct() %>% tally(name = "Calibration bacterial samples")

newsite_n = hindcast_201511_201801 %>% filter(new_site==T & !is.na(truth)) %>% select(siteID) %>% distinct %>% mutate("Calibration fungal samples" = 0,
																																																													"Calibration bacterial samples" = 0)

fieldsites_samples <- merge(bac_site_n, fun_site_n, all=T)
fieldsites_samples <- merge(fieldsites_samples, newsite_n, all=T)
fieldsites_samples <- merge(fieldsites_samples, fieldsites_raw, all=T)

# Clean and aggregate the raw DAAC soil-moisture (SMV) files into per-site
# monthly means. The raw per-site SMV files are not included in this repo
# (see Methods/Zenodo); set the directory below to their location.
library(dplyr)
library(tidyr)


#### 1. GET/PREP SMV DATA ####
# Directory holding the raw DAAC SMV files
smv_dir <- here::here("data/raw/DAAC_SMV/")
files <- list.files(smv_dir, full.names = T)

# Get list of NEON fieldsites/locations
fieldsites <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")

# Loop through files and determine the site each one came from.
df.list <- list()
for (f in files){
  loc <- readLines(f, n = 3)[[3]]
  loc <- strsplit(loc, split = ":|,")
  lat <- gsub(" ", "", loc[[1]][[2]])
  lon <- gsub(" ", "", loc[[1]][[3]])
  latsite <- fieldsites[match(lat, fieldsites$Latitude),]$Site.ID
  lonsite <- fieldsites[match(lon, fieldsites$Longitude),]$Site.ID
  site <- intersect(latsite, lonsite)
  df <- read.csv(f, skip = 4)
  df$siteID <- site
  df.list[[site]] <- df
  print(site)
}
smv_all <- do.call(rbind, df.list)

# Fix weird columns
smv_all <- smv_all %>% mutate(month = substr(time, 1, 7)) %>%
  separate(SMAP_surface, into = c("min_SMAP_s", "max_SMAP_s","mean_SMAP_s"), sep = ";", convert = T) %>%
  separate(SMAP_rootzone, into = c("min_SMAP_r", "max_SMAP_r","mean_SMAP_r"), sep = ";", convert = T) %>%
  separate(SCAN_surface, into = c("min_SCAN_s", "max_SCAN_s","mean_SCAN_s"), sep = ";", convert = T) %>%
  separate(SCAN_rootzone, into = c("min_SCAN_r", "max_SCAN_r","mean_SCAN_r"), sep = ";", convert = T) %>%
  separate(GRACE_surface_pctl, into = c("min_GRACE_s", "mean_GRACE_s","max_GRACE_s"), sep = ";", convert = T)

# Mean per site/month
smv_month <- smv_all %>%
  filter(time > "2013-06-01") %>%
  group_by(siteID, month) %>%
  dplyr::summarise(across(starts_with(c("mean", "min", "max")), ~mean(.x, na.rm = TRUE)), .groups="keep")
smv_month_long <- smv_month %>% pivot_longer(cols = 3:17)
saveRDS(smv_month, here::here("data/raw/DAAC_SMV/monthly_SMV_allsites.rds"))

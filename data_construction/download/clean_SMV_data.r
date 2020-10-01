# download and fix up SMV data from the raw files Jack downloaded.
library(dplyr)
library(tidyr)


#### 1. GET/PREP SMV DATA ####
# Directory where I put all of Jack's files
smv_dir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/DAAC_SMV/"
files <- list.files(smv_dir, full.names = T)

# Get list of NEON fieldsites/locations
fieldsites <- read.csv("https://www.neonscience.org/science-design/field-sites/export")

# Loop through and figure out which site each one came from.
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
saveRDS(smv_month, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/DAAC_SMV/monthly_SMV_allsites.rds")


####################################
## unused/abandoned code below ##
####################################


# Get data from nearby SCAN sites when possible?
library("metScanR")
sm <- getVars("soil moisture", startVarsDate = "2013-06-01", endVarsDate = "2014-12-01")
dsny <- getNearby(siteID = "NEON:DSNY", radius = 150)
intersect(sm, dsny)


library("soilDB")
dsny <- fetchSCAN(site.code = "2012", year = c(2013,2014))
dsny_mois <- dsny$SMS %>% filter(depth < 6) %>% mutate(month = substr(Date, 1, 7)) %>% group_by(month) %>%  summarize(mean = mean(value, na.rm=T))



## Trying the SMOS dataset, netCDF
library(ncdf4)
library(raster)

# Function to get nearest date within netcdf
wherenearest <- function(myPoint, allPoints){
  d <- abs(allPoints-myPoint[1])
  index <- which.min(d)
  return( index )
}

sites <- sort(unique(merged$siteID))
months <- gsub("-","", sort(unique(merged$month)))
months <- months[which(months < "201512")]
# Use descending orbit bc has more data points
allfiles <- list.files("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/SMOS/DES/", full.names = T)
out.list <- list()
for (m in months){
    file <- grep(m, allfiles)
    fname <- allfiles[file]
# Read in netcdf file and get nearest lat/lon indices
nc <- nc_open(fname)
nc_lat <- ncvar_get(nc, attributes(nc$dim)$names[2])
nc_lon <- ncvar_get(nc, attributes(nc$dim)$names[3])
month.df <- cbind.data.frame(month = m, siteID = sites, SM = NA)
for (s in 1:length(sites)){
  # Get lat/lon for site
  site_lon <- fieldsites[fieldsites$Site.ID==sites[s],]$Longitude
  site_lat <- fieldsites[fieldsites$Site.ID==sites[s],]$Latitude
# Get nearest SM value
SM <- ncvar_get(nc, "SM")
SM_val <- SM[wherenearest(site_lon, nc_lon), wherenearest(site_lat, nc_lat)]
month.df[s,]$SM <- SM_val
}
nc_close(nc)
out.list[[m]] <- month.df
}
out <- do.call(rbind, out.list)

merged$month <-  gsub("-","", merged$month)
smos <- merge(merged, out, by=c("siteID", "month"))
smos$SMOS <- smos$SM * 100
smos_need <- smos[smos$siteID %in% need,]

plot(smos$neon_mean, smos$SMOS, col = as.factor(smos$siteID))
plot(smos_need$neon_mean, smos_need$SMOS, col = as.factor(smos_need$siteID))
plot(smos_need$neon_mean, smos_need$SMOS, col = smos_need$siteID)
plot(smos_need$neon_mean, smos_need$SMOS, col = as.numeric(as.factor(smos_need$siteID)))

ggplot(smos) + geom_point(aes(x = neon_mean, y = SMOS, color = as.factor(siteID))) + geom_abline()
ggplot(smos) + geom_point(aes(x = neon_mean, y = mean_SMAP_s, color = as.factor(siteID))) + geom_abline() + facet_grid(~siteID)

## Cleaning the SMOS dataset, converting from netCDF files to data frame
library(ncdf4)
library(raster)

# Function to get nearest date within netcdf
wherenearest <- function(myPoint, allPoints){
  d <- abs(allPoints-myPoint[1])
  index <- which.min(d)
  return( index )
}

# Get NEON data for relevant sites.
neon <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilMois_monthly_allsites.rds")
# Pull out sites/months of interest.
sites <- sort(unique(neon$siteID))
months <- gsub("-","", sort(unique(neon$month)))
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
saveRDS(out, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/SMOS/monthly_SMOS_allsites.rds")

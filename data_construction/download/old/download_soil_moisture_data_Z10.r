# download soil moisture, using modified Z10 functions. Should request multiple cores to run this.

# Setup required for all downloads

library(tidyverse)
library(Z10)
library(parallel)
library(doParallel)
library(foreach)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctionsBasic.r")
outdir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw"

#Get full site list
avail = Z10::dp.avail("DP1.00094.001")
mp.avail = Z10::dp.avail("DP1.00096.001")
all_sites <- unlist(unique(avail$site))

# Create cluster
cl <- makeCluster(10, outfile="")
clusterExport(cl=cl, c("avail","all_sites","getSoilData", "getTISdata","tic","toc","outdir"), envir=environment())

run.raw.soil.moist <- function(s){
  # Get list of all NEON sites
  tic()
  raw.soil.moist.site <- getSoilData(all_sites[s], dp = "DP1.00094.001", avail = avail)
  saveRDS(raw.soil.moist.site, paste0(outdir, "/moist_z10_approach_", all_sites[[s]], ".rds"))
  toc()
  return(raw.soil.moist.site)
}
output.list <- parLapply(cl, 1:length(all_sites), run.raw.soil.moist)
lapply(output.list, dim)

stopCluster(cl)
closeAllConnections()



# 2. CHECK WHAT DID/DIDN'T SUCCESSFULLY DOWNLOAD
# Read in data and see what's successfully downloaded.
existing_files <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw/", pattern = "z10_approach", full.names = T)
df.list <- list()
for (i in existing_files){
  df <- readRDS(i)
  df.list[[i]] <- df %>% distinct(siteID, months = substr(df$startDateTime, 1, 7)) # get unique site/months
  rm(df) # since the original files are so big
}
master.df <- do.call(plyr::rbind.fill, df.list) # combine into df of unique site/months

# Get list of available site/months for data product
avail.site.months <- avail %>% 
  unnest(months) %>% unnest(site) %>% rename(siteID = site) %>% as.data.frame()

# Which ones still need to be downloaded?
to.be.downloaded <- setdiff(avail.site.months, master.df)
sites.to.be.downloaded <- unique(to.be.downloaded$siteID)

clusterExport(cl=cl, c("sites.to.be.downloaded","avail","all_sites","getSoilData", "getTISdata","tic","toc","outdir"), envir=environment())

# Download just those!
run.raw.soil.moist.site <- function(site){
  # Get list of all NEON sites
  tic()
  raw.soil.moist.site <- getSoilData(site, dp = "DP1.00094.001", avail = avail)
  saveRDS(raw.soil.moist.site, paste0(outdir, "/moist_z10_approach_", site, ".rds"))
  toc()
  return(raw.soil.moist.site)
}
output.list <- parLapply(cl = cl, X = sites.to.be.downloaded, fun = run.raw.soil.moist.site)


# size <- list()
# for (obj in ls()) { size[[obj]] <- pryr::object_size(get(obj))}; size

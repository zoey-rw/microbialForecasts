# download soil temp, using modified Z10 functions

# Setup required for all downloads

library(tidyverse)
library(Z10)
library(parallel)
library(doParallel)
library(foreach)
outdir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw"
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/development/Z10_modified_functions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctionsBasic.r")

# Get full site list
avail = Z10::dp.avail("DP1.00041.001")
mp.avail = Z10::dp.avail("DP1.00096.001")
all_sites <- unlist(unique(avail$site))#[11:15]

# Create cluster
cl <- makeCluster(16, outfile="")
clusterExport(cl=cl, c("avail","mp.avail","all_sites","raw.soil.temp2","outdir","getSoilData", "getTISdata"), envir=environment())

#### 1. INITIAL DOWNLOAD (kind of wonky) 
run.raw.soil.temp <- function(s){
  source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctionsBasic.r")
  # Get list of all NEON sites
  raw.soil.temp.site <- raw.soil.temp2(site = all_sites[s], avail, mp.avail)
  saveRDS(raw.soil.temp.site, paste0(outdir, "/z10_approach_sites", s, ".rds"))
  return(raw.soil.temp.site)
}
output.list <- parLapply(cl, 1:length(all_sites), run.raw.soil.temp)
output.list <- parLapply(cl, 40:length(all_sites), run.raw.soil.temp)
names(output.list) <- all_sites
out <- rbind.named.dfs(output.list)

saveRDS(out, paste0(outdir, "/z10_approach_sites16-20.rds"))
stopCluster(cl)
closeAllConnections()

# 2. CHECK WHAT DID/DIDN'T SUCCESSFULLY DOWNLOAD
# Read in data and see what's successfully downloaded.
existing_files <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw/", pattern = "z10_approach", full.names = T)
df.list <- list()
for (i in 1:length(existing_files)){
  print(i)
  df <- readRDS(existing_files[[i]])
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

#### 3. DOWNLOAD ONLY THOSE SITES (should later become site/dates)
#cl <- makeCluster(10, outfile="")
run.raw.soil.temp <- function(site){
  source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctionsBasic.r")
  # Get list of all NEON sites
  raw.soil.temp.site <- getSoilData(site = site, avail = avail)
  saveRDS(raw.soil.temp.site, paste0(outdir, "/z10_approach_", site, ".rds"))
  return(raw.soil.temp.site)
}
output.list <- parLapply(cl = cl, X = sites.to.be.downloaded, fun = run.raw.soil.temp)




registerDoParallel(cl)
output.list <- foreach(s = 1:length(sites.to.be.downloaded), 
                                              .verbose = T, 
                       .packages=c("Z10","dplyr","rjson"), 
                       .errorhandling = "pass") %dopar% {
                         tic()
                         raw.soil.temp.site <- getSoilData(site = sites.to.be.downloaded[[s]], avail=avail)
                         toc()
                         saveRDS(raw.soil.temp.site, paste0(outdir, "/z10_approach_", sites.to.be.downloaded[[s]], ".rds"))
                         return(raw.soil.temp.site)
                       }






# PREV APPROACH - COULDN'T GET TO WORK RIGHT.
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/development/Z10_modified_functions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
outdir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw"
library(foreach)
library(doParallel)
cl<-makeCluster(8, outfile="")
registerDoParallel(cl)
# Get list of all NEON sites
avail = Z10::dp.avail("DP1.00041.001")
mp.avail = Z10::dp.avail("DP1.00096.001")
all_sites <- unlist(unique(avail$site))
output.list <- foreach(s = 41:length(all_sites), 
                       .verbose = T, 
                       .packages=c("Z10","dplyr","rjson"), 
                       .errorhandling = "pass") %dopar% {
  tic()
                         raw.soil.temp.site <- raw.soil.temp2(site = all_sites[s], avail, mp.avail, bgn.date = "2019-01", end.date = "2019-02")
                         raw.soil.temp.site <- getSoilData(site = all_sites[s], avail = avail, bgn.date = "2019-01", end.date = "2019-02")
                         toc()
  saveRDS(raw.soil.temp.site, paste0(outdir, "/z10_approach_sites", s, ".rds"))
  return(raw.soil.temp.site)
}
lapply(output.list, dim)
stopCluster(cl)


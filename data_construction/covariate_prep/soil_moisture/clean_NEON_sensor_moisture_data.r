library(data.table)
library(dplyr)
library(neonstore)

# 12 million rows. oy vey.
#NEON_moist_30m_orig <- fread("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_daily/rawMoistStacked.rds/stackedFiles/SWS_30_minute.csv")


# Read in data and see what's successfully downloaded.
# existing_files <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw/", pattern = "z10_approach", full.names = T)
# df.list <- lapply(existing_files, readRDS)
# NEON_moist_30m_orig <- data.table(do.call(plyr::rbind.fill, df.list))) # combine into df of unique site/months

NEON_moist_30m_orig <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw_allsites.rds")
# subset to second shallowest depth, because this has the closest correlation with sampled moisture values.
#NEON_moist_30m <- NEON_moist_30m_orig[NEON_moist_30m_orig$verticalPosition %in% c(502,2),] # JK something is wrong with the loading of these files so we have to use the shallowest sensors.
NEON_moist_30m <- NEON_moist_30m_orig[NEON_moist_30m_orig$verticalPosition %in% c(501,1),]

# Take mean sensor array
moist_by_site <- NEON_moist_30m %>% 
  group_by(siteID, startDateTime) %>% 
  mutate(by_site = mean(VSWCMean, na.rm=T)) %>% distinct(siteID, startDateTime, .keep_all = T)
#sites <- c("CPER","OSBS","DSNY","STER","HARV")

sites <- unique(moist_by_site$siteID)
# Calculating uncertainties by month
NEON_moist_monthly_allsite <- list()
for (s in 1:length(sites)) { #loop through all sites
  
  siteID <- sites[s]
  print(siteID)
  #### 2. For NEON moist data, get monthly mean and variance ####
  NEON_moist_30m_site <- moist_by_site[which(moist_by_site$siteID==siteID),]
  allmonths <- sort(unique(NEON_moist_30m_site$month))
  # set up output df
  NEON_moist_monthly <- setNames(data.frame(matrix(ncol = 6, nrow = length(months))), 
                                c("siteID", "month", "NEON_moist_mean", "NEON_moist_sd", "simpleMean", "simpleSD"))
  
  monthMeans <- list()
  for (m in 1:length(allmonths)) { 
    
    # Subset 30-min data from one month - get means and uncertainties
    month <- allmonths[m]
    NEON_moist_monthly[m,]$month <- month
    NEON_moist_monthly[m,]$siteID <- siteID
    
    # need to na.omit since we're sampling from these values.s
    moist_month <- NEON_moist_30m_site[which(NEON_moist_30m_site$month == month),] %>% ungroup() %>%  
      mutate(siteID = forcats::fct_explicit_na(siteID)) %>% stats::na.omit()
    
    # discard if we have measurements for less than 15 days.
    if (length(unique(moist_month$day)) < 15) next()
    #if (nrow(moist_month) < 1000) next() # arbitrary cutoff
    print(m)
    tsY <- moist_month$VSWCMean  # 30-min means
    uc <- moist_month$VSWCExpUncert/2 # expanded measurement uncertainty at 95% confidence
    uNAT <- moist_month$VSWCStdErMean
    ub.sq <- (uc^2) - (uNAT^2)
    ub.sq[which(ub.sq < 0)] <- .0001 # negative values from rounded uncertainties...
    ub <- sqrt(ub.sq)
    
    # Estimate monthly mean/SD
    nm = 5000
    tsamp = rep(NA,nm) 
    tsamp_sd = rep(NA,nm) 
    for(i in 1:nm){
      ub.samp <- sample(ub, 1) # sample 1 bias per month
      bias = rnorm(1,0,1) * ub.samp # assume that the bias is perfectly correlated from obs to obs
      
      uNAT.samp <- sample(uNAT, 1) 
      tsim = rnorm(nm,tsY,uNAT.samp) ## sample random error per half hour
      
      tsamp[i] = mean(tsim,na.rm=TRUE)+bias # output monthly mean
      tsamp_sd[i] = sd(tsim,na.rm=TRUE)  # output monthly sd
    }
    NEON_moist_monthly[m,]$NEON_moist_mean <- mean(tsamp)
    NEON_moist_monthly[m,]$NEON_moist_sd <- mean(tsamp_sd)
    
    # just for sanity checking.
    NEON_moist_monthly[m,]$simpleMean <- mean(tsY, na.rm=T)
    NEON_moist_monthly[m,]$simpleSD <- sd(tsY, na.rm = T)
  }
  NEON_moist_monthly_allsite[[s]] <- NEON_moist_monthly
}
NEON_moist_monthly_all <- do.call(rbind, NEON_moist_monthly_allsite)  
saveRDS(NEON_moist_monthly_all, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilMois_monthly_allsites.rds")





## NOW WEEKLY moistS
# 
# library(lubridate)
# library(zoo)
# moist_by_site$week <- lubridate::week(ymd(as.Date(moist_by_site$startDateTime)))
# moist_by_site$year <- lubridate::year(ymd(as.Date(moist_by_site$date_time)))
# moist_by_site$year_week <- paste0(moist_by_site$year, "_", moist_by_site$week)
# 
# NEON_moist_weekly_allsite <- list()
# 
# for (s in 1:length(sites)) { #loop through all sites
#   
#   siteID <- sites[s]
#   print(siteID)
#   #### 2. For NEON moist data, get weekly mean and variance ####
#   
#   #subset by site
#   NEON_moist_30m_site <- moist_by_site[which(moist_by_site$siteID==siteID),]
#   NEON_moist_30m_site <- NEON_moist_30m_site[order(sort(as.Date(NEON_moist_30m_site$startDateTime))),]
#   # get all potential weeks
#   allweeks <- unique(NEON_moist_30m_site$year_week)
#   # set up output df
#   NEON_moist_weekly <- setNames(data.frame(matrix(ncol = 6, nrow = length(allweeks))), 
#                                c("siteID", "year_week", "NEON_moist_mean", "NEON_moist_sd", "simpleMean", "simpleSD"))
#   
#   weekMeans <- list()
#   for (m in 1:length(allweeks)) { 
#     
#     # Subset 30-min data from one week - get means and uncertainties
#     week <- allweeks[m]
#     NEON_moist_weekly[m,]$year_week <- week
#     NEON_moist_weekly[m,]$siteID <- siteID
#     
#     # need to na.omit since we're sampling from these values.s
#     moist_week <- NEON_moist_30m_site[which(NEON_moist_30m_site$year_week == week),] %>% ungroup() %>%  
#       mutate(siteID = forcats::fct_explicit_na(siteID)) %>% stats::na.omit()
#     
#     if (nrow(moist_week) < 500) next() # arbitrary cutoff
#     print(m)
#     tsY <- moist_week$VSWCMean  # 30-min means
#     uc <- moist_week$VSWCExpUncert/2 # expanded measurement uncertainty at 95% confidence
#     uNAT <- moist_week$VSWCStdErMean
#     ub.sq <- (uc^2) - (uNAT^2)
#     ub.sq[which(ub.sq < 0)] <- .0001 # negative values from rounded uncertainties...
#     ub <- sqrt(ub.sq)
#     
#     # Estimate weekly mean/SD
#     nm = 2000
#     tsamp = rep(NA,nm) 
#     tsamp_sd = rep(NA,nm) 
#     for(i in 1:nm){
#       ub.samp <- sample(ub, 1) # sample 1 bias per month
#       bias = rnorm(1,0,1) * ub.samp # assume that the bias is perfectly correlated from obs to obs
#       
#       uNAT.samp <- sample(uNAT, 1) 
#       tsim = rnorm(nm,tsY,uNAT.samp) ## sample random error per half hour
#       
#       tsamp[i] = mean(tsim,na.rm=TRUE)+bias # output weekly mean
#       tsamp_sd[i] = sd(tsim,na.rm=TRUE)  # output weekly sd
#     }
#     NEON_moist_weekly[m,]$NEON_moist_mean <- mean(tsamp)
#     NEON_moist_weekly[m,]$NEON_moist_sd <- mean(tsamp_sd)
#     
#     # just for sanity checking.
#     NEON_moist_weekly[m,]$simpleMean <- mean(tsY, na.rm=T)
#     NEON_moist_weekly[m,]$simpleSD <- sd(tsY, na.rm = T)
#   }
#   NEON_moist_weekly_allsite[[s]] <- NEON_moist_weekly
# }
# NEON_moist_weekly_all <- do.call(rbind, NEON_moist_weekly_allsite)  
# saveRDS(NEON_moist_weekly_all, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilMoist_weekly_5sites.rds")


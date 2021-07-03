library(data.table)
library(dplyr)

# Load in massive dataset!
NEON_temp_30m_orig <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw_allsites.rds")

# subset to shallowest depth, because this has the closest correlation with sampled temperature values.
NEON_temp_30m <- NEON_temp_30m_orig[NEON_temp_30m_orig$verticalPosition %in% c(501,1),]


# Mean of sensor array
temp_by_site <- NEON_temp_30m %>% 
  group_by(siteID, startDateTime) %>% 
  mutate(by_site = mean(soilTempMean, na.rm=T))  %>% 
  distinct(siteID, startDateTime, .keep_all = T)

sites <- unique(temp_by_site$siteID)
# Calculating uncertainties by month
NEON_temp_monthly_allsite <- list()
for (s in 1:length(sites)) { #loop through all sites
  
  siteID <- sites[s]
  print(siteID)
  #### 2. For NEON temp data, get monthly mean and variance ####
  NEON_temp_30m_site <- temp_by_site[which(temp_by_site$siteID==siteID),]
  allmonths <- sort(unique(NEON_temp_30m_site$month))
  # set up output df
  NEON_temp_monthly <- setNames(data.frame(matrix(ncol = 6, nrow = length(months))), 
                                c("siteID", "month", "NEON_temp_mean", "NEON_temp_sd", "simpleMean", "simpleSD"))
  
  monthMeans <- list()
  for (m in 1:length(allmonths)) { 
    
    # Subset 30-min data from one month - get means and uncertainties
    month <- allmonths[m]
    NEON_temp_monthly[m,]$month <- month
    NEON_temp_monthly[m,]$siteID <- siteID
    
    # need to na.omit since we're sampling from these values.s
    temp_month <- NEON_temp_30m_site[which(NEON_temp_30m_site$month == month),] %>% ungroup() %>%  
      mutate(siteID = forcats::fct_explicit_na(siteID)) %>% stats::na.omit()
    
    # discard if we have measurements for less than 15 days.
    if (length(unique(temp_month$day)) < 15) next()
    #if (nrow(temp_month) < 1000) next() # arbitrary cutoff
    print(m)
    tsY <- temp_month$soilTempMean  # 30-min means
    uc <- temp_month$soilTempExpUncert/2 # expanded measurement uncertainty at 95% confidence
    uNAT <- temp_month$soilTempStdErMean
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
    NEON_temp_monthly[m,]$NEON_temp_mean <- mean(tsamp)
    NEON_temp_monthly[m,]$NEON_temp_sd <- mean(tsamp_sd)
    
    # just for sanity checking.
    NEON_temp_monthly[m,]$simpleMean <- mean(tsY, na.rm=T)
    NEON_temp_monthly[m,]$simpleSD <- sd(tsY, na.rm = T)
  }
  NEON_temp_monthly_allsite[[s]] <- NEON_temp_monthly
}
NEON_temp_monthly_all <- do.call(rbind, NEON_temp_monthly_allsite)  
saveRDS(NEON_temp_monthly_all, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilTemp_monthly_allsites.rds")







# ## NOW WEEKLY TEMPS
# 
# library(lubridate)
# library(zoo)
# temp_by_site$week <- lubridate::week(ymd(as.Date(temp_by_site$startDateTime)))
# temp_by_site$year <- lubridate::year(ymd(as.Date(temp_by_site$date_time)))
# temp_by_site$year_week <- paste0(temp_by_site$year, "_", temp_by_site$week)
# 
# NEON_temp_weekly_allsite <- list()
# 
# for (s in 1:length(sites)) { #loop through all sites
# 
#   siteID <- sites[s]
#   print(siteID)
#   #### 2. For NEON temp data, get weekly mean and variance ####
# 
#   #subset by site
#   NEON_temp_30m_site <- temp_by_site[which(temp_by_site$siteID==siteID),]
#   NEON_temp_30m_site <- NEON_temp_30m_site[order(sort(as.Date(NEON_temp_30m_site$startDateTime))),]
#   # get all potential weeks
#   allweeks <- unique(NEON_temp_30m_site$year_week)
#   # set up output df
#   NEON_temp_weekly <- setNames(data.frame(matrix(ncol = 6, nrow = length(allweeks))),
#                                 c("siteID", "year_week", "NEON_temp_mean", "NEON_temp_sd", "simpleMean", "simpleSD"))
# 
#   weekMeans <- list()
#   for (m in 1:length(allweeks)) {
# 
#     # Subset 30-min data from one month - get means and uncertainties
#     week <- allweeks[m]
#     NEON_temp_weekly[m,]$year_week <- week
#     NEON_temp_weekly[m,]$siteID <- siteID
# 
#     # need to na.omit since we're sampling from these values.s
#     temp_week <- NEON_temp_30m_site[which(NEON_temp_30m_site$year_week == week),] %>% ungroup() %>%
#       mutate(siteID = forcats::fct_explicit_na(siteID)) %>% stats::na.omit()
# 
#     if (nrow(temp_week) < 500) next() # arbitrary cutoff
#     print(m)
#     tsY <- temp_week$soilTempMean  # 30-min means
#     uc <- temp_week$soilTempExpUncert/2 # expanded measurement uncertainty at 95% confidence
#     uNAT <- temp_week$soilTempStdErMean
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
#     NEON_temp_weekly[m,]$NEON_temp_mean <- mean(tsamp, na.rm=T)
#     NEON_temp_weekly[m,]$NEON_temp_sd <- mean(tsamp_sd, na.rm=T)
# 
#     # just for sanity checking.
#     NEON_temp_weekly[m,]$simpleMean <- mean(tsY, na.rm=T)
#     NEON_temp_weekly[m,]$simpleSD <- sd(tsY, na.rm = T)
#   }
#   NEON_temp_weekly_allsite[[s]] <- NEON_temp_weekly
# }
# NEON_temp_weekly_all <- do.call(rbind, NEON_temp_weekly_allsite)
# saveRDS(NEON_temp_weekly_all, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilTemp_weekly_5sites.rds")

  
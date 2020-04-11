# Setup temperature time-series for microbial forecasts at NEON sites
# Using NEON temperature data (Single Aspirated Air Temperature) and ERA5 ensembles (n=10)

#### 1. Load NEON data and ERA5 ensemble data. ####
#### 2. For NEON data, get monthly mean and variance ####
#### 3. Use NEON-ERA5 calibration to adjust ERA5 ensemble means.
#### 4. Merge dataframes together, take NEON values when possible, use adjusted ERA5 otherwise.
rm(list=ls())
library(neonUtilities)
library(zoo)
library(rjags)
library(coda)

# 200mb file, didn't wanna move it
load("/usr3/graduate/zrwerbin/temporal_forecasts/data/temp_time_series_input.Rdata")

output.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/temp_time_series_output.rds"

sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")
NEON_temp_monthly_all.sites <- data.frame(NULL)

for (s in 1:5) { #loop through all sites
  
  siteID <- sites[s]
  
  #### 2. For NEON temp data, get monthly mean and variance ####
  NEON_temp_30m_site <- NEON_temp_30m[NEON_temp_30m$siteID==siteID,]
  allmonths <- sort(unique(NEON_temp_30m_site$month))
  # set up output df
  NEON_temp_monthly <- setNames(data.frame(matrix(ncol = 4, nrow = length(months))), 
                                c("siteID", "month", "NEON_temp_mean", "NEON_temp_sd"))
  
  for (m in 1:length(allmonths)) { 
    
    # Subset 30-min data from one month - get means and uncertainties
    month <- allmonths[m]
    NEON_temp_monthly[m,]$month <- month
    NEON_temp_monthly[m,]$siteID <- siteID
    temp_month <- NEON_temp_30m_site[NEON_temp_30m_site$month == month,]
    if (nrow(temp_month) < 1500) next() # arbitrary cutoff
    tsY <- temp_month$tempSingleMean # 30-min means
    tsSE <- temp_month$tempSingleStdErMean
    uc <- temp_month$tempSingleExpUncert/2
    uNAT <- temp_month$tempSingleStdErMean
    ub.sq <- (uc^2) - (uNAT^2)
    ub.sq[which(ub.sq < 0)] <- .0001 # negative values from rounded uncertainties...
    ub <- sqrt(ub.sq)
    
    # Estimate monthly mean/SD
    nm = 1000
    tsamp = rep(NA,nm) 
    monthMean <- mean(temp_month$tempSingleMean, na.rm=T)
    for(i in 1:nm){
      ub.samp <- sample(ub, 1)
      uNAT.samp <- sample(uNAT, 1)
      bias = rnorm(1,0,1)*ub.samp ## sample 1 bias per month
      tsim = rnorm(nm,monthMean,uNAT.samp) ## sample random error per half hour
      tsamp[i] = mean(tsim,na.rm=TRUE)+bias
    }
    NEON_temp_monthly[m,]$NEON_temp_mean <- mean(tsamp)
    NEON_temp_monthly[m,]$NEON_temp_sd <- sd(tsamp)
  }
  
  #### 3. Use NEON-ERA5 calibration to adjust ERA5 ensemble means. Does not incorporate NEON variances.
  
  # Subset ERA5 data by site/month
  ensemble_temp <- clim_ensembles[[siteID]][[1]]
  ERA5 <- data.frame(ERA5_mean = apply(ensemble_temp,2,mean), ERA5_sd = apply(ensemble_temp,2,sd))
  ERA5$Date <- seq.POSIXt(as.POSIXlt("2013-01-01 0:00"), as.POSIXlt("2018-12-01 0:00"), by="month")
  ERA5$month <- substr(ERA5$Date, 1, 7)
  
  temp.model <- " model {
  ### ERA5 calibration priors
  alpha[1] ~ dnorm(0,0.001)   # intercept should be ~0, based on ERA5/NEON plots
  alpha[2] ~ dnorm(0,0.01)    # uninformative prior on coefficient
  sigma ~ dgamma(6,50)        # derived from the distribution of SDs for ERA5 ensemble means
  tau <- 1/(sigma*sigma)      
  
  for(i in 1:N){
  mu[i] <-  alpha[1] + alpha[2]*x[i]      ## process model
  y[i] ~ dnorm(mu[i],tau)               ## data model
  }
}"
merged <- merge(NEON_temp_monthly, ERA5, by="month") # get overlapping months, for calibration
  data.list <- list(x = merged$ERA5_mean, y=merged$NEON_temp_mean, N = nrow(merged))
  # Fit model
  j.model   <- jags.model(file = textConnection(temp.model),
                          data = data.list,
                          n.chains = 3)
  var.out   <- coda.samples (model = j.model,
                             variable.names = c("tau", "alpha[1]", "alpha[2]"),
                             n.iter = 2000)
  
  # Use EIV model to adjust ERA5 values
  mcmc <- do.call(rbind, var.out)
  ERA5[,c("adjusted_mean", "adjusted_sd")] <- NA
  for (e in 1:nrow(ERA5)) {
    convert.out <- list()
    for(i in 1:1000){
      par <- mcmc[sample(nrow(mcmc),1),]
      b <- par[1] + par[2]*ERA5$ERA5_mean[e]
      convert.out[[i]] <- b
    }
    convert.out <- do.call(cbind, convert.out)
    ERA5$adjusted_mean[e] <- rowMeans(convert.out, na.rm = T)
    ERA5$adjusted_sd[e] <- apply(convert.out, 1, sd, na.rm = T)
  }
  
  #### 4. Merge dataframes together, take NEON values when possible, use adjusted ERA5 otherwise.
  output <- merge(ERA5, NEON_temp_monthly, all=T)
  output$monthly_temp <- output$NEON_temp_mean
  output$monthly_temp_sd <- output$NEON_temp_sd 
  output$source <- "NEON"
  output[is.na(output$NEON_temp_mean),]$monthly_temp_sd <- output[is.na(output$monthly_temp),]$adjusted_sd
  output[is.na(output$NEON_temp_mean),]$monthly_temp <- output[is.na(output$monthly_temp),]$adjusted_mean
  output[is.na(output$NEON_temp_mean),]$source <- "ERA5_estimate"
  output$siteID <- siteID
  temp_time_series <- output[,c("siteID", "month", "monthly_temp", "monthly_temp_sd", "source")]
  print(output)
  NEON_temp_monthly_all.sites <- rbind(NEON_temp_monthly_all.sites, temp_time_series)
  }

saveRDS(NEON_temp_monthly_all.sites, output.path)
#load("temp_time_series.Rdata")

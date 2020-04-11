# CODE TO CREATE FIGURE 1 from EE509 preliminary analyses.

library(ecoforecastR)
library(runjags)
library(arm)
library(dplyr)
library(tidyr)
library(RColorBrewer)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepModelData.r")

## Read in model
mod.list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")
#mod.list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS.rds")

## Read in covariate data 
cov <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
weather <- cov[[1]]
plot_pH <- cov[[2]]

# Read in the microbial abundance data
d <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/groupAbundances_EE509.rds")

# Decide which site/plot to look at 
site <- 1 
site.name <- unique(as.factor(weather$siteID))[site]
plot <- 1
plot.name <- unique(as.factor(plot_pH$plotID))[plot]

# Set characteristics for all plots/forecasts
time.t = 1:45    ## total time
time1 = 1:25       ## calibration period
time2 = 26:45   ## forecast period
N.cols <- brewer.pal(4, "Set2") ## set colors
Nmc = 1000
NT <- 20
ngibbs <- 1000
npred <- 25

pl.out <- list()

rank.df <- d[[k]]
rank.df$dates <- gsub("2016", "2015", rank.df$dates)

pl.rank <- list()
#for(j in 1:length(mod.list)) {
  for(j in 1:5) {
  pmod <- mod.list[[j]]

  # Get prediction data
  model.dat <- prepModelData(weather, plot_pH, rank.df, j=j) # Combine using custom function
  plot.truth <- model.dat[[7]]
  site.truth <- model.dat[[8]]
  pH.cov <- model.dat[[5]]
  plot_out <- pmod$summary
  # Get observed values
  obs <- plot_out[grep(paste0("plot_mean[",plot,","), rownames(plot_out),fixed=T),]$truth
  
  # Extract parameter estimates from model
  mfit = pmod$JAGS$mcmc[[1]]
  pred.cols = grep("plot_mean[",colnames(mfit),fixed=TRUE)
  params   = mfit[,-pred.cols]
  prow = sample.int(nrow(params),Nmc,replace=TRUE)
  
  # Format covariate data for calibration/forecast period
  weather.fcast <- weather[!grepl("2015",weather$dateID),]
  weather.fcast$dateID <- as.integer(as.factor(weather.fcast$dateID))
  temp.cper <- weather.fcast[weather.fcast$siteID==site.name,]$monthly_temp
  precip.cper <- weather.fcast[weather.fcast$siteID==site.name,]$monthly_precip
  temp.cal <- temp.cper[time1]
  temp.fcast <- temp.cper[time2]
  precip.cal <- precip.cper[time1]
  precip.fcast <- precip.cper[time2]
  pH <- pH.cov[plot]
  
  # Pull out posterior samples
  site_effect <- params[,grep(paste0("site_effect[",site,"]"), colnames(params),fixed=T)]
  beta_IC           <- params[,grep("beta_IC", colnames(params),fixed=T)]
  plot_effect       <- params[,grep(paste0("plot_effect[",plot,"]"), colnames(params),fixed=T)]
  beta_precip       <- params[,grep("beta_precip",colnames(params))]
  beta_temp         <- params[,grep("beta_temp",colnames(params))]
  beta_pH           <- params[,grep("beta_pH",colnames(params))]
  plot_var          <- sqrt(1/params[,grep("plot_var",colnames(params))])
  tau               <- params[,grep("tau",colnames(params))]
  time_effect       <- params[,grep("time_effect",colnames(params))]
  time_var       <- sqrt(1/params[,grep("time_var",colnames(params))])
  random_time_effect <- rnorm(1500, 0, time_var) 

  # Create and plot confidence/predictive intervals for the calibration period.
  plotpreds <- mfit[,grep(paste0("plot_mean[", plot,","), colnames(mfit),fixed=TRUE)]
  plot.cal.ci <- apply(plotpreds,2,quantile,c(0.025,0.5,0.975))
  plot.cal <- sapply(plot.cal.ci[2,], function (x) rbeta(500, mean(tau) * x, mean(tau) * (1 - x)))
  plot.cal.pi <- apply(plot.cal,2,quantile,c(0.025,0.5,0.975))
  
  plot(time.t,time.t,type='n',ylim=c(0,1),ylab="Relative Abundance", xlab = "Time", main = paste(colnames(rank.df)[j], "at plot:", plot.name))
  ecoforecastR::ciEnvelope(time1,plot.cal.pi[1,],plot.cal.pi[3,],col=col.alpha("lightGreen",0.6))
  ecoforecastR::ciEnvelope(time1,plot.cal.ci[1,],plot.cal.ci[3,],col=col.alpha("lightBlue",0.8))
  lines(time1,plot.cal.ci[2,],col="blue")

  # Plot core-level observations
  core_obs <- rank.df[grep("CPER_001", rownames(rank.df)),]
  core_obs$dateID <- substr(core_obs$dates, 1, 7)
  core_obs$dateCol <- as.Date(paste0(core_obs$dateID, "-01"), "%Y-%m-%d")
  new_obs <- plot_out[grep(paste0("plot_mean[",plot,","), rownames(plot_out),fixed=T),]
  new_obs$time1 <- 1:25
  new_obs <- merge(new_obs, core_obs, all=T, by="dateCol")
  points(new_obs$time1, new_obs[,colnames(rank.df)[j]], pch=16)
  
  
  ## FORECAST!! ###
  ypred.fcast <- matrix(NA,nrow=Nmc,ncol=20)
  ycred.fcast <- matrix(NA,nrow=Nmc,ncol=20)
  
  ICprev <- plot.cal.ci[2,25] # initialize
  for(g in 1:ngibbs){
  for(time in 1:20){
    inv.plot_mean <-  beta_IC[g]*ICprev + beta_precip[g]*precip.fcast[time] + beta_temp[g]*temp.fcast[time] + beta_pH[g]*pH + site_effect[g] + plot_effect[g] + random_time_effect[g]
    ycred.fcast[g,time] <- invlogit(inv.plot_mean)
    
    alpha <-  ycred.fcast[g,time]*tau[g]
    beta <- (1- ycred.fcast[g,time])*tau[g]
    
    ypred.fcast[g,time] <- rbeta(1, alpha, beta)
    ICprev <- ycred.fcast[g,time]
  }
  }
  
  fcast.ci <- apply(ycred.fcast,2,quantile,c(0.025,0.5,0.975))
  fcast.pi <- apply(ypred.fcast,2,quantile,c(0.025,0.5,0.975))
  
  ecoforecastR::ciEnvelope(time2,fcast.pi[1,],fcast.pi[3,],col=col.alpha("lightGreen",0.6))
  ecoforecastR::ciEnvelope(time2,fcast.ci[1,],fcast.ci[3,],col=col.alpha("lightBlue",0.8))
  lines(time2,fcast.ci[2,],col="blue")
  
}

# Add labels if necessary before saving plot. 
mtext("a)", side = 3, cex = 2, at = c(0), line = 1)
mtext("b)", side = 3, cex = 2, at = c(0), line = 1)
mtext("c)", side = 3, cex = 2, at = c(0), line = 1)

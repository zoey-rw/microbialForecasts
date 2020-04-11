library(ecoforecastR)
library(runjags)
library(arm)
library(dplyr)
library(tidyr)
library(RColorBrewer)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepModelData.r")


j = 1
## Read in model
mod.list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")
pmod <- mod.list[[j]]

## Read in covariate data 
cov <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
weather <- cov[[1]]
plot_pH <- cov[[2]]

# Read in the microbial abundance data
d <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/groupAbundances_16S_EE509.rds")
rank.df <- d[[1]]
rank.df$dates <- gsub("2016", "2015", rank.df$dates)
#rank.df <- rank.df[grepl("2013", rownames(rank.df)),]


core_obs <- rank.df[grep("CPER_001", rownames(rank.df)),]
core_obs$dateID <- substr(core_obs$dates, 1, 7)
#core_obs$dateID <- gsub("-", "", core_obs$dateID)
core_obs$dateCol <- as.Date(paste0(core_obs$dateID, "-01"), "%Y-%m-%d")

# Create forecasting function
forecastN <- function(IC, beta_IC, plot_effect, site_effect, time_effect,
                      beta_precip, beta_temp, beta_pH, Nmc, tau = NULL, NT=20){
  N <- matrix(NA,Nmc,NT)  ## storage
  Nprev <- IC           ## initialize
  for(time in 1:NT){
    inv.plot_mean <-  beta_IC*Nprev + beta_precip*precip.fcast[time] + beta_temp*temp.fcast[time] + beta_pH*pH + 
      site_effect + plot_effect + time_effect
    plot_mean <- invlogit(inv.plot_mean)
    if (!is.null(tau)) { # add process error
      alpha <- plot_mean*tau
      beta <- (1-plot_mean)*tau
      N[,time] <- rbeta(Nmc, alpha, beta)
    } else {                                  ## update IC
      N[,time] <- plot_mean
    }
    Nprev <- N[,time]
  }
  return(N)
}

plot.run <- function(){
  plot(time.t,time.t,type='n',ylim=c(0,.6),ylab="Relative Abundance", main = paste(colnames(rank.df)[j], "at plot:", plot.name))
  ecoforecastR::ciEnvelope(time1,plot.cal.pi[1,],plot.cal.pi[3,],col=col.alpha("lightGreen",0.6))
  ecoforecastR::ciEnvelope(time1,plot.cal.ci[1,],plot.cal.ci[3,],col=col.alpha("lightBlue",0.8))
  lines(time1,plot.cal.ci[2,],col="blue")
  lines(time1,plot.cal.ci[2,],col="blue")
  #points(time1, obs, pch = 16)
}


# Decide which site/mean to look at 
site <- 1 
site.name <- unique(as.factor(weather$siteID))[site]
plot <- 1
plot.name <- unique(as.factor(plot_pH$plotID))[plot]

# Get prediction data
model.dat <- prepModelData(weather, plot_pH, rank.df, j=1) # Combine using custom function
plot.truth <- model.dat[[7]]
site.truth <- model.dat[[8]]
pH.cov <- model.dat[[5]]
plot_out <- pmod$summary
# Get observed values
obs <- plot_out[grep(paste0("plot_mean[",plot,","), rownames(plot_out),fixed=T),]$truth
new_obs <- plot_out[grep(paste0("plot_mean[",plot,","), rownames(plot_out),fixed=T),]
new_obs$time1 <- 1:25
new_obs <- merge(new_obs, core_obs, all=T, by="dateCol")

# Extract parameter estimates from model
mfit = pmod$JAGS$mcmc[[1]]
pred.cols = grep("plot_mean[",colnames(mfit),fixed=TRUE)
params   = mfit[,-pred.cols]

# Get estimated means for our calibration period
plotpreds <- mfit[,grep(paste0("plot_mean[", plot,","), colnames(mfit),fixed=TRUE)]
plot.cal.ci <- apply(plotpreds,2,quantile,c(0.025,0.5,0.975))
plot.cal <- sapply(plot.cal.ci[2,], function (x) rbeta(500, mean(tau) * x, mean(tau) * (1 - x)))
plot.cal.pi <- apply(plot.cal,2,quantile,c(0.025,0.5,0.975))


# Set characteristics for all plots/forecasts
time.t = 1:45    ## total time
time1 = 1:25       ## calibration period
time2 = 26:45   ## forecast period
N.cols <- brewer.pal(4, "Set2") ## set colors
Nmc = 1000
NT <- 20

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

# Calculate all parameter means
param.mean <- apply(params,2,mean)
prow = sample.int(nrow(params),Nmc,replace=TRUE)

site_effect <- params[prow,grep(paste0("site_effect[",site,"]"), colnames(params),fixed=T)]
beta_IC           <- params[prow,grep("beta_IC", colnames(params),fixed=T)]
plot_effect       <- params[prow,grep(paste0("plot_effect[",plot,"]"), colnames(params),fixed=T)]
beta_precip       <- params[prow,grep("beta_precip",colnames(params))]
beta_temp         <- params[prow,grep("beta_temp",colnames(params))]
beta_pH           <- params[prow,grep("beta_pH",colnames(params))]
plot_var          <- params[prow,grep("plot_var",colnames(params))]
tau               <- params[prow,grep("tau",colnames(params))]
time_var          <- sqrt(1/params[prow,grep("time_var",colnames(params))])
random_time_effect <- rnorm(1, 0, mean(time_var))


##### Deterministic forecast #####
N.det <- forecastN(site_effect = mean(site_effect),
                   IC = mean(plotpreds[,25]),
                   beta_IC = mean(beta_IC),
                   plot_effect = mean(plot_effect),
                   beta_precip = mean(beta_precip),
                   beta_temp = mean(beta_temp),
                   beta_pH = mean(plot_effect),
                   time_effect = 0,
                   NT = 20,
                   Nmc = 1)

## Plot run
plot.run()  
lines(time2,N.det,col="purple",lwd=3)


##### ADD INITIAL CONDITIONS #####

N.I <- forecastN(site_effect = mean(site_effect),
                 IC = plotpreds[prow,25],
                 beta_IC = beta_IC,
                 plot_effect = mean(plot_effect),
                 beta_precip = mean(beta_precip),
                 beta_temp = mean(beta_temp),
                 beta_pH = mean(beta_pH),
                 time_effect=0,
                 NT=20,
                 Nmc = 1000)
## Plot run
plot.run()  
N.I.ci = apply(N.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],.8))
lines(time2,N.I.ci[2,],lwd=0.5)


##### ADD PARAMETER ERROR #####

N.IP <- forecastN(site_effect = site_effect,
                  IC = plotpreds[prow,25],
                  beta_IC = beta_IC,
                  plot_effect = plot_effect,
                  beta_precip = beta_precip,
                  beta_temp = beta_temp,
                  beta_pH = beta_pH,
                  time_effect = 0,
                  Nmc = 1000)

## Plot run
plot.run()  
N.IP.ci = apply(N.IP,2,quantile,c(0.025,0.5,0.975))
ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],.8))
ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],.8))
lines(time2,N.IP.ci[2,],lwd=0.5)


##### ADD PROCESS ERROR #####
N.IPE <- forecastN(site_effect = site_effect,
                   IC = plotpreds[prow,25],
                   beta_IC = beta_IC,
                   plot_effect = plot_effect,
                   beta_precip = beta_precip,
                   beta_temp = beta_temp,
                   beta_pH = beta_pH,
                   time_effect = 0,
                   tau = tau,
                   Nmc = 1000)
## Plot run
plot.run()  
N.IPE.ci = apply(N.IPE,2,quantile,c(0.025,0.5,0.975))
ciEnvelope(time2,N.IPE.ci[1,],N.IPE.ci[3,],col=col.alpha(N.cols[3],.8))
ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],.8))
ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],.8))
lines(time2,N.IP.ci[2,],lwd=0.5)


##### ADD TIME RANDOM-EFFECT UNCERTAINTY #####
N.IPEA <- forecastN(site_effect = site_effect,
                    IC = plotpreds[prow,25],
                    beta_IC = beta_IC,
                    plot_effect = plot_effect,
                    beta_precip = beta_precip,
                    beta_temp = beta_temp,
                    beta_pH = beta_pH,
                    time_effect = random_time_effect,
                    tau = tau,
                    Nmc = 1000)
## Plot run
plot.run()  
N.IPEA.ci = apply(N.IPEA,2,quantile,c(0.025,0.5,0.975))
ciEnvelope(time2,N.IPEA.ci[1,],N.IPEA.ci[3,],col=col.alpha(N.cols[4],.8))
ciEnvelope(time2,N.IPE.ci[1,],N.IPE.ci[3,],col=col.alpha(N.cols[3],.8))
ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],.8))
ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],.8))
lines(time2,N.IP.ci[2,],lwd=0.5)
legend("topleft",legend=c("RandomTimeEffect","Process","Parameter","InitCond"),col=rev(N.cols[-5]),lty=1,lwd=5)

mtext("a)", side = 3, cex = 2, at = c(0), line = 1)
##### Partition uncertainty #####

### calculation of variances
varI     <- apply(N.I,2,var)
varIP    <- apply(N.IP,2,var)
varIPE  <- apply(N.IPE,2,var)
varIPEA <- apply(N.IPEA,2,var)
varMat   <- rbind(varI,varIP,varIPE,varIPEA)

## in-sample stacked area plot
V.pred.rel.in <- apply(varMat,2,function(x) {x/max(x)})
plot(time2,V.pred.rel.in[1,],ylim=c(0,1),type='n',main="Relative Variance: In-Plot",ylab="Proportion of Variance",xlab="time")
ciEnvelope(time2,rep(0,ncol(V.pred.rel.in)),V.pred.rel.in[1,],col=N.cols[1])
ciEnvelope(time2,V.pred.rel.in[1,],V.pred.rel.in[2,],col=N.cols[2])
ciEnvelope(time2,V.pred.rel.in[2,],V.pred.rel.in[3,],col=N.cols[3])
ciEnvelope(time2,V.pred.rel.in[3,],V.pred.rel.in[4,],col=N.cols[4])
legend("topleft",legend=c("RandomTimeEffect","Process","Parameter","InitCond"),col=rev(N.cols[-5]),lty=1,lwd=5)

mtext("b)", side = 3, cex = 2, line = 1, at = c(26,1.6))

# Create representative forecast figures
library(ecoforecastR)
library(egg)
library(runjags)
library(arm)
library(dplyr)
library(tidyr)
library(RColorBrewer)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepModelData.r")


## Read in model
mod.list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")

#p.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")
p.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS1-5.rds")
# names(p.bac)[6:10] <- names(p.bac)[1:5]
# p.bac <- p.bac[6:10)
c.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/class_bacJAGS.rds")
o.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/order_bacJAGS.rds")
f.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/family_bacJAGS.rds")
g.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/genus_bacJAGS.rds")
all.ranks.bac <- list(p.bac, c.bac, o.bac, f.bac, g.bac)



#p.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS.rds")
p.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS6-10.rds")
names(p.fun)[6:10] <- names(p.fun)[1:5]
p.fun <- p.fun[6:10]
c.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/class_funJAGS.rds")
o.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/order_funJAGS.rds")
f.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/family_funJAGS.rds")
g.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/genus_funJAGS.rds")
all.ranks.fun <- list(p.fun, c.fun, o.fun, f.fun, g.fun)

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


# Set characteristics for all plots/forecasts
time.t = 1:45    ## total time
time1 = 1:25       ## calibration period
time2 = 26:45   ## forecast period
N.cols <- brewer.pal(4, "Set2") ## set colors
Nmc = 500
NT <- 20
ngibbs <- 1000
npred <- 25


pl.out <- list()
for (k in 1:5) {
  
  mod.list <-  all.ranks.bac[[k]]
  #mod.list <-  all.ranks.fun[[k]]
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

# Get estimated means for our calibration period
plotpreds <- mfit[,grep(paste0("plot_mean[", plot,","), colnames(mfit),fixed=TRUE)]
plot.cal.ci <- apply(plotpreds,2,quantile,c(0.025,0.5,0.975))

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
random_time_effect <- dnorm(1, 0, plot_var)

ypred <- matrix(NA,nrow=ngibbs,ncol=npred)
ycred <- matrix(NA,nrow=ngibbs,ncol=npred)

ICprev <- plot.cal.ci[2,]
ICprev <- c(.25, ICprev[1:24])
for(g in 1:ngibbs){
  inv.plot_mean <-  beta_IC[g]*ICprev + beta_precip[g]*precip.cal + beta_temp[g]*temp.cal + beta_pH[g]*pH + 
    site_effect[g] + plot_effect[g] + time_effect[g,]
  ycred[g,] <- invlogit(inv.plot_mean)
    alpha <-  ycred[g,]*tau[g]
    beta <- (1- ycred[g,])*tau[g]
    for (t in 1:25){
      ypred[g,t] <- rbeta(1, alpha[t], beta[t])
    }
}

yobs <- obs
keep <- which(!is.na(yobs))
npred.real <- length(keep)
## Residual variance
ybar <- apply(ycred[,keep],2,mean)
G <- sum((yobs[keep]-ybar)^2, na.rm=T)/npred.real
## Predictive variance
P <- sum(apply(ypred,2,var))/npred
Dpl <- G + P
PL <- c(G,P,Dpl)
names(PL) <- c("Residual variance","Predictive variance","Total variance")
pl.rank[[j]] <- PL
}
pl.out[[k]] <- pl.rank
}

p <- cbind.data.frame(do.call(rbind,pl.out[[1]]), rank = rep("phylum", 10))
c <- cbind.data.frame(do.call(rbind,pl.out[[2]]), rank = rep("class", 5))
o <- cbind.data.frame(do.call(rbind,pl.out[[3]]), rank = rep("order", 5))
f <- cbind.data.frame(do.call(rbind,pl.out[[4]]), rank = rep("family", 5))
g <- cbind.data.frame(do.call(rbind,pl.out[[5]]), rank = rep("genus",  5))

pl.sum.bac <- rbind(p,c,o,f,g)
pl.sum.bac$cat <- "bacteria"

pl.sum.fun <- rbind(p,c,o,f,g)
pl.sum.fun$cat <- "fungi"
#pl.sum <- pl.sum %>% group_by(rank) %>% summarize_all(mean)

pl.sum <- rbind(pl.sum.bac, pl.sum.fun)

n_fun <- function(x){
  return(data.frame(y = -2.302585,
                    label = length(x)))
}

m <- ggplot(pl.sum, aes(x=rank, y=`Total variance`)) + 
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=.2) + 
  facet_grid(. ~ cat) + stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  scale_y_continuous(trans='log', limits = c(NA,0.1)) + xlab("Taxonomic rank") + ylab("log(Total predicive loss)")

dat_text <- data.frame(label = c("a)","b)"), 
                       cat = c("bacteria","fungi"))
m <- m + geom_text(
  data    = dat_text,
  size = 5,
  mapping = aes(x = 1, y = .1, label = label)
)
m
ggsave(filename = "Figure4.png", device="png")

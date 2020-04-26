library(rjags)
library(runjags)
library(dplyr)
library(arm) #for invlogit function
library(tidyr)
library(ggplot2)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepModelData2.r")

# Read in covariate data 
cov <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
weather <- cov[[1]]
plot.preds <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")[[1]]

# Read in the microbial abundance data
d <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/groupAbundances_EE509.rds")

out.allrank <- list()

#for (k in 1:length(d)){
  for (k in c(1,6)){
    
rank.df <- d[[k]]
rank.df$dates <- gsub("2016", "2015", rank.df$dates)
#rank.df <- rank.df[grepl("2013", rownames(rank.df)),]
#n.groups <- ifelse(grepl("fun|bac", names(d)[[k]]), 10, 1)
out <- list()

for (j in 1:10){ # j specifies the group (e.g. Acidobacteria)

model.dat <- prepModelData(weather, plot.preds, rank.df, j=j) # Combine using custom function
dat <- model.dat[[1]]
y <- model.dat[[2]]
time.cov <- model.dat[[3]]
space.cov <- model.dat[[4]]
plot_site <- model.dat[[5]]
plot.truth <- model.dat[[6]]

# Create data object.
data <- list(y = as.matrix(y), 
             N.core = length(dat$y), 
             N.plot = length(unique(dat$plotID)),
             N.site = length(unique(dat$siteID)), 
             N.date = length(unique(dat$dateID)),
             plot_site = droplevels(as.factor(plot_site)),
             plotID = droplevels(as.factor(dat$plotID)),
             time.cov = time.cov,
             space.cov = space.cov
)


# FULL MODEL
model_string <- " model{

#### Core-level observations ####
for (t in 1:N.date){
  for (i in 1:N.core){
   y[i,t] ~ dbeta(shape1[plotID[i],t], shape2[plotID[i],t]) ## Data model: core observations are a function of alpha and beta
  }
}

# Plot means - Process model
for (i in 1:N.plot){
  for (t in 2:N.date){
    logit(plot_mean[i,t]) <-  beta[1]*logit(plot_mean[i,t-1]) + 
                              beta[2]*time.cov[i,t,2] + beta[3]*time.cov[i,t,1] + 
                              inprod(beta[4:6], space.cov[i,1:3]) +
                              plot_effect[i] + time_effect[plot_site[i],t]

    # Convert plot mean into parameters for beta distribution
    shape1[i,t] <- plot_mean[i,t] * tau   
    shape2[i,t]  <- (1-plot_mean[i,t]) * tau 
  }
}

#### Priors ####

for (i in 1:N.site){
  site_effect[i] ~ dnorm(0,site_var)  # Prior on site random effects
  for (t in 1:N.date){
    time_effect[i,t] ~ dnorm(0, time_var) # Prior on time random effects
  }
}

for (i in 1:N.plot){
  plot_effect[i] ~ dnorm(site_effect[plot_site[i]], plot_var) # Prior on plot random effects
  plot_mean[i,1] ~ dbeta(1,3)             # Prior on plot means for first date 
  shape1[i,1] <- plot_mean[i,1] * tau      # alpha for first date
  shape2[i,1]  <- (1-plot_mean[i,1]) * tau  # beta for first date
}

# Coefficient priors
for(i in 1:6){
    beta[i] ~ dnorm(0,0.0001)
}

# Variance priors
tau ~ dgamma(1,.1) 
plot_var ~ dgamma(1,.1) 
site_var ~ dgamma(1,.1) 
time_var ~ dgamma(1,.1) 
}"

monitor <- c("beta","tau","site_var","time_var",
             "plot_var","plot_mean","site_effect","plot_effect","time_effect")

tic()
pmod <- run.jags(model_string,
                 data = data,
                 adapt = 2000,
                 burnin = 1000,
                 sample = 1500,
                 n.chains = 3,
                 thin = 5,
                 #method = "rjparallel",
                 method = "parallel",
                 jags = "/share/pkg.7/jags/4.3.0/install/bin/jags",
                 monitor = monitor, silent.jags = F)
toc()

tic()
# # Pull out plot values
plot_mean <- summary(pmod, vars="plot_mean")
toc()
plot_out <- data.frame(plot_int_mean = plot_mean[,4],
                       plot_int_lo95 = plot_mean[,1],
                       plot_int_hi95 = plot_mean[,3])
plot_out$truth <- plot.truth$plot_mean
plot_out$plotID <- plot.truth$plotID
plot_out$dateCol <- as.Date(paste0(as.character(plot.truth$dateID), '01'), format='%Y%m%d')
plot_out$siteID <- plot.truth$siteID
plot_out$psrf <- plot_mean[,11]
plot_out$SSeff <- plot_mean[,9]

out[[j]] <- list()
out[[j]][[1]] <- plot_out
out[[j]][[2]] <- pmod
saveRDS(out[[j]], paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/",names(d[[k]])[j],"JAGS.rds"))
cat(paste("Model fit for group:", names(d[[k]])[j]))

tic()
out[[j]][[3]] <- "DIC placeholder" #extract(pmod, "dic", n.iter = 500)
names(out[[j]]) <- c("summary", "JAGS", "DIC")
toc()
}
#names(out) <- names(rank.df)[1:n.groups]
names(out) <- names(rank.df)[1:10]
saveRDS(out, paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/",names(d)[[k]],"JAGS1-10.rds"))
out.allrank[[k]] <- out
}
# 
# names(out.allrank) <- names(d)
# saveRDS(out.allrank, "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/allGroupsJAGS.rds")




# # combine two phylum-level runs
# bac1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS1-5.rds")
# bac2 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS6-10.rds")
# names(bac2)[6:10] <- names(bac2)[1:5]
# bac2 <- bac2[6:10]
# bac <- c(bac1, bac2)
# saveRDS(bac, paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds"))
# 
# 
# fun1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS1-5.rds")
# fun2 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS6-10.rds")
# names(fun2)[6:10] <- names(fun2)[1:5]
# fun2 <- fun2[6:10]
# fun <- c(fun1, fun2)
# saveRDS(fun, paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS.rds"))


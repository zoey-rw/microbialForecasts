library(tidyr)
library(nimble)
library(dplyr)
library(reshape2)
library(ggpubr)
library(coda)
library(runjags)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepModelData_clr.r")
weather <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/temporal_covariates.rds")
plot.preds <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilSpatialCov.rds")[[1]]

fixDate <- function(datesToFix){
  as.Date(paste0(datesToFix, "01"), format = "%Y%m%d")
}

plot_model_new <- function(samples_var = samples, model.dat = model.dat,
                           val_model.dat = NULL,
                           group.name = "Acidobacteria", plot.no = 1){
  
  cal.plot.truth <- model.dat$plot.truth
  alldates <- fixDate(colnames(cal.plot.truth))
  plot_key <- cbind.data.frame(plotID = rownames(cal.plot.truth),
                               plot_num = 1:nrow(cal.plot.truth))
  plot.name <- plot_key[plot_key$plot_num==plot.no, ]$plotID
  
  # Extract parameter estimates from model
  if ("WAIC" %in% names(samples_var)) {
    samples_var <- samples_var$samples
  }
  
  all_samp <- do.call(rbind, samples_var)
  ypred <- all_samp[,grep("plot_mean_hat[",colnames(all_samp),fixed=TRUE)]
  ypred_plot <- ypred[,grep(paste0("plot_mean_hat[", plot.no,","), colnames(ypred), fixed=TRUE)]
  ci <- apply(ypred_plot,2,quantile,c(0.025,0.5,0.975))
  cal.observed <- melt(cal.plot.truth[grep(plot.name, rownames(cal.plot.truth)),])
  cal.observed$time <- fixDate(rownames(cal.observed))
  to_plot <- cbind.data.frame(alldates, t(ci), group.name)
  
  p <- ggplot(to_plot, aes(x = alldates)) + 
    geom_ribbon(data = to_plot, aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2) +
    geom_point(data = cal.observed, aes(y = value,  x = time)) + 
    labs(x = "", y = "Abundance (CLR-transformed)", title = paste0(group.name, " at plot: ", plot.name)) +
    theme_classic(base_size = 26)
  
  if (!is.null(val_model.dat)){
    val.plot.truth <- val_model.dat$plot.truth
    val.observed <- melt(val.plot.truth[grep(plot.name, rownames(val.plot.truth)),])
    val.observed$time <- fixDate(rownames(val.observed))
    p <- p + geom_point(data = val.observed, aes(y = value,  x = time), color = "red")
  }
  p
}


create_plotting_data <- function(samples_var = samples, model.dat = model.dat,
                           val_model.dat = NULL,
                           group.name = "Acidobacteria", plot.no = 1){
  
  cal.plot.truth <- model.dat$plot.truth
  plot_key <- cbind.data.frame(plotID = rownames(cal.plot.truth),
                               plot_num = 1:nrow(cal.plot.truth))
  plot.name <- plot_key[plot_key$plot_num==plot.no, ]$plotID
  
  # subset to the correct start date for the site
  site.name <- substr(plot.name, 1, 4)
  # site.start <- model.dat$site_start[which(model.dat$site_start$siteID==site.name),]$start
  # cal.plot.truth <- cal.plot.truth[, site.start:ncol(cal.plot.truth)]
  alldates <- fixDate(colnames(cal.plot.truth))
  
  
  # Extract parameter estimates from model
  if ("WAIC" %in% names(samples_var)) {
    samples_var <- samples_var$samples
  }
  
  all_samp <- do.call(rbind, samples_var)
  # ypred <- all_samp[,grep("plot_mean_hat[",colnames(all_samp),fixed=TRUE)]
  # ypred_plot <- ypred[,grep(paste0("plot_mean_hat[", plot.no,","), colnames(ypred), fixed=TRUE)]
  ypred <- all_samp[,grep("plot_mean[",colnames(all_samp),fixed=TRUE)]
  ypred_plot <- ypred[,grep(paste0("plot_mean[", plot.no,","), colnames(ypred), fixed=TRUE)]
  ci <- apply(ypred_plot,2,quantile,c(0.025,0.5,0.975), na.rm = TRUE)
  cal.observed <- melt(cal.plot.truth[grep(plot.name, rownames(cal.plot.truth)),])
  cal.observed$time <- fixDate(rownames(cal.observed))
  to_plot <- cbind.data.frame(alldates, t(ci), cal.observed, data.type = "Calibration")

  if (!is.null(val_model.dat)){
    val.plot.truth <- val_model.dat$plot.truth
    val.observed <- melt(val.plot.truth[grep(plot.name, rownames(val.plot.truth)),])
    if (nrow(val.observed) == 0) {cat(paste("No validation points for plot:", plot.name)); return(NA)}
    val.observed$time <- fixDate(rownames(val.observed))
    val.observed <- cbind.data.frame(val.observed,data.type = "Validation")
    out <- merge(to_plot, val.observed, all=T, by= c("time"))
    out$data.type <- ifelse(is.na(out$data.type.y), as.character(out$data.type.x), as.character(out$data.type.y))
    out$value <- ifelse(!is.na(out$value.y), out$value.y, out$value.x)
    out <- out[,!colnames(out) %in% c("value.x","value.y", "data.type.x", "data.type.y")]
  } else {
    out <- to_plot
  }
  out$group.name <- group.name
    out$plot.name <- plot.name
  return(out)
}




create_plotting_data_y <- function(samples_var = samples, model.dat = model.dat,
                                 val_model.dat = NULL,
                                 group.name = "Acidobacteria", plot.no = 1){
  
  cal.plot.truth <- model.dat$plot.truth
  alldates <- fixDate(colnames(cal.plot.truth))
  plot_key <- cbind.data.frame(plotID = rownames(cal.plot.truth),
                               plot_num = 1:nrow(cal.plot.truth))
  plot.name <- plot_key[plot_key$plot_num==plot.no, ]$plotID
  date_key <- cbind.data.frame(dateID = colnames(cal.plot.truth),
                               date_num = 1:ncol(cal.plot.truth))
  
  # Extract parameter estimates from model
  if ("WAIC" %in% names(samples_var)) {
    samples_var <- samples_var$samples
  }
  all_samp <- do.call(rbind, samples_var)
  
  y.obs <- model.dat$y
  dimnames(y.obs) <- NULL
  cal.observed.y <- melt(y.obs)
  cal.observed.y$coreID <- paste0("y[",cal.observed.y[,1], ", ", cal.observed.y[,2], ", ", cal.observed.y[,3], "]")
  cal.observed.y <- cal.observed.y[which(cal.observed.y$Var3 == plot.no),]
  cal.observed.y$time <- fixDate(date_key[match(cal.observed.y$Var2, date_key$date_num),]$dateID)
  cal.melt <- cal.observed.y  %>% select(-coreID) %>% pivot_wider( values_from = value, names_from = Var1)
  
  all_y <- all_samp[,grep("y[",colnames(all_samp),fixed=TRUE)]
  ypred <- all_y[,grep(paste0(", ", plot.no, "]"),colnames(all_y),fixed=TRUE)]
  ycore1 <- ypred[,grep("y[1, ", colnames(ypred),fixed=TRUE)]
  ycore2 <- ypred[,grep("y[2, ", colnames(ypred),fixed=TRUE)]
  ycore3 <- ypred[,grep("y[3, ", colnames(ypred),fixed=TRUE)]
  ycore1_ci <- apply(ycore1,2,quantile,c(0.025,0.5,0.975))
  ycore2_ci <- apply(ycore2,2,quantile,c(0.025,0.5,0.975))
  ycore3_ci <- apply(ycore3,2,quantile,c(0.025,0.5,0.975))
  est_y <- cbind.data.frame(t(ycore1_ci),  t(ycore2_ci),  t(ycore3_ci))
  colnames(est_y) <- c("ycore1_2.5", "ycore1_50", "ycore1_97.5", "ycore2_2.5", "ycore2_50", "ycore2_97.5", "ycore3_2.5", "ycore3_50", "ycore3_97.5")
  est_y$coreID <- rownames(est_y)
  est_y$time <- fixDate(date_key$dateID) 
  y_to_plot <- merge(est_y, cal.melt, all.x=T)
  y_to_plot$time <- fixDate(date_key[match(y_to_plot$Var2, date_key$date_num),]$dateID)
  y_to_plot$data.type = "Calibration"
  
  cal.observed.y <- melt(y.obs)
  cal.observed.y$coreID <- paste0("mu[",cal.observed.y[,1], ", ", cal.observed.y[,2], ", ", cal.observed.y[,3], "]")
  cal.observed.y <- cal.observed.y[which(cal.observed.y$Var3 == plot.no),]
  cal.observed.y$time <- fixDate(date_key[match(cal.observed.y$Var2, date_key$date_num),]$dateID)
  cal.melt <- cal.observed.y  %>% select(-coreID) %>% pivot_wider( values_from = value, names_from = Var1)
  all_y <- all_samp[,grep("mu[",colnames(all_samp),fixed=TRUE)]
  ypred <- all_y[,grep(paste0(", ", plot.no, "]"),colnames(all_y),fixed=TRUE)]
  core1 <- ypred[,grep("mu[1, ", colnames(ypred),fixed=TRUE)]
  core2 <- ypred[,grep("mu[2, ", colnames(ypred),fixed=TRUE)]
  core3 <- ypred[,grep("mu[3, ", colnames(ypred),fixed=TRUE)]
  core1_ci <- apply(core1,2,quantile,c(0.025,0.5,0.975))
  core2_ci <- apply(core2,2,quantile,c(0.025,0.5,0.975))
  core3_ci <- apply(core3,2,quantile,c(0.025,0.5,0.975))
  est <- cbind.data.frame(t(core1_ci),  t(core2_ci),  t(core3_ci))
  colnames(est) <- c("core1_2.5", "core1_50", "core1_97.5", "core2_2.5", "core2_50", "core2_97.5", "core3_2.5", "core3_50", "core3_97.5")
  est$time <- fixDate(date_key$dateID) 
  to_plot <- merge(est, cal.melt, all.x=T)
  to_plot$data.type = "Calibration"
  
  
  if (!is.null(val_model.dat)){
    date_key <- cbind.data.frame(dateID = colnames(val_model.dat$plot.truth),
                                 date_num = 1:ncol(val_model.dat$plot.truth))
    val.y <- val_model.dat$y
    dimnames(val.y) <- NULL
    val.observed.y <- melt(val.y)
    val.observed.y$coreID <- paste0("y[",val.observed.y[,1], ", ", val.observed.y[,2], ", ", val.observed.y[,3], "]")
    val.observed.y <- val.observed.y[which(val.observed.y$Var3 == plot.no),]
    val.observed.y$time <- fixDate(date_key[match(val.observed.y$Var2, date_key$date_num),]$dateID)
    val.melt <- val.observed.y  %>% select(-coreID) %>% pivot_wider( values_from = value, names_from = Var1)
    if (nrow(val.observed) == 0) {cat(paste("No validation points for plot:", plot.name)); return(NA)}
    val.observed$time <- fixDate(rownames(val.observed))
    val.observed <- cbind.data.frame(val.observed,data.type = "Validation")
    val_out <- merge(to_plot, val.observed, all=T, by= c("time"))

  
  } else {
    out <- to_plot
  }
  out$group.name <- group.name
  out$plot.name <- plot.name
  return(out)
}



gel <- function(samples){
  if ("WAIC" %in% names(samples)) {
    samples <- samples$samples
  }
  allvars <- c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","alpha","plot_var","time_var","site_var","tau_proc","tau_obs")
  #time_effs <- grep("time_effect", colnames(samples[[1]]), fixed=T)
  vars <- allvars[allvars %in% colnames(samples[[1]])]
  coda::gelman.diag(samples[,vars])
}

var_list <- c("tau_obs","tau_proc","plot_var","site_var", "plot_effect", "time_var", "beta", "site_effect","plot_mean_hat","alpha")
var_listSimple <- c("tau_obs","tau_proc","plot_var", "plot_effect",#"glob_mean",
                    #"time_var", 
                    "plot_mean",
                    "beta","plot_mean_hat","alpha")

# should finish developing this function sometime
oldConstantsFunction <- function(model.dat, #covariates = c("pH","cn","avgprecip","mintemp"), 
                              ...){
  y <- model.dat$y
  constants <-  list(N.plot = dim(y)[3],
                     N.site = length(unique(model.dat$siteID)),
                     N.date = dim(y)[2],
                     siteID = as.numeric(model.dat$siteID))
  # time_covs <- covariates[covariates %in% c("avgtemp","avgprecip","mintemp","maxtemp","precip")]
  #space_covs <-  covariates[covariates %in% c("pH","cn","pC")]                    
  constants$precip = model.dat$time.cov[,colnames(y[,,1]),"mois"]
  constants$temp = model.dat$time.cov[,colnames(y[,,1]),"temp"]
  constants$pH = model.dat$space.cov[,"pH"]
  constants$pC = model.dat$space.cov[,"pC"]
  return(constants)
}

constantsFunction <- function(model.dat, #covariates = c("pH","cn","avgprecip","mintemp"), 
															...){
	y <- model.dat$y
	constants <-  list(N.plot = dim(y)[3],
										 N.site = length(unique(model.dat$siteID)),
										 N.date = dim(y)[2],
										 siteID = as.numeric(model.dat$siteID))
	#constants$nspp = model.dat[["nspp"]]
	constants$site_start = model.dat[["site_start"]][,"start"]
	constants$site_index = model.dat[["site_start"]][,"index"] 
	return(constants)
}

initsFunction <- function(constants = NULL, nbeta = 5){    
  plot_init <- matrix(rnorm(constants$N.plot * constants$N.date, 1, .3), 
                      nrow = constants$N.plot, ncol = constants$N.date)
  list(beta = rnorm(nbeta,0,1), 
       # site_effect =  rnorm(constants$N.site,0,2),
       # plot_effect = rnorm(constants$N.plot),
       # time_effect =  matrix(rnorm(constants$N.site * constants$N.date), 
       #                       nrow = constants$N.site, ncol = constants$N.date),
       tau_obs = rgamma(1, 1, .1), 
       tau_proc = rgamma(1, 1, .1), 
       # plot_var = rgamma(1, 1, .1), 
       # time_var = rgamma(1, 1, .1),
       # site_var = rgamma(1, 1, .1),
       # glob_mean = rnorm(1,0,1),
  		 y = array(rnorm(constants$N.plot * constants$N.date * 3, 1, .3), 
  		 			 dim = c(3, constants$N.date, constants$N.plot)),
       plot_mean = plot_init,
       plot_mean_hat = plot_init,
       init_tau = rgamma(1, 1, .1),
       init_plot = rnorm(1,1,.3),
       alpha = rnorm(1,0,1))
}

### CONVERT PRECISION TO SD
prec_to_sd <- function(x) 1/sqrt(x)


### CONVERT PRECISION TO SD for all samples 
var_to_sd <- function(samples, var.list = c("tau_proc", "tau_obs"), mean = TRUE){
  if (!is.mcmc(samples) && !is.mcmc.list(samples)) {
    cat("Samples must be MCMC or MCMC.list objects"); stop()
  }
  out.list <- list()
  for (i in 1:length(var.list)){
    varname <- var.list[[i]]
    tau_samples <- samples[,varname]
    samples_sd <- unlist(lapply(tau_samples, function(y) lapply(y, function(x) 1/sqrt(x))))
    samples_mean_sd <- mean(samples_sd)
    out.list[[i]] <- samples_mean_sd
  }
  names(out.list) <- var.list
  # 
  # tau_samples_mean <-   summary(tau_samples)[[1]][1]
  # samples_mean_sd <- 1/sqrt(tau_samples_mean)
  return(out.list)
}


plot_betas <- function(samples_var = samples){
  samples <- samples_var$samples
  beta_names <- colnames(samples[[1]])[grep("beta", colnames(samples[[1]]))]
  par(ask = TRUE)
  plot(samples[,beta_names])
}

traceplots <- function(samples = samples, var = "betas"){
  if ("WAIC" %in% names(samples)) {
    samples <- samples$samples
  }
  if (var == "betas"){
    to_plot <- colnames(samples[[1]])[grep("beta", colnames(samples[[1]]))]
  } else if (var == "plot_effect"){
    to_plot <- colnames(samples[[1]])[grep("plot_effect", colnames(samples[[1]]))]
  } else if (var == "alpha"){
    to_plot <- colnames(samples[[1]])[grep("alpha", colnames(samples[[1]]))]
  } else if (var == "tau"){
    to_plot <- colnames(samples[[1]])[grep("tau", colnames(samples[[1]]))]
  } else if (var == "site_effect"){
    to_plot <- colnames(samples[[1]])[grep("site_effect", colnames(samples[[1]]))]
  } else {
    to_plot <- colnames(samples[[1]])[grep(var, colnames(samples[[1]]), fixed=T)]
  }
  par(ask = TRUE)
  plot(samples[,to_plot])
}



# samples_var = samples;
# variable.name = "plot_mean_hat";
# group_name = "Acidobacteria";
# plot.name = "CPER_001";
# plot.no = 1;
# time.cal = 1:63; ## calibration period
# time.val = 64:75; ## forecast period
# start.date = "2013-06-01";
plot_model <- function(samples_var = samples, variable.name = "plot_mean_hat", group_name = "Acidobacteria",
                       plot.name = "CPER_001", plot.no = 1,
                       time.cal = 1:63, ## calibration period
                       time.val = 64:75, ## forecast period
                       start.date = "2013-06-01", ylim=NULL){
  
  ncal <- length(time.cal)
  npred <- length(time.val)
  time.t = c(time.cal, time.val)    ## total time
  
  # Set characteristics for all plots/forecasts
  alldates <- seq(as.Date(start.date), length.out = length(time.t),
                  by = "1 month")
  time1 = alldates[time.cal]
  time2 = alldates[time.val]
  
  # Extract parameter estimates from model
  if ("WAIC" %in% names(samples_var)) {
    samples_var <- samples_var$samples
  }
  sum <- summary(samples_var)
  quartiles <- sum[[2]]
  means = sum[[1]]
  est = quartiles[grep(paste0(variable.name, "[", plot.no,","),rownames(quartiles),fixed=TRUE),]
  obs_error = means[grep("tau_obs", rownames(means)),][1]
  
  
  all_samp <- do.call(rbind, samples_var)
  pred <- all_samp[,grep(paste0(variable.name, "[", plot.no,","),colnames(all_samp),fixed=TRUE)]
  ci <- apply(pred,2,quantile,c(0.025,0.5,0.975))
  
  prow = sample.int(nrow(pred),3000,replace=TRUE)
  
  plot.cal <- matrix(NA, nrow = 3000, ncol = ncal)
  for (i in 1:3000) {
    
    plot.cal[i,] <- sapply(pred[prow[i],], function (x) rnorm(1, x, sqrt(1/obs_error)))
  }
  pi <- apply(plot.cal,2,quantile,c(0.025,0.5,0.975))
  
  if (is.null(ylim)) {
    ylim = c(min(pi[1,]) - 1, max(pi[3,]) + 1)
    # if(abs(ylim[2] - ylim[1]) > 10) ylim = c(-6, 6)
  }
  plot(alldates, time.t, type='n',
       ylim=ylim,
       ylab="Relative Abundance", xlab = NA, main = paste(group_name, "at plot:", plot.name))
  ecoforecastR::ciEnvelope(time1,pi[1,],pi[3,],col=ecoforecastR::col.alpha("lightGreen",0.6))
  ecoforecastR::ciEnvelope(time1,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.8))
  lines(time1,ci[2,],col="blue")
}

# tic/toc functions from Colin Averill:
#' Two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
tic <- function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc <- function() print(Sys.time()-timer)


# create sample information data.frame from NEON sample names
parseNEONsampleIDs <- function(sampleID){
  df <- data.frame(siteID = substr(sampleID, 1, 4), sampleID = sampleID, stringsAsFactors = F) %>% 
  	mutate(sample = sapply(strsplit(sampleID, "-GEN|-gen"),  "[[" , 1)) %>% 
  	mutate(geneticSampleID = sapply(strsplit(sampleID, "-DNA"),  "[[" , 1)) %>% 
  	mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>% 
    mutate(dates = sapply(strsplit(sample, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>% 
    mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>% 
    mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>% 
    mutate(dateID = substr(as.character(dates), 1, 6)) %>% 
    mutate(plotID = substr(sample, 1, 8)) %>% 
    mutate(site_date = paste0(siteID, "-", dateID)) %>% 
    mutate(horizon = ifelse(grepl("-M-", sample), "M", "O")) %>% 
    mutate(without_horizon = gsub("-[M|O]-", "-", sample)) %>% 
    mutate(plot_date = paste0(plotID, "-", dateID)) %>% 
    as.data.frame()
  rownames(df) <- make.unique(sampleID)
  return(df)
}



nimbleModSimple <- nimble::nimbleCode({ 
  
  # Process model
  for (i in 1:N.plot){
    
    # Plot priors
    plot_mean[i,1] ~ dnorm(init_plot, init_tau)  # Plot means for first date
    plot_mean_hat[i,1] ~ dnorm(plot_mean[i,1], tau_proc)  # Plot mean (w/ process error) for first date
    #plot_effect[i] ~ dnorm(site_effect[siteID[i]], plot_var) # Plot random effects, centered on site effects
    #plot_effect[i] ~ dnorm(0, plot_var) # Plot random effects, centered on zero
    plot_effect[i] ~ dnorm(0, plot_var) # Plot random effects, centered on zero
    
    for (t in 2:N.date){ 
      # Dynamic linear model
      plot_mean[i,t] <- alpha + beta[1]*plot_mean[i,t-1] + 
        beta[2]*precip[i,t-1] + beta[3]*temp[i,t] + beta[4]*pH[i] + beta[5]*CN[i] + # Environmental predictors
        plot_effect[i] #+          # Spatial effects 
      #time_effect[siteID[i],t]  # Temporal effects
      # Add process error
      plot_mean_hat[i,t] ~ dnorm(plot_mean[i,t], tau_proc)
    }
  }
  
  # Observation model
  for (i in 1:N.plot){
    for (t in 1:N.date){ 
      for (c in 1:3){ # max 3 cores per plot
        y[c,t,i] ~ dnorm(plot_mean_hat[i,t], tau_obs) 
      }
    }
  }
  
  # Priors: random effects
  # for (i in 1:N.site){
  #   site_effect[i] ~ dnorm(glob_mean, site_var)  # Prior on site random effects
  # for (t in 1:N.date){
  #   time_effect[i,t] ~ dnorm(0, time_var) # Prior on time random effects
  # }
  #}
  
  # Priors: variance
  tau_obs ~ dgamma(1,.1) 
  tau_proc ~ dgamma(1,.1) 
  plot_var ~ dgamma(1,.1) 
  # site_var ~ dgamma(1,.1) 
  #  time_var ~ dgamma(1,.1) 
  init_tau ~ dgamma(1,.1) 
  init_plot ~ dnorm(0,0.1)
  alpha ~ dnorm(0,0.1)
  # glob_mean ~ dnorm(0, .1)
  
  # Priors: coefficients
  for(i in 1:5){
    beta[i] ~ dnorm(0,0.1)
  }
})





nimbleMod <- nimbleCode({ 
  
  # Process model
  for (i in 1:N.plot){
    
    # Plot priors
    plot_mean[i,1] ~ dnorm(init_plot, init_tau)  # Plot means for first date
    plot_mean_hat[i,1] ~ dnorm(plot_mean[i,1], tau_proc)  # Plot mean (w/ process error) for first date
    plot_effect[i] ~ dnorm(site_effect[siteID[i]], plot_var) # Plot random effects, centered on site effects
    
    for (t in 2:N.date){ 
      # Dynamic linear model
      plot_mean[i,t] <- alpha + beta[1]*plot_mean[i,t-1] + 
        beta[2]*precip[i,t] + beta[3]*temp[i,t] + beta[4]*pH[i] + beta[5]*CN[i] + # Environmental predictors
        plot_effect[i] +          # Spatial effects 
        time_effect[siteID[i],t]  # Temporal effects
      # Add process error
      plot_mean_hat[i,t] ~ dnorm(plot_mean[i,t], tau_proc)
    }
  }
  
  # Observation model
  for (i in 1:N.plot){
    for (t in 1:N.date){ 
      for (c in 1:3){ # max 3 cores per plot
        y[c,t,i] ~ dnorm(plot_mean_hat[i,t], tau_obs) 
      }
    }
  }
  
  # Priors: random effects
  for (i in 1:N.site){
    site_effect[i] ~ dnorm(glob_mean, site_var)  # Prior on site random effects
    for (t in 1:N.date){
      time_effect[i,t] ~ dnorm(0, time_var) # Prior on time random effects
    }
  }
  
  # Priors: variance
  tau_obs ~ dgamma(1,.1) 
  tau_proc ~ dgamma(1,.1) 
  plot_var ~ dgamma(1,.1) 
  site_var ~ dgamma(1,.1) 
  time_var ~ dgamma(1,.1) 
  init_tau ~ dgamma(1,.1) 
  init_plot ~ dnorm(0,0.1)
  alpha ~ dnorm(0,0.1)
  glob_mean ~ dnorm(0, .1)
  
  # Priors: coefficients
  for(i in 1:5){
    beta[i] ~ dnorm(0,0.1)
  }
})



spatioTemporalMod <- nimbleCode({ 
  # Process model
  for (i in 1:N.plot){
    
    # Plot priors
    plot_mean[i,1] ~ dnorm(init_plot, init_tau)  # Plot means for first date
    plot_mean_hat[i,1] ~ dnorm(plot_mean[i,1], tau_proc)  # Plot mean (w/ process error) for first date
    plot_effect[i] ~ dnorm(0, plot_var) # Plot random effects, centered on zero
    
    for (t in 2:N.date){ 
      # Dynamic linear model
      plot_mean[i,t] <- alpha + beta[1]*plot_mean[i,t-1] + # Environmental predictors
        plot_effect[i] 
      plot_mean_hat[i,t] ~ dnorm(plot_mean[i,t], tau_proc)
    }
  }
  
  # Observation model
  for (i in 1:N.plot){
    for (t in 1:N.date){ 
      for (c in 1:3){ # max 3 cores per plot
        y[c,t,i] ~ dnorm(plot_mean_hat[i,t], tau_obs) 
      }
    }
  }
  
  # Priors: variance
  tau_obs ~ dgamma(1,.1) 
  tau_proc ~ dgamma(1,.1) 
  plot_var ~ dgamma(1,.1) 
  init_tau ~ dgamma(1,.1) 
  init_plot ~ dnorm(0,0.1)
  alpha ~ dnorm(0,0.1)
# Priors: coefficients
  beta[1] ~ dnorm(0,0.1)

})




rbind.named.dfs <- function(df.list){ 
  # solution from https://stackoverflow.com/questions/15162197/combine-rbind-data-frames-and-create-column-with-name-of-original-data-frames
  dfs <- df.list[sapply(df.list, function(x) !is.null(dim(x)))]
  all.out <- cbind.data.frame(do.call(rbind,dfs), 
                              name = rep(names(dfs), vapply(dfs, nrow, numeric(1))))
  return(all.out)
}


base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}




pretty_breaks <- function(n = 5, ...) {
  function(x) {
    breaks <- pretty(x, n, ...)
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
}

rbind.df.list <- function(pl.out){
  lapply(pl.out, function(x){
    do.call(rbind, x)
  })
}


pretty_names_old <- list("cellulolytic" = "Cellulolytic bacteria",
                     "assim_nitrite_reduction" = "Assimilatory nitrite reducing bacteria",
                     "dissim_nitrite_reduction" = "Dissimilatory nitrite reducing bacteria",
                     "assim_nitrate_reduction" = "Assimilatory nitrate reducing bacteria",
                     "n_fixation" = "Nitrogen fixing bacteria",
                     "dissim_nitrate_reduction" = "Dissimilatory nitrate reducing bacteria",
                     "nitrification" = "Nitrifying bacteria",
                     "denitrification" = "Denitrifying bacteria",
                     "chitinolytic" = "Chitinolytic bacteria",
                     "lignolytic" = "Ligninolytic bacteria",
                     "methanotroph" = "Methanotrophic bacteria",
                     "copiotroph" = "Copiotrophic bacteria",
                     "oligotroph" = "Oligotrophic bacteria",
                     "Arbuscular" = "Arbuscular mycorrhizal fungi",
                     "Animal_Pathogen" = "Animal pathogen fungi",
                     "Plant_Pathogen" = "Plant pathogen fungi",
                     "Saprotroph" = "Saprotrophs fungi",
                     "Wood_Saprotroph" = "Wood saprotrophs (fungi)",
                     "Ectomycorrhizal" = "Ectomycorrhizal fungi",
                     "Bdellovibrio" = "Genus Bdellovibrio",
                     "Pseudogymnoascus" = "Genus Pseudogymnoascus",
                     "Cyanobacteria" = "Phylum Cyanobacteria"
)

pretty_names <- list("cellulolytic" = "Cellulose degraders",
                     "assim_nitrite_reduction" = "Assimilatory nitrite reducers",
                     "dissim_nitrite_reduction" = "Dissimilatory nitrite reducers",
                     "assim_nitrate_reduction" = "Assimilatory nitrate reducers",
                     "n_fixation" = "Nitrogen fixers",
                     "dissim_nitrate_reduction" = "Dissimilatory nitrate reducers",
                     "nitrification" = "Nitrifiers",
                     "denitrification" = "Denitrifiers",
                     "chitinolytic" = "Chitin degraders",
                     "lignolytic" = "Lignin degraders",
                     "methanotroph" = "Methanotrophs",
                     "copiotroph" = "Copiotrophs",
                     "oligotroph" = "Oligotrophs",
                     "Arbuscular" = "Arbuscular mycorrhizae",
                     "Animal_Pathogen" = "Animal pathogens",
                     "Plant_Pathogen" = "Plant pathogens",
                     "Saprotroph" = "Saprotrophs",
                     "Wood_Saprotroph" = "Wood saprotrophs",
                     "Ectomycorrhizal" = "Ectomycorrhizae",
                     "Bdellovibrio" = "Genus Bdellovibrio",
                     "Pseudogymnoascus" = "Genus Pseudogymnoascus",
                     "Cyanobacteria" = "Phylum Cyanobacteria"
)


FG_kingdoms <- list("cellulolytic" = "Bacterial_functional_group",
                     "assim_nitrite_reduction" = "Bacterial_functional_group",
                     "dissim_nitrite_reduction" = "Bacterial_functional_group",
                     "assim_nitrate_reduction" = "Bacterial_functional_group",
                     "n_fixation" = "Bacterial_functional_group",
                     "dissim_nitrate_reduction" = "Bacterial_functional_group",
                     "nitrification" = "Bacterial_functional_group",
                     "denitrification" = "Bacterial_functional_group",
                     "chitinolytic" = "Bacterial_functional_group",
                     "lignolytic" = "Bacterial_functional_group",
                     "methanotroph" = "Bacterial_functional_group",
                     "copiotroph" = "Bacterial_functional_group",
                     "oligotroph" = "Bacterial_functional_group",
                     "Arbuscular" = "Fungal_functional_group",
                     "Animal_Pathogen" = "Fungal_functional_group",
                     "Plant_Pathogen" = "Fungal_functional_group",
                     "Saprotroph" = "Fungal_functional_group",
                     "Wood_Saprotroph" = "Fungal_functional_group",
                     "Ectomycorrhizal" = "Fungal_functional_group")


N_cyclers <- c("assim_nitrite_reduction", "dissim_nitrite_reduction","assim_nitrate_reduction","n_fixation","dissim_nitrate_reduction","nitrification","denitrification")


tukey <- function(df, x, y, extra_info = NULL){
  new.df <- cbind.data.frame(x = df[,x], y = df[,y])
  abs_max <- max(new.df[,"y"])
  maxs <- new.df %>%
    group_by(x) %>%
    summarise(tot=max(y)+ 0.3 * abs_max)
  Tukey_test <- aov(y ~ x, data=new.df) %>%
    agricolae::HSD.test("x", group=TRUE) %>%
    .$groups %>%
    as_tibble(rownames="x") %>%
    rename("Letters_Tukey"="groups") %>% 
    dplyr::select(-y) %>%
    left_join(maxs, by="x") 
  if (!is.null(extra_info)){
    Tukey_test <- cbind.data.frame(Tukey_test, extra_info)
  }
  return(Tukey_test)
}




prepModelData <- function(rank.df, 
                         j=1, 
                          subset_by_prevalence = TRUE, 
                          min.prev = 3,
                          fcast_timesteps = 10,
                          val = F
){
require(tibble)
require(purrr)
rank.df.orig <- rank.df
dates <- as.Date(rank.df.orig$dates)
rank.df.orig$dates <- NULL
rank.df.clr <- as.data.frame(compositions::clr(rank.df.orig + 1))
#rank.df.rel <- as.data.frame(rank.df.orig/rowSums(rank.df.orig)) # make compositional

dat <- cbind.data.frame( # Pull out relevant data
  dateID = as.numeric(as.character(stringr::str_replace_all(substr(dates, 1, 7), "-", ""))),
  siteID = (droplevels(factor(substr(rownames(rank.df.clr), 1, 4)))),
  plotID = (droplevels(factor(substr(rownames(rank.df.clr), 1, 8)))),
  coreID = rownames(rank.df.clr),
  y = rank.df.clr[,j])

# remove organic soil horizons
dat <- dat %>%  mutate(horizon = ifelse(grepl("-M-", coreID), "M", "O")) %>%  filter(!(siteID=="HARV" & horizon == "M"))
dat$horizon  <- NULL
#  dates <- sort(unique(dat$dateID)) # for later
start_date <- paste0(substr(min(dates, na.rm = T), 1, 7), "-01")
poss_dates <- seq.Date(as.Date(start_date), as.Date(max(dates, na.rm = T)), by = "month")
poss_dateID <- as.numeric(as.character(stringr::str_replace_all(substr(poss_dates, 1, 7), "-", "")))

# Expand data frame to include all possible plot-date combinations
all_poss_date_combos <- tidyr::expand(dat, nesting(siteID, plotID), poss_dateID) %>% dplyr::rename(dateID = poss_dateID)
dat <- merge(all_poss_date_combos, dat, all = T) %>% arrange(dateID, siteID, plotID) 

# reformat abundances into separate matrices
plot_dats <- dat %>%  mutate(plot_date = paste0(plotID, "-", dateID)) %>% 
  arrange(dateID, siteID, plotID) %>% mutate(number = 1) %>% group_by(plot_date) %>% 
  dplyr::mutate(core = cumsum(number)) %>% ungroup() %>%  group_split(plotID) %>% 
  map(.f = ~ .x %>% dplyr::select(c(dateID, y, core, plotID)) %>% 
        pivot_wider(names_from = dateID, values_from = y) %>% mutate(plotID = make.unique(as.character(plotID))) %>%  column_to_rownames(var = "plotID") %>%
        dplyr::select(-c(core)) %>% as.matrix()) 


too_small <- which(unlist(lapply(plot_dats, nrow)) < 3)
plot_dats[too_small] <- lapply(plot_dats[too_small], function(x) {
  x <- rbind(x, NA)
  return(x)
})

too_big <- which(unlist(lapply(plot_dats, nrow)) > 3)
plot_dats[too_big] <- lapply(plot_dats[too_big], function(x) {
  x <- x[1:3,] ; return(x)
})

# subset to plots that have been observed for multiple (min.prev) dates 
if (subset_by_prevalence) {
  keep_plots <- which(purrr::map(plot_dats, function(x) sum(colSums(!is.na(x)) > 0)) > min.prev)
  plot_dats <- plot_dats[keep_plots]
}


plotnames <- unlist(lapply(lapply(plot_dats, rownames), '[[', 1))



if(!is.null(fcast_timesteps)){
  names <- substr(seq.Date(max(dates, na.rm=T), length.out = fcast_timesteps+1, by = "month"), 1, 7)
  names <- gsub("-","",names)
  names <- names[-1]
  fcast_NA <- matrix(nrow = 3, ncol = fcast_timesteps) 
  colnames(fcast_NA) <- names
  plot_dats <- lapply(plot_dats, function(x) cbind(x, fcast_NA))
}


# reformat abundances into array
y <- array(unlist(plot_dats), 
           dim=c(3, ncol(plot_dats[[1]]), length(plot_dats)),
           dimnames = list(NULL, colnames(plot_dats[[1]]), plotnames))
plotID <- as.factor(plotnames)



plot.truth <- do.call(rbind, lapply(plot_dats, colMeans, na.rm=T))
rownames(plot.truth) <- plotnames

y.truth <- do.call(rbind, plot_dats)
return(list(y = y,  siteID = as.factor(substr(plotnames, 1, 4)), plotID = as.factor(plotnames), plot.truth = plot.truth, y.truth = y.truth))
}




# clim <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/daymet_monthly.rds")
# sample_dat <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilData.rds")
prep_covariate_data <- function(sample_dat, clim, cov_name = "soilTemp"){
  clim <- clim[clim$dateID > "201305",]
  
  
  new_sample_dat <- sample_dat[,c("siteID","plotID","sampleID","soilTemp","soilInCaClpH","soilMoisture","collectDate")]
  new_sample_dat$dateID <- substr(new_sample_dat$collectDate, 1, 7)
  new_sample_dat$dates <- as.Date(paste0(new_sample_dat$dateID, "-01"))
  new_sample_dat$dateID <- stringr::str_replace_all(substr(new_sample_dat$collectDate, 1, 7), "-", "")
  
#  new_sample_dat <- new_sample_dat %>%  mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%  filter(!(siteID=="HARV" & horizon == "M"))
  new_sample_dat <- new_sample_dat %>%  mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%  filter(!(siteID=="HARV" & horizon == "M")) %>% filter(horizon=="M")
  new_sample_dat$horizon  <- NULL
  
  new_sample_dat <- new_sample_dat[which(!is.na(new_sample_dat[,cov_name])),]
  # remove organic soil horizons
  #  dates <- sort(unique(dat$dateID)) # for later
  start_date <- min(new_sample_dat$dates, na.rm=T)
  start_date <- "2013-06-01"
  poss_dates <- seq.Date(as.Date(start_date), as.Date(max(new_sample_dat$dates, na.rm = T)), by = "month")
  poss_dateID <- as.numeric(as.character(stringr::str_replace_all(substr(poss_dates, 1, 7), "-", "")))
  
  # Expand data frame to include all possible plot-date combinations
  all_poss_date_combos <- tidyr::expand(new_sample_dat, nesting(siteID, plotID), poss_dateID) %>% dplyr::rename(dateID = poss_dateID)
  dat <- merge(all_poss_date_combos, new_sample_dat, all = T) %>% arrange(dateID, siteID, plotID) 
  
  clim$avgtemp <- (clim$maxtemp_mean + clim$mintemp_mean)/2
  dat_clim <- merge(dat, clim, all.x=T)
  
  dat_clim$soilTemp <- ifelse(is.na(dat_clim$soilTemp), dat_clim$avgtemp, dat_clim$soilTemp)
  dat <- dat_clim
  # reformat abundances into separate matrices
  plot_dat_list <- dat %>%  mutate(plot_date = paste0(plotID, "-", dateID)) %>% 
    arrange(dateID, siteID, plotID) %>% mutate(number = 1) %>% group_by(plot_date) %>% 
    dplyr::mutate(core = cumsum(number)) %>% ungroup() %>%  group_split(plotID) 
    
  dateNames <- unique(plot_dat_list[[1]]$dateID)
  
    plot_dats <- plot_dat_list %>%  map(.f = ~ .x %>% dplyr::select(c(dateID, !!!cov_name, core, plotID)) %>% 
          pivot_wider(names_from = dateID, values_from = all_of(cov_name)) %>% mutate(plotID = make.unique(as.character(plotID))) %>%  column_to_rownames(var = "plotID") %>%
          dplyr::select(-c(core)) %>% tidyr::fill(names(.), .direction = "downup") %>% as.data.frame())
  
    too_small <- which(unlist(lapply(plot_dats, nrow)) < 3)
    plot_dats[too_small] <- lapply(plot_dats[too_small], function(x) {
      x <- rbind(x, NA)
      return(x)
    })
    
    too_small <- which(unlist(lapply(plot_dats, nrow)) < 3)
    plot_dats[too_small] <- lapply(plot_dats[too_small], function(x) {
      x <- rbind(x, NA)
      return(x)
    })
    
    too_big <- which(unlist(lapply(plot_dats, nrow)) > 3)
    plot_dats[too_big] <- lapply(plot_dats[too_big], function(x) {
      x <- x[1:3,] ; return(x)
    })
    
    # plot_dats2 <- plot_dats %>% map(.f = ~ .x %>% t() %>% as.data.frame() %>% tidyr::fill(c(1:3), .direction = "downup") %>% t() %>% as.data.frame())
    #plotnames <- unlist(lapply(lapply(plot_dats, rownames), '[[', 1))
    plotnames <- unlist(lapply(plot_dat_list, function(x) unique(x[,"plotID"])))
    
    plot_dats <- array(unlist(plot_dats), 
               dim=c(3, ncol(plot_dats[[1]]), length(plot_dats)),
               dimnames = list(NULL, dateNames, plotnames))
    plot_dats
}



nimbleModNew <- nimbleCode({ 
	# Process model
	for (i in 1:N.plot){
		
		# Plot priors
		plot_mean[i,site_start[siteID[i]]] ~ dnorm(init_plot, init_tau)  # Plot means for first date
		plot_mean_hat[i,site_start[siteID[i]]] ~ dnorm(plot_mean[i,site_start[siteID[i]]], tau_proc) # Plot mean (w/ process error) for first date
		
		for (t in site_index[siteID[i]]:N.date){ 
			# Dynamic linear model
			plot_mean[i,t] <- alpha + beta[1]*plot_mean[i,t-1] + 
				beta[2]*mois[i,t] + 
				beta[3]*temp[i,t] + 
				beta[4]*pH[i,t] + beta[5]*pC[i,t] + 
				beta[6]*nspp[i,t] + beta[7]*rc_grass[i,t] + 
				time_effect[siteID[i],t]  # Temporal effects
			# Add process error
			plot_mean_hat[i,t] ~ dnorm(plot_mean[i,t], tau_proc)
		}
	}
	
	# Observation model
	for (i in 1:N.plot){
		for (t in site_start[siteID[i]]:N.date){ 
			for (c in 1:3){ # max 3 cores per plot
				y[c,t,i] ~ dnorm(plot_mean_hat[i,t], tau_obs) 
			}
		}
	}
	
	# Priors: random effects
	for (i in 1:N.site){
		for (t in site_start[i]:N.date){
			time_effect[i,t] ~ dnorm(0, time_var) # Prior on time random effects
		}
	}
	
	# Priors: variance
	tau_obs ~ dgamma(1,.1) 
	tau_proc ~ dgamma(1,.1) 
	time_var ~ dgamma(1,.1) 
	init_tau ~ dgamma(1,.1) 
	init_plot ~ dnorm(0,0.1)
	alpha ~ dnorm(0,0.1)

	# Priors: coefficients
	for(i in 1:7){
		beta[i] ~ dnorm(0,0.1)
	}
})


nimbleMod_cohesion <- nimbleCode({ 
	# Process model
	for (i in 1:N.plot){
		
		# Plot priors
		plot_mean[i,site_start[siteID[i]]] ~ dnorm(init_plot, init_tau)  # Plot means for first date
		plot_mean_hat[i,site_start[siteID[i]]] ~ dnorm(plot_mean[i,site_start[siteID[i]]], tau_proc) # Plot mean (w/ process error) for first date
		
		for (t in site_index[siteID[i]]:N.date){ 
			# Dynamic linear model
			plot_mean[i,t] <- alpha + beta[1]*plot_mean[i,t-1] + 
				beta[2]*mois[i,t] + 
				beta[3]*temp[i,t] + 
				beta[4]*pH[i,t] + beta[5]*pC[i,t] + 
				beta[6]*nspp[i,t] + beta[7]*rc_grass[i,t] + beta[8]*coh_pos[i,t-1] + beta[9]*coh_neg[i,t-1] +
				time_effect[siteID[i],t]  # Temporal effects
			# Add process error
			plot_mean_hat[i,t] ~ dnorm(plot_mean[i,t], tau_proc)
		}
	}
	
	# Observation model
	for (i in 1:N.plot){
		for (t in site_start[siteID[i]]:N.date){ 
			for (c in 1:3){ # max 3 cores per plot
				y[c,t,i] ~ dnorm(plot_mean_hat[i,t], tau_obs) 
			}
		}
	}
	
	# Priors: random effects
	for (i in 1:N.site){
		for (t in site_start[i]:N.date){
			time_effect[i,t] ~ dnorm(0, time_var) # Prior on time random effects
		}
	}
	
	# Priors: variance
	tau_obs ~ dgamma(1,.1) 
	tau_proc ~ dgamma(1,.1) 
	time_var ~ dgamma(1,.1) 
	init_tau ~ dgamma(1,.1) 
	init_plot ~ dnorm(0,0.1)
	alpha ~ dnorm(0,0.1)
	
	# Priors: coefficients
	for(i in 1:9){
		beta[i] ~ dnorm(0,0.1)
	}
})


#cohesion <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cohesion_all_samples.rds")

prep_cohesion <- function(model.dat, cohesion, param = "PositiveCohesion_network_16S", scale = T){
	keep.plots <- cbind.data.frame(plotID = model.dat$plotID, siteID = substr(model.dat$plotID, 1, 4))
	
	if (scale == T){
	cohesion.scale <- scale(cohesion[,2:11])
	} else cohesion.scale <- cohesion[,2:11]
	
	coh.df <- cbind(cohesion.scale, parseNEONsampleIDs(cohesion$sampleID))
	coh.df[,param] <- as.numeric(gsub("NaN", "NA", coh.df[,param]))
	coh.df$date <- as.Date(gsub("-[0-9][0-9]$", "-01", coh.df$asDate))
	coh_plot <- coh.df %>% #merge(keep.plots, all.y=T) %>% 
		select(plotID, siteID, date, !!as.name(param)) %>% group_by(siteID, plotID, date) %>% 
		mutate(plot_mean = mean(!!as.name(param), na.rm=T)) %>% 
		mutate(plot_mean = as.numeric(gsub("NaN", "NA", plot_mean))) %>% 
		select(plotID, siteID, date, plot_mean) %>% 
		distinct(.keep_all = F) %>% 
		group_by(plotID) %>% 
		padr::pad(start_val = as.Date("2013-06-01"), end_val = as.Date("2018-11-01")) %>% 
		mutate(siteID = substr(plotID, 1, 4)) %>% group_by(plotID) %>% 
		tidyr::fill(plot_mean, .direction = "downup") %>% ungroup() 
	#merge(unique(temp[,c("siteID", "date")]), all=T) %>% 
	#arrange(date, plotID) %>% unique() %>% 
	
	coh_site <- coh.df %>% #merge(keep.plots, all.y=T) %>% 
		select(siteID, date, !!as.name(param)) %>% group_by(siteID, date) %>% 
		mutate(site_mean = mean(!!as.name(param), na.rm=T)) %>% 
		mutate(site_mean = as.numeric(gsub("NaN", "NA", site_mean))) %>% 
		select(siteID, date, site_mean) %>% 
		distinct(.keep_all = F) %>% 
		group_by(siteID) %>% 
		padr::pad(start_val = as.Date("2013-06-01"), end_val = as.Date("2018-11-01"))
	
	coh_out <- merge(coh_plot, coh_site, all=T) %>% 
		mutate(out_val = ifelse(is.na(plot_mean), site_mean, plot_mean)) %>% 
		tidyr::fill(out_val, .direction = "downup") %>% ungroup() %>% select(date, plotID, out_val) %>% 
		pivot_wider(names_from = date, values_from = out_val, values_fn = list)  %>% 
		arrange(plotID)  %>% filter(!is.na(plotID)) %>% 
		column_to_rownames(var = "plotID") %>% as.data.frame()
	
	coh_plot_out <- coh_out[rownames(coh_out) %in% keep.plots$plotID,]
	coh_plot_out <- coh_plot_out[keep.plots$plotID,] %>% data.matrix()
	return(coh_plot_out)
}


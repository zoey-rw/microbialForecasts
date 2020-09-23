# getting plot and site level means of NEON functional groups.
# adapted from Colin Averill's script: https://github.com/colinaverill/NEFI_microbe/blob/master/ITS/data_construction/NEON_ITS/4._sequence_processing_raw_fastq/7._hierarch_group_means_cps.r
rm(list=ls())
library(runjags)
library(foreach)
library(doParallel)

# source hierarch_ddirch_means.r and tic_toc.r from Colin Averill's Github
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/hierarch_ddirch_means.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/tic_toc.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#load data and format.----
d <- readRDS("data/ITS_fg_abundances.rds")

#register parallel environment.----
n.cores <- detectCores()
registerDoParallel(n.cores)

tic()
#loop over levels.----
output <- list()
output <-
  foreach(i = 1:length(d)) %dopar% {
    #Get y multivariate matrix.
    abundances <- d[[i]]$abundances
    seq.depth  <- d[[i]]$seq_total
    #y <- abundances / seq.depth #modified for multinomial dirichlet. Crashes function currently, need a hierarch multinomial dirichlet function.
    y <- as.matrix((abundances + 1) / rowSums(abundances + 1)) #use this for stanadard dirichlet.
    
    #get core_plot and plot_site indexing.
    core_plot <- substr(rownames(y), 1, 8)
    core_site <- substr(rownames(y), 1, 4)
    plot_site <- unique(core_plot)
    plot_site <- substr(plot_site, 1, 4)
    
    #fit the hierarchical means.
    fit <- hierarch_ddirch_means(y=y, core_plot = core_plot, plot_site = plot_site, jags.method = 'parallel')
    
    #add row and column names - plot level (you should put this in the function).
    for(j in 1:length(fit$plot.fit)){
      rownames(fit$plot.fit[[j]]) <- unique(core_plot)
      rownames(fit$site.fit[[j]]) <- unique(plot_site)
      colnames(fit$plot.fit[[j]]) <- colnames(y)
      colnames(fit$site.fit[[j]]) <- colnames(y)
    }
    fit$core.fit <- y #add in the core-level data!
    return(fit)
    cat(names(d)[i],'model fit.',i,'of',length(d),'groups aggregated.\n')
  }
names(output) <- names(d)

saveRDS(output, "data/ITS_fg_hier_means.rds")

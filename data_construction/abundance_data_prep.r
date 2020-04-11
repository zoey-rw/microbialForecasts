# Prep data for EE509 model - fun and bac abundances for 10 most prevalent groups at every taxonomic rank.

library(tidyr)
library(phyloseq)
library(padr)
library(ggplot2)
library(dplyr)
source('/projectnb/talbot-lab-data/zrwerbin/NEFI_microbe/NEFI_functions/crib_fun.r')


#load data and format.----
d.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_groupAbundances_16S.rds")
d.bac$ps_16S <- NULL

# for each rank
out.bac <- list()
#for (i in 1:18) {
  for (i in 1:5) {
    
# Get 15 most prevalent taxa from rank  
prev <- d.bac[[i]]$prevalence[order(-d.bac[[i]]$prevalence$prevalence),]
tax_to_keep <- as.character(prev[1:10,1])

# Convert to relative abundances
rank.df <- d.bac[[i]]$abundances
rank.df <- rank.df[,colnames(rank.df) %in% tax_to_keep]
rank.df <- rank.df/d.bac[[i]]$seq.total
rank.df$other <- 1-rowSums(rank.df)
rank.df <- data.frame(lapply(rank.df, crib_fun, N = nrow(rank.df) * ncol(rank.df)))
rownames(rank.df) <- rownames(d.bac[[i]]$rel.abundances)

core_plot <- unique(substr(rownames(rank.df), 1, 8))
plot_site <- unique(substr(core_plot, 1, 4))

rank.df$dates <- as.Date(stringr::str_extract(rownames(rank.df), "201[0-9][0-9][0-9][0-9][0-9]"), "%Y%m%d")
rank.df <- rank.df[order(rank.df$dates),]
# subset by date (till 2016).
rank.df <- rank.df[rank.df$dates < "2017-01-01",]
out.bac[[i]] <- rank.df
}
names(out.bac) <- paste0(names(d.bac)[1:5],"_bac")



#load data and format.----
d.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_groupAbundances_ITS.rds")
d.fun$ps_ITS <- NULL
    
# for each rank
out.fun <- list()
#for (i in 1:length(d)) {
  for (i in 1:5) {
    
  # Get 15 most prevalent taxa from rank  
  prev <- d.fun[[i]]$prevalence[order(-d.fun[[i]]$prevalence$prevalence),]
  tax_to_keep <- as.character(prev[1:10,1])
  
  # Convert to relative abundances
  rank.df <- d.fun[[i]]$abundances
  rank.df <- rank.df[,colnames(rank.df) %in% tax_to_keep]
  rank.df <- rank.df/d.fun[[i]]$seq.total
  rank.df$other <- NULL
  rank.df <- data.frame(lapply(rank.df, crib_fun, N = nrow(rank.df) * ncol(rank.df)))
  rank.df$other <- 1-rowSums(rank.df)

  rownames(rank.df) <- rownames(d.fun[[i]]$rel.abundances)
  
  core_plot <- unique(substr(rownames(rank.df), 1, 8))
  plot_site <- unique(substr(core_plot, 1, 4))
  
  rank.df$dates <- as.Date(stringr::str_extract(rownames(rank.df), "201[0-9][0-9][0-9][0-9][0-9]"), "%Y%m%d")
  rank.df <- rank.df[order(rank.df$dates),]
  # subset by date (till 2016).
  rank.df <- rank.df[rank.df$dates < "2017-01-01",]
  out.fun[[i]] <- rank.df
}
names(out.fun) <- paste0(names(d.fun)[1:5],"_fun")


out <- c(out.bac, out.fun)
saveRDS(out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/groupAbundances_EE509.rds")

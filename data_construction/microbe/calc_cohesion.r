# script to calculate a few different cohesion metrics for NEON ITS/16S datasets



library(foreach)
library(doParallel)

cl <- parallel::makeForkCluster(28, outfile="")
#cl <- parallel::makeForkCluster(3)
registerDoParallel(cl)


#### Load libraries & scripts ####
library(SpiecEasi)
library(tidyr)
library(phyloseq)
library(padr)
library(ggplot2)
library(dplyr)
library(data.table) 
#library(phylosmith) # devtools::install_github('schuyler-smith/phylosmith')
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
#source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/helperFunctions.r")
#source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/calcCohesion.r")
#source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/calcCohesion_multi.r")

#### 1. Load phyloseq data ####
# # Combine legacy and recent phyloseq objects (16S)
# recent_ps_16S <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_subset.rds")
# legacy_ps_16S <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_legacy.rds")
# new_sample_dat_16S <- parseNEONsampleIDs(as.character(sample_data(recent_ps_16S)$dnaSampleID))
# rownames(new_sample_dat_16S) <- rownames(sample_data(recent_ps_16S))
# sample_data(recent_ps_16S) <- new_sample_dat_16S
# ps_16S <- merge_phyloseq(legacy_ps_16S, recent_ps_16S)
# colnames(tax_table(ps_16S)) <- tolower(colnames(tax_table(ps_16S)))
# 
# # Combine legacy and recent phyloseq objects (ITS)
# recent_ps_ITS <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_subset.rds")
# legacy_ps_ITS <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_legacy_subset.rds")
# new_sample_dat_ITS <- parseNEONsampleIDs(as.character(sample_data(recent_ps_ITS)$dnaSampleID))
# rownames(new_sample_dat_ITS) <- rownames(sample_data(recent_ps_ITS))
# sample_data(recent_ps_ITS) <- new_sample_dat_ITS
# colnames(tax_table(legacy_ps_ITS)) <- tolower(colnames(tax_table(legacy_ps_ITS)))
# ps_ITS <- merge_phyloseq(legacy_ps_ITS, recent_ps_ITS)
# saveRDS(list(ps_16S, ps_ITS), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S_ITS.rds")

# Read in 16S and ITS phyloseq objects
ps_16S_orig <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S_ITS.rds")[[1]]
ps_ITS_orig <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S_ITS.rds")[[2]]

# Align the samples from each phyloseq objects
sample_names(ps_16S_orig) <- make.unique(sample_data(ps_16S_orig)$sample)
sample_names(ps_ITS_orig) <- make.unique(sample_data(ps_ITS_orig)$sample)

ps_16S <- ps_16S_orig
ps_ITS <- ps_ITS_orig


# Remove low-prevalence taxa
ps_16S = filter_taxa(ps_16S_orig, function(x) sum(x > 0) > 3, TRUE)
ps_ITS = filter_taxa(ps_ITS_orig, function(x) sum(x > 0) > 3, TRUE)
# shared_samples <- intersect(sample_data(ps_16S_orig)$sample, sample_data(ps_ITS_orig)$sample)
# ps_16S <- subset_samples(ps_16S_orig, sample_names(ps_16S_orig) %in% shared_samples)
# ps_ITS <- subset_samples(ps_ITS_orig, sample_names(ps_ITS_orig) %in% shared_samples)
# otu_table(ps_16S) <- otu_table(ps_16S)[shared_samples,]
# otu_table(ps_ITS) <- otu_table(ps_ITS)[shared_samples,]
# waldo::compare(sample_names(ps_ITS), sample_names(ps_16S))

# Only keep sites with at least 30 samples total (for bacteria at least..)
site.tally <- parseNEONsampleIDs(sample_names(ps_16S)) %>% 
	group_by(siteID) %>% tally() %>% dplyr::filter(n > 30)
siteID.list <- unique(site.tally$siteID)

# Function for subsetting
subset_by_site <- function(phyloseq_obj, site_name) {
	var_values <- sample_data(phyloseq_obj)[["siteID"]]
	prune_samples(var_values %in% site_name, phyloseq_obj)
}

#cohesion.out.list <- foreach(site=siteID.list, 
connectivity.out.list <- foreach(site=siteID.list, 
														 	.errorhandling='pass', 
														 .export = c("siteID.list","subset_by_site"),
														 .packages = c("phyloseq", "dplyr","SpiecEasi")
														 ) %dopar% {
#for(s in 1:3){
	source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
	source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/helperFunctions.r")
	source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/calcCohesion.r")
	source("/projectnb2/talbot-lab-data/zrwerbin/cohesion/cohesion_git/calcCohesion_multi.r")

# Subset to just one site 
site_16S <- subset_by_site(ps_16S, site)
site_ITS <- subset_by_site(ps_ITS, site)

# Remove low-prevalence taxa
site_16S = filter_taxa(site_16S, function(x) sum(x > 0) > (0.2*length(x)), TRUE)
site_ITS = filter_taxa(site_ITS, function(x) sum(x > 0) > (0.2*length(x)), TRUE)

# Pull out OTU tables (script runs slow if you create null model on phyloseq object)
otu_16S = as(otu_table(site_16S), "matrix") 
otu_ITS = as(otu_table(site_ITS), "matrix") 

#### 2. Calculate various cohesion metrics ####

# A: ORIGINAL: calculate cohesion using the Herren and McMahon 2017 approach (using a null model)
cohesion_orig_16S <- calcCohesion(otu = otu_16S, pers.cutoff = .2, iter = 30, suffix = "_orig_16S")
cohesion_orig_ITS <- calcCohesion(otu = otu_ITS, pers.cutoff = .2, iter = 30, suffix = "_orig_ITS")

# B: SPIEC-EASI, SINGLE DOMAIN: create network via SpiecEasi, and use as the cohesion custom-correlation input
lambda.min.ratio = .1
nlambda = 40
# 16S
network_16S <- spiec.easi(site_16S, method='glasso', lambda.min.ratio=lambda.min.ratio, 
													sel.criterion='bstars', nlambda=nlambda, pulsar.params = list(ncores=3))
evalSpiecEasi(network_16S)
cor_16S <- cov2cor(as.matrix(getOptCov(network_16S))) # create correlation matrix from output
cohesion_network_16S <- calcCohesion(otu = otu_16S, custom.cor = cor_16S, suffix = "_network_16S") 
# ITS
network_ITS <- spiec.easi(site_ITS, method='glasso', lambda.min.ratio=lambda.min.ratio, 
													sel.criterion='bstars', nlambda=nlambda, pulsar.params = list(ncores=3))
evalSpiecEasi(network_ITS)
cor_ITS <- cov2cor(as.matrix(getOptCov(network_ITS))) # create correlation matrix from output
cohesion_network_ITS <- calcCohesion(otu = otu_ITS, custom.cor = cor_ITS, suffix = "_network_ITS")

# C: SPIEC-EASI, MULTI DOMAIN: create a multi-domain network via SPIEC-EASI, get correlation matrix, use for multi-domain cohesion

# check that samples match
shared_samples <- intersect(sample_names(site_16S), sample_names(site_ITS))
otu_table(site_16S) <- otu_table(site_16S)[shared_samples,]
otu_table(site_ITS) <- otu_table(site_ITS)[shared_samples,]
waldo::compare(sample_names(site_ITS), sample_names(site_16S))

network_16S_ITS <- spiec.easi(list(site_16S, site_ITS), method='glasso', lambda.min.ratio = lambda.min.ratio, 
															sel.criterion='bstars',	nlambda=nlambda, pulsar.params = list(ncores=3))
evalSpiecEasi(network_16S_ITS)
cor_16S_ITS <- cov2cor(as.matrix(getOptCov(network_16S_ITS))) # create correlation matrix from output
cohesion_network_16S_ITS <- calcCohesionMulti(otu1 = as(otu_table(site_16S), "matrix") , otu2 = as(otu_table(site_ITS), "matrix") , custom.cor = cor_16S_ITS, suffix = "_multi_network")

coh.merged <- Reduce(
	function(x, y, ...) merge(x, y, all = TRUE, ...),
	list(cohesion_orig_16S$Cohesion, 
			 cohesion_network_16S$Cohesion,
			 cohesion_orig_ITS$Cohesion, 
			 cohesion_network_ITS$Cohesion,
			 cohesion_network_16S_ITS$Cohesion)
)

conn.merged <- data.table::rbindlist(
	list(cbind(cohesion_orig_16S$Connectedness, "coh_orig_16S"), 
			 cbind(cohesion_network_16S$Connectedness,"coh_network_16S"), 
			 cbind(cohesion_orig_ITS$Connectedness, "coh_orig_ITS"), 
			 cbind(cohesion_network_ITS$Connectedness, "coh_network_ITS"), 
			 cbind(cohesion_network_16S_ITS$Connectedness, "coh_multi_network")),  use.names=FALSE)  


#return(coh.merged)
return(conn.merged)
}

#cohesion.out <- do.call(rbind, cohesion.out.list)
connectivity.out <- do.call(rbind, connectivity.out.list)
#cohesion.out.list <- mclapply(1:35, calc_site_cohesion, mc.cores = 28)

#saveRDS(cohesion.out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cohesion_all_samples.rds")
saveRDS(connectivity.out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/connectivity_all_taxa.rds")


stopCluster(cl)




# testing cohesion...
library(tidyverse)
cohesion.out <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cohesion_all_samples.rds")
coh.df <- cbind(cohesion.out, parseNEONsampleIDs(cohesion.out$sampleID)) %>% select(-c(1))

library(ggpubr)
#cohesion.long <- cohesion.out %>% select(sampleID, NegativeCohesion_multi_network, PositiveCohesion_multi_network) %>%  pivot_longer(cols = c(2:3))
p <- ggplot(coh.df[coh.df$siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"),]) + 
	geom_point(aes(x = NegativeCohesion_multi_network, y = PositiveCohesion_multi_network, color = siteID)) + 
	theme_bw(base_size = 24) + xlab("Negative cohesion") + 
	ylab("Positive cohesion") + ggtitle("Fungi-bacteria cohesion, per-sample")
p
sp <- ggscatter(coh.df[coh.df$siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"),], x = "NegativeCohesion_multi_network", y = "PositiveCohesion_multi_network",
								add = "reg.line",  # Add regressin line
								add.params = list(color = "siteID", fill = "siteID"), # Customize reg. line
								conf.int = TRUE # Add confidence interval
)
sp + stat_cor(method = "pearson", label.x = -.05, label.y = .12)



plot(NegativeCohesion_orig_16S ~ PositiveCohesion_orig_16S, cohesion.out)
plot(NegativeCohesion_orig_16S ~ NegativeCohesion_network_16S, cohesion.out)
plot(NegativeCohesion_orig_16S ~ NegativeCohesion_multi_network, cohesion.out)

plot(NegativeCohesion_orig_ITS ~ NegativeCohesion_multi_network, cohesion.out)
plot(NegativeCohesion_orig_ITS ~ NegativeCohesion_orig_16S, cohesion.out)


param <- "NegativeCohesion_orig_ITS"
param <- "NegativeCohesion_orig_16S"

coh.df <- cbind(cohesion, parseNEONsampleIDs(cohesion$sampleID))
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
tidyr::fill(plot_mean, .direction = "downup") %>% ungroup() %>% select(date, plotID, out_val) %>% 
	pivot_wider(names_from = date, values_from = out_val, values_fn = list)  %>% 
	arrange(plotID)  %>% filter(!is.na(plotID)) %>% 
	column_to_rownames(var = "plotID") %>% as.data.frame()

coh_plot_out <- coh_out[rownames(coh_out) %in% keep.plots$plotID,]
return(coh_plot_out)
	





connectivity.out <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/connectivity_all_taxa.rds")
plot(NegativeConnectedness_orig_16S ~ PositiveConnectedness_orig_16S, connectivity.out)

names(connectivity.out) <- c("PositiveConnectedness", "NegativeConnectedness", "Metric")
connectivity.out <- connectivity.out[!grepl("orig", connectivity.out$Metric),]

library(ggplot2)
ggplot(connectivity.out) + geom_point(aes(x = NegativeConnectedness, y = PositiveConnectedness, color = Metric))


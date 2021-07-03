# Script for preparing bacterial abundance data for models
# Input: phyloseq objects. Combines with functional group data and outputs abundances.
# Output is a list of taxonomic rank/functional groups
setwd("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/")

library(tidyr)
library(phyloseq)
library(padr)
library(ggplot2)
library(dplyr)
source('/projectnb/talbot-lab-data/zrwerbin/NEFI_microbe/NEFI_functions/crib_fun.r')
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/createTaxFunction.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/addBacterialFunction.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

#### 16S #####

# Combine legacy and recent phyloseq objects
recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_subset.rds")
legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_legacy.rds")
new_sample_dat <- parseNEONsampleIDs(as.character(sample_data(recent_ps)$dnaSampleID))
rownames(new_sample_dat) <- rownames(sample_data(recent_ps))
sample_data(recent_ps) <- new_sample_dat
master_ps <- merge_phyloseq(legacy_ps, recent_ps)
colnames(tax_table(master_ps)) <- tolower(colnames(tax_table(master_ps)))

# Assign functional groups
tax_ref <- createTaxFunction(ref.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/bacteria_func_groups.csv", 
														 N.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/Npathways_Albright2018.csv", 
														 C.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/cellulolytic_Berlemont.csv",
														 Naylor.path = "/projectnb2/talbot-lab-data/zrwerbin/random_data/Naylor_functional_groups/functional_module_df.rds",
														 out.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/tax_function_ref.csv")
tax_df = as.data.frame(as(tax_table(master_ps), "matrix"))
new_tax <- addBacterialFunction(tax = tax_df[,!colnames(tax_df) %in% "kingdom"], 
																tax_fun_ref = tax_ref)
tax_table(master_ps) <- new_tax

master_ps.filt = filter_taxa(master_ps, function(x) sum(x) > 20, TRUE)


# Make relative
ps.rel <- transform_sample_counts(master_ps, function(x) x/sum(x))

ps.filt = filter_taxa(ps.rel, function(x) sum(x) > .005, TRUE)




out <- get_tax_level_abun(master_ps, 
													tax_rank_list = colnames(tax_table(master_ps)), 
													min_seq_depth = 5000)

# Now go through ranks to get top abundances
tax_rank <- "genus"

n.taxa <- 10

cal.out.bac <- list()
val.out.bac <- list()
for (tax_rank in names(out)){
	
	rank_abun <- out[[tax_rank]]$rel.abundances
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	most_abundant_taxa <- prev_top[,1]
	most_abundant_taxa <- most_abundant_taxa[!most_abundant_taxa=="other"][1:n.taxa]
	most_abundant_taxa <- gsub("\\-|\\(|\\)", "\\.", most_abundant_taxa)
	out_top10 <- rank_abun[,colnames(rank_abun) %in% most_abundant_taxa, drop=F]
	seqDepth <- rowSums(rank_abun)
	out_top10$other <- 1-rowSums(out_top10)
	# Remove samples with a sum of above one (not sure why they exist)
	out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]
	
	
	ps.rank.filt <- prune_samples(sample_names(master_ps) %in% rownames(out_top10), master_ps)
	rank.df <- cbind(sample_data(ps.rank.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)
	
	# organize by date
	rank.df$dates <- as.Date(as.character(rank.df$dates), "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]
	# subset by date (till 2016).
	val <- rank.df[rank.df$dates >= "2017-01-01",]
	cal <- rank.df[rank.df$dates < "2017-01-01",]
	
	cal.out.bac[[tax_rank]]	<- cal
	val.out.bac[[tax_rank]]	<- val
}






ranks <- rank_names(master_ps)
for (tax_rank in ranks){
	
	print(tax_rank)
#for (tax_rank in list("phylum","class","order","family","genus")){
ps.phy <- tax_glom(ps.filt, tax_rank)
glom_melt <- speedyseq::psmelt(ps.phy)
form <- as.formula(paste0("sampleID ~ ", tax_rank))
glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
most_abundant_taxa <- names(sort(colSums(out_abun), decreasing = T)[1:n.taxa])
seqDepth <- rowSums(out_abun)
out_top10 <- out_abun[,colnames(out_abun) %in% most_abundant_taxa]
# Remove samples with a sum of above one (not sure why they exist)
out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]
out_top10$Other <- 1-rowSums(out_top10)
ps.phy.filt <- prune_samples(sample_names(ps.filt) %in% rownames(out_top10), ps.filt)
rank.df <- cbind(sample_data(ps.phy.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)
#saveRDS(rank.df,"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bac_abundances.rds")


# cal.out.bac <- list()
# val.out.bac <- list()
#for (i in 1:length(d.bac)) {
	
	# organize by date
	rank.df$dates <- as.Date(rank.df$dates, "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]
	# subset by date (till 2016).
	val <- rank.df[rank.df$dates >= "2017-01-01",]
	cal <- rank.df[rank.df$dates < "2017-01-01",]
	
	cal.out.bac[[tax_rank]]	<- cal
}
#	val.out.bac[[i]] <- val

names(cal.out.bac)[1:5] <- paste0(names(cal.out.bac)[1:5], "_bac")
names(val.out.bac)[1:5] <- paste0(names(val.out.bac)[1:5], "_bac")

saveRDS(cal.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_10tax.rds")
saveRDS(val.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_10tax.rds")

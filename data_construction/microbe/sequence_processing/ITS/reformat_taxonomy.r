# Converting ITS outputs into formats suitable for modeling.

library(phyloseq)
library(dplyr)
library(tidyr)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/assign_fungal_guilds.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")

# Load combined sequence table and taxonomic table
seqtab_orig <- readRDS("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_subset.rds")
taxa_orig <- read.table("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy/final_tax_table.csv")

# Prep for phyloseq
new_tax <- do.call(rbind, (lapply(taxa_orig$V2, parse_taxonomy_qiime)))
rownames(new_tax) <- taxa_orig$V1
sample_dat <- parseNEONsampleIDs(rownames(seqtab_orig))

# Assign functional groups
tax <- as.data.frame(new_tax)
tax.fun <- assign_fungal_guilds(tax_table = tax, n.cores = 8)

#na <- tax.fun[which(is.na(tax.fun$guild_collapse)),]

ps_its <- phyloseq(otu_table(seqtab_orig, taxa_are_rows = F), tax_table(as.matrix(tax.fun)), sample_data(sample_dat))

out <- get_tax_level_abun(ps_its, 
													tax_rank_list = colnames(tax_table(ps_its)), 
													min_seq_depth = 3000)

master_ps <- ps_its

# Now go through ranks to get top abundances
n.taxa <- 10
cal.out.bac <- list()
val.out.bac <- list()
ranks <- names(out)
ranks <- ranks[ranks != "kingdom"]
for (tax_rank in ranks){
	
	rank_abun <- out[[tax_rank]]$rel.abundances
	
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	
	if (tax_rank=="phylum|genus") prev_top <- prev_top[prev_top$prevalence > .3,]
	most_abundant_taxa <- prev_top[, 1]
	most_abundant_taxa <- most_abundant_taxa[!grepl("unassigned|other|genus|class|family|order|^fungi",most_abundant_taxa)][1:n.taxa]
	
	
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
names(cal.out.bac)[1:5] <- paste0(names(cal.out.bac)[1:5], "_fun")
names(val.out.bac)[1:5] <- paste0(names(val.out.bac)[1:5], "_fun")


saveRDS(cal.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds")
saveRDS(val.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds")



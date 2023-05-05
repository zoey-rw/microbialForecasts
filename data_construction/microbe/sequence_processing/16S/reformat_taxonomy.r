library(phyloseq)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/createTaxFunction.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/addBacterialFunction.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
# RECENT DATA #
# Load combined sequence table and taxonomic table
seqtab_orig <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S.rds")
taxa_orig <- read.csv("/projectnb/microbiome/zrwerbin/NEON_amplicon/16S/NEON_16S_taxonomy/final_tax_table.csv", header=F, sep="\t")

# Prep for phyloseq
new_tax <- do.call(rbind, (lapply(taxa_orig$V2, parse_taxonomy_qiime)))
rownames(new_tax) <- taxa_orig$V1
colnames(new_tax)[1:5] <- c("Kingdom","Phylum","Class","Order","Family")
sample_dat <- parseNEONsampleIDs(rownames(seqtab_orig))

ps_recent <- phyloseq(otu_table(seqtab_orig, taxa_are_rows = F),
							 tax_table(new_tax), sample_data(sample_dat))



# LEGACY #
# Load combined sequence table and taxonomic table
seqtab_legacy <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S_legacy.rds")
taxa_legacy <- read.csv("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/legacy_anomaly_tax_16S/final_tax_table.csv", header=F, sep="\t")

# Prep for phyloseq
new_tax_legacy <- do.call(rbind, (lapply(taxa_legacy$V2, parse_taxonomy_qiime)))
rownames(new_tax_legacy) <- taxa_legacy$V1
sample_dat_legacy <- parseNEONsampleIDs(rownames(seqtab_legacy))

ps_legacy <- phyloseq(otu_table(seqtab_legacy, taxa_are_rows = F),
											tax_table(new_tax_legacy), sample_data(sample_dat_legacy))

ps_legacy_prune <- prune_samples(!sample_names(ps_legacy) %in% c(sample_dat$geneticSampleID, sample_dat$sampleID, sample_dat$sample), ps_legacy)

# Combine legacy and recent!
ps <- merge_phyloseq(ps_recent, ps_legacy_prune)


# Assign functional groups
tax_ref <- createTaxFunction(ref.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/bacteria_func_groups.csv",
														 N.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/Npathways_Albright2018.csv",
														 C.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/cellulolytic_Berlemont.csv",
														 Naylor.path = "/projectnb2/talbot-lab-data/zrwerbin/random_data/Naylor_functional_groups/functional_module_df.rds",
														 out.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/tax_function_ref.csv")
tax_df = as.data.frame(as(tax_table(ps), "matrix"))
new_tax <- addBacterialFunction(tax = tax_df[,!colnames(tax_df) %in% "Kingdom"],
																tax_fun_ref = tax_ref)
tax_table(ps) <- new_tax


# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")


ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")
# Make relative
# ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
# ps.filt = filter_taxa(ps.rel, function(x) sum(x) > .005, TRUE)

# Takes a few minutes to run.
out <- get_tax_level_abun(ps,
													tax_rank_list = colnames(tax_table(ps))[1:6], # appending fg abundances at the end instead
													#tax_rank_list = colnames(tax_table(ps)),
													min_seq_depth = 5000)

# Now go through ranks to get top abundances
n.taxa <- 20
out.bac <- list()
for (tax_rank in names(out)){
	print(tax_rank)
	rank_abun <- out[[tax_rank]]$rel.abundances
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	most_abundant_taxa <- prev_top[,1]
	most_abundant_taxa <- most_abundant_taxa[!grepl("unassigned|other|genus|family|order|^bacteria",most_abundant_taxa)][1:n.taxa]
	most_abundant_taxa <- gsub("\\-|\\(|\\)| ", "\\.", most_abundant_taxa)
	out_top10 <- rank_abun[,colnames(rank_abun) %in% most_abundant_taxa, drop=F]
	seqDepth <- rowSums(rank_abun)
	out_top10$other <- 1-rowSums(out_top10)
	# Remove samples with a sum of above one (not sure why they exist)
	out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]


	ps.rank.filt <- prune_samples(sample_names(ps) %in% rownames(out_top10), ps)
	rank.df <- cbind(sample_data(ps.rank.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)

	# organize by date
	rank.df$dates <- as.Date(as.character(rank.df$dates), "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]
	# saveRDS(rank.df, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/genus_groupAbundances_16S.rds")

	out.bac[[tax_rank]] <- rank.df
}


names(out.bac)[1:6] <- paste0(names(out.bac)[1:6], "_bac")

saveRDS(out.bac_save, here("data/clean/groupAbundances_16S_2023.rds"))


# Append previous fg abundances since they didn't change
old_cal = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds")
old_val = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds")

fg_list = list()
for (i in 7:length(old_cal)){
	fg_list[[i]] <- rbind(old_cal[[i]], old_val[[i]])
}
fg_list = fg_list[7:length(old_cal)]
names(fg_list) <- names(old_cal)[7:length(old_cal)]
out.bac_save <- c(out.bac, fg_list)

saveRDS(out.bac_save, here("data/clean/groupAbundances_16S_2023.rds"))


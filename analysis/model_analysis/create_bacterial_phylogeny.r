#source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

options(scipen=999)

library(phyloseq)
library(DECIPHER) #BiocManager::install("DECIPHER")
library(phangorn)
library(speedyseq)
library(tidyverse)
library(speedyseq)

library(phytools)
library(ggtree)
library(ape)
# library(devtools)
# install_github("YuLab-SMU/ggtree")


# Set output path
output_path <- here("data/clean/bacterial_phylogeny.rds")
# Read in data
ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")
# Remove weird characters from names
ps@tax_table[,c("family")] <- gsub("\\-|\\(|\\)| ", "\\.", ps@tax_table[,c("family")])
ps@tax_table[,c("genus")] <- gsub("\\-|\\(|\\)| ", "\\.", ps@tax_table[,c("genus")])


# Read & reshape model effect estimates
sum.all <- readRDS(here("data/summary/predictor_effects.rds"))
betas_wide = sum.all %>% filter(model_name=="env_cycl" & pretty_group=="Bacteria") %>%
	pivot_wider(id_cols=c("taxon","rank_only","time_period"),
							names_from = "beta", values_from = "Mean")

# Get list of taxa at each rank
summary_unique = sum.all %>% filter(pretty_group=="Bacteria") %>% distinct(only_rank, taxon)
keep_taxa = split(summary_unique, summary_unique$only_rank)

phylum_keep <- keep_taxa[[5]]$taxon
class_keep <- keep_taxa[[4]]$taxon
order_keep <- keep_taxa[[3]]$taxon
family_keep <-keep_taxa[[2]]$taxon
genus_keep <- keep_taxa[[1]]$taxon

rank_list = c("genus","family","order", "class","phylum")

ps_assigned_species = subset_taxa(ps, !grepl("species", ps@tax_table[,"species"]))
ps_assigned_genus = subset_taxa(ps, !grepl("genus|unassigned", ps@tax_table[,"genus"]))
all_rank_taxa_sample <- list()
set.seed(1)
for (k in 1:5) {
	to_keep <- keep_taxa[[k]]$taxon
	rank_name = rank_list[[k]]
rank_taxa_sample <- list()
for (i in 1:length(to_keep)){
	if (to_keep[[i]] %in% unique(ps_assigned_genus@tax_table[,rank_name])){
		# Subset to a single genus/phylum etc.
	temp_ps <- subset_taxa(ps_assigned_genus, eval(as.name(rank_name)) %in% to_keep[[i]])
	# Sample 2 ASVs from this taxon
	rank_taxa_sample[[i]] <- temp_ps %>% taxa_names() %>% sample(size = 5)
	names(rank_taxa_sample)[[i]] <- to_keep[[i]]
	}
}
all_rank_taxa_sample[[k]] <- stack(rank_taxa_sample) %>% mutate(rank=!!rank_name)
}

ASVs_for_phylogeny <- do.call(rbind, all_rank_taxa_sample)
colnames(ASVs_for_phylogeny) <- c("ASV", "taxon", "rank_only")

#####

#sequences <- refseq(ps_gen)
ps_subset <- subset_taxa(ps, taxa_names(ps) %in% ASVs_for_phylogeny$ASV)
sequences <- refseq(ps_subset)

cat("Aligning sequences...")
alignment <- AlignSeqs(sequences, anchor = NA, processors = NULL)
cat("Creating distance matrices...")
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
cat("Neigbour joining...")
treeNJ <- NJ(dm)
fit = pml(treeNJ, data = phang.align)
cat("GTR...")
fitGTR <- update(fit, k = 4, inv = 0.2)
cat("Done.")
tree <- fitGTR$tree
saveRDS(tree, output_path)

# Add tree back to phyloseq objcet
phy_tree(ps_subset) <- tree

# Melt longer, remove all the functional group assignments
melted_ps <- psmelt(ps_subset) %>% select(1:21)


fort <- fortify(tree) %>% dplyr::as_data_frame()
#fort_beta <- merge(fort, ASVs_betas_pred, by.x="label", by.y="ASV", all.x=T)

#ASVs_common <- intersect(ASVs_betas_pred$ASV, fort$label)

fort_tip <- fort %>% filter(isTip)


ps_subset@tax_table[,c("genus")] <- gsub("\\-|\\(|\\)| ", "\\.", ps_subset@tax_table[,c("genus")])
ps_subset@tax_table <- tax_table(cbind(ps_subset@tax_table[,1:6], label = taxa_names(ps_subset)))
tax <- as.data.frame(lapply(data.frame(ps_subset@tax_table), as.factor))[1:7]
tax_full <- as.data.frame(lapply(data.frame(ps@tax_table), as.factor))[1:7]
#tax$label <- rownames(ps_subset@tax_table)

# In case of crash.
save.image(here("data/phylo_workspace.Rdata"))



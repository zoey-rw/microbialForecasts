## Assigning taxonomy to NEON ITS dataset using the UNITE and UTOPIA databases.

library(ranomaly)
library(phyloseq)
library(dplyr)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

#### Prep legacy ASV table (from old processing) ####
library(dada2)

orig.seqtab <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/otuTable_legacy_fwdOnly.rds")


# Remove low-prevalence taxa and low-read count samples.
seqtab_subset <- orig.seqtab[which(rowSums(orig.seqtab) > 3000),]
present <- colSums(seqtab_subset > 0)
keep <- which(present > 2)
seqtab_subset <- seqtab_subset[,keep]
dim(orig.seqtab)
dim(seqtab_subset)
# Collapse ASVs from different runs
seqtab_collapse <- collapseNoMismatch(seqtab_subset, identicalOnly = T)
saveRDS(seqtab_collapse, "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_collapsed_legacy.rds")

out.file <- "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy_legacy"
dada_res <- list()
dada_res$seqtab.export <- seqtab_collapse
dada_res$seqtab.nochim <- seqtab_collapse


tax1 <- "/projectnb/microbiome/ref_db/taxonomy_db/unite_fungi_v8.2_idtaxa.rdata"
tax2 <- "/projectnb/microbiome/ref_db/taxonomy_db/utopia90_201908.rdata"
tax.table = assign_taxo_fun(dada_res = dada_res, id_db = c(tax1,tax2), 
														verbose = 3, output = out.file)


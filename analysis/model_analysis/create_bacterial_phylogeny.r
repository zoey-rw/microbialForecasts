source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

options(scipen=999)

library(phyloseq)
library(DECIPHER)
library(phangorn)
library(speedyseq)
library(tidyverse)
library(speedyseq)

library(phytools)
library(ggtree)
library(ape)


# Set output path
output_path <- here("data/clean/bacterial_phylogeny.rds")
# Read in data
ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")
# Remove weird characters from names
ps@tax_table[,c("family")] <- gsub("\\-|\\(|\\)| ", "\\.", ps@tax_table[,c("family")])
ps@tax_table[,c("genus")] <- gsub("\\-|\\(|\\)| ", "\\.", ps@tax_table[,c("genus")])


# Read & reshape model effect estimates
sum.all <- readRDS("./data/summary/all_fcast_effects.rds")
betas_wide = sum.all %>% filter(model_name=="all_covariates" & pretty_group=="Bacteria") %>%
	pivot_wider(id_cols=c("taxon","rank_only","time_period"),
							names_from = "beta", values_from = "Mean")

# Read in predictability scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))
# Do some fixing that should def have been done already
fcast_info_simple <- scores_list$scoring_metrics_site_lon %>%
	select(fcast_type, pretty_group, model_name, pretty_name, taxon) %>% distinct()
hindcast_rsq_unfilt <- scores_list$scoring_metrics %>%
	mutate(RSQ.1 = ifelse(RSQ.1 < 0, 0, RSQ.1),
				 CRPS_penalty = CRPS_truncated - MAE) %>%
	filter(model_name == "all_covariates") %>%
	distinct()  %>%
	merge(fcast_info_simple, all.x=T, all.y=F)

# Get list of hindcasted taxa
hindcasted_taxa =readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/tax_filter_pass.rds")
#keep_taxa = split(hindcasted_taxa, hindcasted_taxa$pretty_name)

# Get list of taxa at each rank
summary_unique = sum.all %>% distinct(pretty_name, taxon)
keep_taxa = split(summary_unique, summary_unique$pretty_name)

phylum_keep <- keep_taxa[[5]]$taxon
class_keep <- keep_taxa[[4]]$taxon
order_keep <- keep_taxa[[3]]$taxon
family_keep <-keep_taxa[[2]]$taxon
genus_keep <- keep_taxa[[1]]$taxon

rank_list = c("genus","family","order", "class","phylum")

ps_assigned_species = subset_taxa(ps, !grepl("species", ps@tax_table[,"species"]))
ps_assigned_genus = subset_taxa(ps, !grepl("genus|unassigned", ps@tax_table[,"genus"]))
all_rank_taxa_sample <- list()
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

# Merge ASV list with model effect estimates
ASVs_betas <- merge(ASVs_for_phylogeny, betas_wide %>%
											filter(time_period == "2015-11_2018-01"), all=T) %>%
	filter(rank_only !="functional")

# Merge ASV list with predictability scores
hindcast_rsq_unfilt_tax <- hindcast_rsq_unfilt %>% filter(fcast_type=="Taxonomic" & pretty_group=="Bacteria")
ASVs_betas_pred <- merge(ASVs_betas, hindcast_rsq_unfilt_tax, all=T)



fort <- fortify(tree) %>% dplyr::as_data_frame()
fort_beta <- merge(fort, ASVs_betas_pred, by.x="label", by.y="ASV", all.x=T)

#ASVs_common <- intersect(ASVs_betas_pred$ASV, fort$label)

fort_tip <- fort %>% filter(isTip)
fort_tip_beta <- merge(fort_tip, ASVs_betas_pred, by.x="label", by.y="ASV", all=T)
#rownames(fort_tip_beta) <- make.unique(fort_tip_beta$label)

# https://stackoverflow.com/questions/67684167/how-to-use-match-in-dplyr-to-order-a-column-based-on-an-external-vector
fort_tip_beta <- fort_tip_beta %>% summarize(.[match(tree$tip.label, label), ])

identical(fort_tip_beta$label, tree$tip.label)

library(phytools)

phylo_sig_results <- list()
for (beta in c("sin",
							 "cos", "Ectomycorrhizal\ntrees", "LAI", "pC", "pH", "Temperature",
							 "Moisture", "RSQ", "RSQ.1")) {
phylo_sig_results[[beta]] <- phylosig(tree, fort_tip_beta[,beta], test = T, nsim = 10000)

}

sig_for_plot = phylo_sig_results %>% do.call(rbind, .)
sig_for_plot[,"sim.K"] <- 10000
sig_for_plot <- sig_for_plot %>% as.data.frame() %>% rownames_to_column("beta")




library(ape)
ps_subset@tax_table[,c("genus")] <- gsub("\\-|\\(|\\)| ", "\\.", ps_subset@tax_table[,c("genus")])
ps_subset@tax_table <- tax_table(cbind(ps_subset@tax_table[,1:6], label = taxa_names(ps_subset)))
tax <- as.data.frame(lapply(data.frame(ps_subset@tax_table), as.factor))[1:7]
tax_full <- as.data.frame(lapply(data.frame(ps@tax_table), as.factor))[1:7]
#tax$label <- rownames(ps_subset@tax_table)

# In case of crash.
save.image("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/phylo_workspace.Rdata")
load("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/phylo_workspace.Rdata")

# Merge with dendrogram version of phylogeny
frm <- ~phylum/class/order/family/genus/label
tr <- as.phylo(frm, data=tax, collapse = F)
tr_fort <- fortify(tr) %>% dplyr::as_data_frame()


# Merge with model estimates
ASVs_betas_pred$label = ASVs_betas_pred$ASV
tree_toplot <- full_join(tr, ASVs_betas_pred)

# This is the command causing crashes
fort2 <- fortify(tree_toplot) %>% dplyr::as_data_frame()


not_tip <- merge(tr_fort %>% filter(!isTip),
								 ASVs_betas_pred %>% filter(time_period == "2015-11_2018-01"),
								 by.x="label", by.y=c("taxon"), all=T)
not_tip$label.y <- NULL
not_tip$taxon = not_tip$label
fort_tip = fort2 %>% filter(isTip)
out_tree = rbindlist(list(fort_tip, not_tip), fill=T)





moisture_sig = sig_for_plot  %>% filter(beta=="Moisture")

moisture_tree <- ggtree(out_tree, layout="roundrect")
moisture_tree +
	scale_fill_distiller(palette = "Spectral") +
	geom_tiplab(geom = 'label', hjust=.5) +
	geom_nodelab(geom='label', aes(subset = rank_only=="phylum", fill =Moisture)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="class", fill =Moisture)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="order", fill =Moisture)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="family", fill =Moisture)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="genus", fill =Moisture)) +
	ggtitle("Phylogenetic signal of moisture sensitivity", subtitle = paste0("K = ", moisture_sig$K,
																																					 "\nP = ", moisture_sig$P))


RSQ_sig = sig_for_plot %>% filter(beta=="RSQ")
RSQ_tree <- ggtree(out_tree, layout="roundrect")
RSQ_tree +
	#scale_fill_distiller() +
	geom_tiplab(geom = 'label', hjust=.5) +
	geom_nodelab(geom='label', aes(subset = rank_only=="phylum", fill =RSQ)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="class", fill =RSQ)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="order", fill =RSQ)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="family", fill =RSQ)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="genus", fill =RSQ)) +
	ggtitle("Phylogenetic signal of predictability (hindcast RSQ)", subtitle = paste0("K = ", RSQ_sig$K,
																																					 "\nP = ", RSQ_sig$P))

LAI_sig = sig_for_plot %>% filter(beta=="LAI")
LAI_tree <- ggtree(out_tree, layout="roundrect")
LAI_tree +
	scale_fill_distiller(palette = "Spectral") +
	geom_tiplab(geom = 'label', hjust=.5) +
	geom_nodelab(geom='label', aes(subset = rank_only=="phylum", fill =LAI)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="class", fill =LAI)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="order", fill =LAI)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="family", fill =LAI)) +
	geom_nodelab(geom='label', aes(subset = rank_only=="genus", fill =LAI)) +
	ggtitle("Phylogenetic signal of LAI sensitivity", subtitle = paste0("K = ", LAI_sig$K,
																																					 "\nP = ", LAI_sig$P))



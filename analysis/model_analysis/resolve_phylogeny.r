source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr)
library(tidytree)


# Takes a minute to load, huge workspace file
load("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/phylo_workspace.Rdata")


# Create dendrogram version of SILVA phylogeny
frm <- ~phylum/class/order/family/genus/label
tr <- as.phylo(frm, data=tax, collapse = F)
silva_fortify <- fortify(tr) %>% dplyr::as_data_frame()
silva_phylo_tree = as.treedata(tr)

# Create dendrogram version of sampled phylogeny
sampled_phylo_tree = as.treedata(tree)
sampled_phylo_treedata = sampled_phylo_tree %>% as_tibble
# This is removing 1 branch length. Idk why it is essential for this to work.
# silva_phylo_tree@phylo$edge.length <- head(silva_fortify$branch, 1272)

merged_tree <- merge_tree(sampled_phylo_tree, silva_phylo_tree)
merged_tree2 <- merge_tree(silva_phylo_tree, sampled_phylo_tree)

merged_tree2@phylo$tip.label = janitor::make_clean_names(merged_tree2@phylo$tip.label)
merged_tree2@phylo$node.label = janitor::make_clean_names(merged_tree2@phylo$node.label)



# detach("package:speedyseq", unload = TRUE)
# detach("package:phyloseq", unload = TRUE)
# tree_orig = tree
# y = tree_orig %>% as.treedata() %>% as.tibble()
#
# #ASVs_no_dupes = ASVs_for_phylogeny %>% distinct(ASV, .keep_all = T)
# y$rank = NA
# y$genus = NA
# node_no <- 802
#
# for (node_no in y$node){
# 	print(node_no)
# 	# Skip if it's a tip
# 	if (isTip(y, node_no)) {
# 		node_rank = "ASV"; y[y$node == node_no,]$rank = node_rank; next()
# 	}
#
# 	all_offspring <- offspring(y, node_no)
# 	tip_offspring <- all_offspring[isTip(all_offspring, all_offspring$node),]
#
# 	tip_tax_actual_info <- tax[match(tip_offspring$label, tax$label),]
# 	tip_tax_rank_info <- ASVs_no_dupes[match(tip_offspring$label, ASVs_no_dupes$ASV),]
#
#
# 	# If the node contains exactly 2 phyla, it represents a phylum-level split
# 	# if (length(unique(tip_tax_actual_info$phylum)) == 2) {
# 	# 	node_rank = "phylum"
# 	# }
#
# 	# No instead we'll assign if there's any variation at that rank
# 	if (length(unique(tip_tax_actual_info$phylum)) > 1) {
# 		node_rank = "phylum"
# 	} else if (length(unique(tip_tax_actual_info$class)) > 1) {
# 		node_rank = "class"
# 		#node_taxon = unique(tip_tax_actual_info$phylum)
# 	} else if (length(unique(tip_tax_actual_info$order)) > 1) {
# 		node_rank = "order"
# 	} else if (length(unique(tip_tax_actual_info$family)) > 1) {
# 		node_rank = "family"
# 	} else if (length(unique(tip_tax_actual_info$genus)) > 1) {
# 		node_rank = "genus"
# 	} else if (length(unique(tip_tax_actual_info$species)) > 1) {
# 		node_rank = "species"
# 		y[y$node == node_no,]$genus = unique(tip_tax_actual_info$genus)
# 	} else if (length(unique(tip_tax_actual_info$species)) == 1) {
# 		node_rank = "subspecies"
# 		y[y$node == node_no,]$genus = unique(tip_tax_actual_info$genus)
# 	} else {
# 		node_rank = NA
# 	}
# 	print(node_rank)
# 	print(tip_tax_actual_info)
# 	print(tip_tax_rank_info)
#
# 	y[y$node == node_no,]$rank = node_rank
#
# }
#
#
# # Create a numeric variable to represent taxonomic rank
# y$rank_no <- as.numeric(factor(y$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "species","subspecies", "ASV"), ordered =T))

subspecies_nodes = y[y$rank=="subspecies",]
genus_nodes = y[y$rank=="genus",]

genus_offspring_list <- list()
for (i in 1:length(subspecies_nodes$node)){
genus_offspring <- offspring(y, subspecies_nodes$node[[i]])
genus_name = genus_offspring %>% filter(!is.na(genus)) %>% select(genus) %>% unique %>% unlist %>% as.character
if (length(genus_name) ==0) next()
genus_offspring$genus <- genus_name
genus_offspring_list[[i]] <-  genus_offspring
}

genus_tips = rbindlist(genus_offspring_list) %>% select(label,genus_label = genus) %>% na.omit %>% unique

y_genus = merge(y,genus_tips, all.x=T)
back_to_tree_genus <- y_genus %>% as.treedata()

to_keep = y_genus %>% filter(rank=="ASV" & !is.na(genus_label)) %>% group_by(genus_label) %>% dplyr::slice(1)
species<-to_keep$label %>% unlist() %>% unique
genus_tree<-drop.tip(back_to_tree_phylo, setdiff(back_to_tree_phylo$tip.label, species));

genus_tree %>%
	#groupOTU("genus") %>%
	#groupOTU(back_to_tree@data$genus) %>%
	ggtree(aes(color=rank_no)) +
	#geom_tippoint() +
	geom_nodelab(geom = "label", aes(label = rank)) +
	scale_color_viridis_c() +
	theme(legend.position = "right") + layout_dendrogram()





# just checking
unique(y[,c("rank","rank_no")])
table(y$rank)

# DO branch lengths increase or decrease with taxonomic resolution? Increase, except for phyla->class..
plot(y$branch.length ~ y$rank_no)


back_to_tree <- y %>% as.treedata()
back_to_tree_phylo = as.phylo(back_to_tree)




grp = back_to_tree@data$genus



to_drop = y %>% filter(rank %in% c("ASV","subspecies") & is.na(genus)) %>% select(label)
to_drop = y %>% filter(rank %in% c("ASV","subspecies")) %>% select(label) %>% unlist() %>% unique
to_keep = y %>% filter(!rank %in% c("ASV","subspecies"))




grp = split(y$label, y$genus)
grp_clade = split(y$node, y$genus)
#grp_clade = grp_clade[which(leng)]

groupInfo <- split(back_to_tree@data$node, back_to_tree@data$genus)
plotting_tree <- groupOTU(plotting_tree, groupInfo)

plotting_tree %>%
	#groupOTU("genus") %>%
	#groupOTU(back_to_tree@data$genus) %>%
	ggtree(aes(color=rank_no)) +
	#geom_tippoint() +
	geom_nodelab(geom = "label", aes(label = rank)) +
	scale_color_viridis_c() +
	theme(legend.position = "right") + layout_dendrogram()


resolved_fort <- fortify(back_to_tree) %>% dplyr::as_data_frame()
resolved_fort_beta <- merge(resolved_fort, ASVs_betas_pred, by.x="label", by.y="ASV", all=T)
traits3 <- resolved_fort_beta %>% filter(isTip) %>%
	select(name = label, Moisture, Temperature, pH, pC, Ecto = `Ectomycorrhizal\ntrees`, LAI) %>%
	filter(!is.na(Moisture) & name %in% back_to_tree_phylo$tip.label) %>% mutate_at(2:7, as.numeric)
to_drop = back_to_tree_phylo$tip.label[!back_to_tree_phylo$tip.label %in% traits3$name]
back_to_tree_phylo <- drop.tip(back_to_tree_phylo, to_drop)

traits3_reordered <- traits3[match(back_to_tree_phylo$tip.label, traits3$name),]






#### Dec 5 trying again
merged_fort <- fortify(merged_tree2) %>% dplyr::as_data_frame()
merged_fort$label = janitor::make_clean_names(merged_fort$label)
ASVs_betas = ASVs_betas %>% filter(!grepl("other", taxon)) %>% mutate(ASV = janitor::make_clean_names(ASV))
merged_fort_beta <- merge(merged_fort, ASVs_betas, by.x="label", by.y="ASV", all.x=T)

# Working backwards to find the MRCA of each genus, etc
unique_taxa = merged_fort_beta %>% filter(!grepl("other", taxon)) %>% select(rank_only, taxon) %>% distinct(.keep_all = F)
keep_taxa = split(unique_taxa, unique_taxa$rank_only)
phylum_keep <- keep_taxa[[5]]$taxon %>% unique
class_keep <- keep_taxa[[4]]$taxon %>% unique
order_keep <- keep_taxa[[3]]$taxon %>% unique
family_keep <-keep_taxa[[2]]$taxon %>% unique
genus_keep <- keep_taxa[[1]]$taxon %>% unique

MRCA_node_list = list()
rank_list = c("genus","family","order", "class","phylum")
for (rank_no in c(1:5)) {
	taxon_list = keep_taxa[[rank_no]]$taxon
	MRCA_node_list_rank <- list()
	for (i in 1:length(taxon_list)){
		taxon_tips = merged_fort_beta[which(merged_fort_beta$taxon==taxon_list[[i]]),]
		taxon_MRCA <- MRCA(merged_tree2, taxon_tips$label)
		if(length(taxon_MRCA)==0) next()
		MRCA_node_list_rank[[i]] <- cbind.data.frame(taxon = taxon_list[[i]], #rank = rank_no,
																								 node = taxon_MRCA)
	}
	MRCA_node_list[[rank_no]] = rbindlist(MRCA_node_list_rank, fill = T)
}
MRCA_nodes = rbindlist(MRCA_node_list, fill = T)
MRCA_nodes <- merge(MRCA_nodes, unique_taxa)

# Now this can be merged with the aot output!

# to_keep = merged_fort_beta %>% filter(isTip) %>% group_by(genus_label) %>% dplyr::slice(1)
# species<-to_keep$label %>% unlist() %>% unique
# genus_tree<-drop.tip(back_to_tree_phylo, setdiff(back_to_tree_phylo$tip.label, species));

# Subset to the traits we are testing
traits <- merged_fort_beta %>% filter(isTip) %>% select(name = label, Moisture, Temperature) %>% na.omit %>% distinct(name, .keep_all = T)
traits <- apply(traits,2,as.character) %>% as.data.frame()

# Clean up for ph_aot function
#tree_for_aot = reorder(merged_tree2@phylo, "postorder")
tree_for_aot = merged_tree2@phylo
tree_for_aot$node.label = janitor::make_clean_names(tree_for_aot$node.label)
tree_for_aot$tip.label = janitor::make_clean_names(tree_for_aot$tip.label)
traits <- traits[order(match(traits$name, tree_for_aot$tip.label)),]

# Sanity check
identical(traits$name, tree_for_aot$tip.label)

# newick_tree_path = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/bacterial_phylogeny_newick.rds"
# tree_for_aot %>% ape::write.tree(file=newick_tree_path)
# tree_for_aot <- ape::read.tree(newick_tree_path)

# Run analysis of traits
aot_results <- ph_aot(traits = traits, phylo=tree_for_aot)

# Extract df of interest
res_out = aot_results$trait_conservatism
res_out = res_out %>% select(trait.name, node, name, ntaxa, n.nodes, percvaramongnodes, percvaratnode, ssamongnodes, sswithinnodes, contributionindex)
# Add taxon names
#res_out$taxon = MRCA_nodes[match(res_out$node, MRCA_nodes$node),]$taxon
#res$trait_conservatism$label = recode(res$trait_conservatism$node, !!!node_key, .missing = NA)

# Add rank data
tax_long = pivot_longer(tax, cols=1:7, names_to = "rank", values_to = "label") %>%
	mutate(label = janitor::make_clean_names(label))
res_out$name <- recode(res_out$name, "gammaproteobacteriaincertaesedis" = "gammaproteobacteria_incertae_sedis")
res_out$rank = tax_long[match(res_out$name, tax_long$label),]$rank
res_out$rank <- factor(res_out$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)



ggplot(res_out) + geom_jitter(aes(x = rank, y = contributionindex)) +
	facet_grid(~trait.name)



# Visualize phylogenetic signal and contribution index
ggplot(data=res_out %>% filter(!is.na(rank)),
			 aes(x = rank,y = contributionindex))  +
	geom_violin(draw_quantiles = c(.5), alpha = .5) +
	geom_jitter(aes(color = as.factor(trait.name)), size = 4, width=.2, alpha = .5) +
	labs(col = "Site", title = "Phylogenetic contribution index") +
	xlab("Rank")+
	facet_grid(#rows = vars(only_rank),
		rows = vars(trait.name), drop = T,
		scales = "free") +
	ylab(NULL)+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) +
	#geom_text(data = sig_for_plot, aes(x = 1.5, y = 0.01, label = paste0("K: ", K))) +
	stat_compare_means(label.y = 0)

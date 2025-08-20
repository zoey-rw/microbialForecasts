source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr)

# Takes a minute to load, huge workspace file
load("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/phylo_workspace.Rdata")


detach("package:speedyseq", unload = TRUE)
detach("package:phyloseq", unload = TRUE)

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



# edge.lengths = merged_tree@phylo$edge.length
# names(edge.lengths) = merged_tree@phylo$edge

# merged_tree %>%
# 	ggtree() +
# 	geom_nodelab(
# 		mapping = aes(
# 			x = branch),
# 		nudge_y = 0.36
# 	) +
# 	xlim(-.1, 4.5) +
# 	geom_tippoint() +
# 	scale_size_continuous(range = c(3, 10)) +
# 	geom_tiplab(
# 		offset = .14,
# 	) +
# 	geom_nodelab(geom = "label") +
# 	scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu")) +
# 	theme(legend.position = "right")

merged_fort <- fortify(merged_tree2) %>% dplyr::as_data_frame()
merged_fort_beta <- merge(merged_fort, ASVs_betas, by.x="label", by.y="ASV", all=T)
traits2 <- merged_fort_beta %>% filter(isTip) %>% select(name = label, Moisture, Temperature) %>% na.omit

output_path <- here("data/clean/bacterial_phylogeny_newick.rds")

tree2 <- reorder(tree2, "postorder")

tree2$node.label[[1]] <- "unlabeled"
tree2$node.label <- janitor::make_clean_names(tree2$node.label)


#tree2_filt <- tree2 %>% keep.tip(tip = traits2$name)
tree2_filt$node.label[[1]] <- "unlabeled"

# THIS fixes the segmentation fault????
tree2_filt$node.label <- janitor::make_clean_names(tree2_filt$node.label)

traits2 <- traits2 %>% filter(name %in% tree2_filt$tip.label) %>% filter(!duplicated(name))

tree2_filt %>% ape::collapse.singles() %>% ape::ladderize() %>% ape::write.tree(file=output_path)

tree2_filt %>% ape::write.tree(file=output_path)

tree2_newick <- ape::read.tree(output_path)
res2 <- ph_aot(traits2, phylo=tree2_newick)
merged_tree2@phylo$node.label = janitor::make_clean_names(merged_tree2@phylo$node.label)
res2 <- ph_aot(traits2, phylo=merged_tree2@phylo, ebl_unstconst=T)

res2$trait_conservatism$contributionindex


# tree <- readRDS(("data/clean/bacterial_phylogeny.rds"))
# output_path <- here("data/clean/bacterial_phylogeny_newick.rds")
# ape::write.tree(tree, file=output_path)

# traits <- fort_tip_beta %>% select(name = label, Moisture, Temperature)
# res <- ph_aot(traits, phylo=tree)

nodes_labels = merged_fort %>% select(node, label)

#res_out <- merge(res$trait_conservatism, nodes_labels, all=T, by="node")

tax_long = pivot_longer(tax, cols=1:7, names_to = "rank", values_to = "label")

# node_key = unstack(nodes_labels)
# node_key = nodes_labels$node
# names(node_key) = nodes_labels$label
res_out = res2$trait_conservatism
res_out = res_out %>% select(trait.name, node, name, ntaxa, n.nodes, percvaramongnodes, percvaratnode, ssamongnodes, sswithinnodes, contributionindex)


# Add taxon names
res_out$label = nodes_labels[match(res_out$node, nodes_labels$node),]$label
#res$trait_conservatism$label = recode(res$trait_conservatism$node, !!!node_key, .missing = NA)


# Add rank data
res_out$rank = tax_long[match(res_out$label, tax_long$label),]$rank

ggplot(res_out) + geom_jitter(aes(x = rank, y = contributionindex)) + facet_grid(~trait.name)

phylo_file2 <- tempfile()
phylo_str <- readLines(output_path)


res <- ph_aot(traits, phylo=tree)




res$trait_conservatism$contributionindex

fort <- fortify(tree) %>% dplyr::as_data_frame()
fort_beta <- merge(fort, ASVs_betas_pred %>% filter(time_period == "2015-11_2020-01"), by.x="label", by.y="ASV", all.x=T)

fort_tip <- fort %>% filter(isTip)
fort_tip_beta <- merge(fort_tip, ASVs_betas_pred %>% filter(time_period == "2015-11_2020-01"), by.x="label", by.y="ASV", all.x=T)
rownames(fort_tip_beta) <- make.unique(fort_tip_beta$label)
fort_tip_beta <- fort_tip_beta[tree$tip.label, ]
identical(fort_tip_beta$label, tree$tip.label)




silva_phylo_tree, sampled_phylo_tree

## root with random (first) ASV as outgroup
phy1 <- silva_phylo_tree %>% as.phylo() #%>% root(outgroup="ASV25", resolve.root = T, edgelabel = T)
phy1 <-  phy1 %>% root(outgroup="ASV25", resolve.root = T, edgelabel = T)  %>%
	phytools::force.ultrametric(method = "extend")
phy2 <- tree %>%
	as.phylo() %>% root(outgroup="ASV25", resolve.root = T, edgelabel = T)  %>%
	phytools::force.ultrametric(method = "extend")
## convert phylo objects to dendrograms
dnd1 <- as.dendrogram(phy1)
dnd2 <- as.dendrogram(phy2)
## rearrange in ladderized fashion
dnd1 <- ladder(dnd1)
dnd2 <- ladder(dnd2)
## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)



BioGeoBEARS::extend_tips_to_ultrametricize2
phytools::force.ultrametric(method = "extend")

## Okay trying a new approach:
# Go through each node from the sampled tree. For each node, create a subset of the tree. If all of the tips of that subset are in the same phylum, assign that taxon to a node (repeat at all depths)

phy2$edge

tree_orig = tree
#y <- sampled_phylo_treedata
y = tree_orig %>% as.treedata() %>% as.tibble()


ASVs_no_dupes = ASVs_for_phylogeny %>% distinct(ASV, .keep_all = T)

y$rank = NA
y$genus = NA
node_no <- 802

for (node_no in y$node){
	print(node_no)
# Skip if it's a tip
if (isTip(y, node_no)) {
	node_rank = "ASV"; y[y$node == node_no,]$rank = node_rank; next()
}

all_offspring <- offspring(y, node_no)
tip_offspring <- all_offspring[isTip(all_offspring, all_offspring$node),]

tip_tax_actual_info <- tax[match(tip_offspring$label, tax$label),]
tip_tax_rank_info <- ASVs_no_dupes[match(tip_offspring$label, ASVs_no_dupes$ASV),]


# If the node contains exactly 2 phyla, it represents a phylum-level split
# if (length(unique(tip_tax_actual_info$phylum)) == 2) {
# 	node_rank = "phylum"
# }

# No instead we'll assign if there's any variation at that rank
if (length(unique(tip_tax_actual_info$phylum)) > 1) {
	node_rank = "phylum"
} else if (length(unique(tip_tax_actual_info$class)) > 1) {
	node_rank = "class"
	#node_taxon = unique(tip_tax_actual_info$phylum)
} else if (length(unique(tip_tax_actual_info$order)) > 1) {
	node_rank = "order"
} else if (length(unique(tip_tax_actual_info$family)) > 1) {
	node_rank = "family"
} else if (length(unique(tip_tax_actual_info$genus)) > 1) {
	node_rank = "genus"
} else if (length(unique(tip_tax_actual_info$species)) > 1) {
	node_rank = "species"
	y[y$node == node_no,]$genus = unique(tip_tax_actual_info$genus)
} else if (length(unique(tip_tax_actual_info$species)) == 1) {
	node_rank = "subspecies"
	y[y$node == node_no,]$genus = unique(tip_tax_actual_info$genus)
} else {
	node_rank = NA
}
print(node_rank)
print(tip_tax_actual_info)
print(tip_tax_rank_info)

y[y$node == node_no,]$rank = node_rank

}

# Create a numeric variable to represent taxonomic rank
y$rank_no <- as.numeric(factor(y$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "species","subspecies", "ASV"), ordered =T))

subspecies_nodes = y[y$rank=="subspecies",]

# just checking
unique(y[,c("rank","rank_no")])
table(y$rank)

# DO branch lengths increase or decrease with taxonomic resolution? Increase, except for phyla->class..
plot(y$branch.length ~ y$rank_no)


back_to_tree <- y %>% as.treedata()
grp = back_to_tree@data$genus



to_drop = y %>% filter(rank %in% c("ASV","subspecies") & is.na(genus)) %>% select(label)
to_drop = y %>% filter(rank %in% c("ASV","subspecies")) %>% select(label) %>% unlist() %>% unique
to_keep = y %>% filter(!rank %in% c("ASV","subspecies"))

plotting_tree <- drop.tip(back_to_tree, to_drop)

to_keep = y %>% filter(!is.na(genus)) %>% group_by(genus) %>% dplyr::slice(1)
species<-to_keep$label %>% unlist() %>% unique
pruned.tree<-drop.tip(back_to_tree_phylo,back_to_tree_phylo$tip.label[-match(species, back_to_tree_phylo$tip.label)])


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

back_to_tree_phylo = as.phylo(back_to_tree)

resolved_fort <- fortify(back_to_tree) %>% dplyr::as_data_frame()
resolved_fort_beta <- merge(resolved_fort, ASVs_betas_pred, by.x="label", by.y="ASV", all=T)
traits3 <- resolved_fort_beta %>% filter(isTip) %>%
	select(name = label, Moisture, Temperature, pH, pC, Ecto = `Ectomycorrhizal\ntrees`, LAI) %>%
	filter(!is.na(Moisture) & name %in% back_to_tree_phylo$tip.label) %>% mutate_at(2:7, as.numeric)
to_drop = back_to_tree_phylo$tip.label[!back_to_tree_phylo$tip.label %in% traits3$name]
back_to_tree_phylo <- drop.tip(back_to_tree_phylo, to_drop)

traits3_reordered <- traits3[match(back_to_tree_phylo$tip.label, traits3$name),]

# Missing CRPS for certain ASVs such as:
# > ASVs_for_phylogeny[ASVs_for_phylogeny$ASV=="ASV62",]
# ASV             taxon rank_only
# 201 ASV62 hyphomicrobiaceae    family

# Missing env sensitivities for certain ASVs such as:
# resolved_fort_beta[which(resolved_fort_beta$label=="ASV469"),]


back_to_tree_phylo[back_to_tree_phylo$tip.label=="ASV9149"]

library(phylocomr)
res3 <- ph_aot(traits3_reordered, phylo=back_to_tree_phylo)
res3$trait_conservatism
res_out = res3$trait_conservatism
res_out = res_out %>% select(trait.name, node, name, ntaxa, n.nodes, percvaramongnodes, percvaratnode, ssamongnodes, sswithinnodes, contributionindex)

# Add taxon names
res_out$label = resolved_fort[match(res_out$node, resolved_fort$node),]$rank
# Add rank data
#res_out$rank = tax_long[match(res_out$label, tax_long$label),]$rank



#res_out$rank_no = as.numeric(res_out$label)
res_out$rank <- factor(res_out$label, levels = c(NA, "phylum", "class", "order", "family", "genus", "subspecies", "ASV"), ordered =T)

res_out$rank_no <- as.numeric(res_out$rank)

ASVs_betas_long <- ASVs_betas_pred %>% select(c("taxon","rank_only", "time_period","LAI", "pC", "pH", "Temperature",
																								"Moisture","Ectomycorrhizal\ntrees")) %>% distinct(.keep_all = T) %>%
	pivot_longer(cols = c("LAI", "pC", "pH", "Temperature",
																														 "Moisture","Ectomycorrhizal\ntrees"), names_to = "trait.name", values_to = "effect", values_drop_na = T)
ASVs_betas_long$effect = as.numeric(ASVs_betas_long$effect)
ASVs_betas_long$rank <- factor(ASVs_betas_long$rank_only, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)
# Visualize absolute size of site effects
ggplot(data=ASVs_betas_long %>% filter(time_period == "2015-11_2018-01"),
			 aes(x = rank,y = abs(effect))) +
	geom_violin(draw_quantiles = c(.5), alpha = .5) +
	geom_jitter(aes(color = as.factor(trait.name)), size = 4, width=.2, alpha = .5) +
	labs(col = "Site", title = "Effect size") +
	xlab("Rank")+
	facet_grid(#rows = vars(only_rank),
		rows = vars(trait.name), drop = T,
		scales = "free") +
	ylab(NULL)+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) + stat_compare_means(label.y = 0)

beta="Moisture"

library(phytools)
phylo_sig_results <- list()
for (beta in c(#"sin","cos",
							 #"Ectomycorrhizal\ntrees",
							 "Ecto",
							 "LAI", "pC", "pH", "Temperature",
							 "Moisture"
							 #"RSQ", "RSQ.1"
							 )) {
	phylo_sig_results[[beta]] <- phylosig(back_to_tree_phylo, traits3_reordered[,beta], test = T, nsim = 10000)
}

sig_for_plot = phylo_sig_results %>% do.call(rbind, .)
sig_for_plot[,"sim.K"] <- 10000
sig_for_plot <- sig_for_plot %>% as.data.frame() %>% rownames_to_column("trait.name")

# Visualize phylogenetic signal and contribution index
ggplot(data=res_out %>% filter(label !="ASV" & label !="subspecies"),
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
	# geom_smooth(aes(x = rank_no,y = abs(contributionindex))) +
	# stat_cor(aes(x = rank_no,y = abs(contributionindex)))

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr) # remotes::install_github("ropensci/phylocomr")


# Takes a minute to load, huge workspace file
load(here("data/phylo_workspace.Rdata"))

detach("package:speedyseq", unload = TRUE)
detach("package:phyloseq", unload = TRUE)

# Read & reshape model effect estimates
sum.all <- readRDS(here("data/summary/predictor_effects.rds"))
betas_wide = sum.all %>% filter(model_name=="env_cycl" & pretty_group=="Bacteria") %>%
	pivot_wider(id_cols=c("taxon","rank_only","time_period"),
							names_from = "beta", values_from = "Mean")


# Merge ASV list with model effect estimates
ASVs_betas <- merge(ASVs_for_phylogeny, betas_wide %>%
											filter(time_period == "2015-11_2018-01"), all=T) %>%
	filter(rank_only !="functional")

# Create dendrogram version of SILVA phylogeny
frm <- ~phylum/class/order/family/genus/label
tr <- as.phylo(frm, data=tax, collapse = F)
silva_fortify <- fortify(tr) %>% dplyr::as_data_frame()
silva_phylo_tree = as.treedata(tr)

# Create dendrogram version of sampled phylogeny
sampled_phylo_tree = as.treedata(tree)
sampled_phylo_treedata = sampled_phylo_tree %>% as_tibble
merged_tree <- merge_tree(silva_phylo_tree, sampled_phylo_tree)
merged_tree@phylo$tip.label = janitor::make_clean_names(merged_tree@phylo$tip.label)
merged_tree@phylo$node.label = janitor::make_clean_names(merged_tree@phylo$node.label)

merged_fort <- fortify(merged_tree) %>% dplyr::as_data_frame()
ASVs_betas = ASVs_betas %>% filter(!grepl("other", taxon))
ASVs_betas[!is.na(ASVs_betas$ASV), ]$ASV = janitor::make_clean_names(ASVs_betas[!is.na(ASVs_betas$ASV), ]$ASV)

# Add in hindcast scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))
scores_to_merge = scores_list$scoring_metrics %>% filter(#pretty_group="Bacteria" &
																												 	site_prediction=="New time (observed site)" &
																												 		model_name=="env_cycl" &
																												 		!pretty_name %in% c("Functional group", "Diversity")) %>%
	ungroup() %>%
	select(taxon, pretty_name, CRPS = CRPS_truncated, RSQ, RSQ.1)
ASVs_betas_scores <- merge(ASVs_betas, scores_to_merge , all.x=T) %>% distinct()
merged_fort_beta <- left_join(merged_fort, ASVs_betas_scores, by=c("label"="ASV"))
merged_fort_beta$label = janitor::make_clean_names(merged_fort_beta$label)





# ##### Working backwards to find the MRCA of each genus, etc
# ##### This approach didn't work either because the node indices are changed by the AOT command :(  ####
unique_taxa = merged_fort_beta %>% filter(!grepl("other", taxon)) %>% select(rank_only, taxon) %>% distinct(.keep_all = F)
keep_taxa = split(unique_taxa, unique_taxa$rank_only)

MRCA_node_list = list()
rank_list = c("genus","family","order", "class","phylum")
for (rank_no in c(1:5)) {
	taxon_list = keep_taxa[[rank_no]]$taxon
	MRCA_node_list_rank <- list()
	for (i in 1:length(taxon_list)){
		taxon_tips = merged_fort_beta[which(merged_fort_beta$taxon==taxon_list[[i]]),]
		taxon_MRCA <- MRCA(merged_tree, taxon_tips$label)
		if(length(taxon_MRCA)==0) next()
		MRCA_node_list_rank[[i]] <- cbind.data.frame(taxon = taxon_list[[i]], #rank = rank_no,
																								 node = taxon_MRCA)
	}
	MRCA_node_list[[rank_no]] = rbindlist(MRCA_node_list_rank, fill = T)
}
MRCA_nodes = rbindlist(MRCA_node_list, fill = T)
MRCA_nodes <- merge(MRCA_nodes, unique_taxa)

# Now this can be merged with the aot output!

#####

# to_keep = merged_fort_beta %>% filter(isTip) %>% group_by(genus_label) %>% dplyr::slice(1)
# species<-to_keep$label %>% unlist() %>% unique
# genus_tree<-drop.tip(back_to_tree_phylo, setdiff(back_to_tree_phylo$tip.label, species));

# PHYLOCOM ANALYSIS (contribution index)

# Subset to the traits we are testing
traits <- merged_fort_beta %>% filter(isTip) %>%
	select(name = "label", "LAI", "pC", "pH", "Temperature",
				 "Moisture",Ecto = "Ectomycorrhizal\ntrees", "CRPS", "RSQ", "RSQ.1") %>% na.omit %>%
	distinct(name, .keep_all = T)
traits <- apply(traits,2,as.character) %>% as.data.frame()

# Clean up for ph_aot function
#tree_for_aot = reorder(merged_tree2@phylo, "postorder")
tree_for_aot = merged_tree@phylo
tree_for_aot$node.label = janitor::make_clean_names(tree_for_aot$node.label)
tree_for_aot$tip.label = janitor::make_clean_names(tree_for_aot$tip.label)
to_drop <- setdiff(tree_for_aot$tip.label, traits$name)
tree_for_aot <- drop.tip(tree_for_aot, to_drop)

traits <- traits[order(match(traits$name, tree_for_aot$tip.label)),]

# Sanity check
identical(traits$name, tree_for_aot$tip.label)

# Run analysis of traits
aot_results <- ph_aot(traits = traits, phylo=tree_for_aot)

# Extract df of interest
res_out = aot_results$trait_conservatism
res_out = res_out %>% select(trait.name, node, name, ntaxa, n.nodes,
														 percvaramongnodes, percvaratnode, ssamongnodes, sswithinnodes, contributionindex)


# Add rank data
res_out$name <- recode(res_out$name, "gammaproteobacteriaincertaesedis" = "gammaproteobacteria_incertae_sedis")

tax_long = pivot_longer(tax, cols=1:7, names_to = "rank", values_to = "label") %>%
	mutate(label = janitor::make_clean_names(label))
res_out$rank = tax_long[match(res_out$name, tax_long$label),]$rank
res_out$rank <- factor(res_out$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)


# PHYLOSIG ANALYSIS - Blomberg's K (with phylotools)

# Sampled phylo tree
phy2 <- tree %>%
	as.phylo() %>% root(outgroup=1, resolve.root = T, edgelabel = T)  %>%
	phytools::force.ultrametric(method = "extend")
phy2$tip.label = janitor::make_clean_names(phy2$tip.label)


tree_for_k = phy2
# Sanity check
to_drop <- setdiff(tree_for_k$tip.label, traits$name)
tree_for_k <- drop.tip(tree_for_k, to_drop)
traits_reordered = traits[order(match(traits$name, tree_for_k$tip.label)),]
identical(traits_reordered$name, tree_for_k$tip.label)

library(phytools)
phylo_sig_results <- list()
for (beta in c(#"sin","cos",
	#"Ectomycorrhizal\ntrees",
	"Ecto",
	"LAI", "pC", "pH", "Temperature",
	"Moisture"
	#"RSQ", "RSQ.1"
)) {

	# Subset to one trait
	trait_vec = traits_reordered[,beta] %>% as.numeric
	names(trait_vec) <- tree_for_k$tip.label
	trait_vec <- trait_vec[!is.na(trait_vec)]

	# Prune tree
	to_drop = tree_for_k$tip.label[!tree_for_k$tip.label %in% names(trait_vec)]
	tree_for_k.prune <- drop.tip(tree_for_k, to_drop)

	# Sanity check
	identical(names(trait_vec), tree_for_k.prune$tip.label)

	# Run phylosig
	phylo_sig_results[[beta]] <- phylosig(tree = tree_for_k.prune, x = trait_vec, test = T, nsim = 10000)
}

sig_for_plot = phylo_sig_results %>% do.call(rbind, .)
sig_for_plot[,"sim.K"] <- 10000
sig_for_plot <- sig_for_plot %>% as.data.frame() %>% rownames_to_column("trait.name") %>% mutate(P = round(as.numeric(P), 4),
																																																 K = round(as.numeric(K), 3))

phylogenetic_results = res_out %>% filter(!is.na(rank)) %>%
	group_by(trait.name, rank) %>%
	summarize(mean_ContributionIndex = mean(contributionindex, na.rm=T),
						sd_ContributionIndex = sd(contributionindex, na.rm=T))



phylo_res_means = res_out %>% filter(!is.na(rank)) %>%
	group_by(trait.name, rank) %>%
	summarize(mean_ContributionIndex = mean(contributionindex, na.rm=T),
						sd_ContributionIndex = sd(contributionindex, na.rm=T),
						N_N=n(),
						se=sd_ContributionIndex/sqrt(N_N),
						upper_limit=mean_ContributionIndex+se,
						lower_limit=mean_ContributionIndex-se)


results_to_save = list(phylogenetic_results = phylogenetic_results,
											 sig_for_plot = sig_for_plot,
											 aot_results = aot_results,
											 tree_for_k = tree_for_k,
											 trait_data = merged_fort_beta,
											 res_out = res_out,
											 phylo_res_means = phylo_res_means)
											 #genus_treedata = genus_treedata)
saveRDS(results_to_save, here("data/summary/phylo_analysis_results.rds"))




results_to_save = readRDS(here("data/summary/phylo_analysis_results.rds"))
phylogenetic_results = results_to_save$phylogenetic_results
sig_for_plot = results_to_save$sig_for_plot
aot_results = results_to_save$aot_results
tree_for_k = results_to_save$tree_for_k
trait_data = results_to_save$merged_fort_beta




genus_tree <- genus_treedata %>%
	ggtree(aes(color=as.numeric(rank)), show.legend = F) +
	geom_nodelab(geom = "label", aes(label = label), show.legend = F) +
	scale_color_viridis_c() +
	theme(legend.position = "right")




genus_treedata2 <- root(genus_treedata, outgroup = 1, edgelabel = TRUE)

genus_treedata2 <- merged_fort_beta %>% as.treedata()

genus_treedata_pruned <- ggtree:::drop.tip(genus_treedata, "asv1301", trim.internal = T)

genus_treedata


to_keep = merged_fort_beta %>% unnest(cols = c(Temperature, Moisture, pH, pC, `Ectomycorrhizal\ntrees`, LAI,
																										sin, cos)) %>%
	filter(isTip & !is.na(label) & rank_only=="genus") %>% group_by(label) %>% dplyr::slice(1) %>% ungroup
#to_keep = y_genus %>% filter(rank=="ASV" & !is.na(genus_label)) %>% group_by(genus_label) %>% dplyr::slice(1)
species<-to_keep$label %>% unlist() %>% unique
genus_treedata <- merged_fort_beta %>% unnest(cols = c(Temperature, Moisture, pH, pC, `Ectomycorrhizal\ntrees`, LAI,
																											 sin, cos))  %>% filter(!isTip | label %in% species)
genus_treedata$rank = tax_long[match(genus_treedata$label, tax_long$label),]$rank
genus_treedata$rank <- factor(genus_treedata$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)

to_drop =  merged_fort_beta %>% filter(isTip & is.na(pH)) %>% select(label) %>% unlist
#genus_treedata_pruned <- ggtree:::drop.tip(as.treedata(genus_treedata), to_drop)

tree2 = tree_subset(merged_fort_beta %>% as.treedata(), 2, levels_back = 5)
tree2 = tree_subset(merged_fort_beta %>% as.treedata(), "phylum/class/order/family/genus")




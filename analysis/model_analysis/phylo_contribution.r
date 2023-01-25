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

merged_tree <- merge_tree(silva_phylo_tree, sampled_phylo_tree)

merged_tree@phylo$tip.label = janitor::make_clean_names(merged_tree@phylo$tip.label)
merged_tree@phylo$node.label = janitor::make_clean_names(merged_tree@phylo$node.label)

#### Dec 5 trying again
merged_fort <- fortify(merged_tree) %>% dplyr::as_data_frame()
ASVs_betas = ASVs_betas %>% filter(!grepl("other", taxon))
ASVs_betas[!is.na(ASVs_betas$ASV), ]$ASV = janitor::make_clean_names(ASVs_betas[!is.na(ASVs_betas$ASV), ]$ASV)

# Add in hindcast scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))
scores_to_merge = scores_list$scoring_metrics %>% filter(#pretty_group="Bacteria" &
																												 	site_prediction=="New time (observed site)" &
																												 		model_name=="all_covariates" &
																												 		!pretty_name %in% c("Functional group", "Diversity")) %>%
	ungroup() %>%
	select(taxon, pretty_name, CRPS = CRPS_truncated, RSQ, RSQ.1)
#scores_to_merge$taxon = janitor::make_clean_names(scores_to_merge$taxon)
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

# newick_tree_path = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/bacterial_phylogeny_newick.rds"
# tree_for_aot %>% ape::write.tree(file=newick_tree_path)
# tree_for_aot <- ape::read.tree(newick_tree_path)

# Run analysis of traits
aot_results <- ph_aot(traits = traits, phylo=tree_for_aot)

# Extract df of interest
res_out = aot_results$trait_conservatism
res_out = res_out %>% select(trait.name, node, name, ntaxa, n.nodes,
														 percvaramongnodes, percvaratnode, ssamongnodes, sswithinnodes, contributionindex)
# Add taxon names
#res_out$taxon = MRCA_nodes[match(res_out$node, MRCA_nodes$node),]$taxon
#res$trait_conservatism$label = recode(res$trait_conservatism$node, !!!node_key, .missing = NA)

# Add rank data
res_out$name <- recode(res_out$name, "gammaproteobacteriaincertaesedis" = "gammaproteobacteria_incertae_sedis")
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

results_to_save = list(phylogenetic_results = phylogenetic_results,
											 sig_for_plot = sig_for_plot,
											 aot_results = aot_results,
											 tree_for_k = tree_for_k,
											 trait_data = merged_fort_beta,
											 genus_treedata = genus_treedata)
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




phylo_res_means = res_out %>% filter(!is.na(rank)) %>%
	group_by(trait.name, rank) %>%
	summarize(mean_ContributionIndex = mean(contributionindex, na.rm=T),
						sd_ContributionIndex = sd(contributionindex, na.rm=T),
						N_N=n(),
						se=sd_ContributionIndex/sqrt(N_N),
						upper_limit=mean_ContributionIndex+se,
						lower_limit=mean_ContributionIndex-se)



# Visualize phylogenetic signal and contribution index
ggplot(data=res_out %>%  filter(!is.na(rank)),
			 aes(x = rank,y = contributionindex))  +
	geom_point(data=phylo_res_means, aes(y = mean_ContributionIndex, color = as.factor(trait.name)),
						 size = 4, width=.2, alpha = .5, show.legend = F) +
	geom_errorbar(inherit.aes = F, data=phylo_res_means, aes(x = rank, ymin=lower_limit, ymax=upper_limit), show.legend = F) +

	labs(title = "Phylogenetic contribution index") +
	xlab("Rank")+
	facet_grid(#rows = vars(only_rank),
		rows = vars(trait.name), drop = T,
		scales = "free") +
	ylab("Mean phylogenetic contribution of rank")+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) +
#	geom_text(data = sig_for_plot, aes(x = 4, y = 0.03, label = paste0("Blomberg's K: ", K))) +
#	geom_text(data = sig_for_plot, aes(x = 4, y = 0.02, label = paste0("p: ", P))) +
	#stat_compare_means(label.y = 0.001) +
	stat_compare_means(method = "anova", label.y = .003, label.x = 4)+ # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test", label.y = .001) +
	scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))


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
tax_long = pivot_longer(tax, cols=1:7, names_to = "rank", values_to = "label") %>%
	mutate(label = janitor::make_clean_names(label))
genus_treedata$rank = tax_long[match(genus_treedata$label, tax_long$label),]$rank
genus_treedata$rank <- factor(genus_treedata$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)

to_drop =  merged_fort_beta %>% filter(isTip & is.na(pH)) %>% select(label) %>% unlist
#genus_treedata_pruned <- ggtree:::drop.tip(as.treedata(genus_treedata), to_drop)
#genus_treedata_pruned <- drop.tip(genus_treedata2@phylo, to_drop) %>% as.treedata()

genus_tree <- genus_treedata %>% as.treedata %>%
	ggtree(aes(color=as.numeric(rank)), show.legend = F) +
	geom_nodelab(geom = "label", aes(label = label), show.legend = F) +
	scale_color_viridis_c() +
	theme(legend.position = "right")

genus_treedata %>%
	ggtree() +
	geom_nodelab(geom = "label", aes(label = label, fill=as.numeric(unlist(Temperature))), show.legend = F) +
	scale_fill_continuous(low = 'blue', high = 'red',
												na.value = "grey50")



tree2 = tree_subset(merged_fort_beta %>% as.treedata(), 2, levels_back = 5)
tree2 = tree_subset(merged_fort_beta %>% as.treedata(), "phylum/class/order/family/genus")




library(ggridges)
ggplot(res_out %>% filter(!is.na(rank) & trait.name != "RSQ.1"),
			 aes(x = contributionindex,
			 		y = rank,
			 		group= rank, fill = rank))  +
	geom_density_ridges2(aes(alpha=.3), scale = 3) + facet_grid(rows=vars(trait.name), scales = "free") +
	theme_bw(base_size = 18)

ggplot(res_out %>% filter(!is.na(rank) & trait.name != "RSQ.1"),
			 aes(x = contributionindex,
			 		y = rank,
			 		group= rank, fill = stat(x)))  +
	geom_density_ridges_gradient(alpha=.3, scale = 3, show.legend = F) +
	scale_fill_viridis_c(name = "contributionindex", option = "C")  +
	facet_grid(rows=vars(trait.name), scales = "free") +
	theme_ridges()  +
	scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))


betas_to_plot =  sum.all %>% filter(!grepl("other", taxon) & model_name=="all_covariates") %>%
	filter(time_period == "2015-11_2018-01") %>%
	filter(rank_only !="functional")

betas_to_plot_filter <- betas_to_plot %>% filter(!taxon %in% scores_list$unconverged_list$taxon)

ggplot(data=betas_to_plot %>% filter(!beta %in% c("sin","cos")),
			 aes(x = rev(pretty_name),y = effSize, color=beta))  +
	geom_jitter(size = 2, width=.5, alpha = .5, show.legend = F) +
	geom_violin( alpha = .5, show.legend = F) +
	xlab("Rank")+
	facet_grid(cols = vars(pretty_group),
		rows = vars(beta), drop = T,
		scales = "free") +
	ylab("Mean phylogenetic contribution of rank")+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) +
	stat_compare_means(method = "anova")+ # Add global p-value
	stat_compare_means(aes(label = after_stat(p.signif)),
										 method = "t.test")

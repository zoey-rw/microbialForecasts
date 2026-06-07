source("../../source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr) # remotes::install_github("ropensci/phylocomr")


# Load the slim inputs (tree, tax, ASVs_for_phylogeny) -- ~0.2 MB, committable to
# git for reproducibility. Fall back to the full 214 MB workspace only if the slim
# file is absent (e.g. an old checkout before create_bacterial_phylogeny.r was rerun).
slim_path <- here("data/clean/phylo_inputs_slim.rds")
if (file.exists(slim_path)) {
	phylo_inputs <- readRDS(slim_path)
	ASVs_for_phylogeny <- phylo_inputs$ASVs_for_phylogeny
	tree <- phylo_inputs$tree
	tax  <- phylo_inputs$tax
} else {
	message("Slim inputs not found; loading full workspace (slow).")
	load(here("data/phylo_workspace.Rdata"))
}

if ("package:speedyseq" %in% search()) detach("package:speedyseq", unload = TRUE)
if ("package:phyloseq" %in% search()) detach("package:phyloseq", unload = TRUE)

# Read & reshape model effect estimates
sum.all <- readRDS(here("data/summary/predictor_effects.rds"))
betas_wide = sum.all %>% filter(model_name=="env_cycl" & pretty_group=="Bacteria") %>%
	pivot_wider(id_cols=c("taxon","rank_only","time_period"),
							names_from = "beta", values_from = "Mean")


# Merge ASV list with model effect estimates
ASVs_betas <- merge(ASVs_for_phylogeny, betas_wide %>%
											filter(time_period == "2013-06_2018-01"), all=T) %>%
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

# Forecast-score columns (CRPS/RSQ) and their phylogenetic contribution are handled
# separately in phylo_contribution_scores.r, so the heavy scoring_metrics_plsr2.rds
# (26 MB) is not a dependency of the environmental-predictor results below.
merged_fort_beta <- left_join(merged_fort, ASVs_betas, by=c("label"="ASV"))
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
# Select traits for phylogenetic analysis — only require environmental betas (not scores)
traits <- merged_fort_beta %>% filter(isTip) %>%
	select(name = "label", "LAI", "pC", "pH", "Temperature",
				 "Moisture", Ecto = "Ectomycorrhizal\ntrees") %>%
	filter(!is.na(LAI) & !is.na(pH) & !is.na(Temperature) & !is.na(Moisture) & !is.na(pC) & !is.na(Ecto)) %>%
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

# Drop non-finite contribution indices BEFORE averaging. Some AOT nodes are
# degenerate: when a node's entire subtree shares one taxon's broadcasted value,
# the within-node sum-of-squares is zero and ph_aot returns -Inf/NaN. Averaging
# those in corrupts the rank means (previously the pH/order cell came out -Inf).
res_out_finite = res_out %>% filter(!is.na(rank) & is.finite(contributionindex))
n_nonfinite = res_out %>% filter(!is.na(rank)) %>%
	summarize(n = sum(!is.finite(contributionindex))) %>% pull(n)
cat("Dropped", n_nonfinite, "non-finite contributionindex node(s) before averaging\n")

phylogenetic_results = res_out_finite %>%
	group_by(trait.name, rank) %>%
	summarize(mean_ContributionIndex = mean(contributionindex),
						sd_ContributionIndex = sd(contributionindex), .groups = "drop")



phylo_res_means = res_out_finite %>%
	group_by(trait.name, rank) %>%
	summarize(mean_ContributionIndex = mean(contributionindex),
						sd_ContributionIndex = sd(contributionindex),
						N_N=n(),
						se=sd_ContributionIndex/sqrt(N_N),
						upper_limit=mean_ContributionIndex+se,
						lower_limit=mean_ContributionIndex-se, .groups = "drop")


# ============================================================================
# CLADE-SIZE-MATCHED NULL MODEL FOR THE CONTRIBUTION INDEX
# ----------------------------------------------------------------------------
# The contribution index is mechanically larger at deep nodes (they subtend more
# tips), and every taxon's representative ASVs share one broadcasted value, which
# forces among-taxon variance toward deep ranks. To test whether the observed
# deep-rank peak is real or just this geometry, we build a null that holds the
# tree and the within-taxon block structure FIXED and only permutes which signed
# effect-vector each taxon receives (a bijection over taxa). The standardized
# effect size (observed - null mean)/null sd then isolates any deep concentration
# BEYOND what clade-size geometry plus broadcasting produce by chance.
# ============================================================================

beta_cols = c("LAI", "pC", "pH", "Temperature", "Moisture", "Ecto")
rank_levels = c("phylum", "class", "order", "family", "genus")

# One signed effect-vector per taxon (the broadcast value). Mirror the column set
# fed to the real ph_aot so the unpermuted harness reproduces it exactly. Assert
# one row per taxon rather than slicing to dedupe.
taxon_betas = ASVs_betas %>%
	filter(!is.na(taxon) & !grepl("other", taxon)) %>%
	select(taxon, rank_only, "LAI", "pC", "pH", "Temperature", "Moisture",
				 Ecto = "Ectomycorrhizal\ntrees") %>%
	distinct()
stopifnot(!any(duplicated(taxon_betas$taxon)))

# Fixed map from tree tip -> taxon (labels were already clean-named upstream).
tip_taxon_map = merged_fort_beta %>% filter(isTip) %>%
	transmute(name = janitor::make_clean_names(label), taxon) %>%
	filter(!is.na(taxon)) %>% distinct(name, .keep_all = TRUE)

# Rebuild the ph_aot tip table from a taxon -> effect-vector assignment, exactly
# as the real pipeline does (join, NA-filter on the 6 env betas, distinct, order
# to the tree, coerce to character).
build_traits = function(taxon_beta_tbl) {
	tt = tip_taxon_map %>%
		left_join(taxon_beta_tbl, by = "taxon") %>%
		filter(if_all(all_of(beta_cols), ~ !is.na(.))) %>%
		distinct(name, .keep_all = TRUE) %>%
		select(name, all_of(beta_cols))
	tt = tt[tt$name %in% tree_for_aot$tip.label, , drop = FALSE]
	tt = tt[order(match(tt$name, tree_for_aot$tip.label)), , drop = FALSE]
	apply(tt, 2, as.character) %>% as.data.frame()
}

# Run ph_aot on a traits table and return mean contribution index per (trait, rank)
# for the 6 environmental predictors. Re-prunes + re-asserts tip ordering each call
# so any misalignment hard-fails instead of silently corrupting the index.
aot_rank_means = function(traits_tbl) {
	to_drop = setdiff(tree_for_aot$tip.label, traits_tbl$name)
	tr = if (length(to_drop)) drop.tip(tree_for_aot, to_drop) else tree_for_aot
	traits_tbl = traits_tbl[order(match(traits_tbl$name, tr$tip.label)), , drop = FALSE]
	stopifnot(identical(traits_tbl$name, tr$tip.label))
	ao = ph_aot(traits = traits_tbl, phylo = tr)
	ro = ao$trait_conservatism %>% select(trait.name, name, contributionindex)
	ro$name = recode(ro$name, "gammaproteobacteriaincertaesedis" = "gammaproteobacteria_incertae_sedis")
	ro$rank = tax_long[match(ro$name, tax_long$label),]$rank
	ro %>% filter(rank %in% rank_levels &
								is.finite(contributionindex) &
								trait.name %in% beta_cols) %>%
		group_by(trait.name, rank) %>%
		summarize(mean_CI = mean(contributionindex), .groups = "drop")
}

# Correctness gate: the unpermuted harness must reproduce the real pipeline.
stopifnot(setequal(build_traits(taxon_betas)$name, traits$name))
recon = aot_rank_means(build_traits(taxon_betas))
cat("\nCorrectness-gate reconstruction (should match phylo_res_means):\n")
print(recon %>% filter(trait.name %in% c("Temperature", "pH")) %>% arrange(trait.name, rank))
cat("Observed phylo_res_means for comparison:\n")
print(phylo_res_means %>% filter(trait.name %in% c("Temperature", "pH")) %>%
				select(trait.name, rank, mean_ContributionIndex) %>% arrange(trait.name, rank))

observed_means = phylo_res_means %>%
	transmute(trait.name = as.character(trait.name),
						rank = as.character(rank),
						observed = mean_ContributionIndex)

# Null distribution: permute the taxon -> effect-vector assignment NSIM_NULL times.
taxa_vec = taxon_betas$taxon
NSIM_NULL = 999   # finer p-resolution + stable SES (199 was too noisy near p~0.03)
set.seed(1)
null_dist = bind_rows(lapply(seq_len(NSIM_NULL), function(s) {
	perm_tbl = taxon_betas
	perm_tbl$taxon = sample(taxa_vec)
	stopifnot(!any(duplicated(perm_tbl$taxon)))
	aot_rank_means(build_traits(perm_tbl)) %>% mutate(sim = s)
}))

# Standardized effect size + one-sided permutation p per (trait, rank).
null_ses = null_dist %>%
	mutate(trait.name = as.character(trait.name), rank = as.character(rank)) %>%
	left_join(observed_means, by = c("trait.name", "rank")) %>%
	group_by(trait.name, rank) %>%
	summarize(observed = first(observed),
						null_mean = mean(mean_CI),
						null_sd = sd(mean_CI),
						null_lower = quantile(mean_CI, .025),
						null_upper = quantile(mean_CI, .975),
						n_null = n(),
						# One-sided tails plus a two-sided test (deviation from the null mean
						# in either direction). Phylum exceeds the null (upper tail); genus
						# falls below it (lower tail) -- the two-sided p covers both.
						p_greater = (sum(mean_CI >= first(observed)) + 1) / (n() + 1),
						p_less = (sum(mean_CI <= first(observed)) + 1) / (n() + 1),
						p_two_sided = (sum(abs(mean_CI - mean(mean_CI)) >=
															abs(first(observed) - mean(mean_CI))) + 1) / (n() + 1),
						.groups = "drop") %>%
	mutate(SES = (observed - null_mean) / null_sd,
				 rank = factor(rank, levels = rank_levels, ordered = TRUE)) %>%
	arrange(trait.name, rank)

cat("\nClade-size-matched null: SES and permutation p per (trait, rank):\n")
print(as.data.frame(null_ses))

# ============================================================================
# TAXON-LEVEL BLOMBERG'S K (de-broadcast)
# ----------------------------------------------------------------------------
# The ASV-level K above (sig_for_plot) is inflated: each taxon's coefficient is
# copied to ~5 sister ASV tips, creating zero within-taxon variance that reads
# as stronger-than-Brownian signal and shrinks the permutation p. Here we reduce
# to ONE representative tip per taxon (no identical blocks) and recompute K, so
# the signal is driven by closely-related TAXA rather than replicated ASVs.
# Reuses taxon_betas (per-taxon coefficient vectors) and tip_taxon_map from the
# null-model block above.
rep_tips = tip_taxon_map %>%
	filter(name %in% tree_for_k$tip.label) %>%
	arrange(taxon, name) %>%
	group_by(taxon) %>% slice_head(n = 1) %>% ungroup()   # deterministic representative
stopifnot(!any(duplicated(rep_tips$taxon)))               # exactly one tip per taxon
stopifnot(!any(duplicated(rep_tips$name)))
rep_tips = rep_tips %>%
	left_join(taxon_betas %>% select(taxon, all_of(beta_cols)), by = "taxon")

tree_for_k_taxon = drop.tip(tree_for_k, setdiff(tree_for_k$tip.label, rep_tips$name))
cat("\nTaxon-level K tree:", length(tree_for_k_taxon$tip.label),
		"representative tips /", length(unique(rep_tips$taxon)), "taxa\n")

phylo_sig_results_taxon = list()
for (beta in beta_cols) {
	tv = setNames(as.numeric(rep_tips[[beta]]), rep_tips$name)
	tv = tv[!is.na(tv)]
	if (length(tv) < 3) next
	tr = drop.tip(tree_for_k_taxon, setdiff(tree_for_k_taxon$tip.label, names(tv)))
	tv = tv[tr$tip.label]
	stopifnot(identical(names(tv), tr$tip.label))
	phylo_sig_results_taxon[[beta]] = phylosig(tree = tr, x = tv, test = T, nsim = 10000)
}
sig_for_plot_taxon = phylo_sig_results_taxon %>% do.call(rbind, .) %>% as.data.frame() %>%
	rownames_to_column("trait.name") %>%
	mutate(P = round(as.numeric(P), 4), K = round(as.numeric(K), 3), level = "taxon")

sig_for_plot$level = "ASV"
sig_compare = bind_rows(
	sig_for_plot %>% select(trait.name, K, P, level),
	sig_for_plot_taxon %>% select(trait.name, K, P, level)
) %>% arrange(trait.name, level)
cat("\nBlomberg's K: ASV-level (broadcast) vs taxon-level (de-broadcast):\n")
print(as.data.frame(sig_compare))

results_to_save = list(phylogenetic_results = phylogenetic_results,
											 sig_for_plot = sig_for_plot,
											 aot_results = aot_results,
											 tree_for_k = tree_for_k,
											 trait_data = merged_fort_beta,
											 res_out = res_out,
											 phylo_res_means = phylo_res_means,
											 null_dist = null_dist,
											 null_ses = null_ses,
											 sig_for_plot_taxon = sig_for_plot_taxon,
											 sig_compare = sig_compare,
											 # Saved so phylo_contribution_scores.r can reuse the AOT tree and
											 # the tip->taxon / rank lookups without rebuilding them.
											 tree_for_aot = tree_for_aot,
											 tax_long = tax_long,
											 tip_taxon_map = tip_taxon_map)
saveRDS(results_to_save, here("data/summary/phylo_analysis_results.rds"))
saveRDS(results_to_save, here("data/summary/phylo_analysis_results_env_cov.rds"))

# Readable summary of the null-model test (so results can be read without the RDS).
write.csv(null_ses, here("data/summary/phylo_contribution_null_ses.csv"),
					row.names = FALSE)
cat("Saved: data/summary/phylo_contribution_null_ses.csv\n")


results_to_save = readRDS(here("data/summary/phylo_analysis_results_env_cov.rds"))
results_to_save = readRDS(here("data/summary/phylo_analysis_results.rds"))
phylogenetic_results = results_to_save$phylogenetic_results
sig_for_plot = results_to_save$sig_for_plot
aot_results = results_to_save$aot_results
tree_for_k = results_to_save$tree_for_k
trait_data = results_to_save$merged_fort_beta







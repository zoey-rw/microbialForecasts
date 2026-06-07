# Phylogenetic contribution of FORECAST SCORES (CRPS, RSQ, RSQ relative to 1:1).
#
# Split out of phylo_contribution.r so the environmental-predictor results do not
# depend on the 26 MB scoring_metrics_plsr2.rds. This script reuses the AOT tree
# and the tip->taxon / rank lookups that phylo_contribution.r saves into
# phylo_analysis_results.rds, then asks the same question for the per-taxon
# forecast scores: is their variation phylogenetically structured across ranks?
#
# Run from analysis/model_analysis:  Rscript phylogeny/phylo_contribution_scores.r
source("../../source.R")
options(scipen = 999)
suppressMessages({library(dplyr); library(tidyr); library(phylocomr); library(ape)})

score_cols  <- c("CRPS", "RSQ", "RSQ.1")
rank_levels <- c("phylum", "class", "order", "family", "genus")

# Reuse the AOT tree + lookups built by phylo_contribution.r (no workspace needed).
res <- readRDS(here("data/summary/phylo_analysis_results.rds"))
tree_for_aot  <- res$tree_for_aot
tax_long      <- res$tax_long
tip_taxon_map <- res$tip_taxon_map
stopifnot(!is.null(tree_for_aot), !is.null(tax_long), !is.null(tip_taxon_map))

# Per-taxon forecast scores (env_cycl, new-time at observed sites; bacterial taxa).
scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
score_by_taxon <- scores_list$scoring_metrics %>%
	filter(site_prediction == "New time (observed site)" &
				 model_name == "env_cycl" &
				 !pretty_name %in% c("Functional group", "Diversity")) %>%
	ungroup() %>%
	transmute(taxon = species, CRPS = CRPS_truncated, RSQ, RSQ.1) %>%
	filter(if_all(all_of(score_cols), ~ !is.na(.))) %>%
	distinct()
stopifnot(!any(duplicated(score_by_taxon$taxon)))   # one score vector per taxon

# Tip score table from a taxon -> score assignment (mirrors the env pipeline).
build_score_traits <- function(tbl) {
	tt <- tip_taxon_map %>%
		left_join(tbl, by = "taxon") %>%
		filter(if_all(all_of(score_cols), ~ !is.na(.))) %>%
		distinct(name, .keep_all = TRUE) %>%
		select(name, all_of(score_cols))
	tt[tt$name %in% tree_for_aot$tip.label, , drop = FALSE]
}

# Tree pruned once to tips that carry scores; the null keeps this fixed.
obs_traits <- build_score_traits(score_by_taxon)
tree_score <- drop.tip(tree_for_aot, setdiff(tree_for_aot$tip.label, obs_traits$name))
cat("Score AOT tree:", length(tree_score$tip.label), "tips /",
		nrow(score_by_taxon), "scored taxa\n")

aot_score_means <- function(traits_tbl) {
	tt <- traits_tbl[traits_tbl$name %in% tree_score$tip.label, , drop = FALSE]
	tt <- tt[order(match(tt$name, tree_score$tip.label)), , drop = FALSE]
	stopifnot(identical(tt$name, tree_score$tip.label))
	ao <- ph_aot(traits = as.data.frame(apply(tt, 2, as.character)), phylo = tree_score)
	ro <- ao$trait_conservatism %>% select(trait.name, name, contributionindex)
	ro$name <- recode(ro$name, "gammaproteobacteriaincertaesedis" = "gammaproteobacteria_incertae_sedis")
	ro$rank <- tax_long[match(ro$name, tax_long$label), ]$rank
	ro %>% filter(rank %in% rank_levels &
								is.finite(contributionindex) &
								trait.name %in% score_cols) %>%
		group_by(trait.name, rank) %>%
		summarize(mean_CI = mean(contributionindex), .groups = "drop")
}

score_res_means <- aot_score_means(obs_traits)

# Clade-size-matched, two-sided null (same design as the env-predictor analysis):
# permute which score vector each taxon receives, holding tree + tip set fixed.
scored_taxa <- score_by_taxon$taxon
set.seed(1)
score_null <- bind_rows(lapply(seq_len(999), function(s) {
	perm <- score_by_taxon; perm$taxon <- sample(scored_taxa)
	stopifnot(!any(duplicated(perm$taxon)))
	aot_score_means(build_score_traits(perm)) %>% mutate(sim = s)
}))

obs <- score_res_means %>%
	transmute(trait.name = as.character(trait.name), rank = as.character(rank), observed = mean_CI)
score_null_ses <- score_null %>%
	mutate(trait.name = as.character(trait.name), rank = as.character(rank)) %>%
	left_join(obs, by = c("trait.name", "rank")) %>%
	group_by(trait.name, rank) %>%
	summarize(observed = first(observed),
						null_mean = mean(mean_CI), null_sd = sd(mean_CI),
						null_lower = quantile(mean_CI, .025), null_upper = quantile(mean_CI, .975),
						n_null = n(),
						p_two_sided = (sum(abs(mean_CI - mean(mean_CI)) >=
															abs(first(observed) - mean(mean_CI))) + 1) / (n() + 1),
						.groups = "drop") %>%
	mutate(SES = (observed - null_mean) / null_sd,
				 rank = factor(rank, levels = rank_levels, ordered = TRUE)) %>%
	arrange(trait.name, rank)

cat("\nForecast-score phylogenetic contribution (SES vs clade-size null):\n")
print(as.data.frame(score_null_ses))

saveRDS(list(score_res_means = score_res_means,
						 score_null = score_null,
						 score_null_ses = score_null_ses),
				here("data/summary/phylo_score_contribution.rds"))
write.csv(score_null_ses, here("data/summary/phylo_score_contribution_null_ses.csv"),
					row.names = FALSE)
cat("Saved: data/summary/phylo_score_contribution.rds + _null_ses.csv\n")

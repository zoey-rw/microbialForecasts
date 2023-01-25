source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr)

results_to_save = readRDS(here("data/summary/phylo_analysis_results.rds"))
phylogenetic_results = results_to_save$phylogenetic_results
sig_for_plot = results_to_save$sig_for_plot
aot_results = results_to_save$aot_results
tree_for_k = results_to_save$tree_for_k
trait_data = results_to_save$merged_fort_beta

ggplot(phylogenetic_results %>% filter(!is.na(rank) & trait.name != "RSQ.1"),
			 aes(x = mean_ContributionIndex,
			 		y = rank,
			 		color= trait.name))  +
	geom_point(position=position_dodge(width=1), size=4) +
	#facet_grid(rows=vars(trait.name), scales = "free") +
	theme_bw(base_size = 18)


ggplot(phylogenetic_results %>% filter(!is.na(rank) &
																			 	!trait.name %in% c("RSQ.1","CRPS")),
			 aes(x = as.numeric(rank),
			 		y = mean_ContributionIndex,
			 		color= trait.name))  +
	geom_line(size=2, #linetype=2,
						alpha = .5) +
	#facet_grid(rows=vars(trait.name), scales = "free") +
	theme_bw(base_size = 18) + scale_x_continuous(breaks=seq(1,5,1),
																								labels=c("phylum", "class", "order", "family",
																											 "genus")) +xlab("Taxonomic rank") +
	ylab("Mean phylogenetic contribution of rank")+ labs(color="Trait name")



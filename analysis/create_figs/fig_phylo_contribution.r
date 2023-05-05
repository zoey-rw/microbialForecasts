source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(phylocomr)
library(scales)
library(ggrepel)

results_to_save = readRDS(here("data/summary/phylo_analysis_results.rds"))
phylogenetic_results = results_to_save$phylogenetic_results
sig_for_plot = results_to_save$sig_for_plot
aot_results = results_to_save$aot_results
tree_for_k = results_to_save$tree_for_k
trait_data = results_to_save$merged_fort_beta
full_phylo_res = results_to_save$res_out
phylo_res_means = results_to_save$phylo_res_means

ggplot(phylogenetic_results %>% filter(!is.na(rank) & trait.name != "RSQ.1"),
			 aes(x = mean_ContributionIndex,
			 		y = rank,
			 		color= trait.name))  +
	geom_point(position=position_dodge(width=1), size=4) +
	#facet_grid(rows=vars(trait.name), scales = "free") +
	theme_bw(base_size = 18)


phylo_to_plot <- merge(full_phylo_res, phylo_res_means)

tukey_phylo_rank = phylo_to_plot %>%
	filter(!is.na(contributionindex) &
				 	!trait.name %in% c("RSQ.1","CRPS")) %>%
	group_by(trait.name) %>%
	summarize(tukey(x = rank, y = contributionindex, y.offset = 0))
tukey_phylo_rank$rank = tukey_phylo_rank$x %>%
	ordered(levels = rev(c("genus","family","order","class","phylum")))

tukey_phylo_rank <- merge(tukey_phylo_rank, phylo_res_means %>% select(upper_limit, rank, trait.name))
phylo_ci <- ggplot(phylo_res_means %>% filter(!is.na(rank) &
																			 	!trait.name %in% c("RSQ.1","CRPS")),
			 aes(x = as.numeric(rank),
			 		y = mean_ContributionIndex,
			 		color= trait.name))  +
	geom_point(size=2, #linetype=2,
						alpha = .5, show.legend = F) +
	geom_pointrange(aes(as.numeric(rank),
											ymax = upper_limit, ymin = lower_limit, color= trait.name),
									fatten = 3, linewidth=2, alpha=.5, show.legend = F) +
	facet_grid(rows=vars(trait.name), scales = "free") +
	theme_bw(base_size = 18) + scale_x_continuous(breaks=seq(1,5,1),
																								labels=c("phylum", "class", "order", "family",
																											 "genus")) +xlab("Taxonomic rank") +
	ylab("Mean phylogenetic contribution of rank")+ labs(color="Trait name") +
	geom_text(data = tukey_phylo_rank,
						aes(x = as.numeric(rank), y = upper_limit+.01, label = Letters_Tukey),
						show.legend = F, color = 1, size =5, inherit.aes = F)
phylo_ci

png(here("figures","phylo_ci.png"), width = 600, height=1000)
print(phylo_ci)

dev.off()





to_keep = merged_fort_beta %>% unnest(cols = c(Temperature, Moisture, pH, pC, `Ectomycorrhizal\ntrees`, LAI,
																							 sin, cos)) %>%
	filter(isTip & !is.na(label) & rank_only=="genus") %>% group_by(label) %>% dplyr::slice(1) %>% ungroup
#to_keep = y_genus %>% filter(rank=="ASV" & !is.na(genus_label)) %>% group_by(genus_label) %>% dplyr::slice(1)
species<-to_keep$label %>% unlist() %>% unique
genus_treedata <- merged_fort_beta %>% unnest(cols = c(Temperature, Moisture, pH, pC, `Ectomycorrhizal\ntrees`, LAI,
																											 sin, cos))  %>% filter(!isTip | label %in% species)
genus_treedata$rank = tax_long[match(genus_treedata$label, tax_long$label),]$rank
genus_treedata$rank <- factor(genus_treedata$rank, levels = c(NA, "phylum", "class", "order", "family", "genus", "ASV"), ordered =T)

genus_tree <- genus_treedata %>%  filter(!isTip) %>%
	ggtree(aes(color=as.numeric(rank)), show.legend = F) +
	geom_nodelab(geom = "label", aes(label = label), show.legend = F) +
	scale_color_viridis_c() +
	theme(legend.position = "right")

allrank_treedata <- merge( genus_treedata %>%  filter(!isTip) %>% select(-c("taxon","rank_only","time_period","Temperature", "Moisture",
																																						"pH", "pC", "Ectomycorrhizal\ntrees", "LAI", "sin", "cos", "pretty_name",
																																						"CRPS", "RSQ", "RSQ.1")), ASVs_betas_scores %>% select(-ASV) %>% rename("taxon"="label"), all=T)

tree_temperature <- allrank_treedata %>%
	filter(!is.na(Temperature)) %>%
	distinct() %>%
	ggtree(layout = "fan") + geom_label_repel(aes(label=label, fill=Temperature), color="white", size=6, max.overlaps = 100)  +
	scale_fill_gradient(low = muted('blue'),
											#high = muted('red'),
											high = 'red',
											na.value = "grey50") + labs(fill="Temperature effect") + theme(legend.title = element_text(size = 16))
tree_temperature



png(here("figures","phylogeny_temperature.png"), width = 1600, height=1000)
print(tree_temperature)

dev.off()




tree_pH <- allrank_treedata %>%
	filter(!is.na(pH)) %>%
	distinct() %>%
	ggtree(layout = "fan") + geom_label_repel(aes(label=label, fill=pH), color="white", size=6, max.overlaps = 100)  +
	scale_fill_gradient(low = muted('blue'),
											#high = muted('red'),
											high = 'red',
											na.value = "grey50") + labs(fill="pH effect") + theme(legend.title = element_text(size = 12))
tree_pH
png(here("figures","phylogeny_pH.png"), width = 1600, height=1000)
print(tree_pH)

dev.off()

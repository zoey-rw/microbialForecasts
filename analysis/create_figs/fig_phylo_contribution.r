source("source.R")
options(scipen=999)

library(treeio)
library(ggpubr)
library(ggtree)
library(scales)
library(ggrepel)

# Load phylo workspace for tree visualization variables (tax_long, ASVs_betas_scores, etc.)
load(here("data/phylo_workspace.Rdata"))
if ("package:speedyseq" %in% search()) detach("package:speedyseq", unload = TRUE)
if ("package:phyloseq" %in% search()) detach("package:phyloseq", unload = TRUE)

# Load phylo analysis results
results <- readRDS(here("data/summary/phylo_analysis_results.rds"))

phylogenetic_results = results$phylogenetic_results
sig_for_plot = results$sig_for_plot
aot_results = results$aot_results
tree_for_k = results$tree_for_k
trait_data = results$trait_data
full_phylo_res = results$res_out
phylo_res_means = results$phylo_res_means

# Filter to environmental traits only (exclude scoring metrics)
traits_to_keep <- c("Temperature", "Moisture", "pH", "pC", "LAI", "Ecto")
trait_labels <- c(Temperature = "Temperature", Moisture = "Moisture",
									pH = "pH", pC = "Percent C", LAI = "LAI",
									Ecto = "Ectomycorrhizal trees")

# Filter and set factor levels for consistent ordering
phylo_res_means <- phylo_res_means %>%
	filter(!is.na(rank) & trait.name %in% traits_to_keep) %>%
	mutate(trait.name = factor(trait.name, levels = traits_to_keep))

full_phylo_res <- full_phylo_res %>%
	filter(!is.na(rank) & trait.name %in% traits_to_keep) %>%
	mutate(trait.name = factor(trait.name, levels = traits_to_keep))

# Merge full results with means for Tukey test
phylo_to_plot <- merge(full_phylo_res, phylo_res_means)

safe_tukey <- function(x, y, y.offset = 0) {
	# Drop non-finite values and any ranks with fewer than 2 observations
	ok <- is.finite(y)
	x <- x[ok]; y <- y[ok]
	if (length(unique(x)) < 2) {
		maxs <- tibble(x = unique(x), tot = tapply(y, x, max) + y.offset * max(abs(y)))
		return(maxs %>% mutate(Letters_Tukey = "a"))
	}
	keep_levels <- names(which(table(x) >= 2))
	if (length(keep_levels) < 2) {
		maxs <- tibble(x = unique(x), tot = tapply(y, x, max) + y.offset * max(abs(y)))
		return(maxs %>% mutate(Letters_Tukey = "a"))
	}
	tryCatch(tukey(x = x, y = y, y.offset = y.offset),
					 error = function(e) {
					 	maxs <- tibble(x = unique(x), tot = tapply(y, x, max) + y.offset * max(abs(y)))
					 	maxs %>% mutate(Letters_Tukey = "a")
					 })
}

tukey_phylo_rank <- phylo_to_plot %>%
	filter(is.finite(contributionindex)) %>%
	group_by(trait.name) %>%
	reframe(safe_tukey(x = as.character(rank), y = contributionindex, y.offset = 0))
tukey_phylo_rank$rank <- tukey_phylo_rank$x %>%
	ordered(levels = rev(c("genus","family","order","class","phylum")))
tukey_phylo_rank <- merge(tukey_phylo_rank,
													phylo_res_means %>% select(upper_limit, rank, trait.name))
tukey_phylo_rank$trait.name <- factor(tukey_phylo_rank$trait.name, levels = traits_to_keep)

# Contribution index figure
phylo_ci <- ggplot(phylo_res_means,
			 aes(x = as.numeric(rank),
			 		y = mean_ContributionIndex,
			 		color = trait.name))  +
	geom_pointrange(aes(ymax = upper_limit, ymin = lower_limit),
									size = 1.2, linewidth = 1.5, show.legend = FALSE) +
	facet_wrap(~trait.name, scales = "free_y", ncol = 3,
						 labeller = labeller(trait.name = trait_labels)) +
	geom_text(data = tukey_phylo_rank,
						aes(x = as.numeric(rank), y = upper_limit + .005, label = Letters_Tukey),
						show.legend = FALSE, color = "black", size = 4, fontface = "italic",
						inherit.aes = FALSE) +
	scale_x_continuous(breaks = 1:5,
										 labels = c("Phylum", "Class", "Order", "Family", "Genus")) +
	scale_color_manual(values = c(Temperature = "#D55E00", Moisture = "#56B4E9",
																pH = "#0072B2", pC = "#009E73",
																LAI = "#E69F00", Ecto = "#CC79A7")) +
	xlab("Taxonomic rank") +
	ylab("Mean phylogenetic contribution index") +
	theme_bw(base_size = 14) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
				strip.background = element_rect(fill = "white"),
				strip.text = element_text(face = "bold", size = 12),
				panel.grid.minor = element_blank())

ggsave(here("figures","phylo_ci.png"), phylo_ci, width = 9, height = 7, dpi = 300)
cat("Saved: figures/phylo_ci.png\n")

# Tree visualization using the phylo tree object with trait annotations
library(ape)

# Unnest list-cols in trait_data for tips with effect sizes
tip_traits <- trait_data %>% filter(isTip) %>%
	unnest(cols = c(Temperature, Moisture, pH, pC, `Ectomycorrhizal\ntrees`, LAI, sin, cos)) %>%
	distinct(label, .keep_all = TRUE) %>%
	select(label, Temperature, pH, Ecto = `Ectomycorrhizal\ntrees`)

# Map tip labels to phylum using the workspace tax table
tax$label_clean <- janitor::make_clean_names(tax$label)
tip_traits$phylum <- tax$phylum[match(tip_traits$label, tax$label_clean)]
tip_traits$phylum <- tools::toTitleCase(as.character(tip_traits$phylum))

# Use tree_for_k (ultrametric, pruned to tips with traits)
tree_viz <- tree_for_k

# Phylum highlight colors — colorblind-friendly Wong palette
major_phyla <- c("Proteobacteria", "Actinobacteriota", "Acidobacteriota",
								 "Verrucomicrobiota", "Firmicutes", "Bacteroidota",
								 "Gemmatimonadota")
phylum_colors <- c(Proteobacteria = "#E69F00", Actinobacteriota = "#56B4E9",
									 Acidobacteriota = "#009E73", Verrucomicrobiota = "#D55E00",
									 Firmicutes = "#CC79A7", Bacteroidota = "#0072B2",
									 Gemmatimonadota = "#F0E442")

# Helper function to create a fan phylogeny plot with highlighted phylum clades
make_phylo_fan <- function(tree_obj, tip_df, trait_col, trait_label,
													 low_col = "#7b3294", high_col = "#d55e00") {
	# Prune to tips with this trait
	has_trait <- tip_df %>% filter(!is.na(.data[[trait_col]])) %>% pull(label)
	tree_pruned <- drop.tip(tree_obj, setdiff(tree_obj$tip.label, has_trait))
	plot_data <- tip_df %>% filter(label %in% tree_pruned$tip.label)

	# Identify MRCA nodes for major phyla
	phylum_counts <- plot_data %>% count(phylum_group = ifelse(phylum %in% major_phyla, phylum, NA)) %>%
		filter(!is.na(phylum_group) & n >= 5)
	clade_nodes <- list()
	for (i in seq_len(nrow(phylum_counts))) {
		phy <- phylum_counts$phylum_group[i]
		phy_tips <- plot_data %>% filter(phylum == phy) %>% pull(label)
		if (length(phy_tips) >= 2) {
			mrca_node <- getMRCA(tree_pruned, phy_tips)
			if (!is.null(mrca_node)) clade_nodes[[phy]] <- mrca_node
		}
	}

	# Build highlight data for geom_hilight
	highlight_df <- data.frame(
		node = unlist(clade_nodes),
		phylum = names(clade_nodes),
		stringsAsFactors = FALSE
	)

	# Base tree
	p <- ggtree(tree_pruned, layout = "fan", size = 0.3, color = "grey30")

	# Add phylum highlights as subtle background shading
	for (i in seq_len(nrow(highlight_df))) {
		p <- p + geom_hilight(node = highlight_df$node[i],
													fill = phylum_colors[highlight_df$phylum[i]],
													alpha = 0.15, extend = 0.03)
	}

	# Add tip points colored by effect size
	p <- p %<+% plot_data +
		geom_tippoint(aes(color = .data[[trait_col]]), size = 2) +
		scale_color_gradient2(low = muted(low_col), mid = "grey85", high = muted(high_col),
													midpoint = 0, na.value = "grey80",
													name = trait_label)

	# Add phylum legend via annotate trick — invisible points off-canvas
	legend_df <- data.frame(
		phylum = factor(names(clade_nodes), levels = major_phyla)
	)
	p <- p +
		geom_rect(data = legend_df,
							aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf, fill = phylum),
							inherit.aes = FALSE, show.legend = TRUE) +
		scale_fill_manual(name = "Phylum", values = phylum_colors,
											guide = guide_legend(order = 2,
																					 override.aes = list(alpha = 0.4))) +
		theme(legend.position = "right",
					legend.title = element_text(size = 11, face = "bold"),
					legend.text = element_text(size = 9),
					legend.box = "vertical",
					legend.key.size = unit(12, "pt"),
					plot.margin = margin(10, 10, 10, 10))

	p
}

# Temperature phylogeny fan plot
tree_temperature <- make_phylo_fan(tree_viz, tip_traits, "Temperature",
																	 "Temperature\neffect size")
ggsave(here("figures","phylogeny_temperature.png"), tree_temperature,
			 width = 12, height = 12, dpi = 300, bg = "white")
cat("Saved: figures/phylogeny_temperature.png\n")

# pH phylogeny fan plot
tree_pH <- make_phylo_fan(tree_viz, tip_traits, "pH", "pH\neffect size")
ggsave(here("figures","phylogeny_pH.png"), tree_pH,
			 width = 12, height = 12, dpi = 300, bg = "white")
cat("Saved: figures/phylogeny_pH.png\n")

# Ectomycorrhizal trees phylogeny fan plot
tree_ecto <- make_phylo_fan(tree_viz, tip_traits, "Ecto",
														"EM tree effect size")
ggsave(here("figures","phylogeny_ecto.png"), tree_ecto,
			 width = 12, height = 12, dpi = 300, bg = "white")
cat("Saved: figures/phylogeny_ecto.png\n")

# Print Blomberg's K results
cat("\nBlomberg's K results (phylogenetic signal in effect sizes):\n")
print(sig_for_plot %>% select(trait.name, K, P))

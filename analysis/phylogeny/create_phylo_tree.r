# Create phylogenetic tree for NEON bacteria 
# From 29 taxa
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

library(DECIPHER)
library(phangorn)
library(phyloseq)

# Set output path
output_path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phylo_tree.rds"

# Read in data
ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")

# # Phylogeny of genera
# rank.df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/genus_groupAbundances_16S.rds")
# genus_keep <- colnames(rank.df)[7:36]
# # Remove weird characters from names, and subset
# ps@tax_table[,c("genus")] <- gsub("\\-|\\(|\\)| ", "\\.", ps@tax_table[,c("genus")])
# ps_gen <- subset_taxa(ps, genus %in% genus_keep)
# 
# glommed <- tax_glom(ps_gen, taxrank="genus")
# 
# sequences <- refseq(glommed)
# #sequences <- refseq(ps_gen)
# 
# cat("Aligning sequences...")
# alignment <- AlignSeqs(sequences, anchor = NA, processors = NULL)
# cat("Creating distance matrices...")
# phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# dm <- dist.ml(phang.align)
# cat("Neigbour joining...")
# treeNJ <- NJ(dm)
# fit = pml(treeNJ, data = phang.align)
# cat("GTR...")
# fitGTR <- update(fit, k = 4, inv = 0.2)
# cat("Done.")
# tree <- fitGTR$tree
# saveRDS(tree, output_path)
# phy_tree(glommed) <- tree
# plot_tree(glommed, label.tips = "phylum", ladderize = "left", justify = "left" , size = "Abundance")
# plot_tree(glommed, nodelabf=nodeplotblank, label.tips="class", ladderize="left")







harv <- prune_samples(sample_data(ps)$siteID == "HARV",ps)
harv <- prune_taxa(taxa_sums(harv) > 0, harv)
harv <- prune_samples(sample_sums(harv) > 5000, harv)
harv <- prune_samples(grepl("DNA", sample_names(harv)), harv)
# harv_rel <- harv %>% speedyseq::transform_sample_counts(~ . / sum(.))
# 
# most_abundant_taxa <- sort(taxa_sums(harv_rel), TRUE)[1:50]
#harv_asv <- prune_taxa(names(most_abundant_taxa), harv_rel)


n.taxa <- 30
	#for (tax_rank in list("phylum","class","order","family","genus")){
	#ps.phy <- tax_glom(ps.filt, tax_rank)
	glom_melt <- speedyseq::psmelt(harv)
	form <- as.formula(paste0("sampleID ~ OTU"))
	glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
	out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
	
	most_abundant_taxa <- names(sort(colSums(out_abun>0), decreasing = T)[1:n.taxa])
	#most_abundant_taxa <- names(sort(colSums(out_abun), decreasing = T)[1:n.taxa])
	seqDepth <- rowSums(out_abun)
	out_top10 <- out_abun[,colnames(out_abun) %in% most_abundant_taxa]
	# Remove samples with a sum of above one (not sure why they exist)
	#out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]
	out_top10$Other <- seqDepth-rowSums(out_top10)
	out_top10 <- out_top10/rowSums(out_top10)
	ps.filt <- prune_samples(sample_names(harv) %in% rownames(out_top10), harv)
	rank.df <- cbind(sample_data(ps.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)

	
	# organize by date
	rank.df$dates <- as.Date(rank.df$dates, "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]

saveRDS(rank.df,"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/ASV_abundances_16S_HARV.rds")


ps_abun <- prune_taxa(taxa_names(ps.filt) %in% most_abundant_taxa,ps.filt)
# sequences <- refseq(ps_abun)
# 
# cat("Aligning sequences...")
# alignment <- AlignSeqs(sequences, anchor = NA, processors = NULL)
# cat("Creating distance matrices...")
# phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# dm <- dist.ml(phang.align)
# cat("Neigbour joining...")
# treeNJ <- NJ(dm)
# fit = pml(treeNJ, data = phang.align)
# cat("GTR...")
# fitGTR <- update(fit, k = 4, inv = 0.2)
# cat("Done.")
# tree <- fitGTR$tree
# saveRDS(tree, output_path)

tree <- readRDS(output_path)
phy_tree(ps_abun) <- tree



library(phytools)

phylosig(tree, beta_out[beta_out$beta=="Moisture",][1:30,]$Mean, test = T, nsim = 10000)

plot_tree(ps_abun, method = "treeonly", label.tips="class")
plot_tree(ps_abun, nodelabf=nodeplotblank, label.tips="class", ladderize="left", color = "horizon")
library(ggtree)

tax_table(ps_abun)



model.out <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_phylo_min5.rds")
model.sum <- model.out$summary$statistics



# # Prep model inputs/outputs.
# rank.df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/ASV_abundances_16S_HARV.rds")
# model.dat <- prepModelData(rank.df = rank.df, min.prev = 8)
# 
# # taxon name-number key.
# taxon_key <- colnames(model.dat$y)
# names(taxon_key) <- seq(1, length(taxon_key))
# taxon_key2 <- 1:length(colnames(model.dat$y))
# names(taxon_key2) <- colnames(model.dat$y)

beta_out <-  model.sum %>% as.data.frame() %>% 
	rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
	separate(rowname, sep=", ", into=c("taxon_num","beta_num")) %>% 
	mutate(taxon_num = as.numeric(gsub("beta\\[", "", taxon_num)),
				 beta_num = as.numeric(gsub("\\]", "", beta_num))) %>% 
	mutate(beta = recode(beta_num,
											 "1" = "Temperature",
											 "2" = "Moisture",
											 "3" = "pH",
											 "4" = "pC",
											 "5" = "Plant species richness",
											 "6" = "% grasses",
											 "7" = "% invasive species"))
beta_wide <- beta_out %>% pivot_wider(id_cols= taxon_num, names_from=beta, values_from=Mean)

ggtree(ps_abun) + geom_text2(aes(label=label), hjust=-.2, size=4) +
	geom_tiplab(aes(label=class), hjust=-.3) + 
	geom_point(data = beta_wide, aes( color=Temperature),na.rm=TRUE) 

library(metacoder)
x <- parse_phyloseq(ps_abun)
x$data$tax_data$phylum
heat_tree(x, node_label=otu_id, node_color= "green")

x$data$moisture <- mois$Mean[1:30]

x$data$tax_data$moisture <- mois$Mean[1:30]
x$data$taxon_id <- x$data$tax_data$otu_id
heat_tree(x,
					node_color = x$data$tax_data$moisture,
					#node_size = moisture,
					node_label = x$data$tax_data$otu_id
					#tree_label = otu_id
					)




to_plot <- x %>%
	taxa::filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names), reassign_obs = FALSE) 
to_plot$data$taxon_id <- to_plot$data$tax_data$otu_id
to_plot$data$tax_data$taxon_id <- to_plot$data$tax_data$otu_id
to_plot$data$tax_data$supertaxon_id <- to_plot$data$tax_data$otu_id
heat_tree(to_plot, 
					node_label = to_plot$data$taxon_id,
					# node_size = to_plot$data$tax_abund[['Nose']],
					# node_color = to_plot$data$tax_abund[['Nose']],
					layout = "davidson-harel", initial_layout = "reingold-tilford")

heat_tree(to_plot, 
					node_label = taxon_names)

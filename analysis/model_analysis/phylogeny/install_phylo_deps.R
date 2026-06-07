#!/usr/bin/env Rscript
# Install R package dependencies for the bacterial phylogenetic analysis.
#
#   "required" set  -> reproduce the downstream analysis from the committed slim
#                      inputs: phylo_contribution.r, phylo_contribution_scores.r,
#                      fig_phylo_contribution.r (null model, taxon-K, figures).
#   "from-raw" set  -> additionally rebuild the GTR tree from phyloseq_16S.rds
#                      via create_bacterial_phylogeny.r (Bioconductor + GitHub;
#                      heavy, rarely needed since the tree is committed in the
#                      slim inputs).
#
# Usage:
#   Rscript install_phylo_deps.R            # required (downstream) set only
#   Rscript install_phylo_deps.R --fromraw  # also install the from-raw set
#   (or set env MF_PHYLO_FROMRAW=true)

options(Ncpus = max(1L, parallel::detectCores() - 1L))
repos <- "https://cloud.r-project.org"
fromraw <- "--fromraw" %in% commandArgs(TRUE) ||
	tolower(Sys.getenv("MF_PHYLO_FROMRAW")) %in% c("1", "true", "yes")

have <- function(p) p %in% rownames(installed.packages())
need <- function(pkgs) pkgs[!vapply(pkgs, have, logical(1))]

cran_required <- c("phytools", "ape", "ggrepel", "ggpubr", "agricolae",
									 "scales", "janitor", "phylocomr")
bioc_required <- c("ggtree", "treeio")
cran_fromraw  <- c("phangorn")
bioc_fromraw  <- c("DECIPHER", "phyloseq")

cran <- c(cran_required, if (fromraw) cran_fromraw)
bioc <- c(bioc_required, if (fromraw) bioc_fromraw)

if (length(need(cran)))
	install.packages(need(cran), repos = repos, dependencies = NA)

if (length(need(bioc))) {
	if (!have("BiocManager")) install.packages("BiocManager", repos = repos)
	BiocManager::install(need(bioc), update = FALSE, ask = FALSE)
}

# speedyseq is optional (create_bacterial_phylogeny.r loads it only if present).
if (fromraw && !have("speedyseq")) {
	if (!have("remotes")) install.packages("remotes", repos = repos)
	tryCatch(remotes::install_github("mikemc/speedyseq", upgrade = "never"),
					 error = function(e)
					 	message("speedyseq (optional) not installed: ", conditionalMessage(e)))
}

still_missing <- need(c(cran, bioc))
if (length(still_missing))
	stop("Failed to install: ", paste(still_missing, collapse = ", "))
cat("Phylogeny dependencies installed",
		if (fromraw) "(required + from-raw)." else "(required/downstream set).", "\n")

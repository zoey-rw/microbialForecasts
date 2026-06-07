#!/usr/bin/env Rscript
# Download the large input files that are NOT stored in git, needed to reproduce
# the bacterial phylogenetic analysis "from raw". The small inputs and all
# results are committed in the repository, so most users do not need this script:
#
#   * Reproduce the null model, taxon-level K, and figures  -> nothing to download
#     (uses committed data/clean/phylo_inputs_slim.rds + data/summary/*.rds).
#   * Rebuild the GTR tree from scratch (create_bacterial_phylogeny.r)
#                                                          -> need phyloseq_16S.rds
#   * Recompute the forecast-score contribution (phylo_contribution_scores.r)
#                                                          -> need scoring_metrics_plsr2.rds
#
# Setup: after the Zenodo deposit is published, set the record id below (or the
# MF_ZENODO_BASE environment variable) to the record's file base URL, e.g.
#   export MF_ZENODO_BASE="https://zenodo.org/records/1234567/files"

suppressMessages(library(here))

ZENODO_BASE <- Sys.getenv(
	"MF_ZENODO_BASE",
	"https://zenodo.org/records/REPLACE_WITH_RECORD_ID/files"   # <-- edit after upload
)

files <- list(
	list(name = "phyloseq_16S.rds",
			 dest = here("data", "clean", "phyloseq_16S.rds"),
			 md5  = "0bdfa98e04eeec523fe6cfb253319899",
			 note = "16S phyloseq object; input to create_bacterial_phylogeny.r"),
	list(name = "scoring_metrics_plsr2.rds",
			 dest = here("data", "summary", "scoring_metrics_plsr2.rds"),
			 md5  = "2d712c9919f3b22fa15a6a8f6f72d044",
			 note = "forecast scores; input to phylo_contribution_scores.r")
)

if (grepl("REPLACE_WITH_RECORD_ID", ZENODO_BASE))
	stop("Set MF_ZENODO_BASE (or edit ZENODO_BASE) to the published Zenodo record URL first.")

for (f in files) {
	dir.create(dirname(f$dest), showWarnings = FALSE, recursive = TRUE)
	if (file.exists(f$dest) && unname(tools::md5sum(f$dest)) == f$md5) {
		cat("OK (already present):", f$name, "\n"); next
	}
	url <- file.path(ZENODO_BASE, paste0(f$name, "?download=1"))
	cat("Downloading", f$name, "->", f$dest, "\n")
	download.file(url, f$dest, mode = "wb")
	got <- unname(tools::md5sum(f$dest))
	if (got != f$md5)
		stop("md5 mismatch for ", f$name, ": got ", got, ", expected ", f$md5)
	cat("  verified md5", f$md5, "--", f$note, "\n")
}
cat("All Zenodo-hosted input files present and verified.\n")

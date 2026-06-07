# Bacterial phylogenetic analysis — reproducibility

Tests whether the environmental effect sizes estimated by the forecasting models
are phylogenetically structured across taxonomic ranks (Fig S4 in the manuscript).

## What is in git vs. on Zenodo

Everything needed to reproduce the **headline results and figures** is committed
in this repository (all small). The two large *raw* inputs live on Zenodo.

| File | Size | Where | Needed for |
|------|------|-------|-----------|
| `data/clean/phylo_inputs_slim.rds` | 23 KB | **git** | tree + tax + ASV map (downstream analysis) |
| `data/summary/predictor_effects.rds` | 107 KB | **git** | model effect sizes (input) |
| `data/summary/phylo_analysis_results.rds` | 660 KB | **git** | null model + K results |
| `data/summary/phylo_contribution_null_ses.csv` | 4 KB | **git** | per-rank SES + two-sided p |
| `data/summary/phylo_score_contribution_null_ses.csv` | small | **git** | forecast-score contribution |
| `data/clean/phyloseq_16S.rds` | 23 MB | **Zenodo** | rebuild GTR tree from raw |
| `data/summary/scoring_metrics_plsr2.rds` | 26 MB | **Zenodo** | forecast-score script |

The full 214 MB `data/phylo_workspace.Rdata` is **not** needed or archived — it is
a `save.image()` dump dominated by the 3.8 GB in-memory phyloseq object; the slim
file holds the only three objects the analysis uses.

## Reproduce (no downloads needed)

```bash
# 1. install dependencies (downstream/required set)
Rscript analysis/model_analysis/phylogeny/install_phylo_deps.R
# 2. null model + taxon-level Blomberg's K (uses committed slim inputs)
cd analysis/model_analysis && Rscript phylogeny/phylo_contribution.r
# 3. figures (Fig S4 = figures/phylo_ci_ses.png)
cd ../.. && Rscript analysis/create_figs/fig_phylo_contribution.r
```

`set.seed(1)` makes the 999-permutation null and the figures deterministic.

## Optional: from-raw / forecast scores (needs Zenodo files)

```bash
export MF_ZENODO_BASE="https://zenodo.org/records/<RECORD_ID>/files"
Rscript download_data.R                                   # fetch + md5-verify the 2 files
Rscript analysis/model_analysis/phylogeny/install_phylo_deps.R --fromraw
# rebuild the GTR tree from the 16S phyloseq (regenerates phylo_inputs_slim.rds):
cd analysis/model_analysis && Rscript phylogeny/create_bacterial_phylogeny.r
# forecast-score phylogenetic contribution:
Rscript phylogeny/phylo_contribution_scores.r
```

## Zenodo upload manifest

Upload exactly these two files to the Zenodo deposit, then put the published
record id into `download_data.R` (or the `MF_ZENODO_BASE` env var) and into the
manuscript Data Availability statement:

| File to upload | md5 |
|----------------|-----|
| `data/clean/phyloseq_16S.rds` | `0bdfa98e04eeec523fe6cfb253319899` |
| `data/summary/scoring_metrics_plsr2.rds` | `2d712c9919f3b22fa15a6a8f6f72d044` |

`download_data.R` checks these md5s after download, so update them there if the
files are ever regenerated.

## Environment

`install_phylo_deps.R` lists the exact packages (CRAN + Bioconductor `ggtree`/
`treeio`; `phylocomr` from CRAN). The Docker image (`docker/Dockerfile`) runs it
at build time, so `microbial-forecasts:latest` reproduces this analysis. Recorded
package versions are in `SESSIONINFO.txt`.

# microbialForecasts

Hierarchical Bayesian state-space models that forecast soil microbial relative
abundances from [NEON](https://www.neonscience.org/) data, fit with
[NIMBLE](https://r-nimble.org/) (MCMC). This repository holds the analysis code,
the `microbialForecast` R package, and the small inputs/results needed to
reproduce the downstream analyses and figures for the associated manuscript
(*Nature Communications*, in revision).

## Repository layout

| Path | Contents |
|------|----------|
| `source.R` | Environment setup; sourced by every analysis/figure script |
| `microbialForecast/` | R package: data prep, MCMC helpers, summarisation, hindcasting |
| `analysis/model_analysis/` | Numbered pipeline (00–10) + `phylogeny/` + `hpc/` job scripts |
| `analysis/create_figs/` | Figure-generating scripts (manuscript + supplement) |
| `data_construction/` | Raw-data ingestion and covariate preparation |
| `data/` | Inputs (`clean/`), MCMC outputs, summaries; most large files are gitignored |
| `figures/` | Generated figures (output directory) |
| `docker/` | Reproducibility image + container README (see below) |
| `scripts/` | Helper scripts (e.g. `run_all_figures.sh`) |
| `download_data.R` | Fetches the large Zenodo-hosted inputs (md5-verified) |

## Reproducing the analysis

The recommended path is the **single Docker image**, which pins the software
environment and bakes in every git-committed input. See
[`docker/README.md`](docker/README.md) for build and run instructions.

Large inputs that are not in git live on Zenodo and are fetched by
`download_data.R` (the Docker entrypoint runs this automatically when
`MF_ZENODO_BASE` is set):

```bash
export MF_ZENODO_BASE="https://zenodo.org/records/<RECORD_ID>/files"
Rscript download_data.R     # downloads + md5-verifies inputs into data/
```

For a local (non-Docker) run, install the package and its dependencies, then run
the numbered pipeline in `analysis/model_analysis/` and the figure scripts in
`analysis/create_figs/`. The full MCMC fits (step 01, ~100k iterations) require
an HPC cluster; the downstream steps and figures run on a workstation.

## Analysis pipeline (run order)

Input construction (`data_construction/`) builds the cleaned model inputs in
`data/clean/` from NEON amplicon abundances and environmental covariates; the raw
NEON downloads are external (NEON Data Portal) and the large derived inputs are on
Zenodo (`download_data.R`). The model-analysis pipeline then runs in numbered
order from `analysis/model_analysis/` (each script begins with
`source("../../source.R")`):

| Step | Script | Purpose |
|------|--------|---------|
| 00 | `00_createInputDF.r` | Assemble per-group input data frames |
| 01 | `01_fitModels.R` | Fit the hierarchical state-space models (primary cloglog Beta; HPC, ~100k iterations) |
| 02 | `02_combineModelChains.r` | Combine MCMC chains across runs |
| 03 | `03_summarizeModelOutputs.r` | Convergence diagnostics (Gelman–Rubin) and parameter summaries |
| 04 | `04_tidyEffectSizes.r` | Extract and tidy predictor effect sizes |
| 05 | `05_predictSiteEffects.r` | Predict site-level random effects for unobserved sites |
| 06 | `06_createHindcasts_observed.r`, `06_createHindcasts_newsites.r` | Generate hindcasts at observed and new sites |
| 07 | `07_tidyHindcasts.r` | Tidy hindcast outputs |
| 08 | `08_calculateScoringMetrics.r` | Scoring metrics (CRPS, nRMSE, R²) |
| 09 | `09_assignPeakPhenophase.r` | Assign peak phenophase from MODIS land-cover dynamics |
| 10 | `10_calculateFcastHorizon.r` | Estimate per-taxon forecast horizon |
| 11 | `11_siteEffectVariogram.r` | Test residual spatial autocorrelation in site effects |

Scripts with a `_CLR`, `_dirichlet`, or `_truncNorm` suffix repeat a step for the
alternative observation models compared in Appendix S3; the unsuffixed scripts are
the primary (Beta-regression) pipeline. The figure scripts in
`analysis/create_figs/` (run from the repo root) produce the manuscript and
supplement figures once the pipeline outputs exist.

## Data availability

Small inputs and all committed results are in this repository. The larger inputs
(phyloseq objects, forecast scores, tidy summaries, example hindcasts) are
deposited on Zenodo and retrieved via `download_data.R`; the published record id
is set there (or via `MF_ZENODO_BASE`) at release time.

## License

MIT — see [`LICENSE`](LICENSE).

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

## Data availability

Small inputs and all committed results are in this repository. The larger inputs
(phyloseq objects, forecast scores, tidy summaries, example hindcasts) are
deposited on Zenodo and retrieved via `download_data.R`; the published record id
is set there (or via `MF_ZENODO_BASE`) at release time.

## License

MIT — see [`LICENSE`](LICENSE).

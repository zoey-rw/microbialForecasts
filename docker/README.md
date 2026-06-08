# microbialForecasts — Docker reproducibility image

A single image that reproduces the downstream analyses and figures. It pins the
software environment and bakes in every git-committed input; the large
Zenodo-hosted inputs are fetched at runtime.

- **Base / pinning:** `rocker/tidyverse:4.4.2` (R 4.4.2, linux/amd64) with its
  dated Posit Package Manager CRAN snapshot, so package versions are reproducible.
  GitHub dependencies are pinned to commits (`spectratrait` v1.2.6, `speedyseq`
  v0.2.0). The `microbialForecast` package and the extra analysis/phylogeny
  dependencies are installed at build time.
- **Baked-in data:** all git-committed inputs under `data/clean` and
  `data/summary` (see `.dockerignore`). Heavy/uncommitted trees
  (`data/model_outputs`, `data/hindcasts`, the hindcast parquet, megapit dumps)
  and non-runtime dirs (`manuscript/`, `archive3/`, `.claude/`) are excluded.

## Build (from the repository root)

```bash
docker build -f docker/Dockerfile -t microbial-forecasts:latest .
```

The `FROM` line pins `linux/amd64` so Apple Silicon builds without extra flags
(BuildKit may print a harmless `FromPlatformFlagConstDisallowed` warning). For a
native arm64 image, switch the base to `rocker/tidyverse:latest` and drop the
`--platform` flags below.

## Smoke test (no data needed)

```bash
docker run --rm --platform=linux/amd64 microbial-forecasts:latest
```

Loads the project stack and confirms the key packages import.

## Reproduce

The **phylogenetic-signal analysis and figure need no external data** (their
inputs are committed and baked into the image):

```bash
docker run --rm --platform=linux/amd64 \
  -v "$(pwd)/figures:/opt/microbialForecasts/figures" \
  microbial-forecasts:latest \
  Rscript analysis/create_figs/fig_phylo_contribution.r
```

For analyses that need the **large Zenodo-hosted inputs**, set `MF_ZENODO_BASE`;
the entrypoint runs `download_data.R` (md5-verified) before your command:

```bash
docker run --rm --platform=linux/amd64 \
  -e MF_ZENODO_BASE="https://zenodo.org/records/<RECORD_ID>/files" \
  -v "$(pwd)/figures:/opt/microbialForecasts/figures" \
  microbial-forecasts:latest \
  Rscript analysis/create_figs/fig_exampleHindcasts.r
```

## How data reaches the container

1. **Baked in:** the git-committed inputs (copied at build).
2. **Fetched at runtime:** when `MF_ZENODO_BASE` is set, `docker/entrypoint.sh`
   runs `download_data.R`, which downloads each Zenodo file to its `data/...`
   location (and extracts the example-hindcast tarball) with md5 verification.
3. **Mounted:** to use your own full `data/` tree (e.g. MCMC outputs for steps
   06–07), bind-mount it: `-v "$(pwd)/data:/opt/microbialForecasts/data"`.

## Phylogenetic analysis

Tests whether the environmental effect sizes are phylogenetically structured
across taxonomic ranks (Fig S4). Everything for the **downstream** result is
committed (`data/clean/phylo_inputs_slim.rds`, `data/summary/predictor_effects.rds`,
`phylo_analysis_results.rds`, the null-SES csvs); only the *from-raw* tree rebuild
needs the Zenodo `phyloseq_16S.rds`, and the forecast-score contribution needs
`scoring_metrics_plsr2.rds`.

```bash
# downstream (committed inputs only):
Rscript analysis/model_analysis/phylogeny/install_phylo_deps.R          # required deps
cd analysis/model_analysis && Rscript phylogeny/phylo_contribution.r    # null model + Blomberg's K
cd ../.. && Rscript analysis/create_figs/fig_phylo_contribution.r       # Fig S4

# from raw (needs Zenodo files; export MF_ZENODO_BASE first, then download_data.R):
Rscript analysis/model_analysis/phylogeny/install_phylo_deps.R --fromraw
cd analysis/model_analysis && Rscript phylogeny/create_bacterial_phylogeny.r
Rscript phylogeny/phylo_contribution_scores.r
```

`set.seed(1)` makes the 999-permutation null deterministic. Recorded package
versions are in `analysis/model_analysis/phylogeny/SESSIONINFO.txt`.

## Pipeline steps and HPC fits

- Downstream pipeline steps run in the container with a mounted `data/` tree, e.g.
  `... microbial-forecasts:latest Rscript analysis/model_analysis/07_tidyHindcasts.r`.
  Steps 06–07 need mounted `data/model_outputs` and `data/hindcasts`.
- The full MCMC fits (step 01, ~100k iterations) are run on HPC. `run_ec2_batch.sh`
  dispatches warm-start `env_cycl` fits as parallel containers on an EC2/amd64 host
  (`N_PARALLEL=8 DATA_DIR=$HOME/data ./run_ec2_batch.sh`).

## Files in this directory

| File | Purpose |
|------|---------|
| `Dockerfile` | Image definition (pinned env + installs) |
| `entrypoint.sh` | Runtime Zenodo fetch (when `MF_ZENODO_BASE` set), then runs the command |
| `install_extra_pkgs.R` | Extra CRAN/GitHub packages for the analysis scripts |
| `smoke_test.R` | Default `CMD`: loads the stack, checks key packages |
| `run_ec2_batch.sh` | Parallel warm-start MCMC fits across taxa (EC2) |

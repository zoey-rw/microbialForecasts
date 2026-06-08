# microbialForecasts — analysis environment (R + NIMBLE + pipeline extras)
#
# Build (from repository root). FROM pins linux/amd64 so Apple Silicon works without remembering flags.
#   docker build -f docker/Dockerfile -t microbial-forecasts:latest .
# (Optional duplicate: docker build --platform=linux/amd64 ...). BuildKit may warn on the pin; it is safe to ignore.
#
# Smoke test (amd64 image on Apple Silicon needs explicit platform):
#   docker run --rm --platform=linux/amd64 microbial-forecasts:latest
#
# Reproduce from the single image (committed inputs are baked in):
#   # phylogenetic-signal analysis + figures — needs no external data:
#   docker run --rm --platform=linux/amd64 -v "$(pwd)/figures:/opt/microbialForecasts/figures" \
#     microbial-forecasts:latest Rscript analysis/create_figs/fig_phylo_contribution.r
#   # downstream analyses that need the large Zenodo inputs: set MF_ZENODO_BASE and the
#   # entrypoint fetches + md5-verifies them (via download_data.R) before the command runs:
#   docker run --rm --platform=linux/amd64 -e MF_ZENODO_BASE="https://zenodo.org/records/<id>/files" \
#     -v "$(pwd)/figures:/opt/microbialForecasts/figures" microbial-forecasts:latest \
#     Rscript analysis/create_figs/fig_phylo_contribution.r
#
# Run a script with data mounted:
#   docker run --rm --platform=linux/amd64 -v "$(pwd)/data:/opt/microbialForecasts/data" microbial-forecasts:latest \
#     Rscript analysis/model_analysis/07_tidyHindcasts.r
#
# 01_fitModels.R LOCAL_TEST (one chain, first task; mount data/ with clean inputs + model_outputs if using warm-start):
#   docker run --rm --platform=linux/amd64 \
#     -v "$(pwd)/data:/opt/microbialForecasts/data" -w /opt/microbialForecasts \
#     -e LOCAL_TEST=true -e NCHAINS=1 \
#     microbial-forecasts:latest \
#     Rscript analysis/model_analysis/01_fitModels.R --task 1
#
# rocker/tidyverse:4.4.2 is amd64-only — platform required for arm64 builders.
# Native arm64 image: replace with rocker/tidyverse:latest and remove --platform below.

FROM --platform=linux/amd64 rocker/tidyverse:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgit2-dev \
    libxml2-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/microbialForecasts

COPY . /opt/microbialForecasts/

RUN mkdir -p /opt/microbialForecasts/logs

# R CMD INSTALL does not fetch CRAN dependencies; install_local(..., dependencies = TRUE) does.
# All installs use the base image's pinned Posit Package Manager snapshot (rocker sets
# a dated CRAN repo per R release) for reproducible package versions; GitHub deps are
# pinned to commits inside the install scripts.
RUN Rscript -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L)); install.packages("remotes", dependencies = FALSE)' \
 && Rscript -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L)); remotes::install_local("/opt/microbialForecasts/microbialForecast", dependencies = TRUE, upgrade = "never")' \
 && Rscript /opt/microbialForecasts/docker/install_extra_pkgs.R \
 && Rscript /opt/microbialForecasts/analysis/model_analysis/phylogeny/install_phylo_deps.R \
 && Rscript -e 'stopifnot(requireNamespace("microbialForecast", quietly = TRUE))'

ENV MICROBIAL_FORECASTS_ROOT=/opt/microbialForecasts
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV VECLIB_MAXIMUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1

# Entrypoint fetches Zenodo-hosted inputs at runtime when MF_ZENODO_BASE is set,
# then runs the given command (defaults to the smoke test).
RUN chmod +x /opt/microbialForecasts/docker/entrypoint.sh
ENTRYPOINT ["/opt/microbialForecasts/docker/entrypoint.sh"]
CMD ["Rscript", "docker/smoke_test.R"]

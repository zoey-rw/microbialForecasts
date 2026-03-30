# microbialForecasts — analysis environment (R + NIMBLE + pipeline extras)
#
# Build (from repository root). FROM pins linux/amd64 so Apple Silicon works without remembering flags.
#   docker build -f docker/Dockerfile -t microbial-forecasts:latest .
# (Optional duplicate: docker build --platform=linux/amd64 ...). BuildKit may warn on the pin; it is safe to ignore.
#
# Smoke test (amd64 image on Apple Silicon needs explicit platform):
#   docker run --rm --platform=linux/amd64 microbial-forecasts:latest
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

# R CMD INSTALL does not fetch CRAN dependencies; install_local(..., dependencies = TRUE) does.
RUN Rscript -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L)); install.packages("remotes", repos = "https://cloud.r-project.org", dependencies = FALSE)' \
 && Rscript -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L)); remotes::install_local("/opt/microbialForecasts/microbialForecast", dependencies = TRUE, upgrade = "never")' \
 && Rscript /opt/microbialForecasts/docker/install_extra_pkgs.R \
 && Rscript -e 'stopifnot(requireNamespace("microbialForecast", quietly = TRUE))'

ENV MICROBIAL_FORECASTS_ROOT=/opt/microbialForecasts
ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV VECLIB_MAXIMUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1

CMD ["Rscript", "docker/smoke_test.R"]

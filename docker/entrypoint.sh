#!/bin/bash
# Container entrypoint for the microbialForecasts reproducibility image.
#
# The image bakes in the code, pinned R environment, and all git-committed input
# data. The large inputs that live on Zenodo (not in git) are fetched here at
# container start when MF_ZENODO_BASE points at a published record, e.g.:
#
#   docker run --rm -e MF_ZENODO_BASE="https://zenodo.org/records/<id>/files" \
#     microbial-forecasts:latest Rscript analysis/create_figs/fig_phylo_contribution.r
#
# Without MF_ZENODO_BASE, the committed inputs already in the image are used as-is.
set -e
cd "${MICROBIAL_FORECASTS_ROOT:-/opt/microbialForecasts}"

if [ -n "${MF_ZENODO_BASE:-}" ] && [[ "${MF_ZENODO_BASE}" != *REPLACE_WITH_RECORD_ID* ]]; then
  echo "[entrypoint] MF_ZENODO_BASE set -> fetching Zenodo inputs via download_data.R"
  Rscript download_data.R || echo "[entrypoint] WARNING: download_data.R failed; continuing with baked-in data"
fi

exec "$@"

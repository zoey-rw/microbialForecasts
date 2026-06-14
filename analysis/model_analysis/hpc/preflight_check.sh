#!/bin/bash
# Pre-flight check: verify ALL dependencies exist on SCC before submitting pipeline.
# Run from SCC: bash preflight_check.sh
#
# Checks:
#   1. R packages needed by each step
#   2. Data files needed by each step
#   3. Outputs from previous steps
#   4. Disk space

set -euo pipefail

PROJECT=/projectnb/dietzelab/zrwerbin/microbialForecasts
LIB="$PROJECT/R_library"
FAIL=0

red()   { echo -e "\033[31m$1\033[0m"; }
green() { echo -e "\033[32m$1\033[0m"; }
check_file() {
    if [ -f "$1" ]; then
        green "  OK: $(basename $1)"
    else
        red "  MISSING: $1"
        FAIL=$((FAIL + 1))
    fi
}
check_pkg() {
    if Rscript -e ".libPaths(c('$LIB', .libPaths())); if(!requireNamespace('$1',quietly=TRUE)) quit(status=1)" 2>/dev/null; then
        green "  OK: $1"
    else
        red "  MISSING: $1"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== Pre-flight check for pipeline steps 03-10 ==="
echo "Project: $PROJECT"
echo ""

# R packages
echo "--- R packages ---"
for pkg in microbialForecast tidyverse here data.table coda nimble \
    MuMIn pls spectratrait plotrix reshape2 dplyr tidyr stringr \
    doParallel foreach parallel Metrics scoringRules ggpubr \
    lubridate egg nanoparquet DBI duckdb arrow pacman; do
    check_pkg "$pkg"
done

echo ""
echo "--- Data files (input) ---"
check_file "$PROJECT/data/clean/groupAbundances_16S_2023.rds"
check_file "$PROJECT/data/clean/groupAbundances_ITS_2023.rds"
check_file "$PROJECT/data/clean/all_predictor_data.rds"
check_file "$PROJECT/data/clean/dominantHorizonsSite.rds"
check_file "$PROJECT/data/clean/site_effect_predictors.rds"
check_file "$PROJECT/source.R"
check_file "$PROJECT/analysis/model_analysis/logging.R"

echo ""
echo "--- Step 03-04 outputs ---"
check_file "$PROJECT/data/summary/logit_beta_fixed_priors_summaries.rds"
check_file "$PROJECT/data/summary/predictor_effects.rds"
check_file "$PROJECT/data/summary/seasonal_amplitude.rds"
check_file "$PROJECT/data/summary/converged_taxa_list.rds"
check_file "$PROJECT/data/summary/weak_converged_taxa_list.rds"

echo ""
echo "--- Step 05 outputs ---"
check_file "$PROJECT/data/summary/site_effects.rds"
check_file "$PROJECT/data/summary/site_effects_unobserved_env_cycl.rds"
check_file "$PROJECT/data/summary/site_effects_unobserved_env_cov.rds"
check_file "$PROJECT/data/summary/site_effects_unobserved_cycl_only.rds"

echo ""
echo "--- Step 07 outputs (per-model hindcasts, needed by 08-10) ---"
check_file "$PROJECT/data/summary/parquet/hindcasts_env_cycl.parquet"
check_file "$PROJECT/data/summary/parquet/hindcasts_cycl_only.parquet"
check_file "$PROJECT/data/summary/parquet/hindcasts_env_cov.parquet"

echo ""
echo "--- Step 09 dependencies ---"
check_file "$PROJECT/data/summary/plot_estimates.rds"

echo ""
echo "--- HPC task maps ---"
check_file "$PROJECT/data/summary/hpc_model_taskmap.tsv"
check_file "$PROJECT/data/summary/hpc_taxon_taskmap.txt"

echo ""
echo "--- Model sample files ---"
for mtype in env_cycl cycl_only env_cov; do
    n=$(ls "$PROJECT/data/model_outputs/cloglog_beta_driver_uncertainty/$mtype/samples_${mtype}_"*.rds 2>/dev/null | wc -l)
    echo "  $mtype: $n sample files"
done

echo ""
echo "--- Disk space ---"
df -h /projectnb/ | tail -1

echo ""
if [ $FAIL -gt 0 ]; then
    red "=== $FAIL issues found. Fix before submitting. ==="
else
    green "=== All checks passed. Ready to submit. ==="
fi

#!/bin/bash
# Sync, combine, summarize, and clean up env_cycl chain files one taxon at a time.
# This avoids needing 170GB of local disk by processing each taxon individually.
#
# Usage: bash sync_combine_env_cycl.sh [taxon_filter]
# Example: bash sync_combine_env_cycl.sh rhizobiales
#          bash sync_combine_env_cycl.sh   # process all 169 taxa with Nov 11 chains

set -e

REMOTE_USER="zrwerbin"
REMOTE_HOST="scc2.bu.edu"
REMOTE_BASE="/projectnb/dietzelab/zrwerbin/microbialForecasts/data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LOCAL_BASE="$PROJECT_ROOT/data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"
LOGDIR="$PROJECT_ROOT/logs/sync_combine"
mkdir -p "$LOGDIR"

TAXON_FILTER="${1:-}"

# Get list of taxa to process
if [ -n "$TAXON_FILTER" ]; then
    TAXA=("$TAXON_FILTER")
else
    # All taxa with Nov 11 chains (from /tmp/remote_dates.txt)
    if [ -f /tmp/remote_dates.txt ]; then
        TAXA=($(awk '$4>0 {print $1}' /tmp/remote_dates.txt))
    else
        echo "ERROR: /tmp/remote_dates.txt not found. Run the remote comparison first."
        exit 1
    fi
fi

echo "Processing ${#TAXA[@]} taxa"
echo "================================================="

succeeded=0
failed=0

for taxon in "${TAXA[@]}"; do
    echo ""
    echo "--- $taxon ---"

    local_dir="$LOCAL_BASE/$taxon"
    mkdir -p "$local_dir"

    # Step 1: Rsync chain files only (not checkpoints or progress files)
    echo "  [1/4] Syncing chain files from SCC..."
    rsync -avP \
        --include="samples_*chain*.rds" \
        --exclude="*" \
        "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE}/${taxon}/" \
        "$local_dir/" \
        > "$LOGDIR/${taxon}_sync.log" 2>&1

    nchains=$(find "$local_dir" -maxdepth 1 -name "samples_*chain*.rds" | wc -l | tr -d ' ')
    echo "  Synced: $nchains chain files"

    if [ "$nchains" -lt 2 ]; then
        echo "  SKIP: fewer than 2 chains"
        failed=$((failed + 1))
        continue
    fi

    # Step 2+3: Combine chains, summarize, and clean up
    combined_file="$LOCAL_BASE/samples_env_cycl_${taxon}_20130601_20180101_with_legacy_covariate_beta_regression.rds"
    echo "  [2/3] Combining chains + summarizing..."
    Rscript combine_single_taxon.R \
        "$local_dir" "$combined_file" --cleanup \
        > "$LOGDIR/${taxon}_combine.log" 2>&1
    rc=$?

    if [ $rc -ne 0 ]; then
        echo "  FAILED"
        tail -5 "$LOGDIR/${taxon}_combine.log"
        failed=$((failed + 1))
        continue
    fi

    # Show key results from log
    grep -E "^(Selected|Using|Max Rhat|Converged|Saved|Done)" \
        "$LOGDIR/${taxon}_combine.log" 2>/dev/null | sed 's/^/  /'

    succeeded=$((succeeded + 1))
done

echo ""
echo "================================================="
echo "Done: $succeeded succeeded, $failed failed of ${#TAXA[@]} total"

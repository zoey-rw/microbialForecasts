#!/bin/bash
# Sync, combine, summarize, and clean up env_cycl chain files one taxon at a time.
# This avoids needing 170GB of local disk by processing each taxon individually.
#
# Usage: bash sync_combine_env_cycl.sh [taxon_filter]
# Example: bash sync_combine_env_cycl.sh rhizobiales
#          bash sync_combine_env_cycl.sh   # process all 169 taxa with Nov 11 chains

set -e
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REMOTE_MOUNT="$HOME/remote_microbialForecasts/data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"
LOCAL_BASE="$PROJECT_ROOT/data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"
LOGDIR="$PROJECT_ROOT/logs/sync_combine"
mkdir -p "$LOGDIR"

# Verify remote mount is accessible
if [ ! -d "$REMOTE_MOUNT" ]; then
    echo "ERROR: Remote mount not found at $REMOTE_MOUNT"
    echo "Mount with: sshfs zrwerbin@scc2.bu.edu:/projectnb/dietzelab/zrwerbin/microbialForecasts ~/remote_microbialForecasts -o reconnect"
    exit 1
fi

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

    remote_dir="$REMOTE_MOUNT/$taxon"
    local_dir="$LOCAL_BASE/$taxon"
    mkdir -p "$local_dir"

    # Skip if already completed in this run
    if [ -f "$LOGDIR/${taxon}_combine.log" ] && grep -q "^Done" "$LOGDIR/${taxon}_combine.log" 2>/dev/null; then
        echo "  SKIP: already completed"
        rhat=$(grep "^Max Rhat" "$LOGDIR/${taxon}_combine.log" | head -1 | grep -oE "[0-9]+\.[0-9]+")
        echo "  (Rhat=$rhat)"
        succeeded=$((succeeded + 1))
        continue
    fi

    if [ ! -d "$remote_dir" ]; then
        echo "  SKIP: no remote directory"
        failed=$((failed + 1))
        continue
    fi

    # Rsync chain files from remote mount to local (fast: sequential read)
    echo "  [1/2] Copying chain files..."
    rsync -a --include="samples_*chain*.rds" --exclude="*" \
        "$remote_dir/" "$local_dir/" \
        > "$LOGDIR/${taxon}_sync.log" 2>&1

    nchains=$(find "$local_dir" -maxdepth 1 -name "samples_*chain*.rds" 2>/dev/null | wc -l | tr -d ' ')
    echo "  Copied: $nchains chain files"

    if [ "$nchains" -lt 2 ]; then
        echo "  SKIP: fewer than 2 chains"
        failed=$((failed + 1))
        continue
    fi

    # Combine chains, summarize, clean up
    combined_file="$LOCAL_BASE/samples_env_cycl_${taxon}_20130601_20180101_with_legacy_covariate_beta_regression.rds"
    echo "  [2/2] Combining + summarizing..."
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
    grep -E "^(Selected|Using|Keeping|Truncating|Max Rhat|Converged|Saved|Done)" \
        "$LOGDIR/${taxon}_combine.log" 2>/dev/null | sed 's/^/  /'

    succeeded=$((succeeded + 1))
done

echo ""
echo "================================================="
echo "Done: $succeeded succeeded, $failed failed of ${#TAXA[@]} total"

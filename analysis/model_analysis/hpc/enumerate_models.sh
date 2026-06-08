#!/bin/bash
# Scan model_outputs for combined samples files, produce task maps for HPC array jobs.
#
# Outputs:
#   $PROJECT/data/summary/hpc_model_taskmap.tsv  -- TASK_ID MODEL_NAME TAXON MODEL_ID
#   $PROJECT/data/summary/hpc_taxon_taskmap.txt   -- unique taxon names (one per line)
#
# Usage: bash enumerate_models.sh [PROJECT_ROOT]

set -euo pipefail

# Allow PROJECT to be passed as arg or env var; default to SCC path
PROJECT="${1:-${PROJECT:-/projectnb/dietzelab/zrwerbin/microbialForecasts}}"
MODEL_BASE="$PROJECT/data/model_outputs/cloglog_beta_driver_uncertainty"

TASKMAP="$PROJECT/data/summary/hpc_model_taskmap.tsv"
TAXMAP="$PROJECT/data/summary/hpc_taxon_taskmap.txt"

# Ensure output directory exists
mkdir -p "$(dirname "$TASKMAP")"

# Collect all combined samples files across model types
task_id=0
> "$TASKMAP"  # truncate

for model_type in env_cycl env_cov cycl_only; do
    model_dir="$MODEL_BASE/$model_type"
    [ -d "$model_dir" ] || continue

    for f in "$model_dir"/samples_${model_type}_*_beta_regression.rds; do
        [ -f "$f" ] || continue

        base=$(basename "$f")

        # Skip individual chain files and reconstructed checkpoints
        echo "$base" | grep -qE '_chain[0-9]' && continue
        echo "$base" | grep -qi 'reconstructed_from_checkpoints' && continue

        # Parse model_id: strip samples_ prefix and _beta_regression.rds suffix
        model_id="${base#samples_}"
        model_id="${model_id%_beta_regression.rds}"

        # Parse taxon: strip model_name prefix, then strip date suffix
        # e.g. env_cycl_acidibacter_20130601_20180101_with_legacy_covariate -> acidibacter
        taxon="${model_id#${model_type}_}"
        taxon=$(echo "$taxon" | sed 's/_[0-9]\{8\}_.*$//')

        task_id=$((task_id + 1))
        printf '%d\t%s\t%s\t%s\n' "$task_id" "$model_type" "$taxon" "$model_id" >> "$TASKMAP"
    done
done

echo "Wrote $task_id models to $TASKMAP"

# Generate unique taxon list for newsites array job
cut -f3 "$TASKMAP" | sort -u > "$TAXMAP"
n_taxa=$(wc -l < "$TAXMAP")
echo "Wrote $n_taxa unique taxa to $TAXMAP"

echo "N_MODELS=$task_id"
echo "N_TAXA=$n_taxa"

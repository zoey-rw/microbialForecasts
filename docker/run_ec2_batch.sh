#!/bin/bash
set -euo pipefail

# Run warm-start MCMC fits for all env_cycl taxa that have warmstart_inits.rds files.
# Each taxon gets NCHAINS independent chains, dispatched as parallel Docker containers.
#
# Usage:
#   N_PARALLEL=8 DATA_DIR=$HOME/data ./run_ec2_batch.sh
#   N_PARALLEL=4 NCHAINS=2 ./run_ec2_batch.sh          # fewer chains for testing
#   DRY_RUN=true ./run_ec2_batch.sh                     # print commands without running

N_PARALLEL=${N_PARALLEL:-8}
NCHAINS=${NCHAINS:-4}
DATA_DIR=${DATA_DIR:-$HOME/data}
IMAGE=${IMAGE:-microbial-forecasts:latest}
DRY_RUN=${DRY_RUN:-false}
MODEL_NAME=${MODEL_NAME:-env_cycl}
WARMSTART_DIR="$DATA_DIR/model_outputs/cloglog_beta_driver_uncertainty/$MODEL_NAME"

if [ ! -d "$WARMSTART_DIR" ]; then
    echo "ERROR: Warm-start directory not found: $WARMSTART_DIR"
    exit 1
fi

TAXA=()
while IFS= read -r ws_file; do
    taxon=$(basename "$(dirname "$ws_file")")
    TAXA+=("$taxon")
done < <(find "$WARMSTART_DIR" -name "warmstart_inits.rds" -type f | sort)

N_TAXA=${#TAXA[@]}
if [ "$N_TAXA" -eq 0 ]; then
    echo "ERROR: No warmstart_inits.rds files found under $WARMSTART_DIR"
    exit 1
fi

TOTAL_TASKS=$((N_TAXA * NCHAINS))
echo "========================================"
echo "Warm-start MCMC batch runner"
echo "  Model:      $MODEL_NAME"
echo "  Taxa:       $N_TAXA (with warm-start files)"
echo "  Chains:     $NCHAINS per taxon"
echo "  Total tasks: $TOTAL_TASKS"
echo "  Parallel:   $N_PARALLEL workers"
echo "  Data dir:   $DATA_DIR"
echo "  Image:      $IMAGE"
echo "========================================"

TASK_LIST=$(mktemp)
trap 'rm -f "$TASK_LIST"' EXIT

for taxon in "${TAXA[@]}"; do
    for chain in $(seq 1 "$NCHAINS"); do
        echo "$taxon $chain"
    done
done > "$TASK_LIST"

LOG_DIR="$DATA_DIR/logs/batch_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"
echo "Logs: $LOG_DIR"

run_task() {
    local taxon="$1"
    local chain="$2"
    local log_file="$LOG_DIR/${taxon}_chain${chain}.log"

    echo "[$(date +%H:%M:%S)] START $taxon chain $chain"

    docker run --rm \
        -v "$DATA_DIR:/opt/microbialForecasts/data" \
        -e NCHAINS="$NCHAINS" \
        "$IMAGE" \
        Rscript analysis/model_analysis/01_fitModels.R \
            --model-name "$MODEL_NAME" \
            --species "$taxon" \
            --task "$chain" \
        > "$log_file" 2>&1

    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] DONE  $taxon chain $chain"
    else
        echo "[$(date +%H:%M:%S)] FAIL  $taxon chain $chain (exit $exit_code) — see $log_file"
    fi
    return $exit_code
}
export -f run_task
export DATA_DIR IMAGE NCHAINS LOG_DIR MODEL_NAME

if [ "$DRY_RUN" = "true" ]; then
    echo ""
    echo "DRY RUN — commands that would be executed:"
    while read -r taxon chain; do
        echo "  docker run --rm -v $DATA_DIR:/opt/microbialForecasts/data -e NCHAINS=$NCHAINS $IMAGE Rscript analysis/model_analysis/01_fitModels.R --model-name $MODEL_NAME --species $taxon --task $chain"
    done < "$TASK_LIST"
    echo ""
    echo "Total: $TOTAL_TASKS tasks"
    exit 0
fi

echo ""
echo "Starting at $(date)"
FAILED=0
while read -r taxon chain; do
    echo "$taxon $chain"
done < "$TASK_LIST" | xargs -P "$N_PARALLEL" -n 2 bash -c 'run_task "$@" || exit 0' _

echo ""
echo "Batch complete at $(date)"
echo "Logs saved to: $LOG_DIR"

N_DONE=$(find "$LOG_DIR" -name "*.log" -exec grep -l "SUCCESS" {} \; 2>/dev/null | wc -l | tr -d ' ')
N_FAIL=$(find "$LOG_DIR" -name "*.log" -exec grep -l "FAIL\|Error\|error" {} \; 2>/dev/null | wc -l | tr -d ' ')
echo "Results: $N_DONE succeeded, $N_FAIL with errors (of $TOTAL_TASKS total)"

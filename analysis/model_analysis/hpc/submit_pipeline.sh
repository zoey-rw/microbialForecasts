#!/bin/bash -l
# Master submission script for the full downstream pipeline (steps 03-10).
# Enumerates models, then submits all jobs with SGE -hold_jid dependencies.
#
# Usage:
#   bash submit_pipeline.sh                  # submit steps 03-10
#   bash submit_pipeline.sh --from=06        # start from step 06
#   bash submit_pipeline.sh --to=07          # stop after step 07
#   bash submit_pipeline.sh --from=06 --to=07
#   bash submit_pipeline.sh --dry-run        # print commands without submitting
#
# Dependency graph:
#   03 -> 04 -> 05 -> 06_obs + 06_new (parallel) -> 07 -> 08 + 09 + 10 (parallel)

set -euo pipefail

PROJECT=/projectnb/dietzelab/zrwerbin/microbialForecasts
HPC_DIR="$PROJECT/analysis/model_analysis/hpc"

FROM=03
TO=10
DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --from=*) FROM="${arg#*=}" ;;
        --to=*)   TO="${arg#*=}" ;;
        --dry-run) DRY_RUN=true ;;
        *)        echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# Remove leading zeros for arithmetic
FROM=$((10#$FROM))
TO=$((10#$TO))

echo "=== Pipeline submission ==="
echo "Steps: $(printf '%02d' $FROM) through $(printf '%02d' $TO)"
echo "Dry run: $DRY_RUN"
echo ""

# Create log directories
for d in 03_summarize 04_tidy 05_predict 06_obs 06_new 07_tidy 08_scoring 09_pheno 10_horizon; do
    mkdir -p "$PROJECT/logs/hpc/$d"
done

# --- Enumerate models (always runs) ---
echo "=== Enumerating models ==="
bash "$HPC_DIR/enumerate_models.sh" "$PROJECT"
echo ""

N_MODELS=$(wc -l < "$PROJECT/data/summary/hpc_model_taskmap.tsv" | tr -d ' ')
N_TAXA=$(wc -l < "$PROJECT/data/summary/hpc_taxon_taskmap.txt" | tr -d ' ')

# --- Helper: submit and extract job ID ---
submit() {
    local desc="$1"; shift
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] $*"
        echo "DRY_$$_$RANDOM"
        return
    fi
    local output
    output=$("$@" 2>&1) || { echo "SUBMIT FAILED: $output"; exit 1; }
    # Extract numeric job ID from qsub output
    local jid
    jid=$(echo "$output" | grep -oE '[0-9]+' | head -1)
    echo "$desc -> Job $jid"
    echo "$jid"
}

# Track job IDs for dependencies
JID_03="" JID_04="" JID_05="" JID_06_OBS="" JID_06_NEW="" JID_07=""

# --- Step 03 ---
if [ $FROM -le 3 ] && [ $TO -ge 3 ]; then
    echo "--- Step 03: Summarize model outputs ---"
    JID_03=$(submit "Step 03" qsub "$HPC_DIR/submit_03_summarize.qsub")
fi

# --- Step 04 ---
if [ $FROM -le 4 ] && [ $TO -ge 4 ]; then
    HOLD=""
    [ -n "$JID_03" ] && HOLD="-hold_jid $JID_03"
    echo "--- Step 04: Tidy effect sizes ---"
    JID_04=$(submit "Step 04" qsub $HOLD "$HPC_DIR/submit_04_tidy_effects.qsub")
fi

# --- Step 05 ---
if [ $FROM -le 5 ] && [ $TO -ge 5 ]; then
    HOLD=""
    [ -n "$JID_04" ] && HOLD="-hold_jid $JID_04"
    echo "--- Step 05: Predict site effects (8 cores) ---"
    JID_05=$(submit "Step 05" qsub $HOLD "$HPC_DIR/submit_05_predict_site.qsub")
fi

# --- Step 06 (observed + newsites in parallel, both depend on step 05) ---
if [ $FROM -le 6 ] && [ $TO -ge 6 ]; then
    HOLD=""
    [ -n "$JID_05" ] && HOLD="-hold_jid $JID_05"

    echo "--- Step 06_obs: Hindcast observed sites ($N_MODELS array tasks) ---"
    JID_06_OBS=$(submit "Step 06_obs" qsub -t 1-"$N_MODELS" $HOLD "$HPC_DIR/submit_06_obs.qsub")

    echo "--- Step 06_new: Hindcast new sites ($N_TAXA array tasks) ---"
    JID_06_NEW=$(submit "Step 06_new" qsub -t 1-"$N_TAXA" $HOLD "$HPC_DIR/submit_06_new.qsub")
fi

# --- Step 07 (depends on both 06 jobs) ---
if [ $FROM -le 7 ] && [ $TO -ge 7 ]; then
    HOLD=""
    HOLD_PARTS=""
    [ -n "$JID_06_OBS" ] && HOLD_PARTS="$JID_06_OBS"
    [ -n "$JID_06_NEW" ] && HOLD_PARTS="${HOLD_PARTS:+$HOLD_PARTS,}$JID_06_NEW"
    [ -n "$HOLD_PARTS" ] && HOLD="-hold_jid $HOLD_PARTS"
    echo "--- Step 07: Tidy hindcasts ---"
    JID_07=$(submit "Step 07" qsub $HOLD "$HPC_DIR/submit_07_tidy_hindcasts.qsub")
fi

# --- Steps 08, 09, 10 (all in parallel after step 07) ---
if [ $FROM -le 8 ] && [ $TO -ge 8 ]; then
    HOLD=""
    [ -n "$JID_07" ] && HOLD="-hold_jid $JID_07"
    echo "--- Step 08: Scoring metrics ---"
    submit "Step 08" qsub $HOLD "$HPC_DIR/submit_08_scoring.qsub" > /dev/null
fi

if [ $FROM -le 9 ] && [ $TO -ge 9 ]; then
    HOLD=""
    [ -n "$JID_07" ] && HOLD="-hold_jid $JID_07"
    echo "--- Step 09: Phenophases ---"
    submit "Step 09" qsub $HOLD "$HPC_DIR/submit_09_phenophase.qsub" > /dev/null
fi

if [ $FROM -le 10 ] && [ $TO -ge 10 ]; then
    HOLD=""
    [ -n "$JID_07" ] && HOLD="-hold_jid $JID_07"
    echo "--- Step 10: Forecast horizon ---"
    submit "Step 10" qsub $HOLD "$HPC_DIR/submit_10_fcast_horizon.qsub" > /dev/null
fi

echo ""
echo "=== All jobs submitted ==="
echo "Monitor with: qstat -u \$USER"
echo "Logs in: $PROJECT/logs/hpc/"

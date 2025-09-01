#!/bin/bash

# HPC VERSION: Run unconverged models with 28 cores for maximum efficiency
echo "=== LOCAL MODE: UNCONVERGED MODELS WITH 3 CORES ==="
echo "Target: Sequential processing with 3 cores"
echo "Chains per model: 3"
echo "Models per batch: 1 (3 chains × 1 model = 3 cores)"
echo ""

# Configuration
MAX_CORES=3
CHAINS_PER_MODEL=3
MODELS_PER_BATCH=1
LOG_DIR="logs_hpc"

# Create logs directory
mkdir -p "$LOG_DIR"

# Function to run a single model with 4 chains
run_model_with_chains() {
    local model_name=$1
    local model_idx=$2
    local log_dir="$LOG_DIR/${model_name}"
    
    mkdir -p "$log_dir"
    
    echo "=== Starting model: $model_name (index: $model_idx) ==="
    echo "Log directory: $log_dir"
    
    # Function to run a single chain
    run_chain() {
        local chain=$1
        local log_file="$log_dir/chain${chain}.log"
        
        echo "  Starting Chain $chain..."
        echo "  Log file: $log_file"
        
        # Set environment variable for R script (model index, specific chain number)
        export R_ARGS="$model_idx $chain"
        
        # Run the restart script
        R --slave -e "source('analysis/model_analysis/01_restartFitModels.R')" > "$log_file" 2>&1 &
        
        echo "  Chain $chain PID: $!"
    }
    
    # Start all 3 chains simultaneously
    echo "  Starting all 3 chains simultaneously..."
    run_chain "1"
    run_chain "2"
    run_chain "3"
    
    echo "  All chains started for $model_name"
}

# Function to check if a model has converged
check_model_convergence() {
    local model_name=$1
    local output_dir="data/model_outputs/logit_beta_fixed_priors/env_cycl"
    
    # Extract species name from model name (remove env_cycl_ prefix and _date_with_legacy_covariate suffix)
    local species_name=$(echo "$model_name" | sed 's/env_cycl_//g' | sed 's/_20130601_20180101_with_legacy_covariate//g')
    
    # Look for checkpoint files in the species subdirectory
    local species_dir="$output_dir/$species_name"
    local chain1_file="$species_dir/checkpoint_${model_name}_chain1_loop*.rds"
    local chain2_file="$species_dir/checkpoint_${model_name}_chain2_loop*.rds"
    local chain3_file="$species_dir/checkpoint_${model_name}_chain3_loop*.rds"
    
    # Count how many chains have checkpoint files
    local chain_count=0
    ls $chain1_file >/dev/null 2>&1 && ((chain_count++))
    ls $chain2_file >/dev/null 2>&1 && ((chain_count++))
    ls $chain3_file >/dev/null 2>&1 && ((chain_count++))
    
    if [ $chain_count -eq 3 ]; then
        echo "  ✓ All 3 chains completed for $model_name"
        return 0
    elif [ $chain_count -gt 0 ]; then
        echo "  ⚠ $chain_count/3 chains completed for $model_name (looking in $species_dir)"
        return 1
    else
        echo "  ❌ No chains completed for $model_name (looking in $species_dir)"
        return 1
    fi
}

# Function to wait for batch completion
wait_for_batch() {
    local batch_models=("$@")
    local total_chains=$(( ${#batch_models[@]} * CHAINS_PER_MODEL ))
    
    echo "Waiting for batch completion..."
    echo "  Models in batch: ${#batch_models[@]}"
    echo "  Total chains running: $total_chains"
    
    # Wait for all processes to complete
    wait
    
    echo "  ✓ Batch completed!"
    echo ""
}

# Main execution
echo "Loading unconverged taxa list..."
unconverged_models=$(R --slave -e "cat(paste(readRDS('data/summary/unconverged_taxa_list.rds'), collapse='\n'))")

if [ -z "$unconverged_models" ]; then
    echo "Error: No unconverged models found!"
    exit 1
fi

# Convert to array
IFS=$'\n' read -d '' -r -a model_array <<< "$unconverged_models"

echo "Found ${#model_array[@]} unconverged models:"
for i in "${!model_array[@]}"; do
    echo "  $((i+1)). ${model_array[i]}"
done

echo ""
echo "LOCAL MODE: Processing all models sequentially (1 at a time)"
echo "Each model will run with $CHAINS_PER_MODEL chains simultaneously"
echo "Total cores utilized: $MAX_CORES"
echo ""

# Process models in batches
total_models=${#model_array[@]}
total_batches=$(( (total_models + MODELS_PER_BATCH - 1) / MODELS_PER_BATCH ))

echo "Processing $total_models models in $total_batches batches..."
echo ""

for ((batch=0; batch<total_batches; batch++)); do
    start_idx=$((batch * MODELS_PER_BATCH))
    end_idx=$((start_idx + MODELS_PER_BATCH - 1))
    
    # Ensure we don't exceed array bounds
    if [ $end_idx -ge $total_models ]; then
        end_idx=$((total_models - 1))
    fi
    
    echo "=== BATCH $((batch + 1)) of $total_batches ==="
    echo "Processing models $((start_idx + 1)) to $((end_idx + 1))"
    echo "Models in this batch:"
    
    # Start all models in this batch
    batch_models=()
    for ((i=start_idx; i<=end_idx; i++)); do
        model_name="${model_array[i]}"
        model_idx=$((i+1))
        batch_models+=("$model_name")
        
        echo "  $model_idx. $model_name"
        run_model_with_chains "$model_name" "$model_idx"
    done
    
    echo ""
    echo "All models in batch started. Waiting for completion..."
    
    # Wait for this batch to complete
    wait_for_batch "${batch_models[@]}"
    
    # Check completion status for this batch
    echo "Checking batch completion status..."
    for ((i=start_idx; i<=end_idx; i++)); do
        model_name="${model_array[i]}"
        if check_model_convergence "$model_name"; then
            echo "✓ $model_name: COMPLETED"
        else
            echo "❌ $model_name: FAILED"
        fi
    done
    
    echo ""
    echo "Batch $((batch + 1)) completed!"
    echo ""
done

echo "=== ALL MODELS PROCESSED ==="
echo "Check $LOG_DIR/ directory for detailed results"
echo ""
echo "Final status check:"
echo "Looking for completed models..."

# Final status check for all models
completed_count=0
failed_count=0
for i in "${!model_array[@]}"; do
    model_name="${model_array[i]}"
    if check_model_convergence "$model_name"; then
        echo "✓ $model_name: COMPLETED"
        ((completed_count++))
    else
        echo "❌ $model_name: FAILED"
        ((failed_count++))
    fi
done

echo ""
echo "=== FINAL SUMMARY ==="
echo "Total models: $total_models"
echo "Completed: $completed_count"
echo "Failed: $failed_count"
echo "Success rate: $(( (completed_count * 100) / total_models ))%"
echo ""
echo "HPC script completed!"


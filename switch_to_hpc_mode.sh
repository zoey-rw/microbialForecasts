#!/bin/bash

# Script to switch between testing and HPC modes
echo "=== SWITCHING BETWEEN TESTING AND HPC MODES ==="

if [ "$1" = "hpc" ]; then
    echo "Switching to HPC mode (4 chains)..."
    
    # Update restart script to use 4 chains
    sed -i '' 's/nchains <- 3  # Number of chains per model (for testing - can be changed to 4 for production/nchains <- 4  # Number of chains per model (for HPC mode)/' analysis/model_analysis/01_restartFitModels.R
    
    echo "✓ Switched to HPC mode (4 chains)"
    echo "Run with: ./run_unconverged_models_hpc.sh"
    
elif [ "$1" = "test" ]; then
    echo "Switching to testing mode (3 chains)..."
    
    # Update restart script to use 3 chains
    sed -i '' 's/nchains <- 4  # Number of chains per model (for HPC mode)/nchains <- 3  # Number of chains per model (for testing - can be changed to 4 for production/' analysis/model_analysis/01_restartFitModels.R
    
    echo "✓ Switched to testing mode (3 chains)"
    echo "Run with: ./test_unconverged_models.sh"
    
else
    echo "Usage: $0 [test|hpc]"
    echo ""
    echo "  test  - Switch to testing mode (3 chains, 6 cores, 30min timeout)"
    echo "  hpc   - Switch to HPC mode (4 chains, 28 cores, no timeout)"
    echo ""
    echo "Current mode:"
    grep "nchains <-" analysis/model_analysis/01_restartFitModels.R
fi

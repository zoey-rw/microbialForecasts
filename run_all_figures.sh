#!/bin/bash
cd /Users/zoeywerbin/Documents/microbialForecasts

# Discover figure scripts (create_figs only; excludes subdirectories)
# Include both .r and .R so all scripts run regardless of filename case
scripts=( $(find analysis/create_figs -maxdepth 1 \( -name "*.r" -o -name "*.R" \) | sort) )

# Create log file
log_file="logs/figure_generation_$(date +%Y%m%d_%H%M%S).log"
mkdir -p logs
echo "Figure generation log: $log_file" | tee "$log_file"

success_count=0
failed_count=0
skipped_count=0

# Run each script with timeout
for script in "${scripts[@]}"; do
  if [ -f "$script" ]; then
    echo "==========================================" | tee -a "$log_file"
    echo "Running: $script" | tee -a "$log_file"
    echo "==========================================" | tee -a "$log_file"

    # Use timeout to prevent freezing (10 minutes per script)
    timeout 600 Rscript "$script" >> "$log_file" 2>&1
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
      echo "SUCCESS: $script" | tee -a "$log_file"
      ((success_count++))
    elif [ $exit_code -eq 124 ]; then
      echo "TIMEOUT: $script (exceeded 10 minutes)" | tee -a "$log_file"
      ((failed_count++))
    else
      echo "FAILED: $script (exit code: $exit_code)" | tee -a "$log_file"
      ((failed_count++))
    fi
    echo "" | tee -a "$log_file"
  else
    echo "NOT FOUND: $script" | tee -a "$log_file"
    ((skipped_count++))
  fi
done

echo "==========================================" | tee -a "$log_file"
echo "Summary:" | tee -a "$log_file"
echo "  Successful: $success_count" | tee -a "$log_file"
echo "  Failed: $failed_count" | tee -a "$log_file"
echo "  Skipped: $skipped_count" | tee -a "$log_file"
echo "==========================================" | tee -a "$log_file"

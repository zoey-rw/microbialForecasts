# Model Restart Functionality

This directory contains functionality to restart MCMC models using initial values from previous runs, designed to "rescue" models that were fit with loose priors.

## 🎯 Problem Solved

Many MCMC models in the microbial forecasting project were originally fit with loose priors, leading to:
- Extreme parameter values (e.g., rho = -20, site effects = ±40)
- Poor convergence (high R-hat values, low effective sample sizes)
- Numerical instability and model failures

This restart functionality automatically:
1. **Extracts final parameter values** from existing MCMC chains
2. **Detects extreme values** using parameter-specific bounds
3. **Creates fallback values** for extreme/invalid parameters
4. **Restarts models** with improved initial values

## 📁 Files Overview

### Core Functionality
- **`model_restart_functions.R`** - Core restart functions and utilities
- **`01_fitModels_with_restart.R`** - Enhanced model fitting script with restart capability
- **`batch_restart_models.R`** - Batch processing for multiple models

### Testing & Demonstration
- **`test_restart_functionality.R`** - Test script for restart functions
- **`demo_restart_functionality.R`** - Complete demonstration of all features

## 🚀 Quick Start

### 1. Test the Functionality
```bash
# Test restart functions with existing chain files
Rscript test_restart_functionality.R
```

### 2. Run Complete Demonstration
```bash
# See all restart features in action
Rscript demo_restart_functionality.R
```

### 3. Batch Process Multiple Models
```bash
# Process 50 models with 4 cores in parallel
Rscript batch_restart_models.R 50 4
```

## 📊 Key Functions

### `extract_final_values_from_chains()`
Extracts final parameter values from MCMC chain files with diagnostics.

**Parameters:**
- `chain_files`: Vector of RDS files containing MCMC chains
- `min_ess`: Minimum effective sample size required (default: 100)
- `max_rhat`: Maximum R-hat value allowed (default: 1.1)
- `burnin_proportion`: Proportion of samples to discard (default: 0.5)

**Returns:**
- Final parameter values after burnin
- Convergence diagnostics (ESS, R-hat)
- Extreme value flags

### `create_restart_inits()`
Creates restart initial values with fallback strategies for extreme parameters.

**Parameters:**
- `extraction_result`: Result from `extract_final_values_from_chains()`
- `use_fallback_for_extreme`: Whether to use fallback values (default: TRUE)
- `fallback_strategy`: "conservative", "random", or "zero"

**Returns:**
- Initial values suitable for restarting MCMC
- Flags indicating which parameters used fallback values

### `restart_mcmc_model()`
Complete workflow to set up model restart.

**Parameters:**
- `model_name`, `species`, `min_date`, `max_date`: Model identifiers
- `use_legacy_covariate`: Whether model uses legacy covariate
- `chain_files`: Optional specific chain files to use
- `restart_params`: List of restart parameters

**Returns:**
- Setup object ready for MCMC execution

## 🔧 Extreme Value Detection

The system automatically detects extreme values using parameter-specific bounds:

| Parameter Type | Valid Range | Fallback Strategy |
|----------------|-------------|-------------------|
| `precision` | 0.001 - 1000 | 2.0 (conservative) |
| `rho` | 0.001 - 0.999 | 0.5 (center) |
| `beta*` | -50 - 50 | 0.0 (no effect) |
| `site_effect*` | -10 - 10 | 0.0 (no effect) |
| `intercept` | -20 - 20 | 0.0 (no effect) |
| `NA/NaN/Inf` | Any | 0.0 (safe default) |

## 📈 Fallback Strategies

### Conservative (Recommended)
- Uses reasonable default values based on parameter type
- Safe for production use
- Preserves model structure

### Random
- Uses random values within reasonable bounds
- Good for exploring parameter space
- Less predictable results

### Zero
- Sets all extreme parameters to 0
- Simplest approach
- May lose important parameter information

## 🏗️ Integration with Model Fitting

### Enhanced Model Fitting Script
The `01_fitModels_with_restart.R` script extends the original model fitting with:

```r
# Automatic restart detection
RESTART_ENABLED <- TRUE  # Enable/disable restart functionality
RESTART_MIN_ESS <- 10    # Minimum ESS for parameter extraction
RESTART_FALLBACK_STRATEGY <- "conservative"  # Fallback strategy

# Enhanced convergence checking
check_continue_with_restart <- function(samples, model_info = NULL) {
  # Checks both convergence and extreme values
  # Returns restart recommendations
}
```

### Key Improvements
1. **Automatic Chain Detection**: Finds existing chains for each model
2. **Extreme Value Protection**: Prevents starting from extreme parameter values
3. **Enhanced Diagnostics**: Reports both convergence and extreme value issues
4. **Fallback Strategies**: Multiple approaches for handling problematic parameters
5. **Comprehensive Logging**: Detailed reports on restart decisions

## 📊 Batch Processing for 200+ Models

### Create Restart Report
```r
# See which models can be restarted
restart_report <- create_restart_report()
write.csv(restart_report, "restart_candidates.csv")
```

### Batch Restart Models
```r
# Process models in parallel
restart_results <- batch_restart_models(
  max_models = 200,           # Process up to 200 models
  parallel_cores = 8,         # Use 8 CPU cores
  restart_params = list(
    min_ess = 50,            # Higher ESS threshold for production
    max_rhat = 1.2,          # More lenient R-hat for rescue
    fallback_strategy = "conservative",
    niter = 5000,            # Longer chains for better convergence
    nchains = 3              # More chains for better diagnostics
  )
)
```

### Monitor Progress
The batch script provides real-time progress monitoring:
- Parallel task completion status
- Extreme value detection summaries
- Fallback value usage statistics
- Success/failure rates by model type

## 🎯 Expected Results

### For Models with Loose Priors
- **Before**: Extreme parameter values, poor convergence
- **After**: Reasonable parameter ranges, improved convergence

### Typical Improvements
- **Extreme Values**: 20-30 parameters detected and corrected per model
- **Convergence**: R-hat values improve from >2.0 to <1.2
- **Stability**: Eliminates NaN/Inf parameter values
- **Efficiency**: Better starting points reduce burnin time

## 🔍 Troubleshooting

### Common Issues

1. **No Chain Files Found**
   - Ensure model output directories exist
   - Check file naming conventions
   - Verify model completed at least one successful run

2. **All Parameters Extreme**
   - Model may need more fundamental changes
   - Consider different model structure
   - Review original data and model specification

3. **Poor Convergence After Restart**
   - Increase iterations (niter)
   - Try different fallback strategy
   - Review model priors for fundamental issues

### Diagnostic Tools

```r
# Check parameter distributions
hist(extraction_result$final_values, breaks = 50)

# Review extreme value detection
table(extreme_flags)

# Examine convergence diagnostics
summary(extraction_result$diagnostics$ess)
summary(extraction_result$diagnostics$rhat)
```

## 🚀 Production Workflow

### Step 1: Assessment
```bash
# Create inventory of restart candidates
Rscript -e "source('batch_restart_models.R'); create_restart_report()"
```

### Step 2: Test with Single Model
```bash
# Test restart with one model
Rscript demo_restart_functionality.R
```

### Step 3: Batch Processing
```bash
# Process all models in parallel
Rscript batch_restart_models.R 0 8  # 0 = all models, 8 cores
```

### Step 4: Validation
```bash
# Check results and convergence improvements
# Compare before/after parameter distributions
# Validate forecast performance
```

## 📋 File Dependencies

### Required Data Files
- `data/clean/model_input_df.csv` - Model configuration data
- `data/model_outputs/logit_beta_regression/*/samples_*.rds` - Existing chain files

### Required R Packages
- `tidyverse` - Data manipulation
- `coda` - MCMC diagnostics
- `parallel` - Parallel processing
- `doParallel` - Parallel backend
- `microbialForecast` - Project-specific functions

## 🎉 Success Metrics

### Individual Model Success
- ✅ Parameter values within reasonable bounds
- ✅ Improved convergence (R-hat < 1.2)
- ✅ No NaN/Inf values
- ✅ Effective sample sizes > 50 per parameter

### Batch Processing Success
- ✅ 80%+ models successfully restarted
- ✅ 90%+ extreme values detected and corrected
- ✅ Improved convergence across model suite
- ✅ Reduced model failures in production

## 📞 Support

For issues or questions:
1. Run the demonstration script to verify functionality
2. Check the test script output for error patterns
3. Review convergence diagnostics for problematic models
4. Consider adjusting restart parameters for specific model types

---

**🎯 Ready to rescue your 200+ models with improved initial values!**

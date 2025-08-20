# Microbial Forecasts: Spatio-temporal Prediction of Soil Microbiomes

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

This project uses state-space hierarchical Bayesian models to explore spatio-temporal variation and predictability across soil microbiomes of the United States. The models forecast relative abundance patterns for both taxonomic and functional groups using environmental and temporal predictors.

- Relative abundance forecasting for taxonomic groups
- Relative abundance forecasting for functional groups  
- Three different linear model structures for comprehensive analysis
- Separate models for fungi (ITS sequences) and bacteria (16S sequences)
- Plot-level analysis (20m × 20m) with soil core replicates
- ~2000 soil cores from 18 NEON sites across the United States
- Multiple sampling periods per year per site

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Data Sources](#data-sources)
- [Model Types](#model-types)
- [Usage Examples](#usage-examples)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## Installation

### Prerequisites

- R (≥ 4.0.0)
- RStudio (recommended)

### Required R Packages

The project uses the following core packages:

```r
# Core dependencies
install.packages(c(
  "tidyverse", "here", "nimble", "coda", "lubridate", 
  "reshape2", "dplyr", "pacman", "plyr", "tibble",
  "doParallel", "data.table", "Rfast", "moments",
  "scoringRules", "Metrics", "ggpubr"
))
```

### Package Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/microbialForecasts.git
   cd microbialForecasts
   ```

2. **Install the microbialForecast package:**
   ```r
   # From source
   install.packages("microbialForecast", repos = NULL, type = "source")
   ```

3. **Load the environment:**
   ```r
   source("source.R")
   ```

## Quick Start

### Basic Setup

```r
# Load the environment and packages
source("source.R")

# The package uses 'here' for relative file paths
library(here)
```

### Running Analysis

```r
# Example: Prepare taxonomic data for modeling
rank_data <- quick_get_rank_df(k = 1, 
                               min.date = "20151101", 
                               max.date = "20200101")

# Fit a model (example)
# See analysis/workflows/ for complete workflows
```

## Project Structure

```
microbialForecasts/
├── analysis/                    # Main analysis scripts
│   ├── model_analysis/         # Model fitting and analysis (scripts 0-10)
│   ├── create_figs/           # Figure generation scripts
│   └── workflows/             # Organized analysis workflows
│       ├── a_taxonomy/        # Taxonomic model workflows
│       ├── b_functional_groups/ # Functional group workflows
│       ├── c_diversity/       # Diversity model workflows
│       ├── d_taxonomy_dirichlet/ # Dirichlet taxonomic models
│       └── e_taxonomy_beta/   # Beta regression taxonomic models
├── data_construction/          # Data preparation scripts
│   ├── covariate_prep/        # Environmental covariate processing
│   └── microbe/               # Microbial data processing
├── microbialForecast/         # R package source code
├── data/                      # Data directory (see .gitignore)
│   ├── clean/                 # Processed model inputs
│   ├── model_outputs/         # Model results
│   └── summary/               # Summary statistics
├── figures/                   # Output figures
└── shinyapp/                  # Interactive visualization app
```

### Key Directories

- **`analysis/`**: Core analysis scripts for fitting and evaluating Bayesian state-space models
  - Scripts 0-4: Model creation and output processing
  - Scripts 5-8: Forecast creation and evaluation
  
- **`data_construction/`**: Scripts for downloading and preparing data from NEON and other sources

- **`microbialForecast/`**: R package containing functions for model fitting and evaluation

## Data Sources

- **Soil microbiome data**: NEON (National Ecological Observatory Network)
- **Environmental covariates**: 
  - Soil temperature and moisture (NEON sensors, DAYMET, SMOS)
  - Plant diversity and LAI (MODIS, NEON plant sampling)
  - Soil chemistry (NEON soil characterization)

## Model Types

### 1. Taxonomic Models
- **Target**: Relative abundance of taxonomic groups
- **Approach**: Hierarchical Bayesian state-space models
- **Variants**: Dirichlet and Beta regression formulations

### 2. Functional Group Models  
- **Target**: Relative abundance of functional categories
- **Categories**: Based on literature review, genomic pathways, and experimental enrichment
- **Kingdoms**: Separate bacterial and fungal functional classifications

### 3. Diversity Models
- **Target**: Alpha diversity metrics (Shannon diversity)
- **Approach**: Continuous response models with environmental predictors

### Model Structures
- **Environmental predictors**: Soil conditions, plant diversity, climate
- **Seasonality**: Cyclical temporal components  
- **Environmental + Seasonality**: Combined model structure

## Usage Examples

### Fitting a Single Taxon Model

```r
# Load required functions
source("source.R")

# Prepare data for a specific taxonomic rank
model_data <- prepTaxonomicData(rank.df = your_data, 
                                min.prev = 3,
                                min.date = "2015-11-01",
                                max.date = "2020-01-01")

# Fit model (see analysis/workflows/e_taxonomy_beta/ for complete examples)
```

### Creating Forecasts

```r
# Generate hindcast predictions
# See analysis/model_analysis/06_createHindcasts.r for examples
```

### Evaluating Model Performance

```r
# Calculate scoring metrics  
# See analysis/model_analysis/08_calculateScoringMetrics.r
```

## Model Workflow

1. **Data Preparation** (`data_construction/`)
   - Download and clean NEON data
   - Process environmental covariates
   - Prepare microbial abundance matrices

2. **Model Fitting** (`analysis/workflows/`)
   - Fit Bayesian models using NIMBLE
   - Assess convergence
   - Combine MCMC chains

3. **Forecasting** (`analysis/model_analysis/`)
   - Generate hindcast predictions
   - Calculate forecast horizons
   - Evaluate forecast accuracy

4. **Analysis and Visualization** (`analysis/create_figs/`)
   - Create publication figures
   - Analyze model performance
   - Explore spatio-temporal patterns

## Key Functions

The `microbialForecast` package provides:

- `prepTaxonomicData()`: Prepare taxonomic abundance data for modeling
- `prepFunctionalData()`: Prepare functional group data
- `prepDiversityData()`: Prepare diversity data
- `run_MCMC_*()`: Functions for running MCMC sampling
- `summarize_*_model()`: Model summary functions
- `add_scoring_metrics()`: Forecast evaluation metrics

## Computing Requirements

- **Memory**: 8GB+ RAM recommended for full model fitting
- **Storage**: ~50GB for complete data and model outputs
- **Compute**: Models designed for HPC clusters but can run locally
- **Time**: Individual model fits: hours to days depending on complexity

## Reproducibility

- Uses `here` package for portable file paths
- Version-controlled R package for consistent functions
- Documented workflows for each analysis type
- Environment setup via `source.R`

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use this code or approach in your research, please cite:

```
Werbin et al. 2024 preprint [link]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Zoey Werbin**  
Email: zoeywerbin@gmail.com 

## Acknowledgments

- NEON (National Ecological Observatory Network) for providing soil microbiome and environmental data
- NIMBLE development team for Bayesian modeling framework
- Very patient mentors and collaborators
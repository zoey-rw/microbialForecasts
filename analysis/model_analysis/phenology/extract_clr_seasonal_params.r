# Extract Seasonal Parameters from CLR Phenology Models
# This script extracts the sin/cos coefficients from the successfully run CLR models
# to calculate seasonal amplitude and peak timing for phenology analysis

cat("=== EXTRACTING SEASONAL PARAMETERS FROM CLR PHENOLOGY MODELS ===\n\n")

# Load required packages and environment
source("source.R")

cat("✅ Environment loaded successfully\n")

# Function to extract seasonal parameters from CLR model output
extract_seasonal_params <- function(file_path) {
  tryCatch({
    # Load the model output
    model_data <- readRDS(file_path)
    
    # Extract species name from filename
    species <- str_extract(basename(file_path), "CLR_phenology_(.+)_20130601_20151101")
    species <- str_remove(species, "CLR_phenology_")
    species <- str_remove(species, "_20130601_20151101")
    
    # Check if we have samples
    if (!"samples" %in% names(model_data)) {
      cat("   ❌ No samples found in", species, "\n")
      return(NULL)
    }
    
    samples <- model_data$samples
    
    # Check if we have multiple chains
    if (length(samples) < 2) {
      cat("   ⚠️ Only", length(samples), "chain(s) for", species, "\n")
      return(NULL)
    }
    
    # Extract seasonal coefficients (beta parameters)
    # For cycl_only models, we expect:
    # beta[1,1] = sin coefficient for species 1
    # beta[2,1] = cos coefficient for species 1
    # beta[1,2] = sin coefficient for species 2 (other)
    # beta[2,2] = cos coefficient for species 2 (other)
    
    # We want the coefficients for the focal species (first species)
    sin_coef <- mean(samples[[1]][, "beta[1, 1]"])
    cos_coef <- mean(samples[[2]][, "beta[2, 1]"])
    
    # Calculate seasonal amplitude and peak timing
    amplitude <- sqrt(sin_coef^2 + cos_coef^2)
    peak_timing <- atan2(cos_coef, sin_coef) * 12/(2*pi)
    
    # Convert to positive months (0-12)
    if (peak_timing < 0) peak_timing <- peak_timing + 12
    
    # Create result structure
    result <- list(
      species = species,
      sin_coef = sin_coef,
      cos_coef = cos_coef,
      amplitude = amplitude,
      peak_timing = peak_timing,
      peak_month = round(peak_timing, 1),
      n_chains = length(samples),
      n_iterations = nrow(samples[[1]]),
      n_parameters = ncol(samples[[1]])
    )
    
    cat("   ✅", species, "- Amplitude:", round(amplitude, 3), 
        "Peak Month:", round(peak_timing, 1), "\n")
    
    return(result)
    
  }, error = function(e) {
    cat("   ❌ Error processing", basename(file_path), ":", e$message, "\n")
    return(NULL)
  })
}

# Main execution
cat("🔍 Analyzing CLR phenology model outputs...\n")

# Path to CLR phenology models
clr_dir <- here("data/model_outputs/CLR_phenology_2013_2015")

if (!dir.exists(clr_dir)) {
  cat("❌ CLR phenology directory not found:", clr_dir, "\n")
  stop("Please run the phenology models first")
}

# List all CLR model files
clr_files <- list.files(clr_dir, pattern = "CLR_phenology_.*\\.rds$", full.names = TRUE)

if (length(clr_files) == 0) {
  cat("❌ No CLR phenology model files found\n")
  stop("Please run the phenology models first")
}

cat("📁 Found", length(clr_files), "CLR phenology model files\n\n")

# Extract seasonal parameters from each model
cat("📊 Extracting Seasonal Parameters...\n")
seasonal_params <- list()

for (file in clr_files) {
  params <- extract_seasonal_params(file)
  if (!is.null(params)) {
    seasonal_params[[params$species]] <- params
  }
}

# Create summary data frame
if (length(seasonal_params) > 0) {
  cat("\n📋 SEASONAL PARAMETERS SUMMARY\n")
  cat("==============================\n")
  
  # Convert to data frame
  summary_df <- do.call(rbind, lapply(seasonal_params, function(x) {
    data.frame(
      Species = x$species,
      Sin_Coefficient = round(x$sin_coef, 4),
      Cos_Coefficient = round(x$cos_coef, 4),
      Amplitude = round(x$amplitude, 4),
      Peak_Timing = round(x$peak_timing, 2),
      Peak_Month = x$peak_month,
      N_Chains = x$n_chains,
      N_Iterations = x$n_iterations,
      N_Parameters = x$n_parameters,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by amplitude (strongest seasonal patterns first)
  summary_df <- summary_df[order(-summary_df$Amplitude), ]
  
  # Display summary
  print(summary_df)
  
  # Save results
  output_file <- here("figures/phenology_comparison/CLR_seasonal_parameters_2013_2015.csv")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, output_file, row.names = FALSE)
  
  cat("\n💾 Results saved to:", output_file, "\n")
  
  # Summary statistics
  cat("\n📊 SUMMARY STATISTICS\n")
  cat("====================\n")
  cat("Total species analyzed:", nrow(summary_df), "\n")
  cat("Mean amplitude:", round(mean(summary_df$Amplitude), 4), "\n")
  cat("Median amplitude:", round(median(summary_df$Amplitude), 4), "\n")
  cat("Range amplitude:", round(min(summary_df$Amplitude), 4), "to", round(max(summary_df$Amplitude), 4), "\n")
  cat("Mean peak month:", round(mean(summary_df$Peak_Timing), 2), "\n")
  
  # Identify strongest seasonal patterns
  cat("\n🏆 STRONGEST SEASONAL PATTERNS\n")
  cat("==============================\n")
  top_species <- head(summary_df, 3)
  for (i in 1:nrow(top_species)) {
    cat(i, ".", top_species$Species[i], "- Amplitude:", top_species$Amplitude[i], 
        "Peak Month:", top_species$Peak_Month[i], "\n")
  }
  
} else {
  cat("❌ No seasonal parameters could be extracted\n")
}

cat("\n=== SEASONAL PARAMETER EXTRACTION COMPLETE ===\n")

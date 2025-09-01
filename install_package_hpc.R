#!/usr/bin/env Rscript

# Installation script for microbialForecast package on HPC
# This script handles the complete installation process with proper error handling

cat("=== Installing microbialForecast package on HPC ===\n")

# Function to install packages with better error handling
install_with_retry <- function(pkg, retries = 3) {
  for (i in 1:retries) {
    cat("Attempting to install", pkg, "(attempt", i, "of", retries, ")\n")
    result <- try({
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
    }, silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      cat("✓ Successfully installed", pkg, "\n")
      return(TRUE)
    } else {
      cat("✗ Failed to install", pkg, "- attempt", i, "\n")
      if (i < retries) {
        cat("Waiting 10 seconds before retry...\n")
        Sys.sleep(10)
      }
    }
  }
  return(FALSE)
}

# Step 1: Install critical dependencies first
cat("\n1. Installing critical dependencies...\n")
critical_deps <- c("devtools", "roxygen2", "nimble", "dplyr", "tidyverse")

for (dep in critical_deps) {
  if (!require(dep, character.only = TRUE, quietly = TRUE)) {
    if (!install_with_retry(dep)) {
      stop("Failed to install critical dependency: ", dep)
    }
  } else {
    cat("✓", dep, "already installed\n")
  }
}

# Step 2: Load devtools
library(devtools)

# Step 3: Set working directory and build package
cat("\n2. Building package...\n")
if (!file.exists("microbialForecast/DESCRIPTION")) {
  stop("DESCRIPTION file not found. Make sure you're in the correct directory.")
}

setwd("microbialForecast")

# Step 4: Document the package (this will regenerate .Rd files with correct format)
cat("\n3. Generating documentation...\n")
tryCatch({
  document()
  cat("✓ Documentation generated successfully\n")
}, error = function(e) {
  cat("✗ Documentation generation failed:", e$message, "\n")
  stop("Cannot proceed without proper documentation")
})

# Step 5: Build the package
cat("\n4. Building package...\n")
build_result <- tryCatch({
  build()
}, error = function(e) {
  cat("✗ Package build failed:", e$message, "\n")
  return(NULL)
})

if (is.null(build_result)) {
  stop("Package build failed")
} else {
  cat("✓ Package built successfully:", build_result, "\n")
}

# Step 6: Install the built package
cat("\n5. Installing package...\n")
install_result <- tryCatch({
  install.packages(build_result, repos = NULL, type = "source")
  cat("✓ Package installed successfully!\n")
  TRUE
}, error = function(e) {
  cat("✗ Package installation failed:", e$message, "\n")
  FALSE
})

# Step 7: Test the installation
if (install_result) {
  cat("\n6. Testing installation...\n")
  test_result <- tryCatch({
    library(microbialForecast)
    cat("✓ Package loads successfully!\n")
    
    # Test if key functions are available
    if (exists("nimbleModFunctional")) {
      cat("✓ NIMBLE models are available\n")
    } else {
      cat("⚠ Warning: Some NIMBLE models may not be available\n")
    }
    
    TRUE
  }, error = function(e) {
    cat("✗ Package test failed:", e$message, "\n")
    FALSE
  })
  
  if (test_result) {
    cat("\n=== Installation completed successfully! ===\n")
  } else {
    cat("\n=== Installation completed with warnings ===\n")
  }
} else {
  cat("\n=== Installation failed ===\n")
}

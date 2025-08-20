#!/usr/bin/env Rscript

# Post-processing pipeline for microbial forecasting models
# Runs steps 2-8 of the post-processing workflow

cat("Starting post-processing pipeline (steps 2-8)\n")
cat("================================================\n\n")

# Source the main environment setup
source("../../source.R")

cat("✅ Environment setup complete\n\n")

# Step 2: Combine Model Chains
cat("📊 Step 2: Combining model chains...\n")
tryCatch({
  source("02_combineModelChains.r")
  cat("✅ Step 2 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 2 failed:", e$message, "\n\n")
})

# Step 3: Summarize Model Outputs
cat("📊 Step 3: Summarizing model outputs...\n")
tryCatch({
  source("03_summarizeModelOutputs.r")
  cat("✅ Step 3 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 3 failed:", e$message, "\n\n")
})

# Step 4: Tidy Effect Sizes
cat("📊 Step 4: Tidying effect sizes...\n")
tryCatch({
  source("04_tidyEffectSizes.r")
  cat("✅ Step 4 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 4 failed:", e$message, "\n\n")
})

# Step 5: Predict Site Effects
cat("📊 Step 5: Predicting site effects...\n")
tryCatch({
  source("05_predictSiteEffects.r")
  cat("✅ Step 5 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 5 failed:", e$message, "\n\n")
})

# Step 6: Create Hindcasts
cat("📊 Step 6: Creating hindcasts...\n")
tryCatch({
  source("06_createHindcasts.r")
  cat("✅ Step 6 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 6 failed:", e$message, "\n\n")
})

# Step 7: Tidy hindcasts
cat("📊 Step 7: Tidying hindcasts...\n")
tryCatch({
  source("07_tidyHindcasts.r")
  cat("✅ Step 7 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 7 failed:", e$message, "\n\n")
})

# Step 8: Calculate Scoring Metrics
cat("📊 Step 8: Calculating scoring metrics...\n")
tryCatch({
  source("08_calculateScoringMetrics.r")
  cat("✅ Step 8 completed successfully\n\n")
}, error = function(e) {
  cat("❌ Step 8 failed:", e$message, "\n\n")
})

cat("Post-processing pipeline completed!\n")
cat("================================================\n")

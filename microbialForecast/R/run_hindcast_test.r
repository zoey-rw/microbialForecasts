
# Enhanced fcast_logit_beta with comprehensive test mode
fcast_logit_beta_test <- function(plotID,
                                 model.inputs,
                                 param_samples,
                                 truth.plot.long,
                                 plot_summary = NULL,
                                 Nmc = 1000,
                                 predict_site_effects = NULL,
                                 rank.name = NULL,
                                 model_id,
                                 metadata = NULL,
                                 test_mode = FALSE,
                                 test_config = NULL,
                                 ...) {
  
  # Test mode validation
  if (test_mode && !is.null(test_config)) {
    cat("🧪 TEST MODE: Enhanced validation enabled
")
    
    # Validate test configuration
    if (!is.list(test_config)) {
      stop("test_config must be a list")
    }
    
    # Apply test mode limits
    if (test_config$NMC_SAMPLES < Nmc) {
      cat("  ⚠️  Reducing NMC from", Nmc, "to", test_config$NMC_SAMPLES, "for test mode
")
      Nmc <- test_config$NMC_SAMPLES
    }
    
    # Test mode progress reporting
    if (test_config$VERBOSE_OUTPUT) {
      cat("  📊 Test mode parameters:
")
      cat("    - NMC samples:", Nmc, "
")
      cat("    - Model ID:", model_id, "
")
      cat("    - Plot ID:", plotID, "
")
      cat("    - Rank name:", rank.name, "
")
    }
  }
  
  # Call the original function with test mode parameters
  result <- fcast_logit_beta(plotID = plotID,
                            model.inputs = model.inputs,
                            param_samples = param_samples,
                            truth.plot.long = truth.plot.long,
                            plot_summary = plot_summary,
                            Nmc = Nmc,
                            predict_site_effects = predict_site_effects,
                            rank.name = rank.name,
                            model_id = model_id,
                            metadata = metadata,
                            ...)
  
  # Test mode validation of results
  if (test_mode && !is.null(result)) {
    cat("  ✅ Test mode validation:
")
    cat("    - Result dimensions:", nrow(result), "x", ncol(result), "
")
    cat("    - Required columns present:", all(c("plotID", "siteID", "species", "dateID", "dates", "med", "lo", "hi", "fcast_period") %in% colnames(result)), "
")
    cat("    - Truth data rows:", sum(!is.na(result$truth)), "
")
    cat("    - Calibration rows:", sum(result$fcast_period == "calibration"), "
")
    cat("    - Hindcast rows:", sum(result$fcast_period == "hindcast"), "
")
  }
  
  return(result)
}


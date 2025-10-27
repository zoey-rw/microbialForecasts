# Logging setup for microbial forecasts
# Provides structured logging with levels, timestamps, and file output

# Load logger package
library(logger)

# Setup logging configuration
log_setup <- function(logfile = NULL, level = "INFO") {
    # Set log level
    logger::log_threshold(level)
    
    # Set up appenders
    if (!is.null(logfile)) {
        # Write to both console and file
        logger::log_appender(logger::appender_tee(logfile))
    } else {
        # Console only
        logger::log_appender(logger::appender_console)
    }
    
    # Set log format with timestamp and level (use simple formatter to avoid glue issues)
    logger::log_formatter(logger::formatter_sprintf)
    logger::log_layout(logger::layout_simple)
}

# Convenience functions for different log levels
info  <- function(msg, ...) logger::log_info(msg, ...)
warn  <- function(msg, ...) logger::log_warn(msg, ...)
error <- function(msg, ...) logger::log_error(msg, ...)
debug <- function(msg, ...) logger::log_debug(msg, ...)
trace <- function(msg, ...) logger::log_trace(msg, ...)

# Progress logging for long-running operations
log_progress <- function(current, total, operation = "Processing") {
    pct <- round(current / total * 100, 1)
    info("{operation}: {current}/{total} ({pct}%)")
}

# Model-specific logging
log_model_start <- function(model_name, species, rank_name) {
    info("=== Starting model fitting ===")
    info("Model: {model_name}")
    info("Species: {species}")
    info("Rank: {rank_name}")
}

log_model_complete <- function(model_name, species, runtime_minutes) {
    info("=== Model fitting completed ===")
    info("Model: {model_name}")
    info("Species: {species}")
    info("Runtime: {round(runtime_minutes, 2)} minutes")
}

# MCMC-specific logging
log_mcmc_progress <- function(iterations, ess_min, status = "Running") {
    info("MCMC {status}: {iterations} iterations, min ESS: {round(ess_min, 1)}")
}

log_mcmc_convergence <- function(converged, min_ess, target_ess) {
    if (converged) {
        info("✓ MCMC converged: min ESS {round(min_ess, 1)} >= target {target_ess}")
    } else {
        warn("⚠ MCMC not converged: min ESS {round(min_ess, 1)} < target {target_ess}")
    }
}

# Error logging with context
log_error_with_context <- function(error_msg, context = list()) {
    error("Error: {error_msg}")
    if (length(context) > 0) {
        for (name in names(context)) {
            error("  {name}: {context[[name]]}")
        }
    }
}

# Performance logging
log_performance <- function(operation, duration_seconds, details = NULL) {
    info("Performance: {operation} took {round(duration_seconds, 2)}s")
    if (!is.null(details)) {
        info("  Details: {details}")
    }
}

# Data validation logging
log_data_validation <- function(validation_name, passed, details = NULL) {
    if (passed) {
        info("✓ Data validation passed: {validation_name}")
    } else {
        error("✗ Data validation failed: {validation_name}")
    }
    if (!is.null(details)) {
        info("  Details: {details}")
    }
}

# Export the logging functions
cat("Logging setup complete. Available functions:\n")
cat("  - log_setup(logfile, level)\n")
cat("  - info(), warn(), error(), debug(), trace()\n")
cat("  - log_progress(), log_model_start(), log_model_complete()\n")
cat("  - log_mcmc_progress(), log_mcmc_convergence()\n")
cat("  - log_error_with_context(), log_performance(), log_data_validation()\n")

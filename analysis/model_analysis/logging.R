# Logging setup for microbial forecasts
# Provides structured logging with levels, timestamps, and file output

# Try to load logger package, fallback to simple logging if not available
has_logger <- requireNamespace("logger", quietly = TRUE)
if (has_logger) {
  library(logger)
}

# Setup logging configuration
log_setup <- function(logfile = NULL, level = "INFO") {
    # Ensure the log file's directory exists (logger's appender_tee won't create it)
    if (!is.null(logfile)) {
        log_dir <- dirname(logfile)
        if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (has_logger) {
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
    } else {
        # Fallback: simple file connection if logfile provided
        if (!is.null(logfile)) {
            # Ensure directory exists
            log_dir <- dirname(logfile)
            if (!dir.exists(log_dir)) {
                dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
            }
            # Open log file connection (will be used by fallback functions)
            assign(".logfile", logfile, envir = .GlobalEnv)
        }
        cat("Using fallback logging (logger package not available)\n")
    }
}

# Convenience functions for different log levels
if (has_logger) {
    info  <- function(msg, ...) logger::log_info(msg, ...)
    warn  <- function(msg, ...) logger::log_warn(msg, ...)
    error <- function(msg, ...) logger::log_error(msg, ...)
    debug <- function(msg, ...) logger::log_debug(msg, ...)
    trace <- function(msg, ...) logger::log_trace(msg, ...)
} else {
    # Fallback logging functions using cat/print
    .log_write <- function(level, msg, ...) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        formatted_msg <- if (...length() > 0) sprintf(msg, ...) else msg
        output <- sprintf("[%s] [%s] %s\n", timestamp, level, formatted_msg)
        cat(output)
        # Also write to file if logfile was set
        if (exists(".logfile", envir = .GlobalEnv)) {
            logfile <- get(".logfile", envir = .GlobalEnv)
            tryCatch({
                cat(output, file = logfile, append = TRUE)
            }, error = function(e) {
                # Silently fail if file write fails
            })
        }
    }
    info  <- function(msg, ...) .log_write("INFO", msg, ...)
    warn  <- function(msg, ...) .log_write("WARN", msg, ...)
    error <- function(msg, ...) .log_write("ERROR", msg, ...)
    debug <- function(msg, ...) .log_write("DEBUG", msg, ...)
    trace <- function(msg, ...) .log_write("TRACE", msg, ...)
}

# Progress logging for long-running operations
log_progress <- function(current, total, operation = "Processing") {
    pct <- round(current / total * 100, 1)
    if (has_logger) {
        info("{operation}: {current}/{total} ({pct}%)")
    } else {
        info("%s: %d/%d (%.1f%%)", operation, current, total, pct)
    }
}

# Model-specific logging
log_model_start <- function(model_name, species, rank_name) {
    info("=== Starting model fitting ===")
    if (has_logger) {
        info("Model: {model_name}")
        info("Species: {species}")
        info("Rank: {rank_name}")
    } else {
        info("Model: %s", model_name)
        info("Species: %s", species)
        info("Rank: %s", rank_name)
    }
}

log_model_complete <- function(model_name, species, runtime_minutes) {
    info("=== Model fitting completed ===")
    if (has_logger) {
        info("Model: {model_name}")
        info("Species: {species}")
        info("Runtime: {round(runtime_minutes, 2)} minutes")
    } else {
        info("Model: %s", model_name)
        info("Species: %s", species)
        info("Runtime: %.2f minutes", round(runtime_minutes, 2))
    }
}

# MCMC-specific logging
log_mcmc_progress <- function(iterations, ess_min, status = "Running") {
    if (has_logger) {
        info("MCMC {status}: {iterations} iterations, min ESS: {round(ess_min, 1)}")
    } else {
        info("MCMC %s: %d iterations, min ESS: %.1f", status, iterations, round(ess_min, 1))
    }
}

log_mcmc_convergence <- function(converged, min_ess, target_ess) {
    if (converged) {
        if (has_logger) {
            info("✓ MCMC converged: min ESS {round(min_ess, 1)} >= target {target_ess}")
        } else {
            info("✓ MCMC converged: min ESS %.1f >= target %d", round(min_ess, 1), target_ess)
        }
    } else {
        if (has_logger) {
            warn("⚠ MCMC not converged: min ESS {round(min_ess, 1)} < target {target_ess}")
        } else {
            warn("⚠ MCMC not converged: min ESS %.1f < target %d", round(min_ess, 1), target_ess)
        }
    }
}

# Error logging with context
log_error_with_context <- function(error_msg, context = list()) {
    if (has_logger) {
        error("Error: {error_msg}")
    } else {
        error("Error: %s", error_msg)
    }
    if (length(context) > 0) {
        for (name in names(context)) {
            if (has_logger) {
                error("  {name}: {context[[name]]}")
            } else {
                error("  %s: %s", name, context[[name]])
            }
        }
    }
}

# Performance logging
log_performance <- function(operation, duration_seconds, details = NULL) {
    if (has_logger) {
        info("Performance: {operation} took {round(duration_seconds, 2)}s")
    } else {
        info("Performance: %s took %.2fs", operation, round(duration_seconds, 2))
    }
    if (!is.null(details)) {
        if (has_logger) {
            info("  Details: {details}")
        } else {
            info("  Details: %s", details)
        }
    }
}

# Data validation logging
log_data_validation <- function(validation_name, passed, details = NULL) {
    if (passed) {
        if (has_logger) {
            info("✓ Data validation passed: {validation_name}")
        } else {
            info("✓ Data validation passed: %s", validation_name)
        }
    } else {
        if (has_logger) {
            error("✗ Data validation failed: {validation_name}")
        } else {
            error("✗ Data validation failed: %s", validation_name)
        }
    }
    if (!is.null(details)) {
        if (has_logger) {
            info("  Details: {details}")
        } else {
            info("  Details: %s", details)
        }
    }
}

# Export the logging functions
if (has_logger) {
    cat("Logging setup complete (using logger package). Available functions:\n")
} else {
    cat("Logging setup complete (using fallback logging). Available functions:\n")
}
cat("  - log_setup(logfile, level)\n")
cat("  - info(), warn(), error(), debug(), trace()\n")
cat("  - log_progress(), log_model_start(), log_model_complete()\n")
cat("  - log_mcmc_progress(), log_mcmc_convergence()\n")
cat("  - log_error_with_context(), log_performance(), log_data_validation()\n")

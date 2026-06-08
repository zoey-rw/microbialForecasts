#!/usr/bin/env Rscript
# Fit multivariate Dirichlet composition models with driver uncertainty
# - Weak priors for main parameters (Jeffreys, Uniform, wide normal)
# - More informative priors for site effects (dgamma(2, 20))


# Load required packages
if (!require(here)) {
    install.packages("here")
    library(here)
}

# Setup logging first
source(here("analysis/model_analysis/logging.R"))
log_setup(logfile = here("analysis", "model_analysis", "logs", paste0("model_fitting_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")))

# Set project root
here::i_am("analysis/model_analysis/01_fitModels_dirichlet.R")
project_root <- here()

# Validate that the script path matches the expected location
if (!file.exists(here::here("analysis/model_analysis/01_fitModels_dirichlet.R"))) {
    stop("Script path validation failed: analysis/model_analysis/01_fitModels_dirichlet.R not found at here() root")
}

info("here() starts at %s", project_root)
info("Project root set to: %s", getwd())
info("Script path validated: analysis/model_analysis/01_fitModels_dirichlet.R exists")

# Load the microbialForecast package to access helper functions
library(microbialForecast)

# Load nimble explicitly for options
library(nimble)

# Set nimbleOptions for faster compilation
nimbleOptions(buildInterfacesForCompiledNimbleFunctions = FALSE) # faster compile
nimbleOptions(optimize = TRUE)
info("✓ NIMBLE options set for faster compilation")


# One source of truth for driver uncertainty mode
driver_uncertainty_mode <- TRUE  # Set FALSE for fixed drivers

# ---- BEGIN: compat shim to standardize output structure ----
ensure_old_schema <- function(samples, samples2, meta,
                              model_output_root, model_name, species,
                              model_id, chain_no) {

  # 0) Validation guards - prevent empty/NA directory names
  .ensure_nonempty <- function(x, nm) {
    if (is.null(x) || is.na(x) || !nzchar(as.character(x))) {
      stop(sprintf("ensure_old_schema(): '%s' is empty/NA; refusing to save.", nm))
    }
  }
  .ensure_nonempty(model_name, "model_name")
  .ensure_nonempty(species, "species")
  .ensure_nonempty(model_id, "model_id")

  # 1) Ensure matrices
  if (inherits(samples, "mcmc.list")) samples <- as.matrix(do.call(rbind, samples))
  if (inherits(samples, "mcmc"))      samples <- as.matrix(samples)
  if (!is.matrix(samples))            samples <- as.matrix(samples)

  # CRITICAL FIX: Don't fall back to samples if samples2 is missing/empty
  # samples2 should contain plot-level estimates (mu[plot, time]), not parameters
  if (missing(samples2) || is.null(samples2)) {
    # Don't fall back to samples - samples2 should be separate
    warning("samples2 is missing or NULL - expected mu[plot, time] columns from mvSamples2")
    samples2 <- matrix(nrow = 0, ncol = 0)
  }
  if (inherits(samples2, "mcmc.list")) samples2 <- as.matrix(do.call(rbind, samples2))
  if (inherits(samples2, "mcmc"))      samples2 <- as.matrix(samples2)
  if (!is.matrix(samples2))            samples2 <- as.matrix(samples2)
  # CRITICAL FIX: Don't fall back to samples if samples2 is empty - this causes wrong data
  # If mvSamples2 is empty, we should detect this as an error, not silently use parameters
  if (is.matrix(samples2) && nrow(samples2) == 0) {
    warning("samples2 is empty (0 rows) - mvSamples2 may not have been populated. ",
            "This may indicate thin2 is too large or monitors2 is not working correctly.")
    # Keep empty matrix instead of falling back to parameters
    samples2 <- matrix(nrow = 0, ncol = 0)
  }

  # 2) Rename latent columns to plot_mu[*] (drop Ex if present)
  cn2 <- colnames(samples2)
  if (!is.null(cn2)) {
    # CRITICAL FIX: Check for mu/plot_mu/plot_rel columns (Dirichlet uses plot_rel)
    keep <- grepl("^mu\\[", cn2) | grepl("^plot_mu\\[", cn2) | grepl("^plot_rel\\[", cn2)
    
    # CRITICAL FIX: If samples2 contains parameter columns instead of mu columns, this is wrong
    # Check for parameter columns to detect when samples2 was incorrectly populated
    has_param_cols <- any(grepl("^beta\\[", cn2)) || any(grepl("^intercept", cn2)) || 
                      any(grepl("^site_effect\\[", cn2)) || any(grepl("^precision", cn2)) ||
                      any(grepl("^rho", cn2)) || any(grepl("^legacy_effect", cn2))
    
    if (has_param_cols && !any(keep)) {
      # CRITICAL ERROR: samples2 contains parameters, not plot estimates
      # This means mvSamples2 was empty or wrong, and we fell back to samples (parameters)
      stop("CRITICAL: samples2 contains parameter columns (beta, intercept, etc.) instead of mu[plot, time] columns.\n",
           "This indicates mvSamples2 was empty or incorrectly populated.\n",
           "Found columns like: ", paste(head(cn2[grepl("^beta\\[|^intercept|^site_effect", cn2)], 5), collapse=", "), "\n",
           "Expected mu[plot, time] columns from mvSamples2 monitors2.")
    }
    
    if (any(keep)) {
      samples2 <- samples2[, keep, drop = FALSE]
      cn2      <- colnames(samples2)
    } else {
      # No mu/plot_rel columns found - set to empty matrix with warning
      warning("No mu[...], plot_mu[...], or plot_rel[...] columns found in samples2. ",
              "This may indicate mvSamples2 was empty or not properly monitored.")
      samples2 <- matrix(nrow = nrow(samples2), ncol = 0)
      cn2 <- character(0)
    }
    
    # convert mu[...] -> plot_mu[...]; leave plot_rel[...] as is
    if (length(cn2) > 0) {
      cn2 <- sub("^mu\\[", "plot_mu[", cn2)
      colnames(samples2) <- cn2
    }
  }

  # 3) Minimal metadata set with consistent keys/types
  required_meta <- list(
    rank.name           = meta$rank.name              %||% NA_character_,
    species             = meta$species                %||% NA_character_,
    model_name          = meta$model_name             %||% NA_character_,
    model_id            = meta$model_id               %||% NA_character_,
    use_legacy_covariate= isTRUE(meta$use_legacy_covariate),
    has_driver_uncertainty = isTRUE(meta$has_driver_uncertainty),
    scenario            = meta$scenario               %||% NA_character_,
    min.date            = meta$min.date               %||% NA_character_,
    max.date            = meta$max.date               %||% NA_character_,
    niter               = meta$niter                  %||% nrow(samples),
    nburnin             = meta$nburnin                %||% 0L,
    thin                = meta$thin                   %||% 1L,
    thin2               = meta$thin2                  %||% 20L,
    model_structure     = meta$model_structure        %||% "dirichlet_driver_uncertainty"
  )
  # keep any extras too
  for (nm in setdiff(names(meta), names(required_meta))) required_meta[[nm]] <- meta[[nm]]

  # 4) Standard file path (always include species folder)
  species_dir  <- file.path(model_output_root, model_name, species)
  dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
  samples_file <- file.path(species_dir, paste0("samples_", model_id, "_chain", chain_no, ".rds"))

  # 5) Belt-and-suspenders parent directory validation (resilient to symlinks/trailing slashes)
  canon <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  expected_parent <- sub("/+$", "", canon(species_dir))
  actual_parent <- sub("/+$", "", canon(dirname(samples_file)))
  if (expected_parent != actual_parent) {
    stop(sprintf("Path validation failed: refusing to save outside expected directory.\n  expected parent: %s\n  actual parent:   %s",
                 expected_parent, actual_parent))
  }

  # 6) Return final list and path
  list(
    chain_output = list(samples = samples, samples2 = samples2, metadata = required_meta),
    path = samples_file
  )
}

`%||%` <- function(a,b) if(is.null(a)) b else a

# Atomic save wrapper for safe writes on distributed filesystems
atomic_saveRDS <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid(), "_", sample.int(1e6, 1))
  on.exit({ if (file.exists(tmp)) try(unlink(tmp), silent = TRUE) }, add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (!file.rename(tmp, path)) {
    # On some FS, file.rename across devices fails; try copy
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
    if (!ok) stop("Atomic save failed: rename/copy both failed for ", path)
  }
  invisible(path)
}

# Optional safety net: rebuild mu from eta if needed
fix_mu_from_eta <- function(S2, N.plot, N.date) {
  has_mu  <- any(grepl("^mu\\[",  colnames(S2)))
  has_eta <- any(grepl("^eta\\[", colnames(S2)))
  has_plot_rel <- any(grepl("^plot_rel\\[", colnames(S2)))
  if (has_mu || has_plot_rel || !has_eta) return(S2)

  eta_cols <- grep("^eta\\[", colnames(S2), value = TRUE)
  eta_mat  <- as.matrix(S2[, eta_cols, drop = FALSE])
  mu_mat   <- 1 - exp(-exp(eta_mat))
  colnames(mu_mat) <- sub("^eta\\[", "mu[", colnames(eta_mat))
  cbind(S2, mu_mat)
}
# ---- END: compat shim ----

# Load packages and create directories using package functions
load_required_packages()

# Ensure parallel packages are loaded
library(parallel)
library(doParallel)
library(foreach)
library(coda)
create_directories_safe(
    here("data", "model_outputs"), 
    c(if (driver_uncertainty_mode) "dirichlet_driver_uncertainty"
      else "dirichlet_fixed_drivers")
)

# Define output directory based on mode (optional suffix for rho dunif prior test)
base_dir_name <- if (driver_uncertainty_mode) "dirichlet_driver_uncertainty" else "dirichlet_fixed_drivers"
if (nzchar(Sys.getenv("RHO_PRIOR_UNIF", ""))) {
  base_dir_name <- paste0(base_dir_name, "_rho_unif")
  cat("RHO_PRIOR_UNIF set: writing to", base_dir_name, "\n")
}
if (nzchar(Sys.getenv("OUTPUT_SUFFIX", ""))) {
  base_dir_name <- paste0(base_dir_name, "_", Sys.getenv("OUTPUT_SUFFIX"))
  cat("OUTPUT_SUFFIX set: writing to", base_dir_name, "\n")
}
model_output_dir <- here("data", "model_outputs", base_dir_name)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

# Diagnostic: confirm writeability + actual target path
cat("model_output_dir =", model_output_dir, "\n")
cat("file.exists(dir)?", dir.exists(model_output_dir), "\n")
cat("writable? (0 means ok) =>", file.access(model_output_dir, 2), "\n")

info("==================================================")
info("Microbial forecasts environment setup complete!")
info("Ready for analysis.")
info("==================================================")

# ----------------- HPC / CLI ARGUMENTS -----------------
argv <- commandArgs(TRUE)

# helpers
get_flag <- function(name) any(grepl(paste0("^", name, "$"), argv))
get_opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), argv, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", name, "="), "", hit[1]))
  idx <- match(name, argv)
  if (!is.na(idx) && idx < length(argv)) return(argv[idx + 1L])
  return(default)
}

# Chains: default 4 for HPC, 2 for LOCAL_TEST
nchains <- as.integer(Sys.getenv("NCHAINS", ifelse(identical(tolower(Sys.getenv("LOCAL_TEST","false")),"true"), "2", "4")))

# Priority switch stays compatible with your existing logic
use_priority <- Sys.getenv("USE_PRIORITY", ifelse(get_flag("priority") || get_flag("--priority"), "true", "false"))

# Selection methods (only one needs to be provided)
cli_model_id   <- get_opt("--model-id",   NULL)
cli_rank       <- get_opt("--rank",       NULL)
cli_species    <- get_opt("--species",    NULL)
cli_model_name <- get_opt("--model-name", NULL)

# Linear task (1..n_models*nchains) from CLI or scheduler env
cli_task      <- suppressWarnings(as.integer(get_opt("--task", NA)))
sge_task      <- suppressWarnings(as.integer(Sys.getenv("SGE_TASK_ID", NA)))
pbs_task      <- suppressWarnings(as.integer(Sys.getenv("PBS_ARRAYID", NA)))
slurm_task    <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
array_task_id <- suppressWarnings(as.integer(na.omit(c(cli_task, sge_task, pbs_task, slurm_task))[1]))

# Back-compat "k" (ignored previously); we still allow it if provided
k <- suppressWarnings(as.integer(if (length(argv) > 0) argv[1] else NA))
if (!is.na(k) && is.na(array_task_id)) array_task_id <- k
# -------------------------------------------------------

#### Run on all groups ----

# Load data early for filtering
info("Loading data files for filtering...")
bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi)
info("Data loaded successfully for %d ranks", length(all_ranks))

# Function to check if MCMC should continue based on effective sample size
check_continue <- function(samples, min_eff_size = 50) {
    # Validate input
    if (is.null(samples) || nrow(samples) == 0) {
        warn("Empty samples; continuing")
        return(TRUE)
    }
    
    # Calculate effective sample sizes safely (ignoring zero-variance columns)
    eff_sizes <- apply(samples, 2, function(x) {
        if (sum(!is.na(x)) < 2 || is.na(var(x, na.rm = TRUE)) || var(x, na.rm = TRUE) == 0) {
            return(NA)
        }
        tryCatch(coda::effectiveSize(coda::as.mcmc(x)), error = function(e) NA)
    })
    
    # Check if any valid parameter has ESS below threshold
    valid_ess <- eff_sizes[!is.na(eff_sizes)]
    if (length(valid_ess) == 0) {
        warn("  WARNING: No parameters have valid variation for ESS calculation.")
        return(TRUE) # Continue sampling
    }
    
    min_ess <- min(valid_ess)
    continue <- is.finite(min_ess) && (min_ess < min_eff_size)
    
    info("ESS check - Min ESS: %.1f, Target: %d", round(min_ess, 1), min_eff_size)
    info("Continue sampling: %s", continue)
    
    return(continue)
}

params_in = read.csv(here("data/clean/model_input_df.csv"),
                     colClasses = c(rep("character", 4),
                                    rep("logical", 2),
                                    rep("character", 4))) %>%
    # Filter to driver uncertainty models with legacy covariate
    filter(scenario == "Legacy with covariate 2013-2018" &
        min.date == "20130601" &
        max.date == "20180101" &
        temporalDriverUncertainty == TRUE &
        spatialDriverUncertainty == TRUE &
        model_name %in% c("env_cycl", "env_cov", "cycl_only")
    ) %>%
    # CRITICAL FIX: Dirichlet models the whole rank.
    # Change species to a generic folder name and deduplicate tasks.
    mutate(species = "all_taxa") %>%
    distinct(rank.name, model_name, scenario, min.date, max.date, .keep_all = TRUE)

# LOCAL TESTING: Filter early to avoid loading 2,700+ models
if (identical(tolower(Sys.getenv("LOCAL_TEST", "false")), "true")) {
    info("🧪 LOCAL TESTING: Filtering to multiple taxa and model types for testing")
    params_in <- params_in %>% filter(
        # Dirichlet requires >= 3 taxa; use phylum_fun (multiple phyla as columns)
        rank.name == "phylum_fun" &
        # Test all model types
        model_name %in% c("env_cycl", "env_cov", "cycl_only") &
        scenario == "Legacy with covariate 2013-2018" &
        min.date == "20130601" &
        max.date == "20180101" &
        temporalDriverUncertainty == TRUE &
        spatialDriverUncertainty == TRUE
    )
    info("   Filtered to %d models for testing", nrow(params_in))
}

# Priority mode is now handled in CLI parsing above

if (use_priority == "true" || use_priority == "priority") {
    info("🎯 PRIORITY MODE: Using high-priority models with existing progress")
    priority_file <- here("data/summary/priority_rerun_list.rds")
    if (file.exists(priority_file)) {
        rerun_list <- readRDS(priority_file)
        info("   Loaded %d priority models (have chain 1 completed)", length(rerun_list))
    } else {
        warn("   Priority list not found! Falling back to standard unconverged list...")
        rerun_list <- readRDS(here("data/summary/unconverged_taxa_list.rds"))
    }
} else {
    info("📊 STANDARD MODE: Using all unconverged models")
    rerun_list <- readRDS(here("data/summary/unconverged_taxa_list.rds"))
    info("   Tip: Set USE_PRIORITY=true or add 'priority' argument to focus on models with existing progress")
}

converged_list = readRDS(here("data/summary/weak_converged_taxa_list.rds"))

# HPC PRODUCTION CONFIGURATION: Run multiple models with driver uncertainty enabled
# Note: params_in may already be filtered by LOCAL_TEST above

# LOCAL TESTING: For testing, sample from params_in directly and ensure all model types are tested
if (identical(tolower(Sys.getenv("LOCAL_TEST", "false")), "true")) {
    # In testing mode, use all available models regardless of convergence status
    info("🧪 LOCAL TESTING: Skipping converged filter for comprehensive testing")
    info("🧪 Available models before sampling: %d", nrow(params_in))
    set.seed(123)  # For reproducible sampling
    
    # Sample from params_in (already filtered to test taxa/model types)
    # Group by model_name and sample at least one from each model type
    if (nrow(params_in) > 0) {
        params_by_model <- params_in %>%
            group_by(model_name) %>%
            slice_sample(n = 1, replace = FALSE) %>%
            ungroup()
        
        info("🧪 After grouping by model_name and sampling: %d models", nrow(params_by_model))
        info("🧪 Selected models cover: %s", paste(sort(unique(params_by_model$model_name)), collapse=", "))
        
        params <- params_by_model
    } else {
        warn("🧪 No models available after filtering!")
        params <- params_in
    }
} else {
    # PRODUCTION: Apply full filtering
    params <- params_in %>% ungroup %>% filter(
        # Run ALL model types with driver uncertainty enabled
        model_name %in% c("env_cycl", "env_cov", "cycl_only") &
            # Only include models with driver uncertainty enabled
            temporalDriverUncertainty == TRUE &
            spatialDriverUncertainty == TRUE &
            # Focus on 2013-2018 period for legacy analysis
            scenario %in% c("Legacy with covariate 2013-2018")
    ) %>% distinct(.keep_all = TRUE)
    
    # PRODUCTION: Filter out converged models
    params <- params %>% filter(!model_id %in% converged_list)
}

info("LOCAL TESTING: Starting parallel execution for %d models with %d chains", nrow(params), nchains)
info("Expected runtime: Variable (convergence-based sampling)")
info("  - Models to run: %d models (cycl_only, env_cov, env_cycl) with driver uncertainty", nrow(params))
info("  - Chains per model: %d (total %d parallel tasks)", nchains, nrow(params) * nchains)
info("  - Initial iterations: ~ 0.1 minutes per chain")
info("  - Additional iterations: Variable based on convergence")
info("  - Target: ESS >= 10 per parameter")
info("🔧 TESTING MODE: Local testing with 2 cores (not HPC)")

info("LOCAL TESTING: Running %d models with %d chains in parallel", nrow(params), nchains)

# Filter parameters to only include models with available species and ranks
info("Filtering models to only include those with available data...")
original_n_models <- nrow(params)

# Create a function to check if a model has valid data
is_valid_model <- function(rank_name) {
    # Check if rank exists in data
    if (!(rank_name %in% names(all_ranks))) {
        return(FALSE)
    }
    return(TRUE)
}

# Filter the parameters dataframe
valid_indices <- sapply(1:nrow(params), function(i) {
    is_valid_model(params$rank.name[i])
})

params <- params[valid_indices, ]
filtered_n_models <- nrow(params)

# Report filtering results
if (filtered_n_models < original_n_models) {
    n_filtered <- original_n_models - filtered_n_models
    warn("Filtered out %d models with unavailable species/ranks", n_filtered)
    info("  Original models: %d", original_n_models)
    info("  Valid models: %d", filtered_n_models)
} else {
    info("✓ All %d models have valid data", filtered_n_models)
}

# Additional validation: ensure we have at least one valid model
if (filtered_n_models == 0) {
    error("No valid models remaining after filtering!")
    error("Available ranks in data: %s", paste(names(all_ranks), collapse=', '))
    
    # Show some examples of available species for each rank
    for (rank_name in names(all_ranks)) {
        rank_data <- all_ranks[[rank_name]]
        metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
        available_species <- setdiff(colnames(rank_data), metadata_cols)
        info("  %s species (first 5): %s%s", rank_name, paste(head(available_species, 5), collapse=", "),
            if(length(available_species) > 5) paste("... (+", length(available_species) - 5, " more)") else "")
    }
    
    stop("No valid models to run - check species names and rank names in parameters")
}

# Use the filtered parameters from data loading
valid_models <- params

# Apply CLI filters so both single-task and parallel runs respect --species / --model-name / etc.
if (!is.null(cli_model_id)) {
  valid_models <- valid_models %>% dplyr::filter(model_id == cli_model_id)
  info("Filtered to model_id=%s -> %d rows", cli_model_id, nrow(valid_models))
}
if (!is.null(cli_rank)) {
  valid_models <- valid_models %>% dplyr::filter(rank.name == cli_rank)
  info("Filtered to rank.name=%s -> %d rows", cli_rank, nrow(valid_models))
}
if (!is.null(cli_species)) {
  valid_models <- valid_models %>% dplyr::filter(species == cli_species)
  info("Filtered to species=%s -> %d rows", cli_species, nrow(valid_models))
}
if (!is.null(cli_model_name)) {
  valid_models <- valid_models %>% dplyr::filter(model_name == cli_model_name)
  info("Filtered to model_name=%s -> %d rows", cli_model_name, nrow(valid_models))
}
if (nrow(valid_models) == 0) {
  stop("No models remain after CLI filtering; check --species / --model-name / --rank / --model-id.")
}

# Local run: no job array (we drive parallel ourselves); limit to 1 model and 2 workers unless NCORES/LIMIT_ONE_MODEL set
local_run <- is.na(array_task_id)
if (local_run && nrow(valid_models) > 1L && (nzchar(Sys.getenv("LIMIT_ONE_MODEL", "")) || !nzchar(Sys.getenv("NCORES", "")))) {
  valid_models <- valid_models %>% slice(1L)
  info("Local run: limiting to 1 model (%d rows remaining)", nrow(valid_models))
}

# env_cycl_all_taxa: 2 chains; all others: global nchains
valid_models <- valid_models %>%
  mutate(nchains_model = as.integer(nchains))
all_tasks <- bind_rows(lapply(1:nrow(valid_models), function(i) {
  data.frame(model_idx = i, chain_no = 1:valid_models$nchains_model[i])
}))
n_model_tasks <- nrow(all_tasks)
info("After CLI filters: %d model(s), %d total chain(s) (env_cycl+all_taxa=2, others=%d)", nrow(valid_models), n_model_tasks, nchains)
info("Models to run:")
print(valid_models)

# Store single-task mode flag for later use
single_task_mode <- !is.na(array_task_id)

info("=== EXECUTION MODE DEBUGGING ===")
info("array_task_id=%s", ifelse(is.na(array_task_id), "NULL", as.character(array_task_id)))
info("single_task_mode=%s", as.character(single_task_mode))
info("SGE_TASK_ID=%s", Sys.getenv("SGE_TASK_ID", "NULL"))
info("PBS_ARRAYID=%s", Sys.getenv("PBS_ARRAYID", "NULL"))
info("SLURM_ARRAY_TASK_ID=%s", Sys.getenv("SLURM_ARRAY_TASK_ID", "NULL"))
info("Resolved here() root: %s", here::here())
info("model_output_dir: %s", model_output_dir)
info("=================================")

info("HPC PRODUCTION: Running %d models across all ranks with beta regression approach", nrow(valid_models))

# Configuration for parallel execution: default 2 workers locally unless NCORES set
if (identical(tolower(Sys.getenv("LOCAL_TEST","false")), "true")) {
  n_cores <- 2L
  info("🔧 TESTING MODE: Local testing with %d cores", n_cores)
} else {
  n_cores <- as.integer(Sys.getenv("NCORES", "2"))
  info("Workers: %d (set NCORES for more; default 2 for local runs)", n_cores)
}

info("📊 Models to run: %d models, %d total chains (env_cycl+all_taxa=2, others=%d)", nrow(valid_models), n_model_tasks, nchains)
info("⏱️  Expected runtime: Variable (convergence-based sampling)")
info("🎯 Target: ESS >= 10 per parameter")
#

create_nimble_model_with_uncertainty <- function(model_name, use_legacy_covariate = TRUE, 
                                                      temporalDriverUncertainty = TRUE, 
                                                      spatialDriverUncertainty = TRUE) {
    info("Building Nimble model with driver uncertainty: %s", model_name)
    info("  Temporal uncertainty: %s", temporalDriverUncertainty)
    info("  Spatial uncertainty: %s", spatialDriverUncertainty)
    
    if (model_name == "env_cycl" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # 1. PRIORS (Vectorized across N.spp)
            # --- Centered parameterization: site_effect absorbs the intercept ---
            for (s in 1:N.spp) {
                rho[s] ~ dunif(-0.99, 0.99)
                legacy_effect[s] ~ dnorm(0, sd = 10)

                # Centered site effects: intercept is the grand mean (replaces intercept)
                intercept[s] ~ dnorm(0, sd = 10)
                site_effect_sd[s] ~ dgamma(2, 20)
                for (k in 1:N.site) {
                    site_effect[k, s] ~ dnorm(intercept[s], sd = site_effect_sd[s])
                }

                # Multivariate normal prior on beta vector (pools information across covariates)
                beta[s, 1:8] ~ dmnorm(zeros[1:8], omega[1:8, 1:8])

                # Gamma concentration for process error (conjugate-friendly)
                sigma_proc[s] ~ dgamma(2, 0.5)   # concentration; higher = less process noise

                sigma_init[s] ~ dunif(0, 1)
                tau_init[s] <- pow(sigma_init[s], -2)
            }

            # 2. PROCESS MODEL - AR(1) on log-scale, Gamma-distributed alpha
            for (s in 1:N.spp) {
                for (p in 1:N.plot) {
                    # Initial timepoint: log-linear predictor -> Gamma alpha
                    log_mu[p, s, plot_start[p]] <-
                        site_effect[plot_site_num[p], s] +
                        beta[s, 1] * sin_mo[plot_start[p]] + beta[s, 2] * cos_mo[plot_start[p]] +
                        beta[s, 3] * temp_est[plot_site_num[p], plot_start[p]] +
                        beta[s, 4] * mois_est[plot_site_num[p], plot_start[p]] +
                        beta[s, 5] * pH_est[p, plot_start[p]] + beta[s, 6] * pC_est[p, plot_start[p]] +
                        beta[s, 7] * relEM[p, plot_start[p]] +
                        beta[s, 8] * LAI[plot_site_num[p], plot_start[p]] +
                        legacy_effect[s] * legacy[p, plot_start[p]]

                    # Gamma-distributed alpha (conjugate with Dirichlet)
                    # E[alpha] = exp(log_mu), Var[alpha] = exp(log_mu) / sigma_proc
                    alpha[p, s, plot_start[p]] ~ dgamma(
                        shape = exp(log_mu[p, s, plot_start[p]]) * sigma_proc[s],
                        rate = sigma_proc[s])

                    for (t in (plot_start[p] + 1):N.date) {
                        # AR(1) on log scale with covariates
                        log_mu[p, s, t] <-
                            rho[s] * log(max(alpha[p, s, t - 1], 0.001)) +
                            site_effect[plot_site_num[p], s] +
                            beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t] +
                            beta[s, 3] * temp_est[plot_site_num[p], t] +
                            beta[s, 4] * mois_est[plot_site_num[p], t] +
                            beta[s, 5] * pH_est[p, t] + beta[s, 6] * pC_est[p, t] +
                            beta[s, 7] * relEM[p, t] +
                            beta[s, 8] * LAI[plot_site_num[p], t] +
                            legacy_effect[s] * legacy[p, t]

                        alpha[p, s, t] ~ dgamma(
                            shape = exp(log_mu[p, s, t]) * sigma_proc[s],
                            rate = sigma_proc[s])
                    }
                }
            }

            # 3. RELATIVE ABUNDANCE
            for (p in 1:N.plot) {
                for (t in plot_start[p]:N.date) {
                    alpha_sum[p, t] <- sum(alpha[p, 1:N.spp, t])
                    for (s in 1:N.spp) {
                        plot_rel[p, s, t] <- alpha[p, s, t] / alpha_sum[p, t]
                    }
                }
            }

            # 4. DRIVER UNCERTAINTY - Temporal
            if (temporalDriverUncertainty) {
                for (k in 1:N.site) {
                    for (t in site_start[k]:N.date) {
                        mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
                        temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
                    }
                }
            } else {
                for (k in 1:N.site) {
                    for (t in site_start[k]:N.date) {
                        mois_est[k, t] <- mois[k, t]
                        temp_est[k, t] <- temp[k, t]
                    }
                }
            }

            # DRIVER UNCERTAINTY - Spatial
            if (spatialDriverUncertainty) {
                for (p in 1:N.plot) {
                    pH_est_p[p] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
                    pC_est_p[p] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] <- pH_est_p[p]
                        pC_est[p, t] <- pC_est_p[p]
                    }
                }
            } else {
                for (p in 1:N.plot) {
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] <- pH[p, t]
                        pC_est[p, t] <- pC[p, t]
                    }
                }
            }

            # 5. OBSERVATION MODEL (Dirichlet)
            for (i in 1:N.core) {
                y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
            }
        })
    } else if (model_name == "env_cov" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({

            # 1. PRIORS (Vectorized across N.spp)
            for (s in 1:N.spp) {
                rho[s] ~ dunif(-0.99, 0.99)
                intercept[s] ~ dnorm(0, sd = 10)
                legacy_effect[s] ~ dnorm(0, sd = 10)

                site_effect_sd[s] ~ dgamma(2, 20)
                for (k in 1:N.site) {
                    site_effect[k, s] ~ dnorm(0, sd = site_effect_sd[s])
                }

                for (b in 1:6) {
                    beta[s, b] ~ dnorm(0, sd = 10)
                }

                sigma_proc[s] ~ dunif(0, 0.5)
                tau_proc[s] <- pow(sigma_proc[s], -2)
                sigma_init[s] ~ dunif(0, 1)
                tau_init[s] <- pow(sigma_init[s], -2)
            }

            # 2. PROCESS MODEL - AR(1) on eta, map to alpha
            for (s in 1:N.spp) {
                for (p in 1:N.plot) {
                    # Initial state
                    eta[p, s, plot_start[p]] ~ dnorm(
                        intercept[s] +
                        site_effect[plot_site_num[p], s] +
                        beta[s, 1] * temp_est[plot_site_num[p], plot_start[p]] +
                        beta[s, 2] * mois_est[plot_site_num[p], plot_start[p]] +
                        beta[s, 3] * pH_est[p, plot_start[p]] + beta[s, 4] * pC_est[p, plot_start[p]] +
                        beta[s, 5] * relEM[p, plot_start[p]] +
                        beta[s, 6] * LAI[plot_site_num[p], plot_start[p]] +
                        legacy_effect[s] * legacy[p, plot_start[p]],
                        tau_init[s]
                    )
                    alpha[p, s, plot_start[p]] <- exp(eta[p, s, plot_start[p]])

                    # DIRECT AR(1) on eta
                    for (t in (plot_start[p] + 1):N.date) {
                        eta[p, s, t] ~ dnorm(
                            rho[s] * eta[p, s, t - 1] +
                            intercept[s] +
                            site_effect[plot_site_num[p], s] +
                            beta[s, 1] * temp_est[plot_site_num[p], t] +
                            beta[s, 2] * mois_est[plot_site_num[p], t] +
                            beta[s, 3] * pH_est[p, t] + beta[s, 4] * pC_est[p, t] +
                            beta[s, 5] * relEM[p, t] +
                            beta[s, 6] * LAI[plot_site_num[p], t] +
                            legacy_effect[s] * legacy[p, t],
                            tau_proc[s]
                        )
                        alpha[p, s, t] <- exp(eta[p, s, t])
                    }
                }
            }

            # 3. RELATIVE ABUNDANCE (Deterministic extraction for MCMC monitoring)
            for (p in 1:N.plot) {
                for (t in plot_start[p]:N.date) {
                    alpha_sum[p, t] <- sum(alpha[p, 1:N.spp, t])
                    for (s in 1:N.spp) {
                        plot_rel[p, s, t] <- alpha[p, s, t] / alpha_sum[p, t]
                    }
                }
            }

            # 4. DRIVER UNCERTAINTY - Temporal (temperature and moisture)
            if (temporalDriverUncertainty) {
                for (k in 1:N.site) {
                    for (t in site_start[k]:N.date) {
                        mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
                        temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
                    }
                }
            } else {
                for (k in 1:N.site) {
                    for (t in site_start[k]:N.date) {
                        mois_est[k, t] <- mois[k, t]
                        temp_est[k, t] <- temp[k, t]
                    }
                }
            }

            # DRIVER UNCERTAINTY - Spatial (pH and pC - constant over time)
            if (spatialDriverUncertainty) {
                for (p in 1:N.plot) {
                    pH_est_p[p] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
                    pC_est_p[p] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] <- pH_est_p[p]
                        pC_est[p, t] <- pC_est_p[p]
                    }
                }
            } else {
                for (p in 1:N.plot) {
                    for (t in plot_start[p]:N.date) {
                        pH_est[p, t] <- pH[p, t]
                        pC_est[p, t] <- pC[p, t]
                    }
                }
            }

            # 5. OBSERVATION MODEL (Dirichlet)
            for (i in 1:N.core) {
                y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
            }
        })
    } else if (model_name == "cycl_only") {
        modelCode <- nimble::nimbleCode({
            # 1. PRIORS (Vectorized across N.spp)
            for (s in 1:N.spp) {
                rho[s] ~ dunif(-0.99, 0.99)
                intercept[s] ~ dnorm(0, sd = 10)

                site_effect_sd[s] ~ dgamma(2, 20)
                for (k in 1:N.site) {
                    site_effect[k, s] ~ dnorm(0, sd = site_effect_sd[s])
                }

                for (b in 1:2) {
                    beta[s, b] ~ dnorm(0, sd = 10)
                }

                sigma_proc[s] ~ dunif(0, 0.5)
                tau_proc[s] <- pow(sigma_proc[s], -2)
                sigma_init[s] ~ dunif(0, 1)
                tau_init[s] <- pow(sigma_init[s], -2)
            }

            # 2. PROCESS MODEL - AR(1) on eta, map to alpha (seasonal only)
            for (s in 1:N.spp) {
                for (p in 1:N.plot) {
                    eta[p, s, plot_start[p]] ~ dnorm(
                        intercept[s] +
                        site_effect[plot_site_num[p], s] +
                        beta[s, 1] * sin_mo[plot_start[p]] + beta[s, 2] * cos_mo[plot_start[p]],
                        tau_init[s]
                    )
                    alpha[p, s, plot_start[p]] <- exp(eta[p, s, plot_start[p]])

                    for (t in (plot_start[p] + 1):N.date) {
                        eta[p, s, t] ~ dnorm(
                            rho[s] * eta[p, s, t - 1] +
                            intercept[s] +
                            site_effect[plot_site_num[p], s] +
                            beta[s, 1] * sin_mo[t] + beta[s, 2] * cos_mo[t],
                            tau_proc[s]
                        )
                        alpha[p, s, t] <- exp(eta[p, s, t])
                    }
                }
            }

            # 3. RELATIVE ABUNDANCE
            for (p in 1:N.plot) {
                for (t in plot_start[p]:N.date) {
                    alpha_sum[p, t] <- sum(alpha[p, 1:N.spp, t])
                    for (s in 1:N.spp) {
                        plot_rel[p, s, t] <- alpha[p, s, t] / alpha_sum[p, t]
                    }
                }
            }

            # 4. OBSERVATION MODEL (Dirichlet)
            for (i in 1:N.core) {
                y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
            }
        })
    } else {
        stop("Unsupported model combination: ", model_name, " with use_legacy_covariate=", use_legacy_covariate)
    }
    
    return(modelCode)
}

# Function to sanitize driver uncertainty data
sanitize_driver_uncertainty <- function(constants) {
  info("Sanitizing driver uncertainty data (STRICT, window-aware)")

  eps <- 1e-6

  ## helpers
  .summarize_bad <- function(mask, max_show = 10) {
    rc <- which(mask, arr.ind = TRUE)
    if (length(rc) == 0) return("0")
    head_pairs <- apply(head(rc, max_show), 1, \(r) paste0("(", r[1], ",", r[2], ")"))
    paste0(nrow(rc), " cells; e.g., ", paste(head_pairs, collapse = ", "),
           if (nrow(rc) > max_show) ", ..." else "")
  }
  .in_scope_site <- function(Nsite, Ndate, site_start) {
    m <- matrix(FALSE, Nsite, Ndate)
    for (k in seq_len(Nsite)) {
      if (is.finite(site_start[k]) && site_start[k] >= 1 && site_start[k] <= Ndate) {
        m[k, site_start[k]:Ndate] <- TRUE
      }
    }
    m
  }
  .in_scope_plot <- function(Nplot, Ndate, plot_start) {
    m <- matrix(FALSE, Nplot, Ndate)
    for (p in seq_len(Nplot)) {
      if (is.finite(plot_start[p]) && plot_start[p] >= 1 && plot_start[p] <= Ndate) {
        m[p, plot_start[p]:Ndate] <- TRUE
      }
    }
    m
  }

  ## ---------- pC / pC_sd (SPATIAL; should be constant over time) ----------
  bad_pC    <- !is.finite(constants$pC)
    bad_pC_sd <- !is.finite(constants$pC_sd) | (constants$pC_sd <= 0)
    
  # Repair (global medians / epsilon for SD)
    pC_fill <- median(constants$pC[is.finite(constants$pC)], na.rm = TRUE)
    constants$pC[bad_pC] <- pC_fill
  pC_sd_fill <- median(constants$pC_sd[is.finite(constants$pC_sd) & constants$pC_sd > 0], na.rm = TRUE)
  if (!is.finite(pC_sd_fill) || pC_sd_fill <= 0) pC_sd_fill <- eps
    constants$pC_sd[bad_pC_sd] <- pC_sd_fill
    
  # has_pC_plot evaluated at each plot's own start time
  idx_pc <- cbind(seq_len(constants$N.plot), pmax(1, pmin(constants$plot_start, constants$N.date)))
  constants$has_pC_plot <- is.finite(constants$pC[idx_pc]) &
                           is.finite(constants$pC_sd[idx_pc]) &
                           (constants$pC_sd[idx_pc] > 0)

  info("  pC repaired: mean cells=%d, sd cells=%d | valid pC plots at start=%d",
       sum(bad_pC), sum(bad_pC_sd), sum(constants$has_pC_plot))

  ## ---------- pH / pH_sd (SPATIAL; allowed to repair) ----------
  bad_pH    <- !is.finite(constants$pH)
    bad_pH_sd <- !is.finite(constants$pH_sd) | (constants$pH_sd <= 0)
    
    pH_fill <- median(constants$pH[is.finite(constants$pH)], na.rm = TRUE)
    constants$pH[bad_pH] <- pH_fill
  pH_sd_fill <- median(constants$pH_sd[is.finite(constants$pH_sd) & constants$pH_sd > 0], na.rm = TRUE)
  if (!is.finite(pH_sd_fill) || pH_sd_fill <= 0) pH_sd_fill <- eps
    constants$pH_sd[bad_pH_sd] <- pH_sd_fill
    
  idx_ph <- cbind(seq_len(constants$N.plot), pmax(1, pmin(constants$plot_start, constants$N.date)))
  constants$has_pH_plot <- is.finite(constants$pH[idx_ph]) &
                           is.finite(constants$pH_sd[idx_ph]) &
                           (constants$pH_sd[idx_ph] > 0)

  info("  pH repaired: mean cells=%d, sd cells=%d | valid pH plots at start=%d",
       sum(bad_pH), sum(bad_pH_sd), sum(constants$has_pH_plot))

  ## ---------- temp / temp_sd (TEMPORAL; strict inside site windows only) ----------
  if (!is.null(constants$temp)) {
    in_scope <- .in_scope_site(constants$N.site, constants$N.date, constants$site_start)

    bad_temp_all <- !is.finite(constants$temp)
    bad_temp_sd_all <- !is.finite(constants$temp_sd) | (constants$temp_sd <= 0)

    bad_temp_in  <- bad_temp_all  & in_scope
    bad_tsd_in   <- bad_temp_sd_all & in_scope

    denom <- max(1L, sum(in_scope))
    bad_temp_pct <- 100 * sum(bad_temp_in) / denom
    bad_tsd_pct  <- 100 * sum(bad_tsd_in)  / denom

    if (any(bad_temp_in)) {
      stop(sprintf("sanitize_driver_uncertainty: temp has non-finite values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_temp_pct, .summarize_bad(bad_temp_in)))
    }
    if (any(bad_tsd_in)) {
      stop(sprintf("sanitize_driver_uncertainty: temp_sd has non-finite or <=0 values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_tsd_pct, .summarize_bad(bad_tsd_in)))
    }

    info("  temp: 100%% valid within site windows | checked cells=%d | out-of-window ignored=%d",
         sum(in_scope), sum(!in_scope & (bad_temp_all | bad_temp_sd_all)))
    constants$has_temp <- TRUE
  }

  ## ---------- mois / mois_sd (TEMPORAL; strict inside site windows only) ----------
  if (!is.null(constants$mois)) {
    in_scope <- .in_scope_site(constants$N.site, constants$N.date, constants$site_start)

    bad_mois_all <- !is.finite(constants$mois)
    bad_mois_sd_all <- !is.finite(constants$mois_sd) | (constants$mois_sd <= 0)

    bad_mois_in <- bad_mois_all & in_scope
    bad_msd_in  <- bad_mois_sd_all & in_scope

    denom <- max(1L, sum(in_scope))
    bad_mois_pct <- 100 * sum(bad_mois_in) / denom
    bad_msd_pct  <- 100 * sum(bad_msd_in)  / denom

    if (any(bad_mois_in)) {
      stop(sprintf("sanitize_driver_uncertainty: mois has non-finite values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_mois_pct, .summarize_bad(bad_mois_in)))
    }
    if (any(bad_msd_in)) {
      stop(sprintf("sanitize_driver_uncertainty: mois_sd has non-finite or <=0 values within site windows (%.2f%%; %s). No filling allowed.",
                   bad_msd_pct, .summarize_bad(bad_msd_in)))
    }

    info("  mois: 100%% valid within site windows | checked cells=%d | out-of-window ignored=%d",
         sum(in_scope), sum(!in_scope & (bad_mois_all | bad_mois_sd_all)))
    constants$has_mois <- TRUE
  }

  info("✓ Driver uncertainty data sanitized (STRICT, window-aware)")
  constants
}

# Assert helpers
assert_matrix_dims <- function(x, nr, nc, name) {
  if (is.null(dim(x))) {
    stop(sprintf("Expected matrix for %s, got vector of length %d", name, length(x)))
  }
  if (!identical(dim(x), c(nr, nc))) {
    stop(sprintf("Dimension mismatch for %s: got %dx%d; expected %dx%d",
                 name, nrow(x), ncol(x), nr, nc))
  }
  invisible(TRUE)
}

assert_vector_len <- function(x, n, name) {
  if (length(x) != n)
    stop(sprintf("Length mismatch for %s: got %d; expected %d", name, length(x), n))
  invisible(TRUE)
}

create_inits_with_uncertainty <- function(constants, model_name, model_data = NULL, chain_no = 1) {
    info("Creating initial values with driver uncertainty for %s...", model_name)

    # Determine number of beta parameters based on model type
    if (model_name == "env_cycl") {
        n_beta <- 8
    } else if (model_name == "env_cov") {
        n_beta <- 6
    } else {
        n_beta <- 2  # cycl_only
    }

    # Check for warm-start file
    warmstart_file <- Sys.getenv("WARMSTART_FILE", "")
    if (warmstart_file != "" && file.exists(warmstart_file)) {
        info("  Loading warm-start initial values from: %s", warmstart_file)
        ws <- readRDS(warmstart_file)
        pm <- ws$param_means

        # Map flat MCMC parameter names back into structured inits
        # Add small per-chain jitter to break symmetry
        set.seed(chain_no * 42)
        jitter_sd <- 0.05

        inits <- list(
            rho = sapply(1:constants$N.spp, function(s) pm[paste0("rho[", s, "]")] + rnorm(1, 0, jitter_sd)),
            beta = matrix(sapply(1:n_beta, function(b) sapply(1:constants$N.spp, function(s)
                pm[paste0("beta[", s, ", ", b, "]")] + rnorm(1, 0, jitter_sd))),
                nrow = constants$N.spp, ncol = n_beta),
            site_effect_sd = sapply(1:constants$N.spp, function(s)
                abs(pm[paste0("site_effect_sd[", s, "]")] + rnorm(1, 0, jitter_sd * 0.5))),
            site_effect = matrix(sapply(1:constants$N.spp, function(s) sapply(1:constants$N.site, function(k)
                pm[paste0("site_effect[", k, ", ", s, "]")] + rnorm(1, 0, jitter_sd))),
                nrow = constants$N.site, ncol = constants$N.spp),
            intercept = sapply(1:constants$N.spp, function(s) pm[paste0("intercept[", s, "]")] + rnorm(1, 0, jitter_sd)),
            sigma_proc = sapply(1:constants$N.spp, function(s) {
                v <- pm[paste0("sigma_proc[", s, "]")]
                if (is.na(v)) v <- pm[paste0("sigma[", s, ", 1]")]  # try alternative naming
                abs(v + rnorm(1, 0, jitter_sd * 0.2))
            }),
            sigma_init = sapply(1:constants$N.spp, function(s) {
                v <- pm[paste0("sigma_init[", s, "]")]
                if (is.na(v)) v <- pm[paste0("sigma[", s, ", 2]")]
                abs(v + rnorm(1, 0, jitter_sd * 0.2))
            }),
            eta = array(-1, dim = c(constants$N.plot, constants$N.spp, constants$N.date))
        )
        if (model_name %in% c("env_cycl", "env_cov")) {
            inits$legacy_effect <- sapply(1:constants$N.spp, function(s)
                pm[paste0("legacy_effect[", s, "]")] + rnorm(1, 0, jitter_sd))
        }
        info("  Warm-start loaded from %d-iteration run, jitter_sd=%.3f for chain %d", ws$source_iterations, jitter_sd, chain_no)
    } else {
        # Default cold-start initial values (reparameterized model)
        inits <- list(
            rho = rep(0.3, constants$N.spp),
            beta = matrix(0.01, nrow = constants$N.spp, ncol = n_beta),
            site_effect_sd = rep(0.5, constants$N.spp),
            intercept = rep(-1, constants$N.spp),  # centered parameterization (replaces intercept)
            site_effect = matrix(rnorm(constants$N.site * constants$N.spp, -1, 0.1),
                                 nrow = constants$N.site, ncol = constants$N.spp),
            sigma_proc = rep(4, constants$N.spp),  # Gamma concentration (replaces sigma_proc)
            sigma_init = rep(0.5, constants$N.spp),
            alpha = array(1, dim = c(constants$N.plot, constants$N.spp, constants$N.date))
        )
        if (model_name %in% c("env_cycl", "env_cov")) {
            inits$legacy_effect <- rep(0, constants$N.spp)
        }
    }

    info("  ✓ Dirichlet inits (N.spp = %d, n_beta = %d)", constants$N.spp, n_beta)
    
    # Initialize driver uncertainty variables ONLY for environmental models
    if (model_name %in% c("env_cycl", "env_cov")) {
        if (constants$N.site > 0 && constants$N.date > 0) {
            # Initialize near observed values instead of random values
            inits$temp_est <- constants$temp  # Start at observed temperature
            inits$mois_est <- constants$mois  # Start at observed moisture
            info("  ✓ Initialized temp_est and mois_est at observed values")
        }
        
        if (constants$N.plot > 0 && constants$N.date > 0) {
            # Initialize spatial pH/pC variables (one per plot)
            inits$pH_est_p <- constants$pH[cbind(1:constants$N.plot, constants$plot_start)]  # Start at observed pH at plot start
            inits$pC_est_p <- constants$pC[cbind(1:constants$N.plot, constants$plot_start)]  # Start at observed pC at plot start
            # Note: pH_est and pC_est are deterministic from pH_est_p and pC_est_p, so no need to initialize them
            info("  ✓ Initialized pH_est_p, pC_est_p (spatial) at observed values")
        }
    } else {
        info("  ✓ Skipped driver uncertainty variables for %s model (not needed)", model_name)
    }
    
    info("  ✓ Initial values created successfully for %s model", model_name)
    info("    rho: length %d", length(inits$rho))
    info("    beta: %s", paste(dim(inits$beta), collapse = " x "))
    info("    site_effect: %s", paste(dim(inits$site_effect), collapse = " x "))
    info("    eta array dimensions: %s", paste(dim(inits$eta), collapse = " x "))
    
    # Only report driver uncertainty dimensions if they exist
    if ("temp_est" %in% names(inits)) {
        info("    temp_est matrix dimensions: %s", paste(dim(inits$temp_est), collapse = " x "))
    }
    if ("mois_est" %in% names(inits)) {
        info("    mois_est matrix dimensions: %s", paste(dim(inits$mois_est), collapse = " x "))
    }
    if ("pH_est" %in% names(inits)) {
        info("    pH_est matrix dimensions: %s", paste(dim(inits$pH_est), collapse = " x "))
    }
    if ("pC_est" %in% names(inits)) {
        info("    pC_est matrix dimensions: %s", paste(dim(inits$pC_est), collapse = " x "))
    }
    
    return(inits)
}

# Function to create visualization of plot_mu over time
create_plot_mu_visualization <- function(samples, samples2, model_data, metadata, output_dir) {
    
    # Use actual column names from samples
    cn <- colnames(samples)
    param_names <- cn
    
    # Create a basic plot of parameter estimates
    png(file.path(output_dir, "parameter_estimates.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    # Plot key parameters
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(samples[, i], type = "l", main = param_names[i], 
                 xlab = "Iteration", ylab = "Value")
            abline(h = mean(samples[, i]), col = "red", lty = 2)
        }
    }
    
    dev.off()
    
    # Create trace plots
    png(file.path(output_dir, "trace_plots.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(samples[, i], type = "l", main = paste("Trace:", param_names[i]), 
                 xlab = "Iteration", ylab = "Value")
        }
    }
    
    dev.off()
    
    # Create density plots
    png(file.path(output_dir, "density_plots.png"), width = 1200, height = 800)
    par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
    
    for (i in 1:min(8, ncol(samples))) {
        if (i <= ncol(samples)) {
            plot(density(samples[, i]), main = paste("Density:", param_names[i]), 
                 xlab = "Value", ylab = "Density")
            abline(v = mean(samples[, i]), col = "red", lty = 2)
        }
    }
    
    dev.off()
    
    # Create mu over time visualization if samples2 is available
    if (!is.null(samples2) && ncol(samples2) > 0) {
        
        # Extract mu samples and compute summaries
        mu_samples <- samples2
        n_plots <- min(5, sqrt(ncol(mu_samples)))  # Show up to 5 plots
        
        png(file.path(output_dir, "mu_over_time.png"), width = 1200, height = 800)
        par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
        
        for (p in 1:n_plots) {
            # Find mu samples for this plot
            plot_cols <- grep(paste0("mu\\[", p, ","), colnames(mu_samples))
            if (length(plot_cols) > 0) {
                plot_mu <- mu_samples[, plot_cols, drop = FALSE]
                
                # Handle NaN values by removing them
                plot_mu_clean <- plot_mu[complete.cases(plot_mu), , drop = FALSE]
                
                if (nrow(plot_mu_clean) > 0) {
                    # Compute posterior mean and CI for each time point (with na.rm = TRUE)
                    mu_mean <- apply(plot_mu_clean, 2, mean, na.rm = TRUE)
                    mu_ci_lower <- apply(plot_mu_clean, 2, quantile, 0.025, na.rm = TRUE)
                    mu_ci_upper <- apply(plot_mu_clean, 2, quantile, 0.975, na.rm = TRUE)
                    
                    # Check for valid values
                    valid_idx <- !is.na(mu_mean) & !is.na(mu_ci_lower) & !is.na(mu_ci_upper)
                    
                    if (sum(valid_idx) > 0) {
                        # Plot only valid time points
                        time_points <- which(valid_idx)
                        plot(time_points, mu_mean[valid_idx], type = "l", 
                             main = paste("Plot", p, "μ over time"),
                             xlab = "Time", ylab = "μ", ylim = c(0, 1))
                        polygon(c(time_points, rev(time_points)), 
                                c(mu_ci_upper[valid_idx], rev(mu_ci_lower[valid_idx])), 
                                col = "lightblue", border = NA)
                        lines(time_points, mu_mean[valid_idx], col = "blue", lwd = 2)
                    } else {
                        plot(1, 1, type = "n", main = paste("Plot", p, "μ over time (no valid data)"),
                             xlab = "Time", ylab = "μ", ylim = c(0, 1))
                        text(1, 0.5, "No valid data", cex = 1.2)
                    }
                } else {
                    plot(1, 1, type = "n", main = paste("Plot", p, "μ over time (no data)"),
                         xlab = "Time", ylab = "μ", ylim = c(0, 1))
                    text(1, 0.5, "No data", cex = 1.2)
                }
            }
        }
        
        dev.off()
    }
    
    # Create summary statistics using actual column names
    summary_stats <- data.frame(
        Parameter = param_names,
        Mean = apply(samples, 2, mean),
        SD = apply(samples, 2, sd),
        Q025 = apply(samples, 2, quantile, 0.025),
        Q975 = apply(samples, 2, quantile, 0.975),
        ESS = coda::effectiveSize(coda::as.mcmc(samples))
    )
    
    write.csv(summary_stats, file.path(output_dir, "parameter_summary.csv"), row.names = FALSE)
    
    return(summary_stats)
}

# Create function that uses our working approach for each model
run_scenarios_fixed <- function(j, chain_no) {
    # Initialize error tracking and logging
    start_time <- Sys.time()
    error_context <- list()
    
    tryCatch({
        # Load required libraries in each worker using helper function
        load_required_packages()
        info("=== Starting model fitting ===")
        info("Model index: %d Chain: %d", j, chain_no)
        info("Model parameters:")
        print(valid_models[j,])
        info("=============================")
        
        # Debug HPC environment information
        info("=============================")
        
        # Validate input parameters
        if (is.null(valid_models) || nrow(valid_models) < j) {
            stop("Valid_models data frame not available or index out of bounds")
        }
        
        # Extract model parameters
        rank.name <- valid_models$rank.name[[j]]
        species <- valid_models$species[[j]]
        model_id <- valid_models$model_id[[j]]
        model_name <- valid_models$model_name[[j]]
        min.date <- valid_models$min.date[[j]]
        max.date <- valid_models$max.date[[j]]
        scenario <- valid_models$scenario[[j]]
        
        # Validate extracted parameters
        if (is.null(rank.name) || is.na(rank.name) || rank.name == "") {
            stop("Invalid rank.name for model index ", j)
        }
        if (is.null(species) || is.na(species) || species == "") {
            stop("Invalid species for model index ", j)
        }
        if (is.null(model_name) || is.na(model_name) || model_name == "") {
            stop("Invalid model_name for model index ", j)
        }
        
        # Check if this is a legacy covariate model
        use_legacy_covariate <- grepl("Legacy with covariate", scenario)
        
        # Validate data availability
        if (!exists("all_ranks") || is.null(all_ranks)) {
            stop("Data 'all_ranks' not available in worker environment")
        }
        
        # Get the specific group data
        if (!(rank.name %in% names(all_ranks))) {
            stop("Rank name '", rank.name, "' not found in data. Available ranks: ", 
                 paste(names(all_ranks), collapse=", "))
        }
        rank.df <- all_ranks[[rank.name]]
        
        # Validate rank data structure
        if (!is.data.frame(rank.df) || nrow(rank.df) == 0) {
            stop("Rank data for '", rank.name, "' is empty or not a data frame")
        }
        
        info("Preparing model data for %s", rank.name)
        
        # Apply LOCAL_TEST gating to the full matrix
        if (identical(tolower(Sys.getenv("LOCAL_TEST", "false")), "true")) {
            rank.df <- rank.df %>% head(500)
            info("  ✓ LOCAL_TEST mode: Using first 500 rows for faster testing")
        }
        
        # Use prepDirichletData with the full rank data (multi-species)
        info("Calling prepDirichletData...")
        model.dat <- microbialForecast::prepDirichletData(rank.df = rank.df,
                                       min.prev = 3,
                                       min.date = min.date,
                                       max.date = max.date)
        
        info("Data prepared successfully")
        
        # Debug: Check data structure immediately after preparation
        info("DEBUG: Data structure check:")
        info("  model.dat$y dimensions: %s", paste(dim(model.dat$y), collapse = " x "))
        info("  model.dat$y class: %s", class(model.dat$y))
        info("  model.dat$y column names: %s", if(is.matrix(model.dat$y)) paste(colnames(model.dat$y), collapse=", ") else "N/A")
        if (is.matrix(model.dat$y)) {
            info("  model.dat$y first few rows:")
            print(head(model.dat$y, 3))
        }
        info("  N.core calculated: %d", nrow(model.dat$y))
        info("  N.spp calculated: %d", ncol(model.dat$y))
        info("  Model type: %s", model_name)
        
        # Dirichlet model requires composition data (N.spp columns, rows sum to 1)
        if (ncol(model.dat$y) < 1) {
            error("  ❌ ERROR: Dirichlet model requires at least 1 species column in response data")
            stop("Response data must have N.spp >= 1 for Dirichlet model")
        }
        
        # Validate model data structure
        if (!is.list(model.dat) || !("y" %in% names(model.dat))) {
            stop("Invalid model data structure from prepDirichletData")
        }
        if (!is.matrix(model.dat$y) || nrow(model.dat$y) == 0) {
            stop("Model data 'y' is not a valid matrix or is empty")
        }
        
        # Prepare constants with validation
        info("Preparing model constants...")
        required_constants <- c("plotID", "timepoint", "plot_site", "site_start", "plot_start", 
                                "plot_index", "plot_num", "plot_site_num", "N.plot", "N.spp", 
                                "N.core", "N.site", "N.date", "sin_mo", "cos_mo")
        
        missing_constants <- required_constants[!required_constants %in% names(model.dat)]
        if (length(missing_constants) > 0) {
            stop("Missing required constants in model data: ", paste(missing_constants, collapse=", "))
        }
        
        constants <- model.dat[required_constants]
        
        # Add environmental predictors with validation (ONLY for environmental models)
        if (model_name %in% c("env_cycl", "env_cov")) {
            env_predictors <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
            info("Adding environmental predictors with validation...")
            
            for (pred in env_predictors) {
                if (pred %in% names(model.dat)) {
                    constants[[pred]] <- model.dat[[pred]]
                    
                    # Validate predictor dimensions and structure
                    pred_data <- model.dat[[pred]]
                    if (is.matrix(pred_data)) {
                        info("  ✓ Added %s predictor: %s matrix", pred, paste(dim(pred_data), collapse = " x "))
                    } else if (is.vector(pred_data)) {
                        info("  ✓ Added %s predictor: %d vector", pred, length(pred_data))
                    } else {
                        info("  ✓ Added %s predictor: %s object", pred, class(pred_data))
                    }
                    
                    # Check for missing or extreme values
                    if (is.numeric(pred_data)) {
                        missing_pct <- mean(is.na(pred_data)) * 100
                        if (missing_pct > 0) {
                            warn("    WARNING: %s has %.1f%% missing values", pred, missing_pct)
                        }
                        
                        if (is.matrix(pred_data)) {
                            extreme_vals <- sum(abs(pred_data) > 10, na.rm = TRUE)
                            if (extreme_vals > 0) {
                                warn("    WARNING: %s has %d extreme values (>10)", pred, extreme_vals)
                            }
                        }
                    }
                } else {
                    error("  ❌ ERROR: %s predictor not found in model data", pred)
                    stop("Missing required environmental predictor: ", pred)
                }
            }
        } else {
            info("Skipping environmental predictors for %s model", model_name)
        }
        
        # Add explicit dimension checks for all predictors
        if (model_name %in% c("env_cycl", "env_cov")) {
            info("Validating predictor dimensions...")
            assert_matrix_dims(constants$temp,  constants$N.site, constants$N.date, "temp")
            assert_matrix_dims(constants$mois,  constants$N.site, constants$N.date, "mois")
            assert_matrix_dims(constants$LAI,   constants$N.site, constants$N.date, "LAI")
            assert_matrix_dims(constants$pH,    constants$N.plot, constants$N.date, "pH")
            assert_matrix_dims(constants$pC,    constants$N.plot, constants$N.date, "pC")
            assert_matrix_dims(constants$relEM, constants$N.plot, constants$N.date, "relEM")
            info("  ✓ All predictor dimensions validated")
        }
        
        # Add driver uncertainty flags to constants for Nimble model evaluation
        info("Adding driver uncertainty flags to constants...")
        constants$temporalDriverUncertainty <- driver_uncertainty_mode   # Driver uncertainty - temporal uncertainty
        constants$spatialDriverUncertainty <- driver_uncertainty_mode    # Driver uncertainty - spatial uncertainty
        info("  ✓ temporalDriverUncertainty = %s", constants$temporalDriverUncertainty)
        info("  ✓ spatialDriverUncertainty = %s", constants$spatialDriverUncertainty)
        
        # Add driver uncertainty data when flags are enabled
        if (model_name %in% c("env_cycl", "env_cov") &&
            (constants$temporalDriverUncertainty || constants$spatialDriverUncertainty)) {
            info("Adding driver uncertainty data with proper dimension matching...")
            
            # CRITICAL FIX: Extract perfectly aligned arrays directly from model.dat
            if (constants$temporalDriverUncertainty) {
                constants$temp_sd <- model.dat$temp_sd
                constants$mois_sd <- model.dat$mois_sd
                
                assert_matrix_dims(constants$temp_sd, constants$N.site, constants$N.date, "temp_sd")
                assert_matrix_dims(constants$mois_sd, constants$N.site, constants$N.date, "mois_sd")
                info("  ✓ Added temp_sd and mois_sd uncertainty data")
            }
            
            if (constants$spatialDriverUncertainty) {
                constants$pH_sd <- model.dat$pH_sd
                constants$pC_sd <- model.dat$pC_sd
                
                assert_matrix_dims(constants$pH_sd, constants$N.plot, constants$N.date, "pH_sd")
                assert_matrix_dims(constants$pC_sd, constants$N.plot, constants$N.date, "pC_sd")
                info("  ✓ Added pH_sd and pC_sd uncertainty data")
            }
            
            # Strict sanitize: only pH/pC repaired; temp/mois must be clean
            constants <- sanitize_driver_uncertainty(constants)
            
            info("  ✓ Driver uncertainty data added & sanitized")
        } else {
            info("Skipping driver uncertainty data for %s model", model_name)
        }
        
        # Add legacy covariate if needed
        if (use_legacy_covariate) {
            info("Adding legacy covariate using actual dates...")
            legacy_cutoff <- as.Date("2015-01-01")

            # derive a per-timepoint legacy vector of length N.date
            # model.dat should carry a unique ordered list of dates; if not, derive from input df
            if (!"all_dates" %in% names(model.dat)) {
                # fallback from the original df
                all_dates <- sort(unique(as.Date(rank.df$dates)))
                } else {
                all_dates <- as.Date(model.dat$all_dates)
            }

            if (length(all_dates) != constants$N.date) {
                warn("Time axis length mismatch (got %d dates, expected N.date=%d). Rebuilding by index.",
                     length(all_dates), constants$N.date)
                all_dates <- seq_len(constants$N.date) # still deterministic
                legacy_by_time <- all_dates <= floor(0.6 * constants$N.date)
            } else {
                legacy_by_time <- as.numeric(all_dates < legacy_cutoff)
            }

            # expand to plot x time while respecting plot_start windows
            legacy_mat <- matrix(0, nrow = constants$N.plot, ncol = constants$N.date)
            for (p in seq_len(constants$N.plot)) {
                ts <- max(1L, min(constants$plot_start[p], constants$N.date))
                legacy_mat[p, ts:constants$N.date] <- legacy_by_time[ts:constants$N.date]
            }
            constants$legacy <- legacy_mat

            prop_legacy <- mean(constants$legacy)
            info("Legacy covariate set. Proportion legacy cells: %.3f", prop_legacy)
            if (prop_legacy %in% c(0,1)) warn("Legacy covariate is all-%d; model may have identifiability issues.", prop_legacy)
            
            # Validate legacy covariate to prevent extreme values
            legacy_sum <- sum(constants$legacy)
            legacy_total <- length(constants$legacy)
            if (legacy_sum == 0 || legacy_sum == legacy_total) {
                warn("WARNING: Legacy covariate is all 0s or all 1s - this may cause numerical issues")
            }
            
            info("Legacy covariate added: %d legacy observations out of %d", legacy_sum, legacy_total)
            info("Legacy proportion: %.3f", legacy_sum/legacy_total)
            info("Legacy matrix dimensions: %d x %d", nrow(constants$legacy), ncol(constants$legacy))
            
            info("  ✓ Legacy covariate added successfully")
        }
        
        # Model hyperparameters - adjust based on model type
        info("Setting model hyperparameters...")
        if (model_name == "env_cycl") {
            constants$N.beta = 8
            constants$zeros <- rep(0, 8)
            constants$omega <- 0.01 * diag(8)  # precision matrix: SD = 10 per beta
        } else if (model_name == "env_cov") {
            constants$N.beta = 6
            constants$zeros <- rep(0, 6)
            constants$omega <- 0.01 * diag(6)
        } else {
            constants$N.beta = 2
            constants$zeros <- rep(0, 2)
            constants$omega <- 0.01 * diag(2)
        }
        
        info("Constants prepared successfully")
        
        # Keep driver uncertainty flags as logical for NIMBLE (cleaner and safer)
        info("  ✓ Driver uncertainty flags ready for NIMBLE (logical format)")
        
        # Create model with driver uncertainty mode
        info("Building Nimble model with driver uncertainty mode: %s", driver_uncertainty_mode)
        modelCode <- create_nimble_model_with_uncertainty(model_name, use_legacy_covariate, 
                                                         temporalDriverUncertainty = driver_uncertainty_mode, 
                                                         spatialDriverUncertainty = driver_uncertainty_mode)
        
        # Create initial values with drivers
        inits <- create_inits_with_uncertainty(constants, model_name, model_data = model.dat, chain_no = chain_no)
        
        info("Model built successfully")
        
        # STEP: Comprehensive Model Validation
        info("Performing comprehensive model validation...")
        
        # Validate model dimensions
        info("  Validating model dimensions...")
        expected_dims <- list(
            N.plot = constants$N.plot,
            N.date = constants$N.date,
            N.site = constants$N.site,
            N.core = constants$N.core
        )
        
        for (dim_name in names(expected_dims)) {
            if (expected_dims[[dim_name]] <= 0) {
                error("    ❌ ERROR: %s is %d - must be positive", dim_name, expected_dims[[dim_name]])
                stop("Invalid model dimensions")
            }
            info("    ✓ %s = %d", dim_name, expected_dims[[dim_name]])
        }
        
        # Validate environmental predictors for environmental models
        if (model_name %in% c("env_cycl", "env_cov")) {
            info("  Validating environmental predictors...")
            required_env <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
            
            for (pred in required_env) {
                if (pred %in% names(constants)) {
                    pred_data <- constants[[pred]]
                    if (is.matrix(pred_data)) {
                        info("    ✓ %s: %s matrix", pred, paste(dim(pred_data), collapse = " x "))
                        
                        # Check matrix dimensions match expected
                        if (nrow(pred_data) != constants$N.plot && ncol(pred_data) != constants$N.date) {
                            warn("      ⚠️  WARNING: %s dimensions may not match plot/time structure", pred)
                        }
                    } else if (is.vector(pred_data)) {
                        info("    ✓ %s: %d vector", pred, length(pred_data))
                    } else {
                        warn("    ⚠️  WARNING: %s: %s object", pred, class(pred_data))
                    }
                } else {
                    error("    ❌ ERROR: %s missing for environmental model", pred)
                    stop("Missing required environmental predictor: ", pred)
                }
            }
        }
        
        # Validate seasonal predictors
        info("  Validating seasonal predictors...")
        if (length(constants$sin_mo) != constants$N.date) {
            error("    ❌ ERROR: sin_mo length (%d) != N.date (%d)", length(constants$sin_mo), constants$N.date)
            stop("Seasonal predictor dimension mismatch")
        }
        if (length(constants$cos_mo) != constants$N.date) {
            error("    ❌ ERROR: cos_mo length (%d) != N.date (%d)", length(constants$cos_mo), constants$N.date)
            stop("Seasonal predictor dimension mismatch")
        }
        info("    ✓ Seasonal predictors validated")
        
        # Validate response data
        info("  Validating response data...")
        if (nrow(model.dat$y) != constants$N.core) {
            error("    ❌ ERROR: Response data rows (%d) != N.core (%d)", nrow(model.dat$y), constants$N.core)
            stop("Response data dimension mismatch")
        }
        # Dirichlet: response is composition matrix (N.core x N.spp), rows sum to 1
        y_values <- as.vector(model.dat$y)
        y_na <- sum(is.na(y_values))
        y_inf <- sum(is.infinite(y_values))
        y_neg <- sum(y_values < 0, na.rm = TRUE)
        y_gt1 <- sum(y_values > 1, na.rm = TRUE)

        if (y_na > 0) warn("    ⚠️  WARNING: %d missing values in response", y_na)
        if (y_inf > 0) error("    ❌ ERROR: %d infinite values in response", y_inf)
        if (y_neg > 0) error("    ❌ ERROR: %d negative values in response", y_neg)
        if (y_gt1 > 0) error("    ❌ ERROR: %d values > 1 in response", y_gt1)

        if (y_inf > 0 || y_neg > 0 || y_gt1 > 0) {
            stop("Response data contains invalid values")
        }
        
        info("    ✓ Response data validated")
        
        info("  ✓ All model validations passed")
        
        # Build model
        info("Building Nimble model...")
        Rmodel <- nimbleModel(code = modelCode, constants = constants,
                              data = list(y = model.dat$y), inits = inits)
        
        # Debug: Check data dimensions
        info("Data dimensions check:")
        info("  model.dat$y dimensions: %s", paste(dim(model.dat$y), collapse = " x "))
        info("  model.dat$y class: %s", class(model.dat$y))
        info("  constants$N.core: %d", constants$N.core)
        info("  constants$N.spp: %d", constants$N.spp)
        
        # Compile model
        info("Compiling Nimble model...")
        cModel <- compileNimble(Rmodel)
        
        info("Model compiled successfully")
        
        # Configure MCMC with proper sampler management
        info("Configuring MCMC...")
        
        # Enhanced monitoring for comprehensive convergence analysis (matching working approach)
        # Start with core parameters that all models have
        monitored_params <- c(
            "rho", "beta", "site_effect_sd", "site_effect", "intercept",
            "sigma_proc", "sigma_init"
        )
        
        # Add legacy_effect only for models that use it
        if (use_legacy_covariate && model_name %in% c("env_cycl", "env_cov")) {
            monitored_params <- c(monitored_params, "legacy_effect")
        }
        
        # Add driver uncertainty parameters for DRIVER UNCERTAINTY models
        # NOTE: Removed temp_est, mois_est, pH_est, pC_est from monitoring to reduce memory usage
        # These are large arrays that can cause performance issues during testing
        if (model_name %in% c("env_cycl", "env_cov")) {
            # Only monitor the spatial pH/pC variables if using driver uncertainty
            if (constants$spatialDriverUncertainty) {
                monitored_params <- c(monitored_params, "pH_est_p", "pC_est_p")
            }
        }
        
        # Use monitors2 for latent variables (Dirichlet: alpha and plot_rel)
        monitored_latent_params <- c("alpha", "plot_rel")
        
        info("Monitoring parameters for convergence analysis:")
        info("  Core parameters: %s", paste(monitored_params, collapse = ", "))
        info("  Latent variables: %s", paste(monitored_latent_params, collapse = ", "))
        info("  Total beta parameters: %d", constants$N.beta)
        
        # Thinning: configurable via env vars to manage memory for long runs
        thin  <- as.integer(Sys.getenv("THIN",  unset = "5"))
        thin2 <- as.integer(Sys.getenv("THIN2", unset = "10"))
        info("Thinning: thin=%d (parameters), thin2=%d (latent states)", thin, thin2)
        
        mcmcConf <- configureMCMC(
            model = Rmodel,
            monitors = monitored_params,
            monitors2 = monitored_latent_params,  # Use monitors2 for latent variables
            thin = thin,
            thin2 = thin2,
            enableWAIC = FALSE
        )
        
        # Add specialized samplers for better convergence of key parameters (matching working approach)
        info("Adding specialized samplers for convergence improvement...")
        
        # 1. FIRST remove default samplers to prevent conflicts
        info("  Removing default samplers...")
        samplers_to_remove <- c("rho", "beta", "site_effect_sd", "site_effect",
                                "intercept", "sigma_proc", "sigma_init")
        if (use_legacy_covariate && model_name %in% c("env_cycl", "env_cov")) {
            samplers_to_remove <- c(samplers_to_remove, "legacy_effect")
        }
        mcmcConf$removeSamplers(samplers_to_remove)

        # 2. Add specialized samplers for reparameterized Dirichlet
        info("  Adding samplers for reparameterized model...")
        for (s in 1:constants$N.spp) {
            mcmcConf$addSampler(target = paste0("rho[", s, "]"), type = "slice")
            mcmcConf$addSampler(target = paste0("intercept[", s, "]"), type = "slice")
            mcmcConf$addSampler(target = paste0("site_effect_sd[", s, "]"), type = "slice")
            mcmcConf$addSampler(target = paste0("sigma_proc[", s, "]"), type = "slice")
            mcmcConf$addSampler(target = paste0("sigma_init[", s, "]"), type = "slice")
            if (use_legacy_covariate && model_name %in% c("env_cycl", "env_cov")) {
                mcmcConf$addSampler(target = paste0("legacy_effect[", s, "]"), type = "slice")
            }
        }

        # Block sampler for beta vector per species (dmnorm prior enables joint updates)
        for (s in 1:constants$N.spp) {
            tryCatch({
                mcmcConf$addSampler(
                    target = paste0("beta[", s, ", 1:", constants$N.beta, "]"),
                    type = "AF_slice")
                info("  Added AF_slice block sampler for beta[%d, 1:%d]", s, constants$N.beta)
            }, error = function(e) {
                # Fallback to element-wise slice if AF_slice fails
                for (i in 1:constants$N.beta) {
                    mcmcConf$addSampler(target = paste0("beta[", s, ", ", i, "]"), type = "slice")
                }
                info("  Fallback: element-wise slice for beta[%d, *]", s)
            })
        }

        # 3. Site effects sampling (per species) - AF_slice for block updates
        for (s in 1:constants$N.spp) {
            if (constants$N.site > 1) {
                tryCatch({
                    mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, ", ", s, "]"), type = "AF_slice")
                }, error = function(e) {
                    mcmcConf$addSampler(target = paste0("site_effect[1:", constants$N.site, ", ", s, "]"), type = "RW_block")
                })
            } else {
                mcmcConf$addSampler(target = paste0("site_effect[1, ", s, "]"), type = "slice")
            }
        }
        
        info("MCMC configured successfully - ALL advanced features RESTORED!")
        
        # Build and compile MCMC
        info("Building and compiling MCMC...")
        myMCMC <- buildMCMC(mcmcConf)
        compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
        
        info("MCMC configured successfully")
        
        # Run MCMC with convergence-based sampling
        info("Running MCMC with convergence-based sampling...")
        # Conditional parameters: testing vs production
        is_local_test <- identical(tolower(Sys.getenv("LOCAL_TEST", "false")), "true")
        if (is_local_test) {
            # LOCAL TESTING: Reduced values for faster testing
            burnin <- 200
            iter_per_chunk <- 50
            init_iter <- 400
            min_eff_size_perchain <- 10
            max_loops <- 2
            min_total_iterations <- 100
            max_total_iterations <- 1200
            convergence_check_interval <- 1  # Check every loop in testing
            info("  🧪 LOCAL TESTING MODE: Using reduced iteration values")
        } else {
            # PRODUCTION: Values for driver uncertainty models
            # Dirichlet models need ~100k iterations for convergence (based on
            # env_cycl_all_taxa reaching median Rhat 1.18 at 10k iterations)
            burnin <- 5000
            iter_per_chunk <- 5000  # Larger chunks for efficiency (reduces overhead)
            init_iter <- 5000
            min_eff_size_perchain <- 100  # Target ESS for reliable estimates
            max_loops <- 40  # Allow enough loops to reach 200k total
            min_total_iterations <- 50000  # Minimum iterations for reliable estimates
            max_total_iterations <- 200000  # Hard limit to prevent runaway jobs
            convergence_check_interval <- 2  # Check convergence every N loops to reduce overhead
            info("  PRODUCTION MODE: Using iteration values for Dirichlet convergence")
        }
        max_iter_env <- suppressWarnings(as.integer(Sys.getenv("MAX_ITER", NA)))
        if (!is.na(max_iter_env) && max_iter_env > 0) {
            max_total_iterations <- min(max_total_iterations, max_iter_env)
            min_total_iterations <- min(min_total_iterations, max_iter_env)
            info("  MAX_ITER set: capping at %d iterations", max_iter_env)
        }
        max_save_size <- 50000
        
        info("Running MCMC with convergence-based sampling")
        info("  Initial iterations: %d burnin: %d", init_iter, burnin)
        info("  Iterations per chunk: %d max loops: %d", iter_per_chunk, max_loops)
        info("  Target ESS per chain: %d", min_eff_size_perchain)
        info("  Minimum total iterations: %d", min_total_iterations)
        info("  Maximum total iterations: %d", max_total_iterations)
        info("  Convergence check interval: every %d loops", convergence_check_interval)
        
        # Run initial iterations
        info("  Running initial iterations (%d iterations) for adaptation...", init_iter)
        compiled$run(niter = init_iter, thin = thin, thin2 = thin2, nburnin = 0)
        info("  Initial iterations completed")
        
        # Get initial samples and check convergence
        initial_samples <- as.matrix(compiled$mvSamples)
        info("  Initial samples collected, checking convergence...")
        info("  Initial samples dimensions: %s", paste(dim(initial_samples), collapse = " x "))
        
        # Create output directory for checkpoints
        info("  Creating output directory for checkpoints...")
        
        # Create species-specific subdirectory
        species_output_dir <- file.path(model_output_dir, model_name, species)
        
        # Ensure the directory exists
        if (!dir.exists(species_output_dir)) {
            dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
        }
        
        # Verify directory was created
        if (!dir.exists(species_output_dir)) {
            stop("CRITICAL: Failed to create checkpoint directory: ", species_output_dir)
        }
        
        # Create model_id for consistent naming
        model_id <- create_model_id(model_name, species, min.date, max.date, use_legacy_covariate)
        
        # Bulletproof model_id fallback if empty/NA
        if (!is.character(model_id) || !nzchar(model_id)) {
          model_id <- sprintf("mdl_%s_%s_%s",
                             gsub("\\W", "", model_name),
                             gsub("\\W", "", species),
                             format(Sys.time(), "%Y%m%d%H%M%S"))
          warn("create_model_id returned empty; using fallback id: %s", model_id)
        }
        
        # Check if we need to continue sampling for convergence
        continue <- TRUE
        loop_counter <- 0
        total_iterations <- init_iter
        
        # Check convergence
            continue <- check_continue(initial_samples, min_eff_size = min_eff_size_perchain)
        
        # Store all samples as we go
        all_samples <- initial_samples
        
        # Also accumulate samples2 (plot estimates) across loops
        # CRITICAL DEBUG: Check what mvSamples2 actually contains before extraction
        mvSamples2_raw <- compiled$mvSamples2
        info("DEBUG: mvSamples2 raw type: %s", class(mvSamples2_raw))
        if (inherits(mvSamples2_raw, "mcmc.list") && length(mvSamples2_raw) > 0) {
            info("DEBUG: mvSamples2 is mcmc.list with %d chains", length(mvSamples2_raw))
            first_chain_cols <- colnames(mvSamples2_raw[[1]])
            info("DEBUG: First chain columns (first 20): %s", paste(head(first_chain_cols, 20), collapse=", "))
            has_mu_debug <- any(grepl("^mu\\[", first_chain_cols)) || any(grepl("^eta\\[", first_chain_cols)) || any(grepl("^plot_rel\\[", first_chain_cols))
            has_param_debug <- any(grepl("^beta\\[", first_chain_cols)) || any(grepl("^intercept", first_chain_cols))
            info("DEBUG: Has mu/eta columns: %s | Has param columns: %s", has_mu_debug, has_param_debug)
        } else if (is.matrix(mvSamples2_raw)) {
            info("DEBUG: mvSamples2 is matrix with %d rows x %d cols", nrow(mvSamples2_raw), ncol(mvSamples2_raw))
            mv2_cols <- colnames(mvSamples2_raw)
            info("DEBUG: mvSamples2 columns (first 20): %s", paste(head(mv2_cols, 20), collapse=", "))
        }
        
        initial_plot_samples <- as.matrix(compiled$mvSamples2)
        
        # CRITICAL VALIDATION: Check what's actually in mvSamples2 after extraction
        if (ncol(initial_plot_samples) > 0) {
            col_names_mv2 <- colnames(initial_plot_samples)
            has_mu_cols <- any(grepl("^mu\\[", col_names_mv2)) || any(grepl("^eta\\[", col_names_mv2)) || any(grepl("^plot_rel\\[", col_names_mv2))
            has_param_cols <- any(grepl("^beta\\[", col_names_mv2)) || any(grepl("^intercept", col_names_mv2))
            
            if (has_param_cols && !has_mu_cols) {
                error("CRITICAL ERROR: mvSamples2 contains parameter columns instead of mu/eta columns!")
                error("  mvSamples2 dimensions: %d rows x %d cols", nrow(initial_plot_samples), ncol(initial_plot_samples))
                error("  mvSamples2 columns (first 20): %s", paste(head(col_names_mv2, 20), collapse=", "))
                error("  Expected mu[plot, time] or eta[plot, time] columns from monitors2 = c('eta', 'mu')")
                error("  This indicates monitors2 is NOT working - mvSamples2 is getting wrong data!")
                error("  Possible causes:")
                error("    1. NIMBLE monitors2 is not properly configured")
                error("    2. mu/eta nodes don't exist in the model with expected names")
                error("    3. mvSamples2 is accidentally pointing to mvSamples (parameters)")
                stop("mvSamples2 contains parameters instead of mu/eta - cannot continue")
            } else if (has_mu_cols) {
                info("✓ mvSamples2 contains mu/eta columns as expected")
                info("  Found %d mu/eta/plot_rel columns", sum(grepl("^mu\\[|^eta\\[|^plot_rel\\[", col_names_mv2)))
            } else {
                warn("WARNING: mvSamples2 has no mu/eta or parameter columns - may be empty or wrong structure")
                warn("  mvSamples2 dimensions: %d rows x %d cols", nrow(initial_plot_samples), ncol(initial_plot_samples))
            }
        } else {
            warn("WARNING: mvSamples2 is empty (0 columns) - monitors2 may not be working")
            warn("  This may happen if thin2 is too large or monitors2 is not configured correctly")
        }
        
        all_plot_samples <- initial_plot_samples
        cat("  Starting iterative accumulation with", nrow(all_samples), "initial samples\n")
        
        # Save initial samples as checkpoint (including samples2)
        initial_checkpoint_data <- list(
            samples = all_samples,
            samples2 = all_plot_samples,
            iterations = total_iterations,
            loop = 0
        )
        checkpoint_file <- file.path(species_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_initial.rds"))
            saveRDS(initial_checkpoint_data, checkpoint_file)
            cat("  ✓ Initial checkpoint saved with samples and samples2\n")
        
        # Also save a simple progress file
        progress_file <- create_progress_file(species_output_dir, model_id, chain_no, init_iter)
        
        while ((continue || total_iterations < min_total_iterations) && 
               loop_counter < max_loops && 
               total_iterations < max_total_iterations) {
            
            # Continue sampling without resetting
            info("  Loop %d: Running %d iterations (total: %d)", loop_counter + 1, iter_per_chunk, total_iterations)
            compiled$run(niter = iter_per_chunk, thin = thin, thin2 = thin2, nburnin = 0)
            total_iterations <- total_iterations + iter_per_chunk
            
            # Check hard limit first
            if (total_iterations >= max_total_iterations) {
                warn("  Reached maximum iteration limit (%d), stopping", max_total_iterations)
                break
            }
            
            # Get updated samples and accumulate them
            # CRITICAL: NIMBLE's mvSamples resets between runs, so current_samples contains ONLY the latest iteration
            current_samples <- as.matrix(compiled$mvSamples)
            current_sample_count <- nrow(current_samples)
            
            # Always append all current samples since mvSamples resets
            if (current_sample_count > 0) {
                all_samples <- rbind(all_samples, current_samples)
                info("  ✓ Added %d samples, total: %d", current_sample_count, nrow(all_samples))
            } else {
                warn("  WARNING: No samples in current iteration")
            }
            
            # Also accumulate samples2 (plot estimates) - mvSamples2 also resets between runs
            current_plot_samples <- as.matrix(compiled$mvSamples2)
            current_plot_count <- nrow(current_plot_samples)
            
            # CRITICAL VALIDATION: Check what's in mvSamples2 during accumulation
            if (current_plot_count > 0 && ncol(current_plot_samples) > 0) {
                col_names_current <- colnames(current_plot_samples)
                has_mu_cols <- any(grepl("^mu\\[", col_names_current)) || any(grepl("^eta\\[", col_names_current)) || any(grepl("^plot_rel\\[", col_names_current))
                has_param_cols <- any(grepl("^beta\\[", col_names_current)) || any(grepl("^intercept", col_names_current))
                
                if (has_param_cols && !has_mu_cols && loop_counter == 1) {
                    # Only warn once on first loop
                    warn("CRITICAL: mvSamples2 in loop %d contains parameters, not mu/eta!", loop_counter + 1)
                    warn("  Found columns like: %s", paste(head(col_names_current[grepl("^beta\\[|^intercept", col_names_current)], 5), collapse=", "))
                }
                
                all_plot_samples <- rbind(all_plot_samples, current_plot_samples)
                info("  ✓ Added %d plot samples, total: %d", current_plot_count, nrow(all_plot_samples))
            } else {
                warn("  WARNING: No plot samples in current iteration (rows: %d, cols: %d)", 
                     current_plot_count, if(ncol(current_plot_samples) > 0) ncol(current_plot_samples) else 0)
            }
            
            loop_counter <- loop_counter + 1
            
            # Save checkpoint after each loop (including samples2) - but only if not too large
            if (nrow(all_samples) <= max_save_size) {
                loop_checkpoint_data <- list(
                    samples = all_samples,
                    samples2 = all_plot_samples,
                    iterations = total_iterations,
                    loop = loop_counter
                )
                loop_checkpoint_file <- file.path(species_output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_loop", loop_counter, ".rds"))
                tryCatch({
                    saveRDS(loop_checkpoint_data, loop_checkpoint_file)
                    info("  ✓ Checkpoint saved (loop %d)", loop_counter)
                }, error = function(e) {
                    warn("  ✗ Failed to save checkpoint: %s", e$message)
                })
            } else {
                warn("  Skipping checkpoint save - sample size (%d) exceeds limit (%d)", nrow(all_samples), max_save_size)
            }
            
            # Update progress file
            update_progress_file(progress_file, total_iterations, loop_counter)
            
            # Check convergence every N loops to reduce overhead
            # Use recent samples (last 50%) for faster convergence checks on large datasets
            if (loop_counter %% convergence_check_interval == 0 || total_iterations >= min_total_iterations) {
                # For large sample sets, use recent samples for faster convergence checks
                samples_for_check <- if (nrow(all_samples) > 20000) {
                    # Use last 50% of samples for convergence check (keeps it faster)
                    recent_start <- floor(nrow(all_samples) * 0.5)
                    all_samples[recent_start:nrow(all_samples), , drop = FALSE]
                } else {
                    all_samples
                }
                info("  Checking convergence (using %d samples)", nrow(samples_for_check))
                continue <- check_continue(samples_for_check, min_eff_size = min_eff_size_perchain)
            }
            
            info("  Total iterations: %d | Loops: %d/%d | Converged: %s", 
                 total_iterations, loop_counter, max_loops, ifelse(!continue, "YES", "NO"))
        }
        
        if (total_iterations >= max_total_iterations) {
            warn("  Reached hard iteration limit (%d), stopping sampling", max_total_iterations)
        } else if (loop_counter >= max_loops) {
            warn("  Exceeded maximum loops (%d), stopping sampling", max_loops)
        } else if (total_iterations >= min_total_iterations) {
            if (continue) {
                warn("  Minimum iterations reached (%d) but convergence not achieved (ESS < %d)", 
                     total_iterations, min_eff_size_perchain)
            } else {
                info("  SUCCESS: Convergence reached after %d total iterations", total_iterations)
            }
        } else {
            warn("  Stopped before minimum iterations (%d) due to max loops or convergence", min_total_iterations)
        }
        
        # Update final progress status
        tryCatch({
            final_status <- if(total_iterations >= max_total_iterations) "Stopped (max iterations)" else
                if(loop_counter >= max_loops) "Completed (max loops)" else 
                if(total_iterations >= min_total_iterations && !continue) "Converged" else 
                    "Completed (min iterations)"
            writeLines(paste("Completed at:", Sys.time(), "\nTotal iterations:", total_iterations, "\nFinal loop:", loop_counter, "\nStatus:", final_status), progress_file)
            info("  ✓ Final progress status updated: %s", final_status)
        }, error = function(e) {
            warn("  ✗ Failed to update final progress status: %s", e$message)
        })
        
        # Get final samples (use accumulated samples)
        samples <- all_samples
        
        # Use accumulated plot samples instead of just final mvSamples2
        plot_samples <- all_plot_samples
        
        # Actually discard burn-in samples
        if (nrow(samples) > burnin) {
            samples <- samples[(burnin+1):nrow(samples), , drop=FALSE]
        }
        if (nrow(plot_samples) > burnin) {
            plot_samples <- plot_samples[(burnin+1):nrow(plot_samples), , drop=FALSE]
        }
        
        # Validate plot samples structure
        if (nrow(plot_samples) == 0) {
            warn("WARNING: No plot samples collected - mvSamples2 is empty")
            warn("This indicates thin2 (%d) may be too large for total iterations (%d)", thin2, total_iterations - burnin)
        } else {
            # Validate that samples2 has proper column names for combine_chains
            col_names <- colnames(plot_samples)
            if (is.null(col_names) || !any(grepl("^mu\\[|^plot_rel\\[", col_names))) {
                warn("WARNING: samples2 does not contain mu columns")
            } else {
                info("✓ samples2 contains %d mu/plot_rel columns", sum(grepl("^mu\\[|^plot_rel\\[", col_names)))
            }
        }
        
        
        # Final convergence check (safely ignore zero-variance columns)
        final_ess <- apply(samples, 2, function(x) {
            if (sum(!is.na(x)) < 2 || is.na(var(x, na.rm = TRUE)) || var(x, na.rm = TRUE) == 0) return(NA)
            tryCatch(coda::effectiveSize(coda::as.mcmc(x)), error = function(e) NA)
        })
        valid_final_ess <- final_ess[!is.na(final_ess)]
        min_final_ess <- if (length(valid_final_ess) > 0) min(valid_final_ess) else NA
        
        cat("=== ITERATIVE SAVING SUMMARY ===\n")
        cat("  Initial samples:", nrow(initial_samples), "iterations\n")
        cat("  Additional loops:", loop_counter, "iterations\n")
        cat("  Total accumulated samples:", nrow(all_samples), "iterations\n")
        cat("  Total accumulated plot samples:", nrow(all_plot_samples), "iterations\n")
        cat("  Checkpoints saved:", loop_counter + 1, "files\n")
        
        # Optional safety net: rebuild mu from eta if needed
        plot_samples <- fix_mu_from_eta(plot_samples, constants$N.plot, constants$N.date)
        
        # Log save keys for debugging
        info("SAVE KEYS: model_name=%s species=%s model_id=%s", model_name, species, model_id)
        
        # Use compatibility shim to standardize output structure
        compat <- ensure_old_schema(
            samples = all_samples,
            samples2 = plot_samples,
            meta = list(
                rank.name = rank.name,
                species = species,
                model_name = model_name,
                model_id = model_id,
                use_legacy_covariate = use_legacy_covariate,
                has_driver_uncertainty = driver_uncertainty_mode,
                scenario = scenario,
                min.date = min.date,
                max.date = max.date,
                niter = total_iterations,
                nburnin = burnin,
                thin = thin,
                thin2 = thin2,
                model_data = model.dat,
                nimble_code = modelCode,
                model_structure = if (driver_uncertainty_mode)
                    "dirichlet_with_driver_uncertainty"
                  else
                    "dirichlet_fixed_drivers"
            ),
            model_output_root = model_output_dir,
            model_name = model_name,
            species = species,
            model_id = model_id,
            chain_no = chain_no
        )
        
        # Log what we're about to save
        cat("About to save to:", compat$path, "\n")
        cat("Parent dir exists:", dir.exists(dirname(compat$path)), "\n")
        cat("Writable parent? (0 ok):", file.access(dirname(compat$path), 2), "\n")
        
        # Save final samples with standardized structure (atomic save)
        tryCatch({
          atomic_saveRDS(compat$chain_output, compat$path)
          cat("✓ SUCCESS: Saved MCMC samples to:", compat$path, "\n")
        }, error = function(e) {
          cat("SAVE ERROR:", conditionMessage(e), "\n")
          stop(e)
        })
        
        cat("Sample dimensions:", dim(all_samples), "\n")
        cat("=== Model fitting completed (cloglog version) ===\n")
        return(list(
            status = "SUCCESS", 
            samples = all_samples,
            samples2 = plot_samples,  # Include plot-level estimates for consistency
            file = compat$path,
            model_data = model.dat,
            nimble_code = modelCode,
            metadata = list(
                rank.name = rank.name,
                species = species,
                model_name = model_name,
                model_id = model_id,
                use_legacy_covariate = use_legacy_covariate,
                scenario = scenario,
                min.date = min.date,
                max.date = max.date,
                niter = total_iterations,
                nburnin = burnin,
                thin = thin,
                thin2 = thin2,  # Include thin2 for samples2
                model_data = model.dat,
                nimble_code = modelCode,
                model_structure = if (driver_uncertainty_mode)
                    "dirichlet_with_driver_uncertainty"
                  else
                    "dirichlet_fixed_drivers"
            )
        ))
        
    }, error = function(e) {
        # Capture comprehensive error information
        error_time <- Sys.time()
        error_context <- list(
            timestamp = error_time,
            task_idx = j,
            chain_no = chain_no,
            error_message = if(!is.null(e$message) && e$message != "") e$message else "No error message available",
            error_call = if(!is.null(e$call)) paste(deparse(e$call), collapse=" ") else "No call information",
            error_class = class(e)[1],
            system_info = list(
                r_version = R.version.string,
                working_dir = getwd(),
                available_packages = installed.packages()[,"Package"],
                memory_usage = if(exists("gc")) gc() else "GC not available"
            ),
            runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA
        )
        
        
        # Also log to console with detailed information
        cat("ERROR in Model", j, "Chain", error_context$chain_no, ":\n")
        cat("  Message:", error_context$error_message, "\n")
        cat("  Call:", error_context$error_call, "\n")
        cat("  Class:", error_context$error_class, "\n")
        cat("  Runtime:", round(error_context$runtime, 2), "seconds\n")
        
        # Return detailed error information
        return(list(
            status = "ERROR", 
            error = error_context$error_message,
            error_details = error_context
        ))
    })
}

# ---------- SINGLE-TASK SELECTION (HPC) ----------
# Now that run_scenarios_fixed is defined, handle single-task mode
if (single_task_mode) {
  # Apply semantic filters if provided
  if (!is.null(cli_model_id)) {
    valid_models <- valid_models %>% dplyr::filter(model_id == cli_model_id)
    info("Filtered to model_id=%s -> %d rows", cli_model_id, nrow(valid_models))
  }
  if (!is.null(cli_rank)) {
    valid_models <- valid_models %>% dplyr::filter(rank.name == cli_rank)
    info("Filtered to rank.name=%s -> %d rows", cli_rank, nrow(valid_models))
  }
  if (!is.null(cli_species)) {
    valid_models <- valid_models %>% dplyr::filter(species == cli_species)
    info("Filtered to species=%s -> %d rows", cli_species, nrow(valid_models))
  }
  if (!is.null(cli_model_name)) {
    valid_models <- valid_models %>% dplyr::filter(model_name == cli_model_name)
    info("Filtered to model_name=%s -> %d rows", cli_model_name, nrow(valid_models))
  }

  if (nrow(valid_models) == 0) {
    stop("No models remain after filtering; nothing to run.")
  }

  # Rebuild per-model chain count and task list (same logic as main path)
  valid_models <- valid_models %>%
    mutate(nchains_model = as.integer(nchains))
  all_tasks <- bind_rows(lapply(1:nrow(valid_models), function(i) {
    data.frame(model_idx = i, chain_no = 1:valid_models$nchains_model[i])
  }))
  total_tasks <- nrow(all_tasks)
  if (array_task_id < 1 || array_task_id > total_tasks) {
    stop(sprintf("array_task_id=%d out of bounds (1..%d)", array_task_id, total_tasks))
  }
  model_idx <- all_tasks$model_idx[array_task_id]
  chain_no   <- all_tasks$chain_no[array_task_id]
  n_models   <- nrow(valid_models)

  info("HPC single-task mode: task=%d of %d -> model_idx=%d / %d, chain_no=%d",
       array_task_id, total_tasks, model_idx, n_models, chain_no)

  # Set deterministic seed per (model_idx, chain_no) for reproducibility
  set.seed(100000 + model_idx*100 + chain_no)
  info("Set seed to %d for model_idx=%d, chain_no=%d", 100000 + model_idx*100 + chain_no, model_idx, chain_no)

  # run just one (model, chain) and exit
  res <- run_scenarios_fixed(j = model_idx, chain_no = chain_no)

  # save status summary (optional)
  out_dir <- file.path(model_output_dir, "job_summaries")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(out_dir, sprintf("status_task%06d.txt", array_task_id))
  writeLines(c(
    sprintf("task_id=%d", array_task_id),
    sprintf("model_idx=%d", model_idx),
    sprintf("chain_no=%d", chain_no),
    sprintf("status=%s", ifelse(is.list(res) && !is.null(res$status), res$status, "UNKNOWN")),
    sprintf("timestamp=%s", as.character(Sys.time()))
  ), status_file)

  q(save = "no")  # end this job cleanly
}
# ---------- END SINGLE-TASK SELECTION ----------

# Task list already built above with per-model chain counts (env_cycl+all_taxa=2, others=nchains)
cat("Total parallel tasks:", nrow(all_tasks), "(", nrow(valid_models), "models, env_cycl+all_taxa=2 chains, others=", nchains, ")\n")
cat("Task details:\n")
print(all_tasks)
cat("Cluster size:", n_cores, "workers\n")

# Function to monitor progress in real-time


# Run everything in parallel with incremental saving
cat("Preparing parallel execution functions...\n")

# Create a function that saves results as they complete
runAndSave_task <- function(task_idx) {
    
    # Test if we can access the function
    if (!exists("run_scenarios_fixed")) {
        cat("ERROR: run_scenarios_fixed function not found in worker\n")
        return(list(error = "run_scenarios_fixed function not found"))
    }
    
    
    # Initialize error tracking
    error_details <- list()
    start_time <- Sys.time()
    
    tryCatch({
        # Get task details
        task <- all_tasks[task_idx, ]
        model_idx <- task$model_idx
        chain_no <- task$chain_no
        
        
        # Log system information for debugging
        cat("Worker: System info - R version:", R.version.string, "\n")
        cat("Worker: Working directory:", getwd(), "\n")
        cat("Worker: Available packages:", paste(installed.packages()[,"Package"], collapse=", "), "\n")
        
        # Check if required packages are loaded
        required_packages <- c("nimble", "microbialForecast", "here", "tidyverse", "coda")
        missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
        if (length(missing_packages) > 0) {
            stop("Missing required packages: ", paste(missing_packages, collapse=", "))
        }
        
        # Check if data is available
        if (!exists("all_ranks") || is.null(all_ranks)) {
            stop("Data 'all_ranks' not available in worker environment")
        }
        
        # Check if valid_models is available
        if (!exists("valid_models") || is.null(valid_models) || nrow(valid_models) == 0) {
            stop("Parameters 'valid_models' not available or empty in worker environment")
        }
        
        # Validate model index
        if (model_idx > nrow(valid_models)) {
            stop("Model index ", model_idx, " exceeds available models (", nrow(valid_models), ")")
        }
        
        cat("Worker: All checks passed, calling run_scenarios_fixed...\n")
        
        # Run the model
        cat("Worker: About to call run_scenarios_fixed with j =", model_idx, "chain_no =", chain_no, "\n")
        result <- run_scenarios_fixed(j = model_idx, chain_no = chain_no)
        
        cat("Worker: run_scenarios_fixed completed successfully\n")
        cat("Worker: Result class:", class(result), "\n")
        cat("Worker: Result names:", paste(names(result), collapse=", "), "\n")
        if ("status" %in% names(result)) {
            cat("Worker: Result status:", result$status, "\n")
        }
        
        # Validate result structure
        if (!is.list(result) || !("status" %in% names(result))) {
            stop("Invalid result structure from run_scenarios_fixed")
        }
        
        # Skip re-saving if run_scenarios_fixed already saved the file
        if (is.list(result) && identical(result$status, "SUCCESS")) {
            if (!is.null(result$file) && file.exists(result$file)) {
                # File already saved by run_scenarios_fixed - skip duplicate save
                if (exists("info")) info("Worker: chain already saved at %s; skipping duplicate save", result$file)
                cat("Worker: chain already saved by run_scenarios_fixed, skipping duplicate save\n")
                return(list(model_idx = model_idx, chain_no = chain_no, result = result))
            }
            
            # Fallback: If for some reason the file wasn't saved, save it now
            # Create model_id for consistent naming using helper function
            model_id <- create_model_id(
                valid_models$model_name[model_idx], 
                valid_models$species[model_idx],
                valid_models$min.date[model_idx], 
                valid_models$max.date[model_idx],
                grepl("Legacy with covariate", valid_models$scenario[model_idx])
            )
            
            # Bulletproof model_id fallback if empty/NA
            if (!is.character(model_id) || !nzchar(model_id)) {
              model_id <- sprintf("mdl_%s_%s_%s",
                                 gsub("\\W", "", valid_models$model_name[model_idx]),
                                 gsub("\\W", "", valid_models$species[model_idx]),
                                 format(Sys.time(), "%Y%m%d%H%M%S"))
              if (exists("warn")) warn("create_model_id returned empty; using fallback id: %s", model_id)
              cat("Worker fallback: create_model_id returned empty; using fallback id:", model_id, "\n")
            }
            
            # Create the complete chain structure with metadata
            # Use the metadata from the result if available, otherwise create it
            if ("metadata" %in% names(result) && !is.null(result$metadata)) {
                # Use the complete metadata from the result
                metadata <- result$metadata
                # Add parallel execution specific fields
                metadata$task_idx <- task_idx
                metadata$completed_at <- Sys.time()
            } else {
                # Fallback: create metadata if not available in result
                metadata <- list(
                    rank.name = valid_models$rank.name[model_idx],
                    species = valid_models$species[model_idx],
                    model_name = valid_models$model_name[model_idx],
                    model_id = model_id,
                    use_legacy_covariate = grepl("Legacy with covariate", valid_models$scenario[model_idx]),
                    has_driver_uncertainty = driver_uncertainty_mode,  # Flag to identify driver uncertainty mode
                    scenario = valid_models$scenario[model_idx],
                    min.date = valid_models$min.date[model_idx],
                    max.date = valid_models$max.date[model_idx],
                    niter = nrow(result$samples),
                    nburnin = 500,  # Default burnin
                    thin = 1,       # Default thin
                    task_idx = task_idx,
                    completed_at = Sys.time(),
                    model_data = result$model_data,  # Include model_data from result
                    nimble_code = result$nimble_code,  # Include nimble_code if available
                    model_structure = if (driver_uncertainty_mode)
                    "dirichlet_with_driver_uncertainty"
                  else
                    "dirichlet_fixed_drivers"
                )
            }
            
            # Optional safety net: rebuild mu from eta if needed
            samples2_to_use <- if("samples2" %in% names(result)) result$samples2 else result$samples
            # Note: We don't have access to constants here, so we'll skip the fix_mu_from_eta call
            # The main fix is in the run_scenarios_fixed function
            
            # Use compatibility shim to standardize output structure
            compat <- ensure_old_schema(
                samples = result$samples,
                samples2 = samples2_to_use,
                meta = metadata,
                model_output_root = model_output_dir,
                model_name = valid_models$model_name[model_idx],
                species = valid_models$species[model_idx],
                model_id = model_id,
                chain_no = chain_no
            )
            
            # Log fallback save details
            cat("Worker fallback save path:", compat$path, "\n")
            cat("Parent dir exists:", dir.exists(dirname(compat$path)), "\n")
            cat("Writable parent? (0 ok):", file.access(dirname(compat$path), 2), "\n")
            
            # Save MCMC samples with standardized structure (atomic save)
            tryCatch({
              atomic_saveRDS(compat$chain_output, compat$path)
              cat("SAVED (fallback): Chain", chain_no, "for model", model_idx, "to", compat$path, "\n")
            }, error = function(e) {
              cat("FALLBACK SAVE ERROR:", conditionMessage(e), "\n")
              stop(e)
            })
            
        }
        
        cat("Worker: Model", model_idx, "Chain", chain_no, "completed at", format(Sys.time()), "\n")
        return(list(model_idx = model_idx, chain_no = chain_no, result = result))
        
    }, error = function(e) {
        # Capture comprehensive error information
        error_time <- Sys.time()
        error_details <- list(
            timestamp = error_time,
            task_idx = task_idx,
            model_idx = if(exists("model_idx")) model_idx else NA,
            chain_no = if(exists("chain_no")) chain_no else NA,
            error_message = if(!is.null(e$message) && e$message != "") e$message else "No error message available",
            error_call = if(!is.null(e$call)) paste(deparse(e$call), collapse=" ") else "No call information",
            error_class = class(e)[1],
            system_info = list(
                r_version = R.version.string,
                working_dir = getwd(),
                available_packages = installed.packages()[,"Package"],
                memory_usage = if(exists("gc")) gc() else "GC not available"
            ),
            runtime = if(exists("start_time")) difftime(error_time, start_time, units="secs") else NA
        )
        
        
        # Also log to console with detailed information
        cat("ERROR in Model", error_details$model_idx, "Chain", error_details$chain_no, ":\n")
        cat("  Message:", error_details$error_message, "\n")
        cat("  Call:", error_details$error_call, "\n")
        cat("  Class:", error_details$error_class, "\n")
        cat("  Runtime:", round(error_details$runtime, 2), "seconds\n")
        
        # Return detailed error information
        return(list(
            model_idx = error_details$model_idx, 
            chain_no = error_details$chain_no, 
            result = list(
                status = "ERROR", 
                error = error_details$error_message,
                error_details = error_details,
            )
        ))
    })
}


# Set start time for runtime calculation
start_time <- Sys.time()

# Create cluster for parallel execution - LOCAL TESTING with 4 cores
cat("Creating cluster with", n_cores, "cores for", nrow(all_tasks), "tasks (", nrow(valid_models), "models, env_cycl+all_taxa=2 chains, others=", nchains, ")\n")
cat("LOCAL TESTING: 2 cores allocated for local testing\n")

# Worker initialization function to set up worker environment
worker_init <- function(project_root, driver_uncertainty_mode) {
  ## — working directory / here() root —
  setwd(project_root)
  if (requireNamespace("here", quietly = TRUE)) {
    ## anchor here() root properly for workers (explicit, no sentinel needed)
    try(here::set_here(project_root), silent = TRUE)
    assign("here", get("here", asNamespace("here")), envir = .GlobalEnv)
  }

  ## — logging shims if logging.R wasn't sourced —
  if (!exists("info"))  info  <- function(fmt, ...) cat(sprintf(paste0(fmt, "\n"), ...))
  if (!exists("warn"))  warn  <- function(fmt, ...) cat(sprintf(paste0("WARNING: ", fmt, "\n"), ...))
  if (!exists("error")) error <- function(fmt, ...) { msg <- sprintf(fmt, ...); cat("ERROR: ", msg, "\n"); }

  ## if available, wire up your real logger (optional but nice)
  try({
    source(here::here("analysis/model_analysis/logging.R"))
    log_setup(logfile = here::here("analysis", "model_analysis", "logs", sprintf("worker_%s.log", Sys.getpid())))
  }, silent = TRUE)

  ## — packages and nimble options (must run on workers too) —
  suppressPackageStartupMessages({
    library(nimble)
    library(coda)
    library(microbialForecast)
    library(tidyverse)
    library(here)
  })
  nimbleOptions(buildInterfacesForCompiledNimbleFunctions = FALSE)
  nimbleOptions(optimize = TRUE)

  ## prevent thread oversubscription on HPC
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

  ## keep a copy of flags workers need
  assign("driver_uncertainty_mode", driver_uncertainty_mode, envir = .GlobalEnv)

  invisible(TRUE)
}

# Create cluster with explicit error handling
tryCatch({
    cl <- makeCluster(n_cores, type = "PSOCK")
    cat("✓ Cluster created successfully with", length(cl), "workers\n")
    
    
}, error = function(e) {
    cat("✗ ERROR creating cluster:", e$message, "\n")
    cat("Falling back to sequential execution...\n")
    cl <- NULL
})

# Register parallel backend with explicit fallback
cat("Registering parallel backend...\n")
if (is.null(cl)) {
    doParallel::registerDoSEQ()
    cat("Cluster creation failed, using sequential backend\n")
} else {
    doParallel::registerDoParallel(cl)
    cat("Parallel backend registered successfully\n")
}
cat("Backend registration complete. Checking registration...\n")
cat("getDoParWorkers():", getDoParWorkers(), "\n")
cat("getDoParName():", getDoParName(), "\n")

# Only run parallel section if not in single-task mode
if (!single_task_mode) {
  # Export ALL necessary variables and functions to workers
  cat("Exporting all necessary variables to workers...\n")
must_export <- c(
  # task plumbing
  "runAndSave_task", "run_scenarios_fixed", "all_tasks", "valid_models",
  "all_ranks", "model_output_dir", "nchains",

  # config/flags/roots
  "project_root", "driver_uncertainty_mode",

  # helpers used inside run_scenarios_fixed
  "create_nimble_model_with_uncertainty",
  "sanitize_driver_uncertainty",
  "assert_matrix_dims", "assert_vector_len",
  "create_inits_with_uncertainty",
  "check_continue",
  "atomic_saveRDS", "ensure_old_schema",  # Save utilities

  # any other functions you call (from your script)
  "load_required_packages", "create_directories_safe",
  "create_model_id", "save_checkpoint_safe",
  "create_progress_file", "update_progress_file"
)

if (!is.null(cl)) {
  clusterExport(cl, must_export)
  
  # Initialize every worker once
  ## RNG across workers (reproducible & independent)
  parallel::clusterSetRNGStream(cl, iseed = 12345)
  
  ## Make workers look like master
  clusterEvalQ(cl, {
    ## We'll define a trivial shim; we'll replace it next line
    TRUE
  })
  clusterCall(cl, worker_init, project_root = project_root,
              driver_uncertainty_mode = driver_uncertainty_mode)
  
  ## (Optional) sanity checks
  clusterEvalQ(cl, list(
    wd = getwd(),
    has_nimble = "nimble" %in% loadedNamespaces(),
    has_info = exists("info"),
    here_root = try(here::here(), silent = TRUE)
  ))
} else {
  cat("No cluster available; running sequentially with foreach/doSEQ.\n")
}

# Run everything in parallel with incremental saving
cat("Starting foreach loop with", nrow(all_tasks), "tasks\n")

all_results_parallel = foreach(task_idx = 1:nrow(all_tasks),
                                .packages = c("nimble","coda","tidyverse","here","microbialForecast"),
                                .export   = character(0)) %dopar% {
  ## if someone runs this block without worker_init, still be safe:
  if (!exists("info")) info <- function(fmt, ...) cat(sprintf(paste0(fmt,"\n"), ...))
  info("DEBUG: Worker %s starting task %d at %s", Sys.getpid(), task_idx, Sys.time())
  runAndSave_task(task_idx)
                               }

cat("Foreach loop completed. Results length:", length(all_results_parallel), "\n")
cat("Parallel execution completed at:", format(Sys.time()), "\n")

# Stop cluster
if (!is.null(cl)) stopCluster(cl)

} # End of single_task_mode guard

# Only show progress summary if not in single-task mode
if (!single_task_mode) {
# Show progress summary
cat("\n=== PROGRESS SUMMARY ===\n")
cat("Checking which chains have been completed...\n")

# Count completed chains from parallel results
completed_chains <- 0
error_chains <- 0
for (i in 1:length(all_results_parallel)) {
    result <- all_results_parallel[[i]]
    if ("error" %in% names(result)) {
        error_chains <- error_chains + 1
        cat("✗ Task", i, "failed with error:", result$error, "\n")
    } else if ("result" %in% names(result) && "status" %in% names(result$result)) {
        if (result$result$status == "SUCCESS") {
            completed_chains <- completed_chains + 1
            cat("✓ Model", result$model_idx, "Chain", result$chain_no, "completed\n")
        } else {
            error_chains <- error_chains + 1
            cat("✗ Model", result$model_idx, "Chain", result$chain_no, "failed with status:", result$result$status, "\n")
        }
    } else if ("status" %in% names(result)) {
        if (result$status == "SUCCESS") {
            completed_chains <- completed_chains + 1
            cat("✓ Model", result$model_idx, "Chain", result$chain_no, "completed\n")
        } else {
            error_chains <- error_chains + 1
            cat("✗ Model", result$model_idx, "Chain", result$chain_no, "failed with status:", result$status, "\n")
        }
    } else {
        error_chains <- error_chains + 1
        cat("? Task", i, "has unknown result structure\n")
        cat("  Available names:", paste(names(result), collapse=", "), "\n")
        if ("result" %in% names(result)) {
            cat("  Nested result names:", paste(names(result$result), collapse=", "), "\n")
        }
    }
}

total_chain_tasks <- nrow(all_tasks)
cat("  Completed chains:", completed_chains, "/", total_chain_tasks, "\n")
cat("  Failed chains:", error_chains, "/", total_chain_tasks, "\n")
cat("  Success rate:", round(completed_chains / total_chain_tasks * 100, 1), "%\n")

# Reorganize results by model (per-model chain count)
all_results <- list()
for (model_idx in 1:nrow(valid_models)) {
    all_results[[model_idx]] <- list()
    nch <- valid_models$nchains_model[model_idx]
    for (chain_no in 1:nch) {
        # Find the result for this model/chain combination
        task_idx <- which(all_tasks$model_idx == model_idx & all_tasks$chain_no == chain_no)
        if (length(task_idx) > 0) {
            task_result <- all_results_parallel[[task_idx]]
            if ("result" %in% names(task_result)) {
                all_results[[model_idx]][[chain_no]] <- task_result$result
            } else {
                all_results[[model_idx]][[chain_no]] <- task_result
            }
        } else {
            all_results[[model_idx]][[chain_no]] <- list(status = "NOT_FOUND")
        }
    }
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("ALL MODELS COMPLETED\n")
cat("Total runtime:", round(runtime, 1), "minutes\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Summary of all models
cat("\nSummary of All Models:\n")
for (model_idx in 1:nrow(valid_models)) {
    cat("\nModel", model_idx, ":", valid_models$species[model_idx], "(", valid_models$model_name[model_idx], ")\n")
    
    output.list <- all_results[[model_idx]]
    
    # Status summary for this model
    status_summary <- sapply(output.list, function(x) {
        if (is.list(x) && "status" %in% names(x)) {
            x$status
        } else {
            "ERROR"
        }
    })
    
    cat("  Results:", paste(status_summary, collapse = ", "), "\n")
    
    # Detailed status for this model
    for (i in 1:length(output.list)) {
        if (is.list(output.list[[i]]) && "status" %in% names(output.list[[i]])) {
            if (output.list[[i]]$status == "SUCCESS") {
                if ("samples" %in% names(output.list[[i]]) && !is.null(output.list[[i]]$samples)) {
                    sample_dim <- dim(output.list[[i]]$samples)
                    if (length(sample_dim) >= 1) {
                        cat("    Chain", i, ": SUCCESS - Samples:", sample_dim[1], "iterations\n")
                    } else {
                        cat("    Chain", i, ": SUCCESS - Samples: unknown dimensions\n")
                    }
                } else {
                    cat("    Chain", i, ": SUCCESS - No samples data\n")
                }
            } else {
                error_msg <- if ("error" %in% names(output.list[[i]])) output.list[[i]]$error else "Unknown error"
                cat("    Chain", i, ": ERROR -", error_msg, "\n")
            }
        } else {
            cat("    Chain", i, ": ERROR - Unexpected output format\n")
        }
    }
}

} # End of single_task_mode guard for progress summary


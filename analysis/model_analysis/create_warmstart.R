# Create warm-start file from the completed 20k Dirichlet run (3 chains)
library(here)  # resolve project-root paths regardless of working directory
chain_dir <- here("data/model_outputs/dirichlet_driver_uncertainty_reparam_20k/env_cycl/all_taxa")
chain_files <- list.files(chain_dir, pattern="^samples_.*chain[0-9][.]rds$", full.names=TRUE)
cat("Found", length(chain_files), "chain files\n")

# Average parameter means across all chains
all_means <- list()
for (f in chain_files) {
  s <- readRDS(f)
  all_means[[length(all_means)+1]] <- colMeans(s$samples)
}

# Average across chains
param_names <- names(all_means[[1]])
param_means <- setNames(
  sapply(param_names, function(p) mean(sapply(all_means, function(m) m[p]))),
  param_names
)

cat("Parameters:", length(param_means), "\n")
cat("First 10:", paste(head(names(param_means), 10), collapse=", "), "\n")

ws <- list(
  param_means = param_means,
  source_iterations = 20000,
  source_dir = chain_dir,
  created = Sys.time()
)

out_file <- file.path(chain_dir, "warmstart_inits.rds")
saveRDS(ws, out_file)
cat("Saved:", out_file, "\n")

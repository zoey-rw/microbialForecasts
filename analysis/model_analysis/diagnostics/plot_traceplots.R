library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)  # resolve project-root paths regardless of working directory

# Identify worst parameters from new chain 1
new_files <- sort(list.files(here("data/model_outputs/dirichlet_driver_uncertainty_reparam_75k/env_cycl/all_taxa"),
                              pattern="checkpoint.*chain1.*loop", full.names=TRUE))
new1 <- readRDS(tail(new_files, 1))
ess1 <- apply(new1$samples, 2, function(x) tryCatch(effectiveSize(as.mcmc(x)), error=function(e) NA))
worst_params <- names(sort(ess1[!is.na(ess1)]))[1:6]
cat("Worst ESS params:", paste(worst_params, collapse=", "), "\n")

# Load old 20k chains
old_dir <- here("data/model_outputs/dirichlet_driver_uncertainty_reparam_20k/env_cycl/all_taxa")
old_chain_files <- sort(list.files(old_dir, pattern="^samples_.*chain[0-9].rds$", full.names=TRUE))

# Load new chain checkpoints (latest for each chain)
new_dir <- here("data/model_outputs/dirichlet_driver_uncertainty_reparam_75k/env_cycl/all_taxa")

# Build trace data
trace_list <- list()

for (i in seq_along(old_chain_files)) {
  s <- readRDS(old_chain_files[i])
  for (p in worst_params) {
    if (p %in% colnames(s$samples)) {
      trace_list[[length(trace_list)+1]] <- data.frame(
        iteration = seq_len(nrow(s$samples)),
        value = s$samples[, p],
        param = p,
        chain = paste0("Old 20k chain ", i),
        run = "Old (20k cold)"
      )
    }
  }
  rm(s); gc(verbose=FALSE)
}

# New chains - load all checkpoints to get full trace
for (ch in 1:3) {
  ch_files <- sort(list.files(new_dir, pattern=paste0("checkpoint.*chain", ch, ".*loop"), full.names=TRUE))
  if (length(ch_files) == 0) next

  # Use latest checkpoint (has all accumulated samples)
  s <- readRDS(tail(ch_files, 1))
  samp <- s$samples

  for (p in worst_params) {
    if (p %in% colnames(samp)) {
      trace_list[[length(trace_list)+1]] <- data.frame(
        iteration = seq_len(nrow(samp)),
        value = samp[, p],
        param = p,
        chain = paste0("New 75k chain ", ch),
        run = "New (75k warm)"
      )
    }
  }
  rm(s, samp); gc(verbose=FALSE)
}

trace_df <- bind_rows(trace_list)

# Create traceplot
p <- ggplot(trace_df, aes(x = iteration, y = value, color = chain)) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~param, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(
    "Old 20k chain 1" = "#E69F00",
    "Old 20k chain 2" = "#D55E00",
    "Old 20k chain 3" = "#CC79A7",
    "New 75k chain 1" = "#0072B2",
    "New 75k chain 2" = "#009E73",
    "New 75k chain 3" = "#56B4E9"
  )) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95")) +
  labs(title = "Traceplots: lowest-ESS parameters",
       subtitle = "Old = 20k cold start, New = 75k warm start",
       x = "Sample index (thinned)", y = "Value", color = "Chain")

ggsave(here("figures", "dirichlet_traceplots_comparison.png"), p, width = 10, height = 8, dpi = 150)
cat("Saved: figures/dirichlet_traceplots_comparison.png\n")

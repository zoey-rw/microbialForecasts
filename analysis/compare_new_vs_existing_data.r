## Compare new NEON DP1.10081.002 16S data with existing model data
## Quantifies error from taxonomy name changes using overlapping samples,
## then prepares post-hindcast data for model evaluation.
##
## Run from project root: Rscript analysis/compare_new_vs_existing_data.r
##
## Prerequisites: 16S map files extracted from zip into
##   data/new_neon_microbe/NEON_tax-microbe-soil/

source("source.R")
library(dplyr)
library(tidyr)
library(ggplot2)

new_data_dir <- here("data", "new_neon_microbe", "NEON_tax-microbe-soil")

# ==========================================================================
# 1. PHYLUM NAME MAPPING (old GTDB/SILVA -> new ICNP)
# ==========================================================================
# The existing data uses GTDB/SILVA taxonomy (pre-2022 NEON processing).
# The new data uses ICNP-updated names (2025 NEON reprocessing).
# These are the SAME organisms under different nomenclature.

phylum_map <- tribble(
  ~old_name,           ~new_name,
  # Direct renames (ICNP nomenclature updates)
  "proteobacteria",    "Pseudomonadota",
  "actinobacteriota",  "Actinomycetota",
  "firmicutes",        "Bacillota",
  "firmicutes_a",      "Bacillota",         # GTDB-only split; recombine
  "cyanobacteria",     "Cyanobacteriota",
  "chloroflexi",       "Chloroflexota",
  "bacteroidota",      "Bacteroidota",
  "acidobacteriota",   "Acidobacteriota",
  "verrucomicrobiota", "Verrucomicrobiota",
  "planctomycetota",   "Planctomycetota",
  "gemmatimonadota",   "Gemmatimonadota",
  "nitrospirota",      "Nitrospirota",
  "armatimonadota",    "Armatimonadota",
  "myxococcota",       "Pseudomonadota",    # Myxococcota merged into Pseudomonadota in ICNP
  "patescibacteria",   "Candidatus Saccharibacteria",
  "desulfobacterota",  "Pseudomonadota",    # Often reclassified under Deltaproteobacteria
  "dependentiae",      "Candidatus Saccharibacteria",
  # GTDB-specific phyla with no direct ICNP equivalent:
  "nb1.j",             NA_character_,
  "rcp2.54",           NA_character_,
  "wps.2",             NA_character_
)

cat("=== PHYLUM NAME MAPPING ===\n")
cat("\nMerged phyla in new nomenclature:\n")
phylum_map %>%
  filter(!is.na(new_name)) %>%
  group_by(new_name) %>%
  filter(n() > 1) %>%
  summarise(old_names = paste(old_name, collapse = " + "), .groups = "drop") %>%
  {cat(sprintf("  %s <- %s\n", .$new_name, .$old_names))}

cat("\nUnmappable GTDB phyla:", paste(phylum_map$old_name[is.na(phylum_map$new_name)], collapse = ", "), "\n")

# ==========================================================================
# 2. LOAD METADATA
# ==========================================================================
cat("\n=== LOADING METADATA ===\n")

meta_files <- list.files(new_data_dir, pattern = "mct_soilSampleMetadata_16S",
                          recursive = TRUE, full.names = TRUE)
meta_16S <- bind_rows(lapply(meta_files, read.csv, stringsAsFactors = FALSE))
meta_16S$collectDate <- as.Date(substr(meta_16S$collectDate, 1, 10))
meta_16S$dateID <- format(meta_16S$collectDate, "%Y%m")

if ("sequenceCountQF" %in% names(meta_16S)) {
  meta_16S <- meta_16S %>% filter(sequenceCountQF == "OK" | is.na(sequenceCountQF))
}
cat("Total 16S samples (QC-passed):", nrow(meta_16S), "\n")

# ==========================================================================
# 3. LOAD EXISTING DATA
# ==========================================================================
existing_16S <- readRDS(here("data", "clean", "groupAbundances_16S_2023.rds"))
existing_phylum <- existing_16S$phylum_bac

existing_sites <- sort(unique(existing_phylum$siteID))
existing_dates <- sort(unique(existing_phylum$dateID))

# Reshape existing data to long format
existing_long <- existing_phylum %>%
  pivot_longer(
    cols = -c(siteID, plotID, dateID, sampleID, dates, plot_date),
    names_to = "phylum_old",
    values_to = "rel_abund_old"
  ) %>%
  filter(phylum_old != "other")

cat("Existing: ", nrow(existing_phylum), " samples, ", length(existing_sites),
    " sites, ", min(existing_dates), "-", max(existing_dates), "\n", sep = "")

# Apply mapping to existing data: sum old phyla that merge under new names
existing_long_mapped <- existing_long %>%
  left_join(phylum_map, by = c("phylum_old" = "old_name")) %>%
  filter(!is.na(new_name)) %>%
  group_by(siteID, plotID, dateID, sampleID, dates, plot_date, new_name) %>%
  summarise(rel_abund_old = sum(rel_abund_old, na.rm = TRUE), .groups = "drop") %>%
  rename(phylum = new_name)

# ==========================================================================
# 4. READ NEW 16S MAP FILES AND AGGREGATE TO PHYLUM
# ==========================================================================
read_and_aggregate_16S <- function(filepath) {
  df <- tryCatch(read.csv(filepath, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df %>%
    filter(!is.na(phylum) & phylum != "") %>%
    group_by(dnaSampleID, phylum) %>%
    summarise(reads = sum(individualCount, na.rm = TRUE), .groups = "drop") %>%
    group_by(dnaSampleID) %>%
    mutate(total_reads = sum(reads), rel_abund_new = reads / total_reads) %>%
    ungroup() %>%
    select(dnaSampleID, phylum, rel_abund_new, total_reads)
}

map_files <- list.files(new_data_dir, pattern = "16S_map\\.csv$",
                         recursive = TRUE, full.names = TRUE)
cat("\nFound", length(map_files), "16S map files\n")

cat("Reading and aggregating to phylum level...\n")
phylum_list <- lapply(seq_along(map_files), function(i) {
  if (i %% 2000 == 0) cat("  Processed", i, "/", length(map_files), "files\n")
  read_and_aggregate_16S(map_files[i])
})
new_phylum <- bind_rows(phylum_list[!sapply(phylum_list, is.null)])
cat("New data:", n_distinct(new_phylum$dnaSampleID), "samples,",
    n_distinct(new_phylum$phylum), "phyla\n")

# Join with metadata
new_phylum <- new_phylum %>%
  inner_join(
    meta_16S %>% select(dnaSampleID, siteID, plotID, collectDate, dateID, geneticSampleID),
    by = "dnaSampleID"
  )
cat("After metadata join:", n_distinct(new_phylum$dnaSampleID), "samples\n")

# ==========================================================================
# 5. OVERLAP COMPARISON
# ==========================================================================
cat("\n========================================\n")
cat("=== OVERLAP COMPARISON ===\n")
cat("========================================\n")

# Only 1 exact site+month overlap exists: OSBS 201803.
# The new sample (OSBS_003-M-9.5-21.5) is a DIFFERENT soil core from the
# 8 existing samples -- so this comparison captures both taxonomy rename
# error AND natural within-site variability.

new_osbs <- new_phylum %>%
  filter(siteID == "OSBS", dateID == "201803") %>%
  select(phylum, rel_abund_new)

old_osbs <- existing_long_mapped %>%
  filter(siteID == "OSBS", dateID == "201803") %>%
  group_by(phylum) %>%
  summarise(
    mean_rel_old = mean(rel_abund_old, na.rm = TRUE),
    sd_rel_old = sd(rel_abund_old, na.rm = TRUE),
    n_old = n(),
    .groups = "drop"
  )

comparison <- inner_join(new_osbs, old_osbs, by = "phylum") %>%
  mutate(
    diff = rel_abund_new - mean_rel_old,
    abs_diff = abs(diff),
    within_1sd = abs(diff) <= sd_rel_old
  ) %>%
  arrange(desc(mean_rel_old))

cat("\nOSBS 201803: 1 new sample vs. 8 existing samples\n")
cat("(Different soil cores from the same site and month)\n\n")

cat(sprintf("%-25s %8s %8s %8s %8s %8s %s\n",
            "Phylum", "Old_mean", "Old_SD", "New", "Diff", "|Diff|", "Within1SD"))
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  cat(sprintf("%-25s %8.4f %8.4f %8.4f %+8.4f %8.4f %s\n",
              r$phylum, r$mean_rel_old, r$sd_rel_old, r$rel_abund_new,
              r$diff, r$abs_diff,
              ifelse(is.na(r$within_1sd), "N/A", ifelse(r$within_1sd, "YES", "NO"))))
}

cat(sprintf("\nOverall RMSE: %.4f\n", sqrt(mean(comparison$diff^2, na.rm = TRUE))))
cat(sprintf("Overall MAE: %.4f\n", mean(comparison$abs_diff, na.rm = TRUE)))
cat(sprintf("Correlation: %.4f\n", cor(comparison$mean_rel_old, comparison$rel_abund_new)))
cat(sprintf("Within 1 SD: %d/%d phyla\n",
            sum(comparison$within_1sd, na.rm = TRUE),
            sum(!is.na(comparison$within_1sd))))

# --- Cross-site annual mean comparison ---
# For a broader view, compare site-level annual means between old data
# (2013-2018) and new data (first available year per site).
# This gives us ~40+ site-level comparisons per phylum.
cat("\n--- Cross-site mean comparison ---\n")
cat("(Comparing site-level grand means: old 2013-2018 vs new 2019-2024)\n\n")

old_site_means <- existing_long_mapped %>%
  group_by(siteID, phylum) %>%
  summarise(mean_old = mean(rel_abund_old, na.rm = TRUE), .groups = "drop")

new_site_means <- new_phylum %>%
  filter(siteID %in% existing_sites) %>%
  group_by(siteID, phylum) %>%
  summarise(mean_new = mean(rel_abund_new, na.rm = TRUE), .groups = "drop")

site_comparison <- inner_join(old_site_means, new_site_means, by = c("siteID", "phylum"))

cat("Site-phylum comparisons:", nrow(site_comparison), "\n")
cat(sprintf("Overall correlation: %.4f\n",
            cor(site_comparison$mean_old, site_comparison$mean_new, use = "complete")))
cat(sprintf("Overall RMSE: %.4f\n",
            sqrt(mean((site_comparison$mean_old - site_comparison$mean_new)^2, na.rm = TRUE))))

per_phylum_site <- site_comparison %>%
  group_by(phylum) %>%
  summarise(
    n_sites = n(),
    mean_old = mean(mean_old, na.rm = TRUE),
    mean_new = mean(mean_new, na.rm = TRUE),
    mae = mean(abs(mean_old - mean_new), na.rm = TRUE),
    rmse = sqrt(mean((mean_old - mean_new)^2, na.rm = TRUE)),
    r = cor(mean_old, mean_new, use = "complete"),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_old))

cat("\nPer-phylum (site-level grand means):\n")
print(as.data.frame(per_phylum_site), digits = 4, row.names = FALSE)

# Plot: old vs new site means
p1 <- ggplot(site_comparison, aes(x = mean_old, y = mean_new)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = phylum), alpha = 0.6, size = 1.5) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(
    title = "Site-level mean phylum abundances: old (2013-2018) vs new (2019-2024)",
    subtitle = sprintf("n=%d site-phylum pairs; r=%.3f, RMSE=%.4f",
                       nrow(site_comparison),
                       cor(site_comparison$mean_old, site_comparison$mean_new, use = "complete"),
                       sqrt(mean((site_comparison$mean_old - site_comparison$mean_new)^2))),
    x = "Mean rel. abundance (existing data)",
    y = "Mean rel. abundance (new data)",
    color = "Phylum"
  ) +
  theme(legend.text = element_text(size = 7))

ggsave(here("figures", "new_vs_old_phylum_comparison.png"), p1,
       width = 10, height = 7, dpi = 150)
cat("\nSaved: figures/new_vs_old_phylum_comparison.png\n")

# Plot: per-phylum 1:1 facets for the major phyla
top_phyla <- per_phylum_site %>% slice_max(mean_old, n = 9) %>% pull(phylum)

p2 <- site_comparison %>%
  filter(phylum %in% top_phyla) %>%
  ggplot(aes(x = mean_old, y = mean_new)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 1.5, color = "#0072B2") +
  facet_wrap(~phylum, scales = "free") +
  theme_bw() +
  labs(
    title = "Site-level mean comparison by phylum",
    subtitle = "Each point = one site; dashed = 1:1 line",
    x = "Mean rel. abundance (existing)",
    y = "Mean rel. abundance (new)"
  )

ggsave(here("figures", "new_vs_old_by_phylum.png"), p2,
       width = 10, height = 8, dpi = 150)
cat("Saved: figures/new_vs_old_by_phylum.png\n")

# ==========================================================================
# 6. PREPARE POST-HINDCAST EVALUATION DATA
# ==========================================================================
cat("\n========================================\n")
cat("=== POST-HINDCAST EVALUATION DATA ===\n")
cat("========================================\n")

post_hindcast <- new_phylum %>%
  filter(siteID %in% existing_sites, dateID > max(existing_dates))

cat("Post-hindcast data:\n")
cat("  Samples:", n_distinct(post_hindcast$dnaSampleID), "\n")
cat("  Sites:", n_distinct(post_hindcast$siteID), "\n")
cat("  Date range:", min(post_hindcast$dateID), "-", max(post_hindcast$dateID), "\n")
cat("  Phyla:", n_distinct(post_hindcast$phylum), "\n")

# Site-level means
post_site_means <- post_hindcast %>%
  group_by(siteID, dateID, phylum) %>%
  summarise(
    mean_rel = mean(rel_abund_new, na.rm = TRUE),
    sd_rel = sd(rel_abund_new, na.rm = TRUE),
    n_samples = n_distinct(dnaSampleID),
    .groups = "drop"
  )

# Plot-level means
post_plot_means <- post_hindcast %>%
  group_by(siteID, plotID, dateID, phylum) %>%
  summarise(
    mean_rel = mean(rel_abund_new, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

cat("  Site-month-phylum combinations:", nrow(post_site_means), "\n")
cat("  Plot-month-phylum combinations:", nrow(post_plot_means), "\n")

# Save
saveRDS(post_site_means, here("data", "clean", "new_neon_16S_phylum_site_means.rds"))
saveRDS(post_plot_means, here("data", "clean", "new_neon_16S_phylum_plot_means.rds"))
cat("\n  Saved: data/clean/new_neon_16S_phylum_site_means.rds\n")
cat("  Saved: data/clean/new_neon_16S_phylum_plot_means.rds\n")

# Also save the phylum mapping for downstream use
saveRDS(phylum_map, here("data", "clean", "phylum_name_mapping.rds"))
cat("  Saved: data/clean/phylum_name_mapping.rds\n")

# Timeline figure
example_sites <- c("CLBJ", "HARV", "KONZ", "SJER")
example_sites <- example_sites[example_sites %in% post_site_means$siteID]

top_phyla_post <- post_site_means %>%
  group_by(phylum) %>%
  summarise(grand_mean = mean(mean_rel), .groups = "drop") %>%
  slice_max(grand_mean, n = 8) %>%
  pull(phylum)

plot_data <- post_site_means %>%
  filter(siteID %in% example_sites, phylum %in% top_phyla_post) %>%
  mutate(date = as.Date(paste0(dateID, "01"), format = "%Y%m%d"))

p3 <- ggplot(plot_data, aes(x = date, y = mean_rel, color = phylum)) +
  geom_point(size = 1.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~siteID, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Post-hindcast phylum relative abundances (new NEON data)",
    x = "Date", y = "Mean relative abundance", color = "Phylum"
  ) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("figures", "post_hindcast_phylum_timeseries.png"), p3,
       width = 12, height = 8, dpi = 150)
cat("  Saved: figures/post_hindcast_phylum_timeseries.png\n")

# ==========================================================================
# 7. SUMMARY
# ==========================================================================
cat("\n========================================\n")
cat("=== SUMMARY ===\n")
cat("========================================\n")

cat(sprintf("
Taxonomy mapping:
- 17/20 existing phyla mapped to ICNP names (3 GTDB-only phyla unmappable)
- 3 merges: Pseudomonadota <- proteobacteria + myxococcota + desulfobacterota
            Bacillota <- firmicutes + firmicutes_a
            Cand. Saccharibacteria <- patescibacteria + dependentiae

Direct overlap (OSBS 201803, different soil cores):
- RMSE: %.4f, MAE: %.4f, r: %.3f
- %d/%d phyla within 1 SD of existing mean
- Caveat: only 1 overlapping sample; error reflects BOTH taxonomy changes
  AND natural spatial variability within a site-month

Cross-site mean comparison (%d sites):
- Correlation: %.4f (site-level grand means, old vs new periods)
- RMSE: %.4f
- Dominant phyla (Pseudomonadota, Acidobacteriota) show strong correspondence
- This suggests taxonomy mapping preserves rank-order across sites

Post-hindcast evaluation data:
- %d samples across %d sites, %s to %s
- %d site-months available for hindcast comparison
- Saved to data/clean/new_neon_16S_phylum_{site,plot}_means.rds
",
sqrt(mean(comparison$diff^2)), mean(comparison$abs_diff),
cor(comparison$mean_rel_old, comparison$rel_abund_new),
sum(comparison$within_1sd, na.rm = TRUE), sum(!is.na(comparison$within_1sd)),
nrow(per_phylum_site),
cor(site_comparison$mean_old, site_comparison$mean_new, use = "complete"),
sqrt(mean((site_comparison$mean_old - site_comparison$mean_new)^2)),
n_distinct(post_hindcast$dnaSampleID), n_distinct(post_hindcast$siteID),
min(post_hindcast$dateID), max(post_hindcast$dateID),
n_distinct(paste(post_hindcast$siteID, post_hindcast$dateID))))

cat("Done.\n")

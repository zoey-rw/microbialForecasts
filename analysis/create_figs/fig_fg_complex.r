# Fungal functional guild phenology aligned to plant phenophase (standalone)
# Shows modeled abundance of fungal guilds across the plant growing cycle
# (dormancy -> green-up -> peak -> senescence), averaged across all sites.
# Also appears as Panel D in combined_precision_significance_fg_source.

source("source.R")
library(ggrepel)

# ── Data loading ─────────────────────────────────────────────────────────────
scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged <- scores_list$converged_list

# Element [[6]] is the full monthly view (every model x site x month modeled
# estimate with assigned phenophase). Using [[4]] here would average only the
# per-site-year peak month within each phenophase, which is what the model
# placed an annual maximum into — not what abundance actually does across the
# growing season.
pheno_data <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))[[6]]

fungal_guilds <- c("saprotroph", "ectomycorrhizal", "plant_pathogen",
                   "animal_pathogen")

# Restrict to env_cycl so dormancy predictions are anchored to year-round soil
# temperature and moisture sensors, not pure sinusoidal extrapolation from the
# Apr-Oct-dominated NEON sampling window.
guild_data <- pheno_data %>%
  filter(taxon %in% fungal_guilds,
         model_id %in% converged,
         model_name == "env_cycl") %>%
  mutate(pretty_name = recode(taxon, !!!microbialForecast:::pretty_names))

cat("Guild data:", nrow(guild_data), "observations across",
    length(unique(guild_data$siteID)), "sites\n")

# ── Aggregate: mean abundance per guild x phenophase, min-max scaled ─────────
guild_pheno <- guild_data %>%
  group_by(taxon, pretty_name, sampling_season) %>%
  summarise(mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  group_by(taxon) %>%
  mutate(scaled = (mean_abun - min(mean_abun)) / (max(mean_abun) - min(mean_abun))) %>%
  ungroup()

season_labels <- c(dormancy = "Dormancy", greenup = "Green-up",
                   peak = "Peak", greendown = "Senescence")
guild_pheno <- guild_pheno %>%
  mutate(season_label = factor(season_labels[as.character(sampling_season)],
                               levels = season_labels))

# ── Colors, linetypes, and direct labels ─────────────────────────────────────
guild_colors <- c(
  "Saprotrophs"         = "#E69F00",
  "Ectomycorrhizae"     = "#0072B2",
  "Plant pathogens"     = "#009E73",
  "Animal pathogens"    = "#D55E00"
)
guild_linetypes <- c(
  "Saprotrophs"         = "solid",
  "Ectomycorrhizae"     = "dashed",
  "Plant pathogens"     = "dotdash",
  "Animal pathogens"    = "dotted"
)

phenophase_fills <- c(
  "Dormancy"   = "grey85",
  "Green-up"   = "#009E73",
  "Peak"       = "#E69F00",
  "Senescence" = "#D55E00"
)

# Labels at right end of each line
label_data <- guild_pheno %>% filter(season_label == "Senescence")

# ── Figure ───────────────────────────────────────────────────────────────────
fig <- ggplot(guild_pheno, aes(x = season_label, y = scaled,
                               color = pretty_name, linetype = pretty_name,
                               group = pretty_name)) +
  geom_rect(data = data.frame(
    season_label = factor(season_labels, levels = season_labels),
    fill_label = names(phenophase_fills)
  ),
  aes(xmin = as.numeric(season_label) - 0.5,
      xmax = as.numeric(season_label) + 0.5,
      ymin = -Inf, ymax = Inf, fill = fill_label),
  inherit.aes = FALSE, alpha = 0.12) +
  scale_fill_manual(values = phenophase_fills, guide = "none") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text_repel(data = label_data,
                  aes(label = pretty_name),
                  size = 3.5, fontface = "bold",
                  direction = "y", hjust = 0, nudge_x = 0.15,
                  segment.size = 0.3, segment.color = "grey50",
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_linetype_manual(values = guild_linetypes, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.25))) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = "Plant phenophase",
       y = "Relative seasonal abundance\n(min-max scaled within guild)") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_fg_complex.png"), fig,
       width = 8, height = 5.5, dpi = 300)

cat("Saved: figures/fig_fg_complex.png\n")

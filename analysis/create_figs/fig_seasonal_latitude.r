# Latitudinal gradient in forecast accuracy
# Central question: does site latitude predict how well each taxon is forecast,
# and which taxa show the strongest positive/negative latitudinal gradients?
#
# Panel A: 3 representative taxa (most positive, flat, most negative latitude slope)
# Panel B: Spaghetti of per-taxon latitude ~ R² lines (individual + mean per kingdom)
# Panel C: Per-taxon slope vs R² scatter (which taxa drive the gradient)

library(tidyverse)
library(ggpubr)
library(patchwork)
library(broom)

source("source.R")

# ── Data loading ──────────────────────────────────────────────────────────────
scores_list    <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged_base <- gsub("_beta_regression$", "", scores_list$converged_strict_list)

fg_names       <- microbialForecast:::keep_fg_names
rank_spec_names <- microbialForecast:::rank_spec_names
pretty_names_vec <- unlist(microbialForecast:::pretty_names)
pretty_rank_vec  <- unlist(microbialForecast:::pretty_rank_names)

all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
all_fungi    <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])

# Site-level R² scores (env_cycl, new-time observed-site predictions)
site_scores <- scores_list$scoring_metrics_site_long %>%
  filter(!siteID %in% "MLBS",
         model_id %in% converged_base,
         metric == "RSQ") %>%
  # Ensure model_name, rank_name, species columns
  { if (!"model_name" %in% names(.)) {
      left_join(., scores_list$scoring_metrics_long %>%
                     select(model_id, rank_name, species, model_name) %>% distinct(),
                by = "model_id")
    } else . } %>%
  mutate(
    pretty_group = ifelse(species %in% all_bacteria, "Bacteria",
                          ifelse(species %in% all_fungi, "Fungi", NA_character_))
  ) %>%
  filter(!is.na(pretty_group))

# Merge site latitude
site_descr  <- readRDS(here("data/clean/site_effect_predictors.rds"))
site_env_cycl <- site_scores %>%
  filter(model_name == "env_cycl") %>%
  left_join(site_descr %>% select(siteID, latitude), by = "siteID") %>%
  filter(!is.na(latitude))

# ── Per-taxon latitude slopes ─────────────────────────────────────────────────
taxon_slopes <- site_env_cycl %>%
  group_by(pretty_group, model_id, species, rank_name) %>%
  filter(n() >= 10) %>%
  summarise(
    slope   = coef(lm(score ~ latitude))[["latitude"]],
    r2      = summary(lm(score ~ latitude))$r.squared,
    p_value = glance(lm(score ~ latitude))$p.value,
    n_sites = n(),
    .groups = "drop"
  ) %>%
  mutate(
    significant  = p_value < 0.1,
    pretty_species = ifelse(species %in% names(pretty_names_vec),
                            pretty_names_vec[species],
                            tools::toTitleCase(gsub("_", " ", species))),
    pretty_rank = case_when(
      rank_name %in% names(pretty_rank_vec) ~ pretty_rank_vec[rank_name],
      species %in% fg_names ~ "Functional group",
      TRUE ~ tools::toTitleCase(gsub("_", " ", rank_name))
    )
  )

# ── Shared aesthetics ─────────────────────────────────────────────────────────
kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

# ── Panel A: 3 representative taxa ───────────────────────────────────────────
taxon_slopes_good <- taxon_slopes %>% filter(n_sites >= 15)

rep_positive <- taxon_slopes_good %>% filter(significant) %>% slice_max(slope, n = 1)
rep_negative <- taxon_slopes_good %>% filter(significant) %>% slice_min(slope, n = 1)
rep_flat     <- taxon_slopes_good %>% slice_min(abs(slope), n = 1)

representative_taxa <- bind_rows(
  rep_positive %>% mutate(example = "Positive"),
  rep_negative %>% mutate(example = "Negative"),
  rep_flat     %>% mutate(example = "No relationship")
) %>%
  mutate(
    facet_label = paste0(example, "\n", pretty_species, " (", pretty_rank, ")")
  )

# Fix factor order: Positive first, then No relationship, then Negative
representative_taxa$facet_label <- factor(
  representative_taxa$facet_label,
  levels = representative_taxa$facet_label[order(
    match(representative_taxa$example, c("Positive", "No relationship", "Negative"))
  )]
)

rep_site_data <- site_env_cycl %>%
  filter(model_id %in% representative_taxa$model_id) %>%
  left_join(representative_taxa %>% select(model_id, facet_label), by = "model_id") %>%
  mutate(facet_label = factor(facet_label, levels = levels(representative_taxa$facet_label)))

panel_a <- ggplot(rep_site_data, aes(x = latitude, y = score)) +
  geom_point(aes(color = pretty_group), alpha = 0.5, size = 2.5,
             position = position_jitter(height = 0, width = 0.3)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~")),
           size = 3.2, label.x.npc = 0.02, label.y.npc = 0.95, color = "black") +
  facet_wrap(~facet_label, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = seq(30, 70, by = 10)) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(x = "Latitude (\u00b0N)", y = expression(R^2)) +
  base_theme + theme(legend.position = "none")

# ── Panel B: Spaghetti of per-taxon lines ─────────────────────────────────────
panel_b <- ggplot(site_env_cycl, aes(x = latitude, y = score)) +
  geom_smooth(aes(group = model_id, color = pretty_group),
              method = "lm", se = FALSE, linewidth = 0.3, alpha = 0.2) +
  geom_smooth(aes(color = pretty_group, fill = pretty_group),
              method = "lm", linewidth = 1.5, alpha = 0.15) +
  scale_x_continuous(breaks = seq(30, 70, by = 10)) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  scale_fill_manual(values  = kingdom_colors, guide = "none") +
  labs(x = "Latitude (\u00b0N)", y = expression(Site~R^2)) +
  base_theme +
  theme(legend.position = c(0.12, 0.82),
        legend.background = element_rect(fill = "white", color = NA))

# ── Panel C: Slope vs R² scatter ──────────────────────────────────────────────
taxon_slopes <- taxon_slopes %>%
  mutate(is_representative = model_id %in% representative_taxa$model_id)

panel_c <- ggplot(taxon_slopes, aes(x = slope, y = r2, color = pretty_group)) +
  geom_point(aes(shape = significant), size = 1.8, alpha = 0.6) +
  geom_point(data = taxon_slopes %>% filter(is_representative),
             size = 4, shape = 21, fill = NA, color = "black", stroke = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.1", "FALSE" = "n.s."),
                     name = "Significance") +
  labs(x = expression(Slope~of~R^2~vs.~latitude),
       y = expression(R^2~of~latitude~model)) +
  base_theme +
  theme(legend.position = c(0.28, 0.80),
        legend.background = element_rect(fill = "white", color = NA),
        legend.spacing.y  = unit(0.1, "cm")) +
  guides(color = "none")

# ── Combine and save ──────────────────────────────────────────────────────────
fig_lat_grad <- panel_a / (panel_b | panel_c) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  ) +
  plot_layout(heights = c(0.8, 1))

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS9_latitude_vs_accuracy.png"), fig_lat_grad,
       width = 13, height = 9, dpi = 200)
cat("Saved: figures/figS9_latitude_vs_accuracy.png\n")

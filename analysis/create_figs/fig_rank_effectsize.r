# Effect sizes by taxonomic rank
# Shows how environmental predictor effect sizes vary across taxonomic ranks
# (phylum, class, order, family, functional group) for fungi and bacteria.
# Uses Tukey HSD to annotate significant rank differences.

source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, ggpubr, rstatix, agricolae)

# ── Data loading ─────────────────────────────────────────────────────────────
sum.all <- readRDS(here("data", "summary/predictor_effects.rds"))
seasonal_amplitude_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))

beta_names <- c("Ectomycorrhizal\ntrees", "LAI", "pC",
                "pH", "Temperature", "Moisture")

# ── Prepare plotting data ───────────────────────────────────────────────────
plotting_df <- sum.all %>%
  filter(beta %in% beta_names,
         model_name == "env_cycl",
         !grepl("other", taxon),
         time_period == "2013-06_2018-01")

# ── Tukey HSD tests per group x predictor ───────────────────────────────────
tukey_list <- list()

for (pg in c("Bacteria", "Fungi")) {
  group_df <- plotting_df %>% filter(pretty_group == !!pg)
  tukey_group_list <- list()

  for (b in beta_names) {
    df <- group_df %>% filter(beta == b)
    new.df <- data.frame(x = df$rank_only, y = df$effSize) %>%
      filter(!is.na(y))
    abs_max <- max(new.df$y, na.rm = TRUE)
    maxs <- new.df %>%
      group_by(x) %>%
      summarise(tot = max(y, na.rm = TRUE) + 0.2 * abs_max, .groups = "drop")

    unique_groups <- unique(new.df$x)
    if (length(unique_groups) < 2 || any(table(new.df$x) < 2)) {
      Tukey_test <- maxs %>%
        mutate(Letters_Tukey = "a") %>%
        rename("rank_only" = "x")
    } else {
      Tukey_test <- tryCatch({
        aov(y ~ x, data = new.df) %>%
          agricolae::HSD.test("x", group = TRUE) %>%
          .$groups %>%
          as_tibble(rownames = "x") %>%
          rename("Letters_Tukey" = "groups") %>%
          select(-y) %>%
          left_join(maxs, by = "x") %>%
          rename("rank_only" = "x")
      }, error = function(e) {
        maxs %>%
          mutate(Letters_Tukey = "a") %>%
          rename("rank_only" = "x")
      })
    }
    Tukey_test$beta <- b
    Tukey_test$pretty_group <- pg
    tukey_group_list[[b]] <- Tukey_test
  }
  tukey_list[[pg]] <- data.table::rbindlist(tukey_group_list)
}

tukey_df <- data.table::rbindlist(tukey_list)

# ── Colorblind-friendly palette ─────────────────────────────────────────────
# kingdom_colors comes from source.R

# ── Main figure: effect sizes by rank, faceted by predictor x kingdom ────────
p <- ggplot(plotting_df, aes(x = rank_only, y = effSize, color = pretty_group)) +
  geom_jitter(aes(shape = as.factor(significant)),
              width = 0.1, height = 0, size = 2, alpha = 0.2, show.legend = FALSE) +
  geom_violin(draw_quantiles = 0.5, show.legend = FALSE) +
  geom_text(data = tukey_df,
            aes(x = rank_only, y = tot, label = Letters_Tukey),
            show.legend = FALSE, color = "black", size = 3.5) +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "anova", size = 3.5, label.y.npc = 0.5) +
  facet_grid(cols = vars(beta), rows = vars(pretty_group),
             drop = TRUE, scales = "free", space = "free") +
  scale_color_manual(values = kingdom_colors) +
  scale_shape_manual(values = c("0" = 1, "1" = 16),
                     labels = c("Not significant", "Significant")) +
  labs(x = "Taxonomic rank", y = "Absolute effect size") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 320, vjust = 1, hjust = -0.05),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  guides(color = "none")

# ── Save ────────────────────────────────────────────────────────────────────
out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS2_parameter_estimates_rank.png"), p,
       width = 14, height = 8, dpi = 200)

cat("Saved: figures/figS2_parameter_estimates_rank.png\n")

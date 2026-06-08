# Microbial phenology strength figure
# Shows that env_cycl models detect more robust seasonal patterns
# than cycl_only models.
# Panel A: % of taxa with significant seasonality.
# Panel B: Paired amplitude difference (env_cycl minus cycl_only)
#          per taxon, showing env_cycl amplitudes are significantly
#          larger (paired Wilcoxon p < 0.01 in both groups).

library(tidyverse)
library(ggpubr)

source("source.R")

# ── Data loading ───────────────────────────────────────────
phenophase_in <- readRDS(
  here("data/clean/pheno_group_peak_phenophases.rds")
)
mode_data <- phenophase_in[[1]]

# ── Display settings ───────────────────────────────────────
model_labels <- c(
  cycl_only = "Seasonality only",
  env_cycl  = "Seasonality + environment"
)
model_colors <- c(
  "Seasonality + environment" = "#009E73",
  "Seasonality only"          = "#0072B2"
)
fcast_labels <- c(
  Functional = "Functional groups",
  Taxonomic  = "Taxonomic groups"
)

# ── Panel A: % of taxa with significant seasonality ────────
pct_sig <- mode_data %>%
  filter(model_name %in% c("cycl_only", "env_cycl")) %>%
  group_by(model_name, fcast_type, pretty_group) %>%
  summarise(
    total   = n(),
    n_sig   = sum(significant_sin == 1 |
                    significant_cos == 1, na.rm = TRUE),
    pct_sig = n_sig / total * 100,
    .groups = "drop"
  ) %>%
  mutate(
    model_name = recode(model_name, !!!model_labels),
    fcast_type = recode(fcast_type, !!!fcast_labels)
  )

pA <- ggplot(pct_sig,
             aes(x = pretty_group, y = pct_sig,
                 fill = model_name)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65) +
  geom_text(aes(label = paste0("n=", n_sig, "/", total)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2, color = "grey30") +
  facet_wrap(~fcast_type) +
  scale_fill_manual(values = model_colors, name = "Model") +
  scale_y_continuous(
    limits = c(0, 115),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL,
       y = "% of groups with significant seasonality") +
  theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey92",
                                      color = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "top",
    panel.grid.major.x = element_blank()
  )

# ── Panel B: paired amplitude difference ───────────────────
# Each point is one taxon: env_cycl amplitude minus cycl_only
amp_wide <- mode_data %>%
  filter(model_name %in% c("cycl_only", "env_cycl"),
         !is.na(amplitude)) %>%
  select(taxon, fcast_type, pretty_group,
         model_name, amplitude) %>%
  pivot_wider(names_from = model_name,
              values_from = amplitude) %>%
  filter(!is.na(cycl_only) & !is.na(env_cycl)) %>%
  mutate(
    diff = env_cycl - cycl_only,
    fcast_type = recode(fcast_type, !!!fcast_labels)
  )

# Per-facet Wilcoxon p-values
p_labels <- amp_wide %>%
  group_by(fcast_type) %>%
  summarise(
    p_val = wilcox.test(env_cycl, cycl_only,
                        paired = TRUE)$p.value,
    n_higher = sum(diff > 0),
    n_total  = n(),
    med_diff = median(diff),
    .groups  = "drop"
  ) %>%
  mutate(
    p_text = ifelse(p_val < 0.001,
                    format(p_val, digits = 2,
                           scientific = TRUE),
                    paste0("p = ", round(p_val, 3))),
    label = paste0(n_higher, "/", n_total,
                   " taxa higher\n", p_text)
  )

pB <- ggplot(amp_wide,
             aes(x = pretty_group, y = diff,
                 shape = pretty_group)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50") +
  geom_boxplot(aes(group = pretty_group),
               outlier.shape = NA, width = 0.5,
               fill = "grey90", alpha = 0.5) +
  geom_jitter(aes(color = diff > 0),
              width = 0.15, alpha = 0.5, size = 2) +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.3, hjust = 0.5, size = 3.3,
    color = "grey30"
  ) +
  facet_wrap(~fcast_type) +
  scale_color_manual(
    values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"),
    labels = c("TRUE" = "Higher with env.",
               "FALSE" = "Lower with env."),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(Bacteria = 16, Fungi = 17),
    guide = "none"
  ) +
  labs(x = NULL,
       y = expression(
         Delta ~ "seasonal amplitude"
         ~ "(env_cycl" ~ "-" ~ "cycl_only)"
       )) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92",
                                    color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    panel.grid.major.x = element_blank()
  )

# ── Combine and save ──────────────────────────────────────
fig_pheno <- ggarrange(
  pA, pB,
  ncol = 2, labels = c("A", "B"),
  widths = c(1, 1)
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS8_abundance_error_amplitude.png"),
       fig_pheno, width = 12, height = 5.5, dpi = 200)

cat("Saved: figures/figS8_abundance_error_amplitude.png\n")

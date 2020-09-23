## code to create sample plots for N-cycling gene abundance and activity ##
rm(list = ls())
library(ggplot2)
library(gridExtra)

# let's say we have 6 sites, each with 10 plots
set.seed(5)
sites <- LETTERS[1:6]

site_genes_activity <- data.frame()
for (i in 1:length(sites)) {
  site <- sites[i]
  
  # simulating data for nitrification
  site_amoA_avg <- rnorm(1, 10, 2) # pick some number that the site will center on
  plot_amoA_counts <- rnorm(10, site_amoA_avg, 1) # generate plot-level gene counts
  nitr_rates <- sapply(plot_amoA_counts, function(x)
                         rnorm(1, 2 * x + 3, 3)) # create activity rates based on gene counts
  
  # simulating data for n-mineralization
  site_chiA_avg <- rnorm(1, 10, 2) # pick some number that the site will center on
  plot_chiA_counts <- rnorm(10, site_chiA_avg, 1) # generate plot-level gene counts
  nmin_rates <- sapply(plot_chiA_counts, function(x)
                         rnorm(1, 4 * x + 2, 5)) # create activity rates based on gene counts
  
  site_vals <- data.frame(
    site = factor(site),
    amoA_count = plot_amoA_counts,
    chiA_count = plot_chiA_counts,
    nmin_rate = nmin_rates,
    nitr_rate = nitr_rates
  )  # wrap into data.frame!
  site_genes_activity <-
    rbind(site_genes_activity, site_vals) # add to main data.frame
}

fit_nitr <- summary(lm(nitr_rate ~ amoA_count, site_genes_activity)) # linear regression
r2_nitr <- round(fit_nitr$r.squared, 3) # pull out R2 value and round to 3 decimal places
p_nitr <- round(fit_nitr$coefficients[, 4][2], 3) # same for p-value
nitr_plot <-  ggplot(data = site_genes_activity,
                      aes(x = amoA_count, y = nitr_rate, col = site)) +
                      geom_point(size = 3) +
                      ylab("Nitrification rate") + xlab("amoA gene count") +
                      geom_text(aes(x = -Inf, y = Inf, hjust = -.1, vjust = 2, label = paste0("R2 = ", r2_nitr, ", p = ", p_nitr)), size = 8, show.legend = FALSE) +
                      ggtitle("Nitrification rates vs functional genes") +
                      theme(text = element_text(size = 18))

fit_nmin <-
  summary(lm(nmin_rate ~ chiA_count, site_genes_activity)) # linear regression
r2_nmin <-
  round(fit_nmin$r.squared, 3) # pull out R2 value and round to 3 decimal places
p_nmin <- round(fit_nmin$coefficients[, 4][2], 3) # same for p-value
nmin_plot <-  ggplot(data = site_genes_activity, 
                      aes(x = chiA_count, y = nmin_rate, col = site, show.legend = FALSE)) +
                      geom_point(size = 3) +
                      ylab("N mineralization rate") + xlab("chiA gene count") +
                      geom_text(aes(x = -Inf, y = Inf, hjust = -.1, vjust = 2, label = paste0("R2 = ", r2_nmin, ", p = ", p_nmin)), size = 8, show.legend = FALSE) +
                      ggtitle("N-mineralization rates vs functional genes") +
                      theme(text = element_text(size = 18))

grid.arrange(
  nitr_plot + theme(legend.position = "none"),
  nmin_plot,
  ncol = 2, widths = c(1.7, 2)
)

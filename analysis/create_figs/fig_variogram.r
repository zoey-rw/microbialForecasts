source("source.R")

library(scales)
library(spgs) # for uniformity test

# run siteEffectVariogram.r to generate
variograms = readRDS(here("data/summary/site_effect_variograms.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

# look at example
variograms[[2]][12]

sig_results = variograms[[1]]

sig_results_long = sig_results %>%
	pivot_longer(cols=5:6) %>%
	mutate(value = as.numeric(value)) %>%
	filter(model_id %in% converged_strict)

var_plot = ggplot(sig_results_long) +
	#geom_point(position = position_jitterdodge(jitter.width = .05, jitter.height = 0)) +
	geom_density(aes(x = value), show.legend = F) +
	facet_grid(rows=vars(name), scales="free") +
	geom_vline(xintercept = .05, color=2, linetype=2) +
	theme_classic(base_size = 20)+
	ylab("Variogram p-value") +
	xlab(NULL) +
	ggtitle("Site effect autocorrelation is mostly captured by PLSR") +
	#annotate("text", x = .1, y = 1.5, label = "p.value < .05 indicates \nspatial autocorrelation")  +
	geom_rug(aes(x = value), alpha=.3,length= unit(0.06, "npc"), show.legend = F) + xlab("variogram p-value")

var_plot <- tag_facet(var_plot)

# Neither set of values follows a uniform distribution
chisq.unif.test(sig_results_long[sig_results_long$name=="site effect",]$value)
# X-squared = 387.23, df = 2, a = 0, b = 1, p-value < 2.2e-16
chisq.unif.test(sig_results_long[sig_results_long$name=="site effect residuals",]$value)
# X-squared = 39.581, df = 17, a = 0, b = 1, p-value = 0.001482

# Non-uniform residuals are driven by env_cov model (without seasonality predictors)
sig_results_long %>% filter(name=="site effect residuals" & model_name=="cycl_only") %>%
	select(value) %>% unlist() %>% chisq.unif.test()
#X-squared = 7.1852, df = 4, a = 0, b = 1, p-value = 0.1264
sig_results_long %>% filter(name=="site effect residuals" & model_name=="env_cov") %>%
	select(value) %>% unlist() %>% chisq.unif.test()
# X-squared = 11.27, df = 4, a = 0, b = 1, p-value = 0.02369
sig_results_long %>% filter(name=="site effect residuals" & model_name=="env_cycl") %>%
	select(value) %>% unlist() %>% chisq.unif.test()
# X-squared = 6.9259, df = 5, a = 0, b = 1, p-value = 0.2262

png(here("figures","variogram_p_val.png"), width = 800, height=1000)
print(var_plot)
dev.off()

png(here("figures","variogram_example.png"), width = 600, height=800)
variograms[[2]][12]
dev.off()


by_model = ggplot(sig_results_long) +
	#geom_point(position = position_jitterdodge(jitter.width = .05, jitter.height = 0)) +
	geom_density(aes(x = value), show.legend = F) +
	facet_grid(rows=vars(name), cols=vars(model_name), scales="free") +
	geom_vline(xintercept = .05, color=2, linetype=2) +
	theme_classic(base_size = 20)+
	ylab("Variogram p-value") +
	xlab(NULL) +
	ggtitle("Site effect autocorrelation is mostly captured by PLSR") +
	#annotate("text", x = .1, y = 1.5, label = "p.value < .05 indicates \nspatial autocorrelation")  +
	geom_rug(aes(x = value), alpha=.3,length= unit(0.06, "npc"), show.legend = F) + xlab("variogram p-value")

png(here("figures","variogram_p_val_by_model.png"), width = 1200, height=1000)
print(by_model)
dev.off()

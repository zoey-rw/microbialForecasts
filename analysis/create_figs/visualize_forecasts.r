library(ggforce)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data_construction/microbe/prep_model_data.r")
# view output of functional group forecasts.

dat1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_data_coh_multi_its.rds")
dat2 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_data_no_coh.rds")
dat3 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_data_coh_16s.rds")
dat4 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_multi_no.rds")
dat4 <- dat4 %>% filter(rank == "No_cohesion")
plot.allrank <- data.table::rbindlist(list(dat1, dat2, dat3, dat4))

# temp file
#plot.allrank <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_multi_no.rds")
#plot.allrank <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/forecast_plotting_data.rds")
plot.subset <- plot.allrank

plot.name <- "CPER_002"
plot.name <- "OSBS_001"

# acetogen_anaerobic, alkaline_stress, denitrification, erythromycin_antibiotic, heat_stress, nitrification
# Ectomycorrhizal - somethign wrong with validation?
plot.allrank$pretty_rank <- ifelse(plot.allrank$rank == "No_cohesion", "No cohesion",
																	 ifelse(plot.allrank$rank == "coh_multi", "Fungal-bacterial cohesion", NA))
# Subset to group that cohesion mattered
#plot.subset <- plot.allrank[plot.allrank$group.name %in% c("propionate_simple","rhamnose_simple"),]
#plot.subset <- plot.allrank[plot.allrank$group.name %in% c("streptomycin_antibiotic"),]
#plot.subset <- plot.allrank[plot.allrank$rank %in% c("No_cohesion","coh_multi") & plot.allrank$group.name %in% c("acetogen_anaerobic", "alkaline_stress", "denitrification", "erythromycin_antibiotic", "heat_stress", "nitrification"),]

plot.subset <- plot.allrank[plot.allrank$rank %in% c("No_cohesion","coh_multi") & plot.allrank$group.name %in% c("erythromycin_antibiotic", "nitrification"),]
plot.subset[plot.subset$group.name=="erythromycin_antibiotic",]$pretty_name <- "Erythromycin-resistant"
#plot.subset <- plot.allrank[plot.allrank$group.name %in% c("light_stress"),]


# WITH validation points
ggplot(plot.subset, aes(x = time, y = value)) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2) +
  geom_point(aes(color = data.type), size = 4) + 
  labs(x = "", y = "Abundance (CLR-transformed)", title = paste0("Microbial abundance at NEON plot: ", plot.name)) +
  #  facet_wrap(~group.name, scales="free") +
  #facet_wrap(~rank, scales="free") +
  #theme_classic(base_size = 28) + 
  theme_minimal(base_size = 26) + 
  theme(legend.title = element_blank()) + theme(panel.grid.major = element_blank(), 
  																							panel.spacing=unit(0, "lines"))  +
	facet_wrap_paginate(pretty_name ~ pretty_rank, ncol = 2, nrow = 2, #scales = "free",
											page = 1)


### without validation points
plot.subset[plot.subset$data.type=="Validation",]$value <- NaN
ggplot(plot.subset, aes(x = time, y = value)) + 
	geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2) +
	geom_point(aes(color = data.type), size = 4) + 
	labs(x = "", y = "Abundance (CLR-transformed)", title = paste0("Microbial abundance at NEON plot: ", plot.name)) +
	#  facet_wrap(~group.name, scales="free") +
	#facet_wrap(~rank, scales="free") +
	#theme_classic(base_size = 28) + 
	theme_minimal(base_size = 26) + 
	theme(legend.title = element_blank()) + theme(panel.grid.major = element_blank(), 
																								panel.spacing=unit(0, "lines"))  +
	facet_wrap_paginate(pretty_name ~ pretty_rank, ncol = 2, nrow = 2, #scales = "free",
											page = 1)




fg_plots <- plot.allrank[[1]]
bac_phy_plots <- plot.allrank[[2]]
fun_phy_plots <- plot.allrank[[3]]


genus_plots <- c(plot.allrank[[4]], plot.allrank[[5]])


ggarrange(plotlist = fg_plots, nrow = 4, ncol = 2)
ggarrange(plotlist = bac_phy_plots, ncol = 2, nrow = 4)
ggarrange(plotlist = fun_phy_plots, ncol = 2, nrow = 4)
ggarrange(plotlist = genus_plots, ncol = 2, nrow = 4)

fg_plots[[1]] <- fg_plots[[1]] + scale_y_continuous(breaks=c(-1,-2))+ labs(y="")
fg_plots[[8]] <- fg_plots[[8]] + scale_y_continuous(breaks=c(-1,-2), limits = c( -2.2, -.8))+ labs(y="")
fg_plots[[16]] <- fg_plots[[16]] + labs(title = paste0("Plant pathogens at plot: CPER_001")) + labs(y="")
genus_plots[[1]] <- genus_plots[[1]] + labs(y="")
genus_plots[[17]] <- genus_plots[[17]] + labs(y="") + 
plots_for_mike <- list(fg_plots[[1]], fg_plots[[8]], fg_plots[[16]], genus_plots[[1]], genus_plots[[17]])
figure <- ggarrange(plotlist = plots_for_mike, nrow=5, ncol = 1)
annotate_figure(figure, left = text_grob("Abundance (CLR-transformed)", size = 29, , rot = 90), "common.legend = TRUE")

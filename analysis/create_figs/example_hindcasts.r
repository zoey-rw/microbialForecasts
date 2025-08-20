# Create example figure with one hindcast from each type
library(tidyverse)
source("source.R")

tax_hindcast <- readRDS(here("data/summary/hindcast_tax.rds"))
tax_hindcast$dates <- fixDate(tax_hindcast$dateID)
div_hindcast <- readRDS(here("data/summary/hindcast_div.rds"))
fg_hindcast <- readRDS(here("data/summary/hindcast_fg.rds"))

summaries_all <- readRDS(here("data/summary/all_fcast_effects.rds"))
summaries_examples <- summaries_all %>% filter(fcast_type == "Taxonomic" & taxon == "actinobacteriota" |
																							 	fcast_type == "Diversity" &	scenario == "full_uncertainty_16S" |
																							 	fcast_type == "Functional group" & taxon == "plant_pathogen") %>%
	filter(!is.na(beta))
summaries_examples$example <- recode(summaries_examples$fcast_type,
																		 Taxonomic = "Phylum: Actinobacteriota",
																		 Diversity = "Bacterial diversity",
																		 `Functional group` = "Plant pathogens")

tax_hindcast_example <- tax_hindcast %>% filter(taxon_name == "actinobacteriota") %>% 
	mutate(example = "Phylum: Actinobacteriota", truth = as.numeric(truth)) 
div_hindcast_example <- div_hindcast %>% filter(scenario == "full_uncertainty_16S") %>% 
	mutate(example = "Bacterial diversity")  
fg_hindcast_example <- fg_hindcast %>% filter(scenario == "full_uncertainty" & 
																								taxon_name == "plant_pathogen") %>%
	mutate(example = "Plant pathogens") 

example_fcasts_allplots  <- data.table::rbindlist(list(tax_hindcast_example, div_hindcast_example, fg_hindcast_example), fill = T)

example_fcasts_oneplot <- example_fcasts_allplots %>% filter(plotID %in% c("CPER_003"))

#example_fcasts_oneplot <- example_fcasts_oneplot %>% filter(!is.na(example_fcasts_oneplot$))

example_fcasts_oneplot$lo <- ifelse(!is.na(example_fcasts_oneplot$`2.5%`),
																		example_fcasts_oneplot$`2.5%`, example_fcasts_oneplot$lo)
example_fcasts_oneplot$med <- ifelse(!is.na(example_fcasts_oneplot$`50%`),
																		example_fcasts_oneplot$`50%`, example_fcasts_oneplot$med)
example_fcasts_oneplot$hi <- ifelse(!is.na(example_fcasts_oneplot$`97.5%`),
																		example_fcasts_oneplot$`97.5%`, example_fcasts_oneplot$hi)

fcasts <-  ggplot(example_fcasts_oneplot) +
	facet_grid(rows=vars(example), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6,  na.rm = F) +
	labs(title = "Example hindcasts at CPER_003 (calibration: 2013-2016)") + 
	theme_bw()+
	scale_fill_brewer(palette = "Paired") + 
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"), 
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) + 
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') 
fcasts

effects <- ggplot(data=summaries_examples,
			 aes(x = beta,y = Mean, color = fcast_period)) +
	geom_point(size = 5, 	alpha=0.6,
		show.legend = T) +
	labs(title = "Effect size") + 
	xlab("Predictor")+ 
	ylab(NULL) + scale_y_continuous(breaks = seq(-.5, .5, by = .1)) +
	scale_color_brewer(palette = "Paired", labels=c("calibration","full dataset")) + 
	geom_hline(yintercept = 0, linetype=2) +
	facet_grid(rows = vars(example),# cols = vars(group), 
						 scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(
		text = element_text(size = 14),
		axis.text.x=element_text(
			angle = 320, vjust=1, hjust = -0.05),
			plot.margin = unit(c(.2, .2, .2, .2), "cm"),
		axis.title=element_text(size=16,face="bold"))  + labs(color='')  
effects

library(ggpubr)
fig <- ggarrange(fcasts, effects, ncol=2, widths = c(2,1))
fig






not_na <- example_fcasts_allplots[!is.na(example_fcasts_allplots$truth),]
top_plots <- names(tail(sort(table(not_na$plotID)), 15))
fg_top_plots <- fg_hindcast_example[fg_hindcast_example$plotID %in% top_plots,]

ggplot(fg_top_plots) +
	#	ggplot(out.data[out.data$group=="16S",]) +
	
	facet_grid(rows=vars(plotID), #cols=vars(scenario), 
						 drop=T, scales="free"#, space="free"
	) +
	geom_line(aes(x = dates, y = mean), show.legend = F) + 
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), 
							#fill = "darkblue", 
							alpha=0.4#, show.legend = F
	) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Example hindcasts at CPER_003 (calibration: 2013-2016)")) + #scale_x_date() +	
	#theme_minimal(base_size=20) +
	theme_bw()+
	#scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab(NULL) + 
	geom_point(aes(x = dates, y = as.numeric(truth)))

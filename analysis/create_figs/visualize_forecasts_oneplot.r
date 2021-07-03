library(coda)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

d1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_10tax.rds")
out.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/allrank_summaries_10tax.rds"

rank.df <- d1[["phylum_bac"]] 

model.dat <- prepModelData(rank.df = rank.df, min.prev = 3)

# # Read in samples for visualization
read_in <- readRDS(out.path)
plot_est <- read_in$plot_est
full_summary <- read_in[[2]]

means <- full_summary[[2]]
plot_mu_out <- means[grep("plot_mu", rownames(means)),]
ypred_out <- means[grep("y", rownames(means)),]
site_eff_out <- means[grep("site_eff", rownames(means)),]
beta_out <- means[grep("beta", rownames(means)),]

# 
phyla <- plot_est[plot_est$rank=="phylum_bac",]
pred.plot1 <- phyla %>% filter(plotID == "HARV_001")
ggplot(pred.plot1, aes(x = dates)) +
	facet_wrap(~taxon, scales="free") +
	geom_line(aes(y = `50%`), show.legend = F) + #facet_wrap(~species, scales = "free") +
	geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = taxon), alpha=0.2, show.legend = F) +
	geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = taxon), alpha=0.4, show.legend = F) +
	geom_point(aes(y = truth)) + theme_ipsum(base_size = 14, strip_text_size = 22) + ggtitle(paste0("Phylum abundances at ", unique(pred.plot1$plotID))) + scale_x_date() +	scale_fill_brewer(palette = "Paired") + 
	theme(panel.spacing = unit(0, "lines"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) 




pred.plot1 <- phyla %>% filter(plotID == "HARV_001")
pred.plot1$observed <- ifelse(is.na(pred.plot1$truth), "Estimated", "Observed")
# Stacked barplot
ggplot(pred.plot1, aes(dates, `50%`, fill = taxon, color = observed)) + 
	geom_bar(stat = "identity", width = 23) + theme_light() + 
	scale_fill_brewer(palette = "Paired") + 
	scale_color_manual(values = c("darkgrey","black")) +
	ylab("Median abundance estimate") + xlab("Date") + 
	guides(colour = guide_legend(title=NULL,override.aes = list(fill="white"))) + ggtitle(paste0("Phylum abundances at ", unique(pred.plot1$plotID))) 





single_tax_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_201801.rds")
pacman::p_load(gridExtra, ggpubr)


gelman.summary <- single_tax_summaries$gelman.summary %>%
	separate(rank.name, sep = "_", into = c("only_rank", "group", NA), remove = F) %>%
	separate(rank.name, sep = "_bac_|_fun_", into = c(NA, "taxon.name"), remove = F) %>%
	mutate(rank_group = paste(only_rank, group, sep = "_"))
by.group <- gelman.summary %>%
	group_by(group, rank_group, only_rank, rank.name, taxon.name) %>%
	summarize(median_gbr = median(`Point est.`))

ggplot(by.group) + geom_jitter(aes(x = only_rank, y = median_gbr, color = group)) + ylim(c(0,10))

keep <- by.group[by.group$median_gbr <= 2,]
keep_list <- split(keep, keep$rank_group)


cycl_only <- single_tax_summaries$summary_df %>% filter(model_name == "cycl_only")
all_covariates <- single_tax_summaries$summary_df %>% filter(model_name == "all_covariates")

ggplot(all_covariates,
			 aes(x = beta,y = effSize)) +
	geom_jitter(aes(shape = as.factor(significant),
									color = rank_only), size = 5, height = 0, width=.1, alpha = .5,
							shape = 16) + ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Absolute effect size (refit)") +
	xlab("Parameter")+
	ylab(NULL) +
	facet_grid(rows = vars(fcast_type), cols = vars(group), drop = F,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL,
												 labels = c("Not significant","Significant"))


cycl_only_est <- single_tax_summaries$plot_est %>% filter(model_name == "cycl_only")
all_covariates_est <- single_tax_summaries$plot_est %>% filter(model_name == "all_covariates")
ggplot(all_covariates_est %>% filter(taxon=="clostridia"),
			 aes(x = dates,y = effSize)) + geom_point()

input_df <- all_covariates_est %>% filter(taxon=="acidobacteriota" & siteID=="HARV")
plot_model(input_df)
p1 <- plot_model(input_df, site_plots = "facet")
p1 + ylim(c(0,.6))






#' @title 			plot_model
#' @description Visualize model fits and forecasts at NEON plots and sites
#'
#' @export
siteID = "HARV"
site_plots = "overlap"
site_plots = NULL
model_type = "all"
model_type = "all_covariates"
input_df = all_covariates
taxon = "actinobacteriota"
plot_model <- function(input_df,
											 site_plots = "overlap",
											 siteID = "HARV",
											 plotID = NULL,
											 taxon = NULL,
											 model_type = "all_covariates") {

	if (model_type != "all") {
		message("Filtering to: model type: ", model_type)
		input_df <- input_df %>% filter(model_name == !!model_type)
	}


	if (!is.null(siteID)) {
		message("Filtering to: sites(s): ", siteID)
		siteID_vector <- ifelse(length(siteID) > 1,
														paste(siteID, collapse = "|"), siteID)
		input_df <- input_df %>% filter(grepl(siteID_vector, siteID))
	}

	if (!is.null(plotID)) {
		message("Filtering to: plot(s): ", plotID)
		plotID_vector <- ifelse(length(plotID) > 1,
														paste(plotID, collapse = "|"), plotID)
		input_df <- input_df %>% filter(grepl(plotID_vector, plotID))
	}

	if (!is.null(taxon)) {
		message("Filtering to: ", taxon)
		input_df <- input_df %>% filter(taxon == !!taxon)
	}

	n.taxa <- length(unique(input_df$species))
	if (n.taxa > 1) {
		message("More than one taxon/group in dataframe. Filtering to: ", input_df$species[[1]])
		taxon = input_df$species[[1]]
		input_df <- input_df %>% filter(species == !!taxon)
	}

	if (!"med" %in% colnames(input_df)) {
		input_df <- input_df %>% mutate(med = `50%`,
											 lo = `2.5%`,
											 hi = `97.5%`)
	} else {
		input_df <- input_df %>% mutate(med = ifelse(is.na(`50%`), med, `50%`),
											 lo = ifelse(is.na(`2.5%`), lo,`2.5%`),
											 hi = ifelse(is.na(`97.5%`), hi, `97.5%`))
	}

	if (!"fcast_period" %in% colnames(input_df)) {
		input_df <- input_df %>% mutate(fcast_period = "calibration")
		}

	input_df <- input_df %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4))

p <- ggplot(input_df, aes(x = dates))  +
	ylab("Abundance") +
	theme_bw() +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ggtitle(label = taxon)  +
	scale_fill_brewer(palette = "Paired")

if (!taxon %in% c("diversity_ITS","diversity_16S")) p <- p +
	scale_y_sqrt()

if (site_plots == "overlap") {
	p <- p +
		geom_ribbon(aes(ymin = lo, ymax = hi, fill = plotID), alpha=0.5) +
		geom_ribbon(aes(ymin = lo_25, ymax = hi_75, fill = plotID), alpha=0.7) +
		geom_line(aes(y = med, color = plotID, group = plotID), show.legend = F) +
		geom_jitter(aes(y = as.numeric(truth)), width = .2, height = 0) +
		xlab(NULL)

} else if (site_plots == "facet") {
	p <- p +
		geom_ribbon(aes(ymin = lo, ymax = hi, fill = plotID), alpha=0.3) +
		geom_ribbon(aes(ymin = lo_25, ymax = hi_75, fill = plotID), alpha=0.6) +
		geom_line(aes(y = med), show.legend = F) +
		facet_grid(rows = vars(plotID), drop=T, scales="free")  +
		geom_point(aes(y = as.numeric(truth))) + 	xlab(NULL) +
		labs(fill='') #+ scale_alpha(range = c(0.4, 0.8))
}

if (model_type == "all") {
	p <- p + facet_grid(rows = vars(plotID),
						 cols = vars(model_name),
						 drop=T, scales="free")
		}

print(p)
}





siteID = "HARV"
site_plots = "overlap"
site_plots = NULL
model_type = "all"
model_type = "all_covariates"
input_df = hindcasts
taxon = "acidobacteriota"
beta_df = beta_effects_min %>% filter(time_period=="2015-11_2018-01")
plot_summary_inset <- function(input_df,
															 beta_df,
											 site_plots = "overlap",
											 siteID = "HARV",
											 plotID = NULL,
											 taxon = NULL,
											 model_type = "all_covariates") {

	if (model_type != "all") {
		message("Filtering to: model type: ", model_type)
		input_df <- input_df %>% filter(model_name == !!model_type)
	}


	if (!is.null(siteID)) {
		message("Filtering to: sites(s): ", siteID)
		siteID_vector <- ifelse(length(siteID) > 1,
														paste(siteID, collapse = "|"), siteID)
		input_df <- input_df %>% filter(grepl(siteID_vector, siteID))
	}

	if (!is.null(plotID)) {
		message("Filtering to: plot(s): ", plotID)
		plotID_vector <- ifelse(length(plotID) > 1,
														paste(plotID, collapse = "|"), plotID)
		input_df <- input_df %>% filter(grepl(plotID_vector, plotID))
	}

	if (!is.null(taxon)) {
		message("Filtering to: ", taxon)
		input_df <- input_df %>% filter(taxon == !!taxon)
	}

	n.taxa <- length(unique(input_df$species))
	if (n.taxa > 1) {
		message("More than one taxon/group in dataframe. Filtering to: ", input_df$species[[1]])
		taxon = input_df$species[[1]]
		input_df <- input_df %>% filter(species == !!taxon)
	}

	if (!"med" %in% colnames(input_df)) {
		input_df <- input_df %>% mutate(med = `50%`,
																		lo = `2.5%`,
																		hi = `97.5%`)
	} else {
		input_df <- input_df %>% mutate(med = ifelse(is.na(`50%`), med, `50%`),
																		lo = ifelse(is.na(`2.5%`), lo,`2.5%`),
																		hi = ifelse(is.na(`97.5%`), hi, `97.5%`))
	}

	if (!"fcast_period" %in% colnames(input_df)) {
		input_df <- input_df %>% mutate(fcast_period = "calibration")
	}

	input_df <- input_df %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4))

	# We need each y = 0 points at both groups ("positive" and "negative")
	input_df2 <- input_df %>% group_by(plotID) %>% filter(fcast_period == "calibration") %>% filter(timepoint == max(timepoint)) %>%
		mutate(fcast_period = "hindcast")

	# Now we adding y = 0 to original dataset
	input_df3 <- bind_rows(input_df, input_df2) %>%
		arrange(plotID, timepoint)

	p <- ggplot(input_df3, aes(x = dates ,  group = plotID))  +
		ylab("Abundance") +
		theme_bw() +
		theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
					legend.position = "bottom",legend.title = element_text(NULL),
					plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ggtitle(label = taxon)  +
		scale_fill_brewer(palette = "Paired") +
		geom_ribbon(data = input_df3 %>% filter(fcast_period == "calibration"),
								aes(ymin = lo, ymax = hi, group = plotID), alpha = .3, fill = "deepskyblue3") +
		geom_ribbon(data = input_df3 %>% filter(fcast_period == "hindcast"),
								aes(ymin = lo, ymax = hi, group = plotID), alpha = .1, fill = "deepskyblue3") +
		geom_line(data = input_df3 %>% filter(fcast_period == "calibration"), aes(y = med), show.legend = F, alpha = .6) +
		geom_line(data = input_df3 %>% filter(fcast_period == "hindcast"), aes(y = med), show.legend = F, alpha = .5, linetype=2) +
		geom_point(aes(y = as.numeric(truth)), position = position_jitter()) + 	xlab(NULL) +
		labs(fill='')  +
		scale_y_sqrt()

	beta_df_taxon <- beta_df %>% filter(taxon == !!taxon & model_name == "all_covariates")
	max_point_size <- ifelse(sum(beta_df_taxon$effSize) < 1, 1, sum(beta_df_taxon$effSize)^2)
	p2 <- ggplot(beta_df_taxon) + geom_point(aes(x = beta, y = 0, size = effSize, color = beta), show.legend = F) + theme_minimal(base_size = 16) + theme(
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.text.y=element_blank(), axis.ticks = element_blank()) +
		xlab(NULL) +
		ylab(NULL) +
		ylim(0,0) +
		ggtitle("Predictor importance") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_size(range = c(.1, max_point_size))


	# p3 <- p +
	# 	inset_element(
	# 		p2,
	# 		left = 0.5,
	# 		bottom = 0.6,
	# 		right = 1,
	# 			top=1,
	# 		align_to = "plot"
	# 	)
	p3 = p
	print(p3)
	return(p3)
}


library(patchwork)



plot_model(hindcast_data, taxon = "animal_pathogen", site_plots="facet", siteID= "CPER") + ylim(c(0,.07))
plot_model(hindcast_data, taxon = "animal_pathogen", site_plots="facet") + ylim(c(0,.07))


hindcast_data %>% filter(taxon == "animal_pathogen") %>% ggplot() + geom_boxplot(aes(x = reorder(siteID, truth), y = truth))


plot_model(hindcast_data, taxon = "animal_pathogen", site_plots="facet", siteID= "ORNL") + ylim(c(0,.1))
plot_model(hindcast_data, taxon = "plant_pathogen", site_plots="facet", siteID= "ORNL") + ylim(c(0,.1))


hindcast_data %>% filter(taxon %in% c("animal_pathogen","plant_pathogen") & siteID=="ORNL") %>%
ggplot() +
	facet_grid(rows=vars(plotID), cols=vars(taxon), drop=T, scales="free") +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`, fill=taxon), alpha=0.2) +
	geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`, fill=taxon), alpha=0.5) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth), fill=plotID)) + xlab(NULL) + labs(fill='') +
	scale_y_log10()

# hindcast_data %>% filter(taxon == "animal_pathogen") %>% ggplot() + geom_point(aes(x = siteID, y = truth))


# Visualize phenological parameters from all models
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(ggrepel, ggforce, gridExtra, ggpubr)


#######
# convert sin and cos to maximum and minimum of amplitude
# from Stolwijk 1999

# Agrees with results from Stolwijk 1999
# getMaxMin(sin = 0.009822, cos = 0.06929, max_only = F)

##### Create season data frame for plotting #####
mo <- 1:12
mo_sin <- sin((2*pi*mo)/12)
mo_cos <- cos((2*pi*mo)/12)
seasons <- cbind.data.frame(mo, mo_sin, mo_cos)
seasons$season <- "winter"
seasons$season[3:5] <- "spring"
seasons$season[6:8] <- "summer"
seasons$season[9:11] <- "fall"

if(model_name == "cycl_only") cov_key <- microbialForecast:::cycl_only_key


#####


df_orig <- readRDS("./data/summary/all_fcast_effects.rds")


df_orig$taxon <- ifelse(is.na(df_orig$taxon),
												paste(df_orig$pretty_group, df_orig$fcast_type),
												df_orig$taxon)

df <- df_orig %>%
	filter(beta %in% c("cos", "sin") &
				 	time_period == "2015-11_2018-01") %>%
	group_by(taxon, time_period, model_name) %>%
	mutate(cyc_significant = ifelse(any(significant==1), 1, 0))


vals <- df %>% filter(cyc_significant == 1) %>%
	pivot_wider(id_cols = c("taxon","model_name","fcast_type","pretty_group"),
							values_from = "Mean", names_from = "beta")

cycl_vals <- vals[vals$model_name=="cycl_only",]


# Couldn't figure out how to vectorize.
vals$max <- NA
vals$amplitude <- NA
for (i in 1:nrow(vals)){
	vals$max[[i]] <- getMaxMin(vals$sin[[i]],
														 vals$cos[[i]])
	vals$amplitude[[i]] <- sqrt(vals$sin[[i]]^2 + vals$cos[[i]]^2)
}
vals <- as.data.frame(vals)
vals$max <- unlist(vals$max)
vals$amplitude <- unlist(vals$amplitude)


# Get plotting coordinates for param estimates


vals$plot_sin <- sin((2*pi*vals$max)/12)
vals$plot_cos <- cos((2*pi*vals$max)/12)
title <- paste0("Peak cyclical abundances")
ggplot(vals)  +
	geom_circle(aes(x0 = 0, y0 = 0, r = 1)) + # base circle
	geom_label_repel(data = seasons[c(1,4,7,10),], # season names
									 aes(mo_sin*.6,mo_cos*.6, label = season), force = 2, size = 7)  +
	geom_label(data = seasons, aes(mo_sin*1.1,mo_cos*1.1, label = mo), size = 5)  + # month numbers
	theme_minimal(base_size = 16) +
	ggtitle(title) +
	geom_polygon(data = seasons, aes(mo_sin, mo_cos, fill = season), # season colors
							 size = 3, alpha = .2, show.legend = F) +
	geom_point(data = vals, aes(plot_sin, plot_cos,
															color = pretty_group, size = (amplitude*1000)^2),
						 shape = 18,  alpha = 1) + # mean estimate +
geom_text_repel(aes(x = plot_sin, y = plot_cos, label = taxon), force = 5) + xlab(NULL) + ylab(NULL)





fg_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")
fg_effects <- fg_in$summary_df %>%
	filter(beta %in% c("cos", "sin") &
				 	taxon != "other" &
				 	model_name == "cycl_only") %>%
	group_by(group, fg_cat, rank, taxon) %>%
	mutate(cyc_significant = ifelse(any(significant==1), 1, 0),
				 rank_only = "functional")



convergence_in <-
	readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxon_convergence_summaries_201511_201801.rds")
tax_effects <- convergence_in$summary_df

convergence_in$gelman.summary

# remove calculated in 02_summary_convergence
tax_effects <- tax_effects %>%
	filter(!rank %in% remove) %>%
	filter(beta %in% c("cos", "sin") &
				 	taxon != "other") %>%
	group_by(group, rank_only, rank) %>%
	mutate(cyc_significant = ifelse(any(significant==1), 1, 0))

vals <- plyr::rbind.fill(fg_effects, tax_effects) %>% filter(cyc_significant == 1) %>%
	pivot_wider(id_cols = c("fg_cat","rank_only","rank","group","taxon"),
							values_from = "Mean", names_from = "beta")

# Couldn't figure out how to vectorize.
vals$max <- NA
vals$amplitude <- NA
for (i in 1:nrow(vals)){
	vals$max[[i]] <- getMaxMin(vals$sin[[i]],
														 vals$cos[[i]])
	vals$amplitude[[i]] <- sqrt(vals$sin[[i]]^2 + vals$cos[[i]]^2)
}
vals <- as.data.frame(vals)
vals$max <- unlist(vals$max)
vals$amplitude <- unlist(vals$amplitude)

vals$plot_sin <- sin((2*pi*vals$max)/12)
vals$plot_cos <- cos((2*pi*vals$max)/12)
vals$only_rank <- ordered(vals$rank_only, levels = c("genus",
																								 "family",
																								 "order",
																								 "class",
																								 "phylum","functional"))
ggplot(vals %>% filter(only_rank != "functional"),
			 aes(x = only_rank, y = amplitude, color = group)) +
	geom_jitter() +
	geom_smooth(aes(x = as.numeric(only_rank)))

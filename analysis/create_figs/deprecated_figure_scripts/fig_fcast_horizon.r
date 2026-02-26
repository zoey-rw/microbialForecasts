# Visualize the forecast horizon
source("source.R")
library(tagger) # For tag_facet function


# Load the plotting data from 10_calculateFcastHorizon.r
# Structure: list(to_plot, to_plot_null, fcast_horizon_null_site, fcast_horizon_model_mean)
fcast_horizon_input <- readRDS("data/summary/fcast_horizon_input.rds")
to_plot = fcast_horizon_input[[1]]
to_plot_null = fcast_horizon_input[[2]] 
fcast_horizon_null_site = fcast_horizon_input[[3]]
fcast_horizon_model_mean = fcast_horizon_input[[4]]

# Load the final results
fcast_horizon_results <- readRDS("data/summary/fcast_horizon_df.rds")
fcast_horizon_df = fcast_horizon_results[[1]]  # Main results dataframe
fcast_horizon_plotting = fcast_horizon_results[[2]]  # Figure list (may be NULL if MAKE_PLOTS=FALSE)
fcast_horizon_long = fcast_horizon_results[[3]]  # Long format data

weak_converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
strict_converged <- readRDS(here("data/summary/converged_taxa_list.rds"))
converged = weak_converged
converged = strict_converged

# Fix model_id matching by removing _beta_regression suffix from converged models
converged_base <- gsub("_beta_regression$", "", converged)

# Fix pretty_group column using proper taxonomy structure from microbialForecast package
# Get the correct taxonomy lists
fg_names <- microbialForecast:::keep_fg_names
rank_spec_names <- microbialForecast:::rank_spec_names

# Create vectors of all bacteria and fungi species
all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
all_fungi <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])

# Assign pretty_group using same logic as 08_calculateScoringMetrics.r and 10_calculateFcastHorizon.r
assign_pretty_group <- function(data) {
  if (!"pretty_group" %in% colnames(data)) {
    data$pretty_group <- NA_character_
  }
  if ("species" %in% colnames(data)) {
    # Functional groups first
    if (any(data$species %in% fg_names)) {
      fg_cat <- microbialForecast:::assign_fg_categories(data$species)
      group <- microbialForecast:::assign_fg_kingdoms(fg_cat)
      data$pretty_group <- ifelse(group == "16S", "Bacteria", ifelse(group == "ITS", "Fungi", data$pretty_group))
    }
    # Taxonomic species
    data$pretty_group <- ifelse(data$species %in% all_bacteria, "Bacteria",
      ifelse(data$species %in% all_fungi, "Fungi", data$pretty_group))
    # Fallback: rank_name pattern (same as 08 and 10)
    still_na <- is.na(data$pretty_group)
    if (any(still_na) && "rank_name" %in% colnames(data)) {
      data$pretty_group[still_na] <- ifelse(
        grepl("_bac$|16S", data$rank_name[still_na]), "Bacteria",
        ifelse(grepl("_fun$|ITS", data$rank_name[still_na]), "Fungi", NA_character_)
      )
    }
  }
  data
}

# Apply pretty_group assignment to all data objects
to_plot <- assign_pretty_group(to_plot)
to_plot_null <- assign_pretty_group(to_plot_null)
fcast_horizon_null_site <- assign_pretty_group(fcast_horizon_null_site)
fcast_horizon_model_mean <- assign_pretty_group(fcast_horizon_model_mean)
fcast_horizon_df <- assign_pretty_group(fcast_horizon_df)
fcast_horizon_long <- assign_pretty_group(fcast_horizon_long)

# Filter by converged models if data exists
if (nrow(fcast_horizon_df) > 0) {
  fcast_horizon_df <- fcast_horizon_df %>% filter(model_id %in% converged_base) 
}
if (nrow(fcast_horizon_long) > 0) {
  fcast_horizon_long <- fcast_horizon_long %>% filter(model_id %in% converged_base)
}
#fcast_horizon_df <- fcast_horizon_df[fcast_horizon_df$model_id %in% converged_strict,]
#


# Remove infinite horizons 
# fcast_horizon_long <- fcast_horizon_long %>%
# 	filter(parameter_type == "horizon" &
# 				 	value !=Inf) 

# fcast_horizon_df <- fcast_horizon_df %>%
# 	filter(!is.infinite(crps_fcast_horizon) & 
# 				 	!is.infinite(rmse_fcast_horizon) & 
# 				 	!is.infinite(crps_fcast_horizon)) 


# # Set infinite horizons to the max value (20)
# fcast_horizon_long <- fcast_horizon_long %>%
# 	filter(model_id %in% converged) %>%
# 	mutate(value = ifelse(parameter_type == "horizon" &
# 																															 	value ==Inf,
# 																															 20, value))
# 
# fcast_horizon_df <- fcast_horizon_df %>%
# 	filter(model_id %in% converged) %>%
# 	mutate(crps_fcast_horizon = ifelse(crps_fcast_horizon ==Inf,
# 												20, crps_fcast_horizon),
# 				 rmse_fcast_horizon = ifelse(rmse_fcast_horizon ==Inf,
# 				 														20, rmse_fcast_horizon),
# 				 rsq_fcast_horizon = ifelse(rsq_fcast_horizon ==Inf,
# 				 														20, rsq_fcast_horizon))
# Filter to only finite months_since_obs and valid mean_crps for proper plotting
to_plot_filtered <- to_plot %>%
  filter(is.finite(months_since_obs) & months_since_obs >= 0 & 
         is.finite(mean_crps) & mean_crps > 0)

if (nrow(to_plot_filtered) > 0) {
  p1 <- ggplot(to_plot_filtered,
  			 aes(x = months_since_obs, y = mean_crps, color = pretty_group)) +
  	#	geom_point(alpha=.3, position=position_jitter(height=0, width=.3), size=3) +
  	facet_grid(rows=vars(model_name),
  						 labeller = labeller(model_name = model.labs)) +
  	geom_smooth(method="loess", span=1, se=F) +
  	stat_cor(
  		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  		p.digits = 1, label.x.npc = .65, label.y.npc = .75
  	) +
  	theme_minimal(base_size = 18) +
  	ylab("Forecast error (CRPS)") +
  	xlab("Months since last observation") +
  	labs(color="Kingdom") + ylim(c(0,.25))  +
  	scale_y_log10()
  
  # Save the plot
  png(here("figures", "fcast_horizon_error_over_time.png"), width = 1200, height = 800)
  print(p1)
  dev.off()
  cat("Figure saved to figures/fcast_horizon_error_over_time.png\n")
} else {
  cat("Warning: No valid data for forecast error plot (all months_since_obs are Inf or mean_crps is invalid)\n")
}


fcast_horizon_df$seasonal_predictors = as.factor(ifelse(fcast_horizon_df$model_name %in% c("env_cycl","cycl_only"), 1, 0))
fcast_horizon_df$environmental_predictors = as.factor(ifelse(fcast_horizon_df$model_name %in% c("env_cycl","env_cov"), 1, 0))
fcast_horizon_long$seasonal_predictors = as.factor(ifelse(fcast_horizon_long$model_name %in% c("env_cycl","cycl_only"), 1, 0))
fcast_horizon_long$environmental_predictors = as.factor(ifelse(fcast_horizon_long$model_name %in% c("env_cycl","env_cov"), 1, 0))

# Check what variables are available
cat("Available variables in fcast_horizon_df:", colnames(fcast_horizon_df), "\n")
cat("Number of rows in fcast_horizon_df:", nrow(fcast_horizon_df), "\n")

# Check if required variables exist before plotting
if (nrow(fcast_horizon_df) > 0 && all(c("model_name", "pretty_group", "crps_fcast_horizon", "rmse_fcast_horizon") %in% colnames(fcast_horizon_df))) {
  mean_fcast_horizon = fcast_horizon_df %>% 
  	group_by(model_id, model_name, pretty_group, rank_only, seasonal_predictors, environmental_predictors) %>% 
  	summarize(mean_horizon=mean(c(crps_fcast_horizon, rmse_fcast_horizon), 
  															#rsq_fcast_horizon, 
  															na.rm=T)) %>% filter(!is.infinite(mean_horizon))
  
  if (nrow(mean_fcast_horizon) > 0 && length(unique(mean_fcast_horizon$pretty_group)) > 0) {
    ggplot(mean_fcast_horizon,# %>% filter(rank_only %in% c("functional","phylum")),
    			 aes(x = model_name, y = mean_horizon, color = model_name)) +
    	facet_grid(~pretty_group, 
    						 labeller = labeller(model_name = model.labs), scales="free") + 
    	coord_flip() +
    	geom_boxplot() +
    	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw()
  } else {
    cat("Warning: No data available for mean_fcast_horizon plot\n")
  }
} else {
  cat("Warning: Required variables not found in fcast_horizon_df\n")
}

# Statistical analysis moved inside the conditional block above


# Statistical analysis with error checking
if (nrow(fcast_horizon_long) > 0) {
  for_stats = fcast_horizon_long %>% 
  	filter(parameter_type=="horizon" & metric=="rsq") %>% 
  	filter(!is.na(value) & !is.infinite(value))
  
  if (nrow(for_stats) > 0 && length(unique(for_stats$model_name)) > 1) {
    for_stats %>% 
    	group_by(pretty_group) %>%
    	filter(n_distinct(model_name) > 1) %>%  # Only run Tukey if multiple model_names in group
    	summarize(tukey(x = model_name, y = value)) %>%
    	rename(model_name = x)
  } else {
    cat("Warning: Insufficient data for statistical analysis\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_long for statistical analysis\n")
}


# Additional statistical analysis with error checking
if (exists("for_stats") && nrow(for_stats) > 0 && length(unique(for_stats$model_name)) > 1) {
  model_stat_pvalue <- for_stats %>% 
  	#group_by(pretty_group) %>% 
  	rstatix::tukey_hsd(value ~ model_name) %>%
  	#filter(p.adj < 0.05) %>% 
  	rstatix::add_y_position(step.increase = .4) %>% 
  	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
} else {
  cat("Warning: Skipping model_stat_pvalue calculation - insufficient data\n")
  model_stat_pvalue <- NULL
}

# View horizon results by model type with error checking
if (nrow(fcast_horizon_long) > 0) {
  horizon_data <- fcast_horizon_long %>% filter(parameter_type=="horizon" & metric=="rsq")  %>% 
  	filter(!is.na(value) & !is.infinite(value))
  
  if (nrow(horizon_data) > 0) {
    horizon_by_model_type <- ggplot(horizon_data,
    			 aes(x = model_name, y = value, color = model_name)) +
    	# facet_grid(~pretty_group, 
    	# 					 labeller = labeller(model_name = model.labs), scales="free") + 
    	coord_flip() +
    	geom_boxplot() +
    	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw()
    
    if (!is.null(model_stat_pvalue)) {
      horizon_by_model_type <- horizon_by_model_type +
        ggpubr::stat_pvalue_manual(model_stat_pvalue, label = "p.adj.signif", #bracket.nudge.y = -.4, 
    														 size=4, hide.ns = T)
    }
    horizon_by_model_type <- horizon_by_model_type + coord_flip()
    horizon_by_model_type
  } else {
    cat("Warning: No data available for horizon_by_model_type plot\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_long for horizon_by_model_type plot\n")
}




# Kingdom statistical analysis with error checking
if (exists("for_stats") && nrow(for_stats) > 0 && length(unique(for_stats$pretty_group)) > 1) {
  kingdom_stat_pvalue <- for_stats %>% 
  	#group_by(pretty_group) %>% 
  	rstatix::tukey_hsd(value ~ pretty_group) %>%
  	#filter(p.adj < 0.05) %>% 
  	rstatix::add_y_position(step.increase = .4) %>% 
  	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
} else {
  cat("Warning: Skipping kingdom_stat_pvalue calculation - insufficient data\n")
  kingdom_stat_pvalue <- NULL
}

# View horizon results by kingdom with error checking
if (nrow(fcast_horizon_long) > 0) {
  kingdom_data <- fcast_horizon_long %>% filter(parameter_type=="horizon" & metric=="rsq")  %>% 
  	filter(!is.na(value) & !is.infinite(value))
  
  if (nrow(kingdom_data) > 0) {
    horizon_by_kingdom <- ggplot(kingdom_data,
    	aes(x = pretty_group, y = value, color = pretty_group)) +
    	# facet_grid(~pretty_group, 
    	# 					 labeller = labeller(model_name = model.labs), scales="free") + 
    	coord_flip() +
    	geom_boxplot() +
    	geom_point(position=position_jitterdodge(jitter.width = .2), alpha=.3) + theme_bw()
    
    if (!is.null(kingdom_stat_pvalue)) {
      horizon_by_kingdom <- horizon_by_kingdom +
        ggpubr::stat_pvalue_manual(kingdom_stat_pvalue, label = "p.adj.signif", #bracket.nudge.y = -.4, 
    														 size=4, hide.ns = T)
    }
    horizon_by_kingdom <- horizon_by_kingdom + coord_flip()
    horizon_by_kingdom
  } else {
    cat("Warning: No data available for horizon_by_kingdom plot\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_long for horizon_by_kingdom plot\n")
}






# Additional statistical analyses with error checking
if (nrow(fcast_horizon_long) > 0) {
  env_data <- fcast_horizon_long %>% 
  	filter(rank_only %in% c("functional","phylum")) %>% 
  	filter(!is.na(value) & !is.infinite(value))
  
  if (nrow(env_data) > 0 && length(unique(env_data$environmental_predictors)) > 1) {
    env_data %>% 
    	group_by(pretty_group, metric) %>% 
    	rstatix::tukey_hsd(value ~ environmental_predictors)
  } else {
    cat("Warning: Insufficient data for environmental predictors analysis\n")
  }
} else {
  cat("Warning: No data for environmental predictors analysis\n")
}

# Linear model analysis with error checking (need 2+ levels for both factors for contrasts)
n_model <- length(unique(na.omit(fcast_horizon_df$model_name)))
n_kingdom <- length(unique(na.omit(fcast_horizon_df$pretty_group)))
if (nrow(fcast_horizon_df) > 0 && n_model > 1 && n_kingdom > 1) {
  model_parameters_lm <- lm(rsq_fcast_horizon ~
  	model_name * pretty_group,
  	data = fcast_horizon_df)
  sjPlot::plot_model(model_parameters_lm, main = "Effects on predictability (RSQ)", intercept = F)
} else {
  cat("Warning: Insufficient data for linear model analysis (need 2+ model_name and 2+ pretty_group, have ", n_model, " and ", n_kingdom, ")\n")
}  



# Environmental predictors statistical analysis with error checking
if (nrow(fcast_horizon_df) > 0) {
  env_df <- fcast_horizon_df %>% 
  	filter(rank_only %in% c("functional","phylum","genus")) %>% 
  	filter(!is.na(rsq_fcast_horizon) & !is.infinite(rsq_fcast_horizon))
  
  if (nrow(env_df) > 0 && length(unique(env_df$environmental_predictors)) > 1) {
    stat_pvalue_env <- env_df %>% 
    	group_by(pretty_group) %>% 
    	rstatix::tukey_hsd(rsq_fcast_horizon ~ environmental_predictors) %>%
    	#filter(p.adj < 0.05) %>% 
    	rstatix::add_y_position(step.increase = .4) %>% 
    	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
  } else {
    cat("Warning: Insufficient data for environmental predictors statistical analysis\n")
    stat_pvalue_env <- NULL
  }
} else {
  cat("Warning: No data for environmental predictors statistical analysis\n")
  stat_pvalue_env <- NULL
}

# Rank-only statistical analysis with error checking
if (nrow(fcast_horizon_df) > 0) {
  rank_df <- fcast_horizon_df %>% 
  	filter(rank_only %in% c("functional","phylum","genus")) %>% 
  	filter(!is.na(rsq_fcast_horizon) & !is.infinite(rsq_fcast_horizon))
  
  if (nrow(rank_df) > 0 && length(unique(rank_df$rank_only)) > 1) {
    stat_pvalue <- rank_df %>% 
    	group_by(pretty_group) %>% 
    	rstatix::tukey_hsd(rsq_fcast_horizon ~ rank_only) %>%
    	#filter(p.adj < 0.05) %>% 
    	rstatix::add_y_position(step.increase = .4) %>% 
    	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))
  } else {
    cat("Warning: Insufficient data for rank-only statistical analysis\n")
    stat_pvalue <- NULL
  }
} else {
  cat("Warning: No data for rank-only statistical analysis\n")
  stat_pvalue <- NULL
}



# Rank-only plotting with error checking
if (nrow(fcast_horizon_df) > 0) {
  rank_plot_data <- fcast_horizon_df %>% filter(rank_only %in% c("functional","phylum","genus"))
  
  if (nrow(rank_plot_data) > 0) {
    rank_plot <- ggplot(rank_plot_data,
    			 aes(x = rank_only, y = rsq_fcast_horizon, color = model_name)) +
    	facet_grid(~pretty_group, 
    						 labeller = labeller(model_name = model.labs), scales="free") + 
    	geom_boxplot() +
    	coord_flip() +
    	theme_minimal(base_size = 20) + 
    	geom_point(position=position_jitterdodge(), alpha=.3) + 
    	scale_y_continuous(n.breaks=10)
    
    if (!is.null(stat_pvalue)) {
      rank_plot <- rank_plot +
        ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.adj.signif", bracket.nudge.y = -.4, size=4, hide.ns = T)
    }
    rank_plot
  } else {
    cat("Warning: No data available for rank-only plot\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_df for rank-only plot\n")
}


# View horizon results by model type with error checking
if (nrow(fcast_horizon_long) > 0) {
  phylum_data <- fcast_horizon_long %>% filter(parameter_type=="horizon") %>% 
  	filter(rank_only %in% c("phylum"))
  
  if (nrow(phylum_data) > 0 && length(unique(phylum_data$pretty_group)) > 0) {
    ggplot(phylum_data,
    			 aes(x = model_name, y = value, color = pretty_group)) +
    	facet_grid(metric~pretty_group,  scales="free") + 
    	coord_flip() +
    	geom_boxplot() +
    	geom_point(position=position_jitter(), alpha=.3)
  } else {
    cat("Warning: No data available for phylum horizon plot\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_long for phylum horizon plot\n")
} 




# View mean values per group, RSQ.1 with error checking
cat("Debug: fcast_horizon_model_mean rows:", nrow(fcast_horizon_model_mean), "\n")
b <- NULL  # Initialize b outside tryCatch
if (nrow(fcast_horizon_model_mean) > 0) {
  model_mean_data <- fcast_horizon_model_mean %>% filter(months_since_obs < 14)
  cat("Debug: model_mean_data rows after filtering:", nrow(model_mean_data), "\n")
  
  if (nrow(model_mean_data) > 0 && length(unique(model_mean_data$pretty_group[!is.na(model_mean_data$pretty_group)])) > 0) {
    cat("Debug: Creating b plot...\n")
    tryCatch({
      b = ggplot(model_mean_data, # %>% filter(rank_only=="phylum"),
      			 aes(x = months_since_obs, y = RSQ.1, color = pretty_group)) +
      	
      	facet_grid(pretty_group ~model_name, 
      						 labeller = labeller(model_name = model.labs), scales="free") +
      	geom_smooth(span=2,  method="loess", show.legend = F, se=F) +
      	stat_cor(
      		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
      		p.digits = 1, label.x.npc = .15, label.y.npc = .15
      	) +
      	theme_minimal(base_size = 20) +
      	ylab("Mean forecast accuracy (RSQ 1:1)") +
      	xlab("Forecast horizon (months since last observation)") +
      	labs(color="Kingdom") +
      	#	geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      	# geom_hline(data =fcast_horizon_null_site, aes(yintercept = mean(null_RSQ.1)), na.rm = T) 
      
      if (!is.null(b)) {
        b <- tag_facet(b, size=7)
      }

      png(here("figures","forecast_horizon.png"), width = 1400, height=1000)
      print(b)
      dev.off()
      cat("Debug: b plot created and saved successfully\n")
    }, error = function(e) {
      cat("Debug: Error creating b plot:", e$message, "\n")
    })
  } else {
    cat("Warning: No data available for RSQ.1 plot\n")
  }
} else {
  cat("Warning: No data in fcast_horizon_model_mean for RSQ.1 plot\n")
}



# Plot CRPS values for all observations
# to_plot is now loaded from fcast_horizon_input.rds above
a <- ggplot(to_plot,
						aes(x = months_since_obs, y = mean_crps, color = pretty_group)) +
#	geom_point(alpha=.3, position=position_jitter(height=0, width=.3), size=3) +
	facet_grid(rows=vars(model_name),
						 labeller = labeller(model_name = model.labs)) +
	geom_smooth(method="loess", span=1, se=F) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .65, label.y.npc = .75
	) +
	theme_minimal(base_size = 18) +
	ylab("Forecast error (CRPS)") +
	xlab("Months since last observation") +
	labs(color="Kingdom") + ylim(c(0,.25))  +
	scale_y_log10()
a

# Site-means plotting removed - data structure changed



# Data is now loaded from fcast_horizon_input.rds above



# Plot CRPS values for all observations

# View mean values per group, RSQ, with fcast horizon with error checking
c_plot <- NULL  # Initialize c_plot outside tryCatch
if (nrow(to_plot) > 0 && all(c("months_since_obs", "RSQ.1", "pretty_group", "model_name") %in% colnames(to_plot)) && length(unique(to_plot$pretty_group[!is.na(to_plot$pretty_group)])) > 0) {
  tryCatch({
    c_plot = ggplot(to_plot, # %>% filter(rank_only=="phylum"),
    					 aes(x = months_since_obs, y = RSQ.1, color = pretty_group)) +
    	
    	facet_grid(pretty_group ~model_name, 
    						 labeller = labeller(model_name = model.labs), scales="free") +
    	geom_smooth(span=2,  method="loess", show.legend = F, se=F) +
    	stat_cor(
    		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    		p.digits = 1, label.x.npc = .15, label.y.npc = .15
    	) +
    	theme_minimal(base_size = 20) +
    	ylab("Mean forecast accuracy (RSQ 1:1)") +
    	xlab("Forecast horizon (months since last observation)") +
    	labs(color="Kingdom") +
    	#	geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    	# geom_hline(data =fcast_horizon_null_site, aes(yintercept = mean(null_RSQ.1)), na.rm = T) 
    c_plot
  }, error = function(e) {
    cat("Error creating c_plot:", e$message, "\n")
    NULL
  })
} else {
  cat("Warning: No data available for c_plot (RSQ plot)\n")
}


# All plotting data is now loaded from fcast_horizon_input.rds above

# Additional plotting code removed - data structure changed 

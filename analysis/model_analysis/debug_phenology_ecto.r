library(NISTunits)
library(ggrepel)
library(ggforce)


#######
# convert sin and cos to maximum and minimum of amplitude
# from Stolwijk 1999
getMaxMin <- function(sin, cos, T = 12, max_only = T) {
	
	print(sin[[1]]);
	print(cos[[1]])
	
	print(length(sin));
	print(length(cos))
	
	sin <- sin[[1]]
	cos <- cos[[1]]
	
	t <- atan(sin/cos) * T/(2*pi)
	
	if ((sin/cos) > 0){
		extreme1 <- t
		extreme2 <- extreme1 + T/2 
	} else if ((sin/cos) <= 0){
		extreme1 <- t + T/2 
		extreme2 <- t + T
	}
	
	if (sin > 0){
		max <- extreme1 
		min <- extreme2 
	} else if (sin <= 0){
		min <- extreme1 
		max <- extreme2 
	}
	if (max_only) {
		return(max)
	} else {
		return(list("min" = min, "max" = max))
	}
}

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

cycl_only_key <- list("1" = "sin",
											"2" = "cos")
if(model_name == "cycl_only") cov_key <- cycl_only_key


#####

read_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups/cycl_only/refit_samples_ectomycorrhizal_full_uncertainty.rds")
model_name = "cycl_only"
time_period = "refit"
rank.name = "Ectomycorrhizal"


read_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups/cycl_only/refit_samples_saprotroph_full_uncertainty.rds")
model_name = "cycl_only"
time_period = "refit"
rank.name = "Saprotroph"


read_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity//cycl_only/refit_samples_div_full_uncertainty_ITS.rds")
model_name = "cycl_only"
time_period = "refit"
rank.name = "Fungal diversity"

samples <- read_in$samples
param_summary <- read_in$param_summary
means <- param_summary[[1]]

# Get beta sizes per rank
beta_out <-  means %>% as.data.frame() %>% 
	rownames_to_column("rowname") %>% filter(grepl("beta|rho", rowname)) %>% 
	mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
	mutate(beta = recode(beta_num, !!!cov_key),
				 taxon = rank.name,)
beta_out[grep("rho", beta_out$rowname),]$beta = "rho"
beta_out[grep("rho", beta_out$rowname),]$beta_num = "0"
# Use quantiles to assign significance to beta parameters.
beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
	rownames_to_column("rowname") %>% filter(grepl("beta|rho", rowname)) %>% 
	mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) 
beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
															 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
															 1, 0)
beta_out$effSize <- abs(beta_out$Mean)

# Combine parameter estimates into summary
summary_df <- plyr::rbind.fill(beta_out)
summary_df$scenario <- scenario
summary_df$time_period <- time_period
summary_df$model_name <- model_name
summary_df$taxon <- rank.name

df <- summary_df

df <- df %>% 
	filter(beta %in% c("cos", "sin") & 
				 	time_period == "refit") %>% 
	group_by(taxon, time_period) %>% 
	mutate(cyc_significant = ifelse(any(significant==1), 1, 0))


vals <- df %>% #filter(cyc_significant == 1) %>% 
	pivot_wider(id_cols = c("taxon","model_name","scenario"),#,"fcast_type"), 
							values_from = "Mean", names_from = "beta")

cycl_vals <- vals[vals$model_name=="cycl_only",]


# Get actual MCMC samples for cos/sin params
beta_samps <- as.data.frame(do.call(rbind, samples[,c("beta[1]", "beta[2]")]))
colnames(beta_samps) <- c("sin","cos")
beta_samps_reduced <- beta_samps[sample(1:nrow(beta_samps), 1000),]

# Couldn't figure out how to vectorize.
beta_samps_reduced$max <- NA
for (i in 1:nrow(beta_samps_reduced)){
	beta_samps_reduced$max[[i]] <- getMaxMin(beta_samps_reduced$sin[[i]],
																					 beta_samps_reduced$cos[[i]])
}
samps_phenology <- as.data.frame(beta_samps_reduced)
samps_phenology$max <- unlist(samps_phenology$max)
									 

vals$max <- getMaxMin(vals$sin, vals$cos)

# Get plotting coordinates for param estimates
samps_phenology$plot_sin <- sin((2*pi*samps_phenology$max)/12)
samps_phenology$plot_cos <- cos((2*pi*samps_phenology$max)/12)

vals$plot_sin <- sin((2*pi*vals$max)/12)
vals$plot_cos <- cos((2*pi*vals$max)/12)
title <- paste0("Peak cyclical abundance: ", rank.name)
ggplot(samps_phenology)  + 
	geom_circle(aes(x0 = 0, y0 = 0, r = 1)) + # base circle
	geom_label_repel(data = seasons[c(1,4,7,10),], # season names
									 aes(mo_sin*.6,mo_cos*.6, label = season), force = 2, size = 7)  + 
	geom_label(data = seasons, aes(mo_sin*1.1,mo_cos*1.1, label = mo), size = 7)  + # month numbers
		theme_minimal(base_size = 16) +
	ggtitle(title) + 
	geom_polygon(data = seasons, aes(mo_sin, mo_cos, fill = season), # season colors
							 size = 3, alpha = .2, show.legend = F) +
	geom_point(data = samps_phenology, aes(plot_sin, plot_cos), # parameter estimates
						 size = 4, show.legend = F, alpha = .1) + 
	geom_point(data = vals, aes(plot_sin, plot_cos), 
						 shape = 18, color = 2, size = 6, alpha = 1) # mean estimate

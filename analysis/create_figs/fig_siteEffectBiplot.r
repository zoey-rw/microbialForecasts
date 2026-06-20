source("source.R")
pacman::p_load(ggforce, ggrepel)

#https://stackoverflow.com/questions/39137287/plotting-partial-least-squares-regression-plsr-biplot-with-ggplot2

pred_cat_key = list("nitrogenTot" = "macronutrient",
										"totalP" = "macronutrient",
										"sulfurTot" = "macronutrient",
										"estimatedOC" = "macronutrient",
										"pMjelm"= "macronutrient",
										"feMjelm" = "metal",
										"alOxalate" = "metal",
										"cecdNh4" = "cations",
										"mgNh4d" = "cations",
										"naNh4d" = "cations",
										"mnMjelm" = "cations",
										"kNh4d"  = "cations",
										"caNh4d" = "cations",
										"MAT"	     = "climate",
										"MAP"      = "climate",
										"latitude_scaled"      = "climate",
										"so4Satx"= "micronutrient",
										"siMjelm"  = "micronutrient")

pred_name_key = list("no2Satx" = "nitrite",
										 "totalP" = "total phosphorus",
										 "nitrogenTot" = "total nitrogen",
										 "estimatedOC" = "organic carbon",
										 "no3Satx" = "nitrate",
										 "pMjelm"= "phosphorus (Mehlich)",
										 "alOxalate"= "aluminum",
										 "feMjelm" = "iron (Mehlich)",
														 "alKcl" = "aluminum",
														 "cecdNh4" = "cation exchange",
														 "mgNh4d" = "magnesium",
														 "naNh4d" = "sodium",
														 "mnMjelm" = "manganese (Mehlich)",
														 "kNh4d"  = "potassium",
														 "caNh4d" = "calcium",
														 "MAT"	     = "mean annual temperature",
														 "MAP"      = "mean annual precipitation",
														 "latitude_scaled" = "latitude",
										 "so4Satx"= "sulfate",
										 "sulfurTot"= "sulfur",
										 "siMjelm"  = "silicon (Mehlich)")


read_in = readRDS(here("data/summary/site_effects_dredged.rds"))
# Check structure: should be list with 5 elements:
# [[1]] = dredged_predictor_importance
# [[2]] = pred_sites
# [[3]] = pred_sites_plsr
# [[4]] = plsr_model_importance
# [[5]] = plsr_model_scores (contains plsr_scores for all models)

component_df <- NULL

# First try: Check if plsr_model_scores is in the 5th element (new structure)
if (length(read_in) >= 5 && !is.null(read_in[[5]]) && nrow(read_in[[5]]) > 0) {
  plsr_scores_all <- read_in[[5]]
  cat("Found plsr_model_scores in read_in[[5]] with", nrow(plsr_scores_all), "rows\n")
  
  # Filter for a specific model (use env_cycl_saprotroph as in the title, or first available)
  target_model <- "env_cycl_saprotroph_20151101_20180101"
  
  # Check if predictor names are in a column or as row names
  has_predictor_col <- "predictor" %in% names(plsr_scores_all)
  
  if (target_model %in% plsr_scores_all$model_id) {
    component_df <- plsr_scores_all %>% 
      filter(model_id == target_model) %>%
      select(-model_id) %>%
      as.data.frame()
    
    # Handle predictor names - either from column or row names
    if (has_predictor_col && "predictor" %in% names(component_df)) {
      rownames(component_df) <- component_df$predictor
      component_df$predictor <- NULL
    } else if (nrow(component_df) > 0 && length(rownames(component_df)) == 0) {
      # If no row names and no predictor column, use row numbers (shouldn't happen)
      rownames(component_df) <- paste0("predictor_", 1:nrow(component_df))
    }
    cat("Using plsr_scores for model:", target_model, "\n")
  } else if (nrow(plsr_scores_all) > 0) {
    # Use first available model
    first_model <- unique(plsr_scores_all$model_id)[1]
    component_df <- plsr_scores_all %>% 
      filter(model_id == first_model) %>%
      select(-model_id) %>%
      as.data.frame()
    
    # Handle predictor names
    if (has_predictor_col && "predictor" %in% names(component_df)) {
      rownames(component_df) <- component_df$predictor
      component_df$predictor <- NULL
    } else if (nrow(component_df) > 0 && length(rownames(component_df)) == 0) {
      rownames(component_df) <- paste0("predictor_", 1:nrow(component_df))
    }
    cat("Using plsr_scores for first available model:", first_model, "\n")
  }
}

# Fallback: Try old structure (read_in[[4]] might contain individual model results)
if (is.null(component_df) && length(read_in) >= 4) {
  cat("Trying fallback: checking read_in[[4]] for plsr_scores\n")
  if (is.list(read_in[[4]]) && length(read_in[[4]]) > 0) {
    # Try to find an element with plsr_scores
    for (i in seq_along(read_in[[4]])) {
      elem <- read_in[[4]][[i]]
      # Check if element is a list (not atomic) and has plsr_scores
      if (is.list(elem) && !is.atomic(elem) && !is.null(elem$plsr_scores)) {
        component_df = elem$plsr_scores
        cat("Using element", i, "from read_in[[4]]\n")
        break
      } else if (is.list(elem) && !is.atomic(elem) && "plsr_scores" %in% names(elem)) {
        component_df = elem[["plsr_scores"]]
        cat("Using element", i, "from read_in[[4]] (accessed via [[)\n")
        break
      }
    }
  }
}

if (is.null(component_df)) {
  cat("Error: No plsr_scores found in site_effects_dredged.rds\n")
  cat("Length of read_in:", length(read_in), "\n")
  if (length(read_in) >= 5) {
    cat("read_in[[5]] structure:\n")
    print(str(read_in[[5]]))
  }
  stop("Cannot proceed without plsr_scores data. Please regenerate site_effects_dredged.rds with updated script.")
}
component_df$predictor_category = recode(rownames(component_df), !!!pred_cat_key)
rownames(component_df) = recode(rownames(component_df), !!!pred_name_key)
p1 <- ggplot(data=component_df, aes(x = `Comp 1`,y = `Comp 2`))+
	ylab("Component 2")+
	xlab("Component 1")+
	ggtitle("Site effects, partial least squares regression", "model_id: env_cycl_saprotroph_20151101_20180101")+
	theme_bw(base_size = 18) +
	theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),panel.grid.major = element_blank(),
										 panel.grid.minor = element_blank(),
										 axis.line = element_line(colour = "black"))+
	#geom_text(aes(label=rownames(component_df), color=predictor_category))+
	geom_text_repel(aes(label=rownames(component_df), color=predictor_category), max.overlaps = 20, size=5)+
	coord_fixed(ylim=c(-1, 1),xlim=c(-1,1))+
	theme(axis.ticks = element_line(colour = "red")) +
	ggforce::geom_circle(data=data.frame(x0=0, y0=0, r=sqrt(1/2)), aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE) +
	ggforce::geom_circle(data=data.frame(x0=0, y0=0, r=1), aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE)  + labs(color = "Predictor type")
# Check if Comp 3 and Comp 4 exist before creating p2
if ("Comp 3" %in% names(component_df) && "Comp 4" %in% names(component_df)) {
  p2 <- ggplot(data=component_df, aes(x = `Comp 3`,y = `Comp 4`))+
  	ylab("Component 4")+
  	xlab("Component 3")+
  	ggtitle("","")+
  	#ggtitle("")+
  	theme_bw(base_size = 18) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
  										 panel.grid.major = element_blank(),
  										 panel.grid.minor = element_blank(),
  										 axis.line = element_line(colour = "black"))+
  #	geom_text(aes(label=rownames(component_df), color=predictor_category))+
  	geom_text_repel(aes(label=rownames(component_df), color=predictor_category), max.overlaps = 20, size=5)+
  
  	coord_fixed(ylim=c(-1, 1),xlim=c(-1,1))+
  	theme(axis.ticks = element_line(colour = "red")) +
  	ggforce::geom_circle(data=data.frame(x0=0, y0=0, r=sqrt(1/2)), mapping=aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE) +
  	ggforce::geom_circle(data=data.frame(x0=0, y0=0, r=1), mapping=aes(x0=x0, y0=y0, r=r), inherit.aes=FALSE) + labs(color = "Predictor type")
  combined_plot <- ggarrange(p1, p2, common.legend = T)
} else {
  cat("Warning: Comp 3 or Comp 4 not found in component_df. Available columns:", paste(names(component_df), collapse=", "), "\n")
  cat("Creating plot with only p1\n")
  combined_plot <- p1
}

# Save the plot
png(here("figures","site_effect_biplot.png"), width = 1600, height=800)
print(combined_plot)
dev.off()
cat("Figure saved to figures/site_effect_biplot.png\n")


# scoreplot and plot functions from pls package - commented out to avoid errors if model not available
# library(pls)
# if (length(read_in[[4]]) >= 660 && !is.null(read_in[[4]][[660]])) {
#   scoreplot(read_in[[4]][[660]], comps = 4, labels="names")
#   plot(read_in[[4]][[660]], "biplot")
# }

# pred_cat_key and pred_name_key are already defined at the top of the file

pred_sites <- read_in[[2]]
pred_sites_plsr <- read_in[[3]]

# Sanity checks: write model summaries to file instead of plotting
sanity_file <- here("figures", "site_effect_biplot_sanity.txt")
if (!is.null(pred_sites_plsr) && nrow(pred_sites_plsr) > 0 && 
    "Median" %in% names(pred_sites_plsr) && "plsr_pred_observed_4comp" %in% names(pred_sites_plsr)) {
  lm_plsr <- lm(pred_sites_plsr$Median ~ pred_sites_plsr$plsr_pred_observed_4comp)
  write("=== PLSR: Median ~ plsr_pred_observed_4comp ===\n", sanity_file)
  capture.output(summary(lm_plsr), file = sanity_file, append = TRUE)
} else {
  write("Warning: pred_sites_plsr data not available or missing required columns\n", sanity_file)
}

if (!is.null(pred_sites) && nrow(pred_sites) > 0 && 
    "Median" %in% names(pred_sites) && "pred" %in% names(pred_sites)) {
  lm_pred <- lm(pred_sites$Median ~ pred_sites$pred)
  write("\n=== Observed: Median ~ pred ===\n", sanity_file, append = TRUE)
  capture.output(summary(lm_pred), file = sanity_file, append = TRUE)
} else {
  write("Warning: pred_sites data not available or missing required columns\n", sanity_file, append = TRUE)
}

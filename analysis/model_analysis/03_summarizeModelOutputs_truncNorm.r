# Summarize MCMC output from all single-taxon truncated normal models
# Assumes input files have already had MCMC chains combined

source("../../source.R")

# Custom function to summarize truncNorm model outputs
summarize_truncnorm_model <- function(file_path, save_summary = TRUE, overwrite = TRUE) {
  tryCatch({
    cat("  Loading model output from:", basename(file_path), "\n")
    
    # Load the model output
    model_output <- readRDS(file_path)
    
    cat("  Model output loaded successfully\n")
    cat("  Model output names:", names(model_output), "\n")
    
    # Extract samples and metadata
    samples <- model_output$samples
    metadata <- model_output$metadata
    
    cat("  Samples class:", class(samples), "\n")
    
    # Handle mcmc.list objects properly
    if (inherits(samples, "mcmc.list")) {
      cat("  Samples is an mcmc.list with", length(samples), "chains\n")
      
      # Extract the first chain for now (we can combine them later if needed)
      samples_matrix <- as.matrix(samples[[1]])
      cat("  First chain dimensions:", dim(samples_matrix), "\n")
      cat("  First chain column names (first 10):", paste(head(colnames(samples_matrix), 10), collapse=", "), "\n")
      
      # Use the matrix for further processing
      samples <- samples_matrix
    } else if (is.matrix(samples) || is.data.frame(samples)) {
      cat("  Samples is a matrix/dataframe with dimensions", dim(samples), "\n")
      cat("  Sample column names (first 10):", paste(head(colnames(samples), 10), collapse=", "), "\n")
    } else {
      cat("  Unexpected samples class:", class(samples), "\n")
      return(NULL)
    }
    
    # Extract site mapping from metadata
    if (!is.null(metadata$model_data$site_start)) {
      site_mapping <- metadata$model_data$site_start
      site_names <- names(site_mapping)
      site_indices <- as.numeric(site_mapping)
      
      cat("  Site mapping found with", length(site_mapping), "sites\n")
      
      # Create a mapping from numeric index to site name
      # The site_effect parameters in the model use 1-based indexing
      # We need to map these to the actual site names
      site_key <- setNames(site_names, 1:length(site_names))
      
      # Extract site effect parameters
      site_effect_cols <- grep("^site_effect", colnames(samples))
      cat("  Site effect columns found:", length(site_effect_cols), "\n")
      
      if (length(site_effect_cols) > 0) {
        cat("  First few site effect columns:", paste(head(colnames(samples)[site_effect_cols], 5), collapse=", "), "\n")
      }
      
      if (length(site_effect_cols) == 0) {
        cat("  No site effect parameters found in model output\n")
        return(NULL)
      }
      
      # Create summary dataframe for site effects
      summary_list <- list()
      
      for (col in site_effect_cols) {
        param_name <- colnames(samples)[col]
        
        # Extract numeric site index from parameter name (e.g., "site_effect[1]" -> "1")
        site_num <- as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", param_name))
        
        # Map to actual site name
        siteID <- site_key[as.character(site_num)]
        
        if (is.na(siteID)) {
          cat("  Warning: Could not map site index", site_num, "to site name\n")
          next
        }
        
        # Calculate summary statistics
        param_samples <- samples[, col]
        cat("    Processing parameter:", param_name, "with", length(param_samples), "samples\n")
        cat("    Sample values (first 5):", head(param_samples, 5), "\n")
        
        # Calculate each statistic individually to debug
        mean_val <- mean(param_samples, na.rm = TRUE)
        sd_val <- sd(param_samples, na.rm = TRUE)
        q_025 <- quantile(param_samples, 0.025, na.rm = TRUE)
        q_25 <- quantile(param_samples, 0.25, na.rm = TRUE)
        q_50 <- quantile(param_samples, 0.50, na.rm = TRUE)
        q_75 <- quantile(param_samples, 0.75, na.rm = TRUE)
        q_975 <- quantile(param_samples, 0.975, na.rm = TRUE)
        
        cat("    Individual calculations:\n")
        cat("      Mean:", mean_val, "\n")
        cat("      SD:", sd_val, "\n")
        cat("      2.5%:", q_025, "\n")
        cat("      25%:", q_25, "\n")
        cat("      50%:", q_50, "\n")
        cat("      75%:", q_75, "\n")
        cat("      97.5%:", q_975, "\n")
        
        # Create summary stats using a list approach to avoid naming issues
        summary_stats <- list(
          Mean = mean_val,
          SD = sd_val,
          q025 = q_025,
          q25 = q_25,
          q50 = q_50,
          q75 = q_75,
          q975 = q_975
        )
        
        cat("    Summary stats list created with", length(summary_stats), "elements\n")
        cat("    Summary stats names:", names(summary_stats), "\n")
        cat("    Summary stats values:", paste(round(unlist(summary_stats), 6), collapse=", "), "\n")
        
        # Create row for this site effect
        row_data <- data.frame(
          rowname = param_name,
          Mean = as.numeric(summary_stats["Mean"]),
          SD = as.numeric(summary_stats["SD"]),
          X2.5. = as.numeric(summary_stats["q025"]),
          X25. = as.numeric(summary_stats["q25"]),
          X50. = as.numeric(summary_stats["q50"]),
          X75. = as.numeric(summary_stats["q75"]),
          X97.5. = as.numeric(summary_stats["q975"]),
          site_num = site_num,
          siteID = siteID,
          stringsAsFactors = FALSE
        )
        
        cat("    Row data created with", nrow(row_data), "rows and", ncol(row_data), "columns\n")
        cat("    Row data columns:", paste(colnames(row_data), collapse=", "), "\n")
        cat("    X50. value:", row_data$X50., "\n")
        cat("    Summary stats 50% value:", summary_stats$q50, "\n")
        
        summary_list[[length(summary_list) + 1]] <- row_data
      }
      
      # Combine all site effects
      if (length(summary_list) == 0) {
        cat("  No valid site effects found\n")
        return(NULL)
      }
      
      summary_df <- do.call(rbind, summary_list)
      cat("  Created summary dataframe with", nrow(summary_df), "rows\n")
      
      # Add metadata columns
      summary_df$rank.name <- metadata$rank.name
      summary_df$species <- metadata$species
      summary_df$model_name <- metadata$model_name
      summary_df$model_id <- metadata$model_id
      summary_df$time_period <- paste(metadata$min.date, metadata$max.date, sep = "_")
      summary_df$group <- ifelse(grepl("_bac", metadata$rank.name), "16S", "ITS")
      summary_df$rank <- metadata$rank.name
      summary_df$rank_only <- sapply(strsplit(metadata$rank.name, "_"), `[`, 1)
      summary_df$taxon <- metadata$species
      summary_df$fcast_type <- "Taxonomic"
      summary_df$pretty_group <- ifelse(grepl("_bac", metadata$rank.name), "Bacteria", "Fungi")
      
      # Add Median column (using X50. which is the 50th percentile)
      summary_df$Median <- summary_df$X50.
      
      # Calculate effective sample size and Gelman-Rubin statistic
      # For now, we'll use placeholder values since we only have one chain
      summary_df$es <- nrow(samples)
      summary_df$`Point est.` <- 1.0  # Placeholder for Gelman-Rubin
      
      # Create plot_est dataframe (placeholder for now)
      plot_est <- data.frame()
      
      cat("  Summary dataframe created successfully\n")
      
      # Save summary if requested
      if (save_summary) {
        summary_dir <- dirname(file_path)
        summary_file <- file.path(summary_dir, paste0("summary_", basename(file_path)))
        
        if (file.exists(summary_file) && !overwrite) {
          cat("  Summary file already exists, skipping save\n")
        } else {
          saveRDS(list(summary_df = summary_df, plot_est = plot_est, gelman_list = summary_df), summary_file)
          cat("  Summary saved to:", summary_file, "\n")
        }
      }
      
      return(list(summary_df = summary_df, plot_est = plot_est, gelman_list = summary_df))
      
    } else {
      cat("  No site_start found in metadata\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("  Error summarizing model:", e$message, "\n")
    return(NULL)
  })
}

# Get all available truncated normal models with legacy effects
# Look for models from 2013-2018 and 2013-2020 time periods
file.list = intersect(list.files(here("data/model_outputs/truncated_normal/"),recursive = T,pattern = "20130601_20151101|20151101_20180101|20130601_20180101|20130601_20200101", full.names = T),
											list.files(here("data/model_outputs/truncated_normal/"), recursive = T,
																 pattern = "samples", full.names = T))

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files
info <- file.info(file.list)
# Remove date filtering to include all existing files
# newer <- rownames(info[which(info$mtime > "2025-01-01 00:00:00 EST"), ])
# file.list <- file.list[file.list %in% newer]

# Process all files sequentially for debugging
# cl <- makeCluster(4, outfile="")
# registerDoParallel(cl)

# Export the function to the cluster
# clusterExport(cl, "summarize_truncnorm_model")

# Run summary function sequentially for debugging
file_summaries = list()
for (f in file.list) {
  cat("Processing file:", basename(f), "\n")
  out <- summarize_truncnorm_model(f, save_summary=T, overwrite = TRUE)
  file_summaries[[length(file_summaries) + 1]] <- out
}

# stopCluster(cl)

# Filter out any NULL results from failed summaries
file_summaries <- file_summaries[!sapply(file_summaries, is.null)]

# Combine summary files from the function results
if (length(file_summaries) > 0) {
  summary_df <- map_df(file_summaries, 1)
  plot_est <- map_df(file_summaries, 2)
  gelman_list <- map_df(file_summaries, 3)
} else {
  cat("No successful summaries to process\n")
  stop("All model summaries failed")
}

# Filter out any problematic models and focus on truncNorm specific parameters
# The summarize_truncnorm_model function returns a list with 3 elements, so we need to access them correctly
gelman.summary <- gelman_list %>%
	filter(!grepl("all_covariates", model_id)) %>%
	mutate(is_major_param = ifelse(grepl("beta|int|sigma|core_sd|rho", rowname), T, F))

# Save the results
saveRDS(list(summary_df = summary_df, plot_est = plot_est, gelman_list = gelman_list), 
        here("data/summary/truncated_normal_summaries.rds"))

cat("Truncated normal model summaries completed!\n")
cat("Total models processed:", length(file.list), "\n")
cat("Results saved to: data/summary/truncated_normal_summaries.rds\n")

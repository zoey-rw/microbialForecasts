# Tidy and format truncated normal hindcast outputs for analysis
# This script processes the raw hindcast outputs and creates clean, analysis-ready datasets

source("../../source.R")
pacman::p_load(dplyr, tidyr, stringr, lubridate)

cat("🧹 Starting truncated normal hindcast tidying process...\n")

# Check if hindcast files exist
hindcast_dir <- here("data/hindcasts/truncated_normal")
if (!dir.exists(hindcast_dir)) {
  cat("❌ Hindcast directory not found:", hindcast_dir, "\n")
  cat("Creating directory for future use...\n")
  dir.create(hindcast_dir, recursive = TRUE, showWarnings = FALSE)
}

# Look for existing hindcast files
hindcast_files <- list.files(hindcast_dir, pattern = "*.rds", full.names = TRUE, recursive = TRUE)

if (length(hindcast_files) == 0) {
  cat("⚠️  No hindcast files found in:", hindcast_dir, "\n")
  cat("This is expected in testing mode.\n")
  
  # NO DUMMY DATA - stop execution if no real hindcast files exist
  cat("❌ CRITICAL ERROR: No hindcast files found and dummy data generation is disabled\n")
  cat("Please run the hindcast generation script first to create real hindcast data\n")
  stop("Cannot proceed without real hindcast data")
} else {
  cat("📁 Found", length(hindcast_files), "hindcast files\n")
}

# Process each hindcast file
cat("🔄 Processing hindcast files...\n")

all_hindcasts <- list()
for (i in seq_along(hindcast_files)) {
  file_path <- hindcast_files[i]
  file_name <- basename(file_path)
  
  cat("  Processing:", file_name, "\n")
  
  tryCatch({
    # Read the hindcast file
    hindcast_data <- readRDS(file_path)
    
    if (is.data.frame(hindcast_data)) {
      # If it's already a dataframe, use it directly
      processed_data <- hindcast_data
    } else if (is.list(hindcast_data)) {
      # If it's a list, try to extract the relevant components
      if ("hindcasts" %in% names(hindcast_data)) {
        processed_data <- hindcast_data$hindcasts
      } else if ("predictions" %in% names(hindcast_data)) {
        processed_data <- hindcast_data$predictions
      } else {
        # Try to find any dataframe in the list
        df_indices <- sapply(hindcast_data, is.data.frame)
        if (any(df_indices)) {
          processed_data <- hindcast_data[[which(df_indices)[1]]]
        } else {
          cat("    ⚠️  No dataframe found in hindcast data\n")
          next
        }
      }
    } else {
      cat("    ⚠️  Unexpected data type:", class(hindcast_data), "\n")
      next
    }
    
    # Ensure we have the required columns
    required_cols <- c("model_id", "siteID", "plotID", "dateID")
    missing_cols <- required_cols[!required_cols %in% colnames(processed_data)]
    
    if (length(missing_cols) > 0) {
      cat("    ⚠️  Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
      cat("    Available columns:", paste(colnames(processed_data), collapse = ", "), "\n")
      next
    }
    
    # Add file identifier
    processed_data$source_file <- file_name
    
    # Store in list
    all_hindcasts[[i]] <- processed_data
    
    cat("    ✅ Successfully processed\n")
    
  }, error = function(e) {
    cat("    ❌ Error processing file:", e$message, "\n")
  })
}

# Combine all hindcasts
if (length(all_hindcasts) > 0) {
  cat("🔗 Combining hindcast data...\n")
  
  # Remove NULL entries
  all_hindcasts <- all_hindcasts[!sapply(all_hindcasts, is.null)]
  
  if (length(all_hindcasts) > 0) {
    # Combine dataframes
    combined_hindcasts <- bind_rows(all_hindcasts, .id = "file_index")
    
    # Clean up the data
    cat("🧹 Cleaning hindcast data...\n")
    
    # Remove duplicate rows
    initial_rows <- nrow(combined_hindcasts)
    combined_hindcasts <- combined_hindcasts %>% distinct()
    final_rows <- nrow(combined_hindcasts)
    cat("  Removed", initial_rows - final_rows, "duplicate rows\n")
    
    # Ensure numeric columns are numeric
    numeric_cols <- c("predicted_abundance", "observed_abundance", 
                     "prediction_interval_lower", "prediction_interval_upper")
    for (col in numeric_cols) {
      if (col %in% colnames(combined_hindcasts)) {
        combined_hindcasts[[col]] <- as.numeric(combined_hindcasts[[col]])
      }
    }
    
    # Add date information if dateID exists
    if ("dateID" %in% colnames(combined_hindcasts)) {
      combined_hindcasts$date <- as.Date(paste0(combined_hindcasts$dateID, "01"), format = "%Y%m%d")
      combined_hindcasts$year <- year(combined_hindcasts$date)
      combined_hindcasts$month <- month(combined_hindcasts$date)
    }
    
    # Save the combined data
    output_file <- here("data/hindcasts/truncated_normal/combined_hindcasts.rds")
    saveRDS(combined_hindcasts, output_file)
    
    cat("✅ Combined hindcasts saved to:", output_file, "\n")
    cat("📊 Final dataset has", nrow(combined_hindcasts), "rows and", ncol(combined_hindcasts), "columns\n")
    
    # Show summary statistics
    cat("\n📈 Hindcast Summary:\n")
    cat("  Models:", length(unique(combined_hindcasts$model_id)), "\n")
    cat("  Sites:", length(unique(combined_hindcasts$siteID)), "\n")
    cat("  Plots:", length(unique(combined_hindcasts$plotID)), "\n")
    if ("dateID" %in% colnames(combined_hindcasts)) {
      cat("  Date range:", min(combined_hindcasts$dateID), "to", max(combined_hindcasts$dateID), "\n")
    }
    
  } else {
    cat("❌ No valid hindcast data to combine\n")
  }
} else {
  cat("❌ No hindcast files were successfully processed\n")
}

cat("\n🎉 Truncated normal hindcast tidying process completed!\n")

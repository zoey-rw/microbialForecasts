# Robust version of add_scoring_metrics function that handles edge cases
robust_add_scoring_metrics = function(observed,
															 median_predicted,
															 mean_predicted,
															 sd_predicted,
															 type=c("RMSE","BIAS","MAE",
															 			 "CRPS", "RSQ", "RSQ.1",
															 			 "RMSE.norm","RMSE.norm.orig",
															 			 "RMSE.iqr",
															 			 "residual_variance",
															 			 "predictive_variance",
															 			 "total_PL", "CRPS_truncated"),
															 			 use_median=TRUE){

	require(Metrics, scoringRules)

	# Input validation and cleaning
	if(length(observed) == 0 || length(mean_predicted) == 0) {
		stop('Error: Empty input vectors.')
	}
	
	if(sum(is.na(observed)) > 0 || sum(is.na(mean_predicted)) > 0) {
		stop('Error: NAs in observed or predicted vectors.')
	}
	
	# Ensure all vectors have the same length
	n <- length(observed)
	if(length(mean_predicted) != n || length(median_predicted) != n || length(sd_predicted) != n) {
		stop('Error: All input vectors must have the same length.')
	}
	
	# Remove any infinite values
	finite_mask <- is.finite(observed) & is.finite(mean_predicted) & is.finite(median_predicted) & is.finite(sd_predicted)
	if(sum(finite_mask) == 0) {
		stop('Error: No finite values in input vectors.')
	}
	
	# Apply mask
	observed_clean <- observed[finite_mask]
	mean_predicted_clean <- mean_predicted[finite_mask]
	median_predicted_clean <- median_predicted[finite_mask]
	sd_predicted_clean <- sd_predicted[finite_mask]
	
	# Check if we have enough data
	if(length(observed_clean) < 2) {
		stop('Error: Need at least 2 finite values for calculations.')
	}
	
	# Ensure sd values are positive
	sd_predicted_clean <- pmax(sd_predicted_clean, 1e-10)  # Small positive value to avoid division by zero

	# These CRPS stats require distributions for each forecast
	# and cannot be calculated from median
	tryCatch({
		CRPS <- mean(scoringRules::crps_norm(observed_clean, mean_predicted_clean, sd_predicted_clean))
	}, error = function(e) {
		cat("Warning: CRPS calculation failed, setting to NA. Error:", e$message, "\n")
		CRPS <<- NA
	})
	
	tryCatch({
		CRPS_truncated <- mean(scoringRules::crps(observed_clean, family = "tnorm", 
												 location = mean_predicted_clean, 
												 scale = sd_predicted_clean, 
												 lower = 0, upper = 1))
	}, error = function(e) {
		cat("Warning: CRPS truncated calculation failed, setting to NA. Error:", e$message, "\n")
		CRPS_truncated <<- NA
	})

	# The rest of these metrics can use forecast median as the best estimate
	if (use_median==TRUE) mean_predicted_clean = median_predicted_clean

	# Calculate other metrics with error handling
	tryCatch({
		RMSE = Metrics::rmse(actual = observed_clean, predicted = mean_predicted_clean)
		RSQ.1 = 1 - (RMSE^2)/var(observed_clean)
		BIAS = Metrics::bias(actual = observed_clean, predicted = mean_predicted_clean)
		MAE = Metrics::mae(actual = observed_clean, predicted = mean_predicted_clean)
		RSQ = summary(lm(observed_clean ~ mean_predicted_clean))$r.squared
		mean_abundance = mean(observed_clean, na.rm=T)
		abundance = ifelse(mean_abundance < .005, .005, mean_abundance)
		q1 = quantile(observed_clean, .25)
		q3 = quantile(observed_clean, .75)
		IQR = q3-q1
		RMSE.iqr = RMSE/IQR
		RMSE.norm = RMSE/abundance
	}, error = function(e) {
		cat("Warning: Some metrics calculation failed, setting to NA. Error:", e$message, "\n")
		RMSE <<- NA; RSQ.1 <<- NA; BIAS <<- NA; MAE <<- NA; RSQ <<- NA
		mean_abundance <<- NA; abundance <<- NA; q1 <<- NA; q3 <<- NA; IQR <<- NA
		RMSE.iqr <<- NA; RMSE.norm <<- NA
	})

	# Create output data frame with all expected columns
	out_df = data.frame(
		RMSE = RMSE,
		BIAS = BIAS,
		MAE = MAE,
		CRPS = CRPS,
		RSQ = RSQ,
		RSQ.1 = RSQ.1,
		RMSE.iqr = RMSE.iqr,
		RMSE.norm = RMSE.norm,
		CRPS_truncated = CRPS_truncated,
		residual_variance = var(observed_clean - mean_predicted_clean, na.rm=TRUE),
		predictive_variance = var(mean_predicted_clean, na.rm=TRUE),
		total_PL = var(observed_clean - mean_predicted_clean, na.rm=TRUE) + var(mean_predicted_clean, na.rm=TRUE)
	)

	# Lower limit if RSQ 1:1 is 0
	out_df$RSQ.1 = ifelse(out_df$RSQ.1 < 0, 0, out_df$RSQ.1)

	# Upper limit of RMSE.normalized is 5
	out_df$RMSE.norm = ifelse(out_df$RMSE.norm > 5, 5, out_df$RMSE.norm)

	return(out_df)
}

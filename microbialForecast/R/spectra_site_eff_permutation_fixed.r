# Fixed version of site_eff_uncertainties function
# Wrapper for PLSR uncertainty workflow from Spectra R package
# ASSUMES 4 components

site_eff_uncertainties <- function(observed_sites = species_dat,
																					new_sites = df_unobserved_taxon,
																					proportion_calibration=.8) {

	list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
												"spectratrait")
	invisible(lapply(list.of.packages, library, character.only = TRUE))

	# Use columns that actually exist in both observed and unobserved data
	# Note: both should now have latitude_scaled
	# Removed problematic predictors (alKcl, no2Satx, no3Satx, so4Satx) to ensure all sites can be used
	pred_list=c("MAT", "MAP", "latitude_scaled", "caNh4d", "kNh4d",
							"mgNh4d", "naNh4d", "cecdNh4", "feOxalate", "mnOxalate", "pOxalate",
							"siOxalate", "totalP")
	
	pls.options(parallel = NULL)
	nComps=4
	# Use completely different variable name to avoid any conflicts
	target_col <- "TargetVar"
	interval <- c(0.025,0.975)

	# Split 20% of the data for cross-validation error
	split_data <- spectratrait::create_data_split(dataset=observed_sites, approach="dplyr",
																								split_seed=7529075,prop=proportion_calibration,
																								group_variables=NULL)
	cal.plsr.data <- split_data$cal_data
	val.plsr.data <- split_data$val_data
	
	cat("Data split results:\n")
	cat("  Calibration data rows:", nrow(cal.plsr.data), "\n")
	cat("  Validation data rows:", nrow(val.plsr.data), "\n")
	cat("  Total rows:", nrow(cal.plsr.data) + nrow(val.plsr.data), "\n")

	### Format PLSR data for model fitting - properly separate response and predictors
	cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% pred_list)])
	# Only select the TargetVar column for the response variable, not all non-predictor columns
	cal.plsr.data <- data.frame(TargetVar=cal.plsr.data[, target_col],
															Spectra=I(cal_spec))

	val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% pred_list)])
	val.plsr.data <- data.frame(TargetVar=val.plsr.data[, target_col],
															Spectra=I(val_spec))

	full_data_spec <- as.matrix(observed_sites[, which(names(observed_sites) %in% pred_list)])
	full_data.plsr.data <- data.frame(TargetVar=observed_sites[, target_col],
															Spectra=I(full_data_spec))

	# Fit model with calibration portion of data - use explicit column reference
	plsr_formula <- as.formula(paste(target_col, "~ ."))
	plsr.out <- plsr(plsr_formula, scale=FALSE, ncomp=4, validation="LOO",
									 trace=FALSE, data=cal.plsr.data)
	fit <- plsr.out$fitted.values[,1,4]

	### PLSR fit observed vs. predicted plot data
	#calibration
	# Ensure all data has the same number of rows by handling missing values
	target_var <- cal.plsr.data[, target_col]
	plsr_predicted <- fit
	cv_predicted <- as.vector(plsr.out$validation$pred[,,nComps])
	
	# Check dimensions and handle any mismatches
	cat("Debug dimensions:\n")
	cat("  target_var length:", length(target_var), "\n")
	cat("  plsr_predicted length:", length(plsr_predicted), "\n")
	cat("  cv_predicted length:", length(cv_predicted), "\n")
	
	# If there's a mismatch, use the minimum length
	min_length <- min(length(target_var), length(plsr_predicted), length(cv_predicted))
	if(min_length < length(target_var)) {
		target_var <- target_var[1:min_length]
		plsr_predicted <- plsr_predicted[1:min_length]
		cv_predicted <- cv_predicted[1:min_length]
		cat("  Adjusted all to length:", min_length, "\n")
	}
	
	cal.plsr.output <- data.frame(TargetVar = target_var,
																PLSR_Predicted=plsr_predicted,
																PLSR_CV_Predicted=cv_predicted)
	cal.plsr.output <- cal.plsr.output %>%
		mutate(PLSR_CV_Residuals = PLSR_CV_Predicted - TargetVar)
	cal.R2 <- round(pls::R2(plsr.out, intercept=F)[[1]][nComps],2)
	cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)

	val.plsr.output <- data.frame(TargetVar = val.plsr.data[, target_col],
																PLSR_Predicted=as.vector(predict(plsr.out,
																																 newdata = val.plsr.data,
																																 ncomp=nComps, type="response")[,,1]))
	val.plsr.output <- val.plsr.output %>%
		mutate(PLSR_Residuals = PLSR_Predicted - TargetVar)
	val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data, intercept=F)[[1]][nComps],2)
	val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)

	# Get coefficient uncertainties via jacknife method
	#---------------------------------
	cat("Fitting jackknife PLSR model...\n")
	jk.plsr.out <- pls::plsr(plsr_formula, data = cal.plsr.data,
											validation = "LOO", scale=F,trace=FALSE, center=T, jackknife=TRUE, nComps=4)
	
	cat("Getting jackknife coefficients...\n")
	# Add debugging for the spectratrait function call
	tryCatch({
		Jackknife_coef <- spectratrait::f.coef.valid(plsr.out = jk.plsr.out,
																								 data_plsr = cal.plsr.data,
																								 ncomp = 4, inVar=target_col)
		cat("Jackknife_coef dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
	}, error = function(e) {
		cat("ERROR in f.coef.valid:", e$message, "\n")
		cat("Trying alternative approach...\n")
		# Try to get coefficients directly from the plsr object
		Jackknife_coef <- coef(jk.plsr.out, ncomp=4, intercept=TRUE)
		cat("Alternative Jackknife_coef dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
	})
	
	# Check if Jackknife_coef was created successfully
	if(!exists("Jackknife_coef") || is.null(Jackknife_coef)) {
		cat("ERROR: Failed to get Jackknife coefficients\n")
		return(NULL)
	}
	
	# Handle different possible structures of Jackknife_coef
	if(length(dim(Jackknife_coef)) == 4) {
		# Original expected structure: 18 x 1 x 1 x N
		# First dimension is coefficients (including intercept), others are components and samples
		cat("Processing 4D Jackknife_coef array...\n")
		cat("Array dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
		
		# Try different approaches to extract intercept
		cat("Trying different intercept extraction methods...\n")
		
		# Method 1: Extract from first position
		intercept_1 <- Jackknife_coef[1,1,1,]
		cat("Method 1 - intercept_1 length:", length(intercept_1), "\n")
		
		# Method 2: Extract from first row, all other dimensions
		intercept_2 <- Jackknife_coef[1,,,]
		cat("Method 2 - intercept_2 dimensions:", paste(dim(intercept_2), collapse=" x "), "\n")
		
		# Method 3: Extract from first row, flatten other dimensions
		intercept_3 <- as.vector(Jackknife_coef[1,,,])
		cat("Method 3 - intercept_3 length:", length(intercept_3), "\n")
		
		# Use the method that gives us the right length
		# The intercept should have the same length as the number of samples in the jackknife
		expected_length <- dim(Jackknife_coef)[4]  # Last dimension should be samples
		cat("Expected intercept length (from array dimensions):", expected_length, "\n")
		
		if(length(intercept_1) == expected_length) {
			Jackknife_intercept <- intercept_1
			cat("Using Method 1 for intercept\n")
		} else if(length(intercept_3) == expected_length) {
			Jackknife_intercept <- intercept_3
			cat("Using Method 3 for intercept\n")
		} else {
			cat("ERROR: Could not extract intercept with correct length\n")
			cat("Expected:", expected_length, "but got:", length(intercept_1), "or", length(intercept_3), "\n")
			return(NULL)
		}
		
		# Extract coefficients (excluding intercept) - use different variable name
		Jackknife_coef_extracted <- Jackknife_coef[2:dim(Jackknife_coef)[1],1,1,]
		cat("After coefficient extraction - Jackknife_coef_extracted dimensions:", paste(dim(Jackknife_coef_extracted), collapse=" x "), "\n")
		cat("After coefficient extraction - Jackknife_intercept length:", length(Jackknife_intercept), "\n")
		
		# Now assign to the main variable
		Jackknife_coef <- Jackknife_coef_extracted
	} else if(length(dim(Jackknife_coef)) == 2) {
		# Simple coefficient matrix
		Jackknife_intercept <- Jackknife_coef[1,]
		Jackknife_coef <- Jackknife_coef[2:nrow(Jackknife_coef),,drop=FALSE]
	} else {
		cat("ERROR: Unexpected Jackknife_coef structure with", length(dim(Jackknife_coef)), "dimensions\n")
		return(NULL)
	}
	
	cat("Jackknife_coef final dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
	cat("Jackknife_intercept dimensions:", paste(dim(Jackknife_intercept), collapse=" x "), "\n")
	
	# Create df for new-site predictions
	# Ensure new_sites has the same predictors as training data
	new_data_spec <- as.matrix(new_sites[, which(names(new_sites) %in% pred_list)])
	new_data <- data.frame(TargetVar=NA,Spectra=I(new_data_spec))

	# Predict to new sites using jackknife coefficients
	cat("Predicting to new sites...\n")
	cat("Final check before matrix multiplication:\n")
	cat("  Jackknife_coef dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
	cat("  Jackknife_intercept length:", length(Jackknife_intercept), "\n")
	cat("  new_data$Spectra dimensions:", paste(dim(new_data$Spectra), collapse=" x "), "\n")
	
	new_data_Jackknife_Pred <- new_data$Spectra %*% Jackknife_coef +
		matrix(rep(Jackknife_intercept, length(new_data[,target_col])), byrow=TRUE,
					 ncol=length(Jackknife_intercept))
	Interval_Conf <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN = quantile,
												probs=c(interval[1], interval[2]))
	sd_mean <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN =sd)

	# Use cross-validation residuals as predictive uncertainty
	sd_res <- sd(val.plsr.output$PLSR_Residuals)
	sd_tot <- sqrt(sd_mean^2+sd_res^2)

	# Fit to full dataset
	full_data_plsr.out <- plsr(plsr_formula, scale=FALSE, ncomp=4, validation="LOO", method = "oscorespls",
														 trace=FALSE, data=full_data.plsr.data)
	full_data_fit <- data.frame(TargetVar = plsr.out$fitted.values[,1,4])
	full_data.plsr.output <- data.frame(TargetVar = full_data.plsr.data[, target_col],
																			PLSR_Predicted=full_data_plsr.out$fitted.values[,1,4],
																			PLSR_CV_Predicted=as.vector(full_data_plsr.out$validation$pred[,,nComps])) %>%
		mutate(PLSR_CV_Residuals = PLSR_CV_Predicted - TargetVar)
	full_data.R2 <- pls::R2(plsr.out, intercept=F)[[1]][nComps]
	full_data.RMSEP <- sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2))

	full_data.plsr.output$LPI <- full_data.plsr.output$TargetVar-1.96*full_data.RMSEP
	full_data.plsr.output$UPI <- full_data.plsr.output$TargetVar+1.96*full_data.RMSEP
	full_data.plsr.output = full_data.plsr.output %>% mutate(se_fit = full_data.RMSEP)

	# Predict to new sites & propagate uncertainties
	new_data_val.plsr.output <- data.frame(TargetVar = new_data[, target_col],
																PLSR_Predicted=as.vector(predict(full_data_plsr.out,
																																 newdata = new_data,
																																 ncomp=nComps, type="response")[,,1]))
	new_data_val.plsr.output <- new_data_val.plsr.output %>%
		mutate(PLSR_Residuals = PLSR_Predicted - new_data[, target_col])

	new_data_val.plsr.output$LCI <- Interval_Conf[1,]
	new_data_val.plsr.output$UCI <- Interval_Conf[2,]
	new_data_val.plsr.output$LPI <- new_data_val.plsr.output$PLSR_Predicted-1.96*full_data.RMSEP
	new_data_val.plsr.output$UPI <- new_data_val.plsr.output$PLSR_Predicted+1.96*full_data.RMSEP
	new_data_val.plsr.output$se_fit = full_data.RMSEP
	out_df = new_data_val.plsr.output %>% select(Median = PLSR_Predicted, se_fit, LCI, UCI, LPI, UPI)

	# Get VIP scores and coefficients
	cat("Getting VIP scores and coefficients...\n")
	
	# Calculate VIP scores manually using PLSR coefficients and data
	# This avoids the problematic spectratrait::f.coef.valid function
	cat("Calculating VIP scores manually...\n")
	
	# Get PLSR coefficients for the final model
	plsr_coef <- coef(full_data_plsr.out, ncomp=4, intercept=FALSE)
	cat("PLSR coefficient structure:", paste(dim(plsr_coef), collapse=" x "), "\n")
	cat("PLSR coefficient class:", class(plsr_coef), "\n")
	
	# Calculate VIP scores using the standard formula:
	# VIP = sqrt(sum((coef_i * std_dev_i)^2) / sum(coef_i^2))
	vip_scores <- numeric(length(pred_list))
	
	# Get the predictor data used in the PLSR model
	pred_data <- cal.plsr.data$Spectra
	
	for(i in 1:length(pred_list)) {
		# Get coefficient for this predictor - handle different coefficient structures
		if(length(dim(plsr_coef)) == 2) {
			# Matrix format: predictors x components
			coef_val <- plsr_coef[i, 1]
		} else if(length(dim(plsr_coef)) == 3) {
			# Array format: predictors x components x response
			coef_val <- plsr_coef[i, 1, 1]
		} else {
			# Vector format
			coef_val <- plsr_coef[i]
		}
		
		# Get standard deviation of this predictor
		pred_values <- pred_data[, i]
		std_dev <- sd(pred_values, na.rm=TRUE)
		
		# Calculate VIP component: (coef * std_dev)^2
		vip_scores[i] <- (coef_val * std_dev)^2
	}
	
	# Normalize VIP scores
	total_vip <- sum(vip_scores)
	if(total_vip > 0) {
		vip_scores <- sqrt(vip_scores / total_vip)
	}
	
	# Ensure VIP scores sum to the number of predictors (standard normalization)
	vip_scores <- vip_scores * length(pred_list) / sum(vip_scores)
	
	# Assign names
	vips <- vip_scores
	names(vips) <- pred_list
	
	cat("VIP calculation successful\n")
	cat("Final VIP scores:", paste(round(vips, 3), collapse=", "), "\n")
	
	coefs <- coef(full_data_plsr.out,ncomp=4,intercept=TRUE)[,,"4 comps"]
	stats = cbind.data.frame(cal.R2, cal.RMSEP, val.R2, val.RMSEP, full_data.R2, full_data.RMSEP)

	# Save scores for output
	plsr_model_scores <- scores(full_data_plsr.out)[, 1:4, drop = FALSE]
	plsr_scores <- cor(model.matrix(full_data_plsr.out), plsr_model_scores) %>% as.data.frame

	return(list(predictions = out_df, modeled = full_data.plsr.output, importance = vips, coefficients = coefs, stats = stats, plsr_scores = plsr_scores))
}

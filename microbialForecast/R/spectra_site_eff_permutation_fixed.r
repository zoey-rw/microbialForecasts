# Fixed version of site_eff_uncertainties function
# Wrapper for PLSR uncertainty workflow from Spectra R package
# Fits up to 3 PLSR components; returns per-component LOO-CV stats for analysis
# Uses scale=TRUE and drops predictors with low training-data variance (SD < 0.3)
# to avoid inflating influence of predictors that are nearly constant across observed sites.
#
# Prediction uncertainty is leverage-aware: sites far from the training data
# (in predictor space) receive wider uncertainty via the hat-matrix diagonal.
# se_fit_i = sigma * sqrt(1 + h_i) where h_i is the leverage of site i and
# sigma is the LOO-CV RMSEP.

site_eff_uncertainties <- function(observed_sites = species_dat,
																					new_sites = df_unobserved_taxon,
																					min_predictor_sd = 0.3) {

	list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
												"spectratrait")
	invisible(lapply(list.of.packages, library, character.only = TRUE))

	# Use columns that actually exist in both observed and unobserved data
	# Note: both should now have latitude_scaled
	# Mehlich III extractable metals (feMjelm, mnMjelm, pMjelm, siMjelm) replace
	# former oxalate columns (feOxalate, mnOxalate, pOxalate, siOxalate) from megapit update
	pred_list=c("MAT", "MAP", "latitude_scaled", "caNh4d", "kNh4d",
							"mgNh4d", "naNh4d", "cecdNh4", "feMjelm", "mnMjelm", "pMjelm",
							"siMjelm", "totalP")

	pls.options(parallel = NULL)
	nComps=3
	# Use completely different variable name to avoid any conflicts
	target_col <- "TargetVar"
	interval <- c(0.025,0.975)

	# Filter to predictors present in both datasets
	pred_list <- pred_list[pred_list %in% names(observed_sites) & pred_list %in% names(new_sites)]

	# Drop predictors with low training-data SD
	# Predictors are pre-scaled to SD=1 across all sites, but the training subset
	# (observed sites) may have much lower variance. When PLSR scales internally,
	# low-SD predictors get inflated and dominate predictions at new sites.
	train_sds <- apply(observed_sites[, pred_list, drop=FALSE], 2, sd, na.rm=TRUE)
	low_sd <- names(train_sds[train_sds < min_predictor_sd])
	if (length(low_sd) > 0) {
		cat("Dropping low-variance predictors (train SD <", min_predictor_sd, "):", paste(low_sd, collapse=", "), "\n")
		pred_list <- setdiff(pred_list, low_sd)
	}

	if (length(pred_list) < 2) {
		cat("ERROR: Too few predictors remaining after filtering\n")
		return(NULL)
	}

	# Use ALL observed sites for PLSR training (no cal/val split)
	# LOO cross-validation provides residual uncertainty without holding out data
	full_data_spec <- as.matrix(observed_sites[, pred_list])
	full_data.plsr.data <- data.frame(TargetVar=observed_sites[, target_col],
															Spectra=I(full_data_spec))

	cat("PLSR training with", nrow(full_data.plsr.data), "sites,", length(pred_list), "predictors\n")

	plsr_formula <- as.formula(paste(target_col, "~ ."))

	# Get coefficient uncertainties via jackknife method
	# scale=TRUE: PLSR internally standardizes predictors using training-data statistics,
	# ensuring each predictor contributes proportionally regardless of its raw variance.
	#---------------------------------
	cat("Fitting jackknife PLSR model (scale=TRUE)...\n")
	jk.plsr.out <- pls::plsr(plsr_formula, data = full_data.plsr.data,
											validation = "LOO", scale=TRUE, trace=FALSE, center=TRUE, jackknife=TRUE, nComps=nComps)

	# Add debugging for the spectratrait function call
	tryCatch({
		Jackknife_coef <- spectratrait::f.coef.valid(plsr.out = jk.plsr.out,
																								 data_plsr = full_data.plsr.data,
																								 ncomp = nComps, inVar=target_col)
		cat("Jackknife_coef dimensions:", paste(dim(Jackknife_coef), collapse=" x "), "\n")
	}, error = function(e) {
		cat("ERROR in f.coef.valid:", e$message, "\n")
		cat("Trying alternative approach...\n")
		# Try to get coefficients directly from the plsr object
		Jackknife_coef <- coef(jk.plsr.out, ncomp=nComps, intercept=TRUE)
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

		# Method 1: Extract from first position
		intercept_1 <- Jackknife_coef[1,1,1,]
		# Method 2: Extract from first row, all other dimensions
		intercept_2 <- Jackknife_coef[1,,,]

		# Method 3: Extract from first row, flatten other dimensions
		intercept_3 <- as.vector(Jackknife_coef[1,,,])

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
	# Use the same filtered pred_list as training data
	new_data_spec <- as.matrix(new_sites[, pred_list])
	new_data <- data.frame(TargetVar=NA,Spectra=I(new_data_spec))

	# Predict to new sites using jackknife coefficients
	cat("Predicting to new sites...\n")
	new_data_Jackknife_Pred <- new_data$Spectra %*% Jackknife_coef +
		matrix(rep(Jackknife_intercept, length(new_data[,target_col])), byrow=TRUE,
					 ncol=length(Jackknife_intercept))
	Interval_Conf <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN = quantile,
												probs=c(interval[1], interval[2]))
	sd_mean <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN =sd)

	# Fit PLSR on all training data (same data used for jackknife above)
	full_data_plsr.out <- plsr(plsr_formula, scale=TRUE, ncomp=nComps, validation="LOO", method = "oscorespls",
														 trace=FALSE, data=full_data.plsr.data)

	# Compute LOO-CV stats for each number of components (for later analysis)
	y_obs <- full_data.plsr.data[, target_col]
	ss_tot <- sum((y_obs - mean(y_obs))^2)
	cv_stats_per_ncomp <- data.frame(
		ncomp = seq_len(nComps),
		loo_cv_r2 = NA_real_,
		loo_cv_rmsep = NA_real_
	)
	for (nc in seq_len(nComps)) {
		cv_pred_nc <- as.vector(full_data_plsr.out$validation$pred[,, nc])
		cv_resid_nc <- cv_pred_nc - y_obs
		cv_stats_per_ncomp$loo_cv_r2[nc] <- 1 - sum(cv_resid_nc^2) / ss_tot
		cv_stats_per_ncomp$loo_cv_rmsep[nc] <- sqrt(mean(cv_resid_nc^2))
	}

	# Use nComps for final predictions (downstream can select fewer if desired)
	full_data.plsr.output <- data.frame(TargetVar = y_obs,
																			PLSR_Predicted=full_data_plsr.out$fitted.values[,1,nComps],
																			PLSR_CV_Predicted=as.vector(full_data_plsr.out$validation$pred[,,nComps])) %>%
		mutate(PLSR_CV_Residuals = PLSR_CV_Predicted - TargetVar)
	full_data.R2 <- pls::R2(full_data_plsr.out, intercept=F)[[1]][nComps]
	full_data.RMSEP <- sqrt(mean(full_data.plsr.output$PLSR_CV_Residuals^2))

	cat("PLSR full-data LOO CV R2:", round(full_data.R2, 3), "RMSEP:", round(full_data.RMSEP, 3), "\n")
	cat("LOO-CV R2 by ncomp:", paste(round(cv_stats_per_ncomp$loo_cv_r2, 3), collapse=", "), "\n")

	# --- Leverage-based per-site prediction uncertainty ---
	# PLSR projects predictors into a score space T = X W* where W* are loading weights.
	# Prediction uncertainty depends on how far each new site is from the training data
	# in this score space (leverage). We compute the hat matrix in score space:
	#   H = T (T'T)^{-1} T'  for training data
	#   h_new_i = t_new_i' (T'T)^{-1} t_new_i  for new sites
	# Then se_fit_i = sigma * sqrt(1 + h_new_i) accounts for both residual variance
	# and extrapolation distance.
	sigma_cv <- full_data.RMSEP  # LOO-CV residual standard error

	# Get PLSR scores for training data
	TT_scores <- scores(full_data_plsr.out)[, 1:nComps, drop = FALSE]
	# Regularized (T'T)^{-1} for numerical stability
	TtT <- t(TT_scores) %*% TT_scores
	TtT_inv <- tryCatch(
		solve(TtT + 1e-6 * diag(nComps)),
		error = function(e) {
			cat("WARNING: Score matrix near-singular, using diagonal approximation\n")
			diag(1 / (diag(TtT) + 1e-6))
		}
	)

	# Training leverages (for diagnostics and modeled output)
	h_train <- rowSums((TT_scores %*% TtT_inv) * TT_scores)
	se_fit_train <- sigma_cv * sqrt(1 + h_train)

	full_data.plsr.output$LPI <- full_data.plsr.output$TargetVar - 1.96 * se_fit_train
	full_data.plsr.output$UPI <- full_data.plsr.output$TargetVar + 1.96 * se_fit_train
	full_data.plsr.output$se_fit <- se_fit_train

	cat("Training leverage: median=", round(median(h_train), 3),
			", max=", round(max(h_train), 3), "\n")

	# Predict to new sites & propagate uncertainties
	new_data_val.plsr.output <- data.frame(TargetVar = new_data[, target_col],
																PLSR_Predicted=as.vector(predict(full_data_plsr.out,
																																 newdata = new_data,
																																 ncomp=nComps, type="response")[,,1]))
	new_data_val.plsr.output <- new_data_val.plsr.output %>%
		mutate(PLSR_Residuals = PLSR_Predicted - new_data[, target_col])

	# New-site leverages in PLSR score space
	# Project new data into PLSR score space using the model
	new_scores <- predict(full_data_plsr.out, newdata = new_data, type = "scores")
	if (length(dim(new_scores)) == 3) new_scores <- new_scores[, 1:nComps, 1]
	else new_scores <- new_scores[, 1:nComps, drop = FALSE]
	h_new <- rowSums((new_scores %*% TtT_inv) * new_scores)
	se_fit_new <- sigma_cv * sqrt(1 + h_new)

	cat("New-site leverage: median=", round(median(h_new), 3),
			", max=", round(max(h_new), 3), "\n")
	cat("New-site se_fit: median=", round(median(se_fit_new), 3),
			", range=[", round(min(se_fit_new), 3), ",", round(max(se_fit_new), 3), "]\n")

	new_data_val.plsr.output$LCI <- Interval_Conf[1,]
	new_data_val.plsr.output$UCI <- Interval_Conf[2,]
	new_data_val.plsr.output$LPI <- new_data_val.plsr.output$PLSR_Predicted - 1.96 * se_fit_new
	new_data_val.plsr.output$UPI <- new_data_val.plsr.output$PLSR_Predicted + 1.96 * se_fit_new
	new_data_val.plsr.output$se_fit <- se_fit_new
	out_df = new_data_val.plsr.output %>% select(Median = PLSR_Predicted, se_fit, LCI, UCI, LPI, UPI)

	# Compute proper PLSR VIP scores using the standard formula
	# VIP_j = sqrt(p * sum_a(w_ja^2 * SS_a) / sum_a(SS_a))
	# where w_ja = loading weight for variable j in component a
	# and SS_a = explained sum of squares by component a
	cat("Computing VIP scores...\n")
	W <- loading.weights(full_data_plsr.out)  # p x nComps loading weights matrix
	TT <- scores(full_data_plsr.out)          # n x nComps scores matrix
	Q <- loadings(full_data_plsr.out)         # not needed directly but used via SS

	p <- ncol(full_data_spec)

	# SS per component: variance explained in Y by each component
	SS <- numeric(nComps)
	for (a in seq_len(nComps)) {
		# t_a' * y / (t_a' * t_a) gives the Y loading for component a
		q_a <- sum(TT[, a] * y_obs) / sum(TT[, a]^2)
		SS[a] <- q_a^2 * sum(TT[, a]^2)
	}

	vip_scores <- numeric(p)
	total_SS <- sum(SS)
	for (j in seq_len(p)) {
		vip_scores[j] <- sqrt(p * sum(W[j, ]^2 * SS) / total_SS)
	}
	names(vip_scores) <- colnames(full_data_spec)

	cat("VIP scores:", paste(round(vip_scores, 3), collapse=", "), "\n")

	coefs <- coef(full_data_plsr.out, ncomp=nComps, intercept=TRUE)[,,paste0(nComps, " comps")]
	stats <- cbind.data.frame(full_data.R2, full_data.RMSEP)

	# Save scores for output
	plsr_model_scores <- scores(full_data_plsr.out)[, 1:nComps, drop = FALSE]
	plsr_scores <- cor(model.matrix(full_data_plsr.out), plsr_model_scores) %>% as.data.frame

	return(list(predictions = out_df, modeled = full_data.plsr.output,
	            importance = vip_scores, coefficients = coefs, stats = stats,
	            plsr_scores = plsr_scores, cv_stats_per_ncomp = cv_stats_per_ncomp))
}

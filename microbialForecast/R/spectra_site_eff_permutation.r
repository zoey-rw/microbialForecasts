# very messy, sorry
# Wrapper for PLSR uncertainty workflow from Spectra R package
# ASSUMES 4 components

# proportion_calibration=.8
# testing = site_eff_uncertainties(observed_sites = species_dat, new_sites = df_unobserved_taxon)
#
# ggplot(testing$predictions) + geom_line(aes(x = 1:16, y = Median)) +
# 	geom_ribbon(aes(x = 1:16, ymin=LCI, ymax=UCI),alpha=.3) +
# 	geom_ribbon(aes(x = 1:16, ymin=LPI, ymax=UPI),alpha=.1)

site_eff_uncertainties <- function(observed_sites = species_dat,
																					new_sites = df_unobserved_taxon,
																					proportion_calibration=.8) {

	list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
												"spectratrait")
	invisible(lapply(list.of.packages, library, character.only = TRUE))


	pred_list=c("MAT", "MAP", "latitude_scaled", "caNh4d", "kNh4d",
							"mgNh4d", "naNh4d", "cecdNh4", "feOxalate", "mnOxalate", "pOxalate",
							"siOxalate", "alKcl", "no2Satx", "no3Satx", "so4Satx",
							"totalP")
	pls.options(parallel = NULL)
	nComps=4
	inVar="Median"
	interval <- c(0.025,0.975)

# Split 20% of the data for cross-validation error
split_data <- spectratrait::create_data_split(dataset=observed_sites, approach="dplyr",
																							split_seed=7529075,prop=proportion_calibration,
																							group_variables=NULL)
cal.plsr.data <- split_data$cal_data
val.plsr.data <- split_data$val_data


### Format PLSR data for model fitting
cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% pred_list)])
cal.plsr.data <- data.frame(Median=cal.plsr.data[, which(names(cal.plsr.data) %notin% pred_list)],
														Spectra=I(cal_spec))

val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% pred_list)])
val.plsr.data <- data.frame(Median=val.plsr.data[, which(names(val.plsr.data) %notin% pred_list)],
														Spectra=I(val_spec))


full_data_spec <- as.matrix(observed_sites[, which(names(observed_sites) %in% pred_list)])
full_data.plsr.data <- data.frame(Median=observed_sites[, which(names(observed_sites) %notin% pred_list)],
														Spectra=I(full_data_spec))

# Fit model with calibration portion of data
plsr.out <- plsr(as.formula(Median ~ .),scale=FALSE,ncomp=4,validation="LOO",
								 trace=FALSE,data=cal.plsr.data)
fit <- plsr.out$fitted.values[,1,4]


### PLSR fit observed vs. predicted plot data
#calibration
cal.plsr.output <- data.frame(Median = cal.plsr.data[, "Median"],
															PLSR_Predicted=fit,
															PLSR_CV_Predicted=as.vector(plsr.out$validation$pred[,,nComps]))
cal.plsr.output <- cal.plsr.output %>%
	mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-Median)
cal.R2 <- round(pls::R2(plsr.out, intercept=F)[[1]][nComps],2)
cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)

val.plsr.output <- data.frame(Median = val.plsr.data[, "Median"],
															PLSR_Predicted=as.vector(predict(plsr.out,
																															 newdata = val.plsr.data,
																															 ncomp=nComps, type="response")[,,1]))
val.plsr.output <- val.plsr.output %>%
	mutate(PLSR_Residuals = PLSR_Predicted-Median)
val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data, intercept=F)[[1]][nComps],2)
val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)

# Get coefficient uncertainties via jacknife method
#---------------------------------
jk.plsr.out <- pls::plsr(Median ~ ., data = cal.plsr.data,
										validation = "LOO", scale=F,trace=FALSE, center=T, jackknife=TRUE, nComps=4)
Jackknife_coef <- spectratrait::f.coef.valid(plsr.out = jk.plsr.out,
																						 data_plsr = cal.plsr.data,
																						 ncomp = 4, inVar="Median")
Jackknife_intercept <- Jackknife_coef[1,,,]
Jackknife_coef <- Jackknife_coef[2:dim(Jackknife_coef)[1],,,]
Jackknife_Pred <- val.plsr.data$Spectra %*% Jackknife_coef +
	matrix(rep(Jackknife_intercept, length(val.plsr.data[,"Median"])), byrow=TRUE,
				 ncol=length(Jackknife_intercept))
Interval_Conf <- apply(X = Jackknife_Pred, MARGIN = 1, FUN = quantile,
											 probs=c(interval[1], interval[2]))
sd_mean <- apply(X = Jackknife_Pred, MARGIN = 1, FUN =sd)
sd_res <- sd(val.plsr.output$PLSR_Residuals)
sd_tot <- sqrt(sd_mean^2+sd_res^2)
val.plsr.output$LCI <- Interval_Conf[1,]
val.plsr.output$UCI <- Interval_Conf[2,]
val.plsr.output$LPI <- val.plsr.output$PLSR_Predicted-1.96*sd_tot
val.plsr.output$UPI <- val.plsr.output$PLSR_Predicted+1.96*sd_tot

# Create df for new-site predictions
new_data_spec <- as.matrix(new_sites[, which(names(new_sites) %in% pred_list)])
new_data <- data.frame(Median=NA,Spectra=I(new_data_spec))

# Predict to new sites using jackknife coefficients
new_data_Jackknife_Pred <- new_data$Spectra %*% Jackknife_coef +
	matrix(rep(Jackknife_intercept, length(new_data[,"Median"])), byrow=TRUE,
				 ncol=length(Jackknife_intercept))
Interval_Conf <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN = quantile,
											 probs=c(interval[1], interval[2]))
sd_mean <- apply(X = new_data_Jackknife_Pred, MARGIN = 1, FUN =sd)

# Use cross-validation residuals as predictive uncertainty
sd_res <- sd(val.plsr.output$PLSR_Residuals)
sd_tot <- sqrt(sd_mean^2+sd_res^2)


# Fit to full dataset
full_data_plsr.out <- plsr(as.formula(Median ~ .),scale=FALSE,ncomp=4,validation="LOO", method = "oscorespls",
													 trace=FALSE,data=full_data.plsr.data)
full_data_fit <- data.frame(Median = plsr.out$fitted.values[,1,4])
full_data.plsr.output <- data.frame(Median = full_data.plsr.data[, "Median"],
																		PLSR_Predicted=full_data_plsr.out$fitted.values[,1,4],
																		PLSR_CV_Predicted=as.vector(full_data_plsr.out$validation$pred[,,nComps])) %>%
	mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-Median)
full_data.R2 <- pls::R2(plsr.out, intercept=F)[[1]][nComps]
full_data.RMSEP <- sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2))

full_data.plsr.output$LPI <- full_data.plsr.output$Median-1.96*full_data.RMSEP
full_data.plsr.output$UPI <- full_data.plsr.output$Median+1.96*full_data.RMSEP
full_data.plsr.output = full_data.plsr.output %>% mutate(se_fit = full_data.RMSEP)


# Predict to new sites & propagate uncertainties
new_data_val.plsr.output <- data.frame(Median = new_data[, "Median"],
															PLSR_Predicted=as.vector(predict(full_data_plsr.out,
																															 newdata = new_data,
																															 ncomp=nComps, type="response")[,,1]))
new_data_val.plsr.output <- new_data_val.plsr.output %>%
	mutate(PLSR_Residuals = PLSR_Predicted-Median)

new_data_val.plsr.output$LCI <- Interval_Conf[1,]
new_data_val.plsr.output$UCI <- Interval_Conf[2,]
new_data_val.plsr.output$LPI <- new_data_val.plsr.output$PLSR_Predicted-1.96*full_data.RMSEP
new_data_val.plsr.output$UPI <- new_data_val.plsr.output$PLSR_Predicted+1.96*full_data.RMSEP
new_data_val.plsr.output$se_fit = full_data.RMSEP
out_df = new_data_val.plsr.output %>% select(Median = PLSR_Predicted, se_fit, LCI, UCI, LPI, UPI)


vips <- spectratrait::f.coef.valid(plsr.out = full_data_plsr.out,
																						 data_plsr = full_data.plsr.data,
																						 ncomp = 4, inVar="Median")
coefs <- coef(full_data_plsr.out,ncomp=4,intercept=TRUE)[,,"4 comps"]
stats = cbind.data.frame(cal.R2, cal.RMSEP, val.R2, val.RMSEP, full_data.R2, full_data.RMSEP)

# Save scores for output
plsr_model_scores <- scores(full_data_plsr.out)[, 1:4, drop = FALSE]
plsr_scores <- cor(model.matrix(full_data_plsr.out), plsr_model_scores) %>% as.data.frame

return(list(predictions = out_df, modeled = full_data.plsr.output, importance = vips, coefficients = coefs, stats = stats, plsr_scores = plsr_scores))
}



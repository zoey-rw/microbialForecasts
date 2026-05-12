
#' @title 	add_scoring_metrics
#' @description Score models and forecasts. Handles edge cases (infinite values,
#'   zero SDs, short vectors) gracefully by filtering and using tryCatch.
#'
#' @export
add_scoring_metrics = function(observed,
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

	# Input validation
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

	observed <- observed[finite_mask]
	mean_predicted <- mean_predicted[finite_mask]
	median_predicted <- median_predicted[finite_mask]
	sd_predicted <- sd_predicted[finite_mask]

	if(length(observed) < 2) {
		stop('Error: Need at least 2 finite values for calculations.')
	}

	# Ensure sd values are positive
	sd_predicted <- pmax(sd_predicted, 1e-10)

	# CRPS stats require distributions and cannot be calculated from median
	tryCatch({
		CRPS <- mean(scoringRules::crps_norm(observed, mean_predicted, sd_predicted))
	}, error = function(e) {
		cat("Warning: CRPS calculation failed, setting to NA. Error:", e$message, "\n")
		CRPS <<- NA
	})

	tryCatch({
		CRPS_truncated <- mean(scoringRules::crps(observed, family = "tnorm",
												 location = mean_predicted,
												 scale = sd_predicted,
												 lower = 0, upper = 1))
	}, error = function(e) {
		cat("Warning: CRPS truncated calculation failed, setting to NA. Error:", e$message, "\n")
		CRPS_truncated <<- NA
	})

	# The rest of these metrics can use forecast median as the best estimate
	if (use_median==TRUE) mean_predicted = median_predicted

	tryCatch({
		RMSE = Metrics::rmse(actual = observed, predicted = mean_predicted)
		RSQ.1 = 1 - (RMSE^2)/var(observed)
		BIAS = Metrics::bias(actual = observed, predicted = mean_predicted)
		MAE = Metrics::mae(actual = observed, predicted = mean_predicted)
		RSQ = summary(lm(observed ~ mean_predicted))$r.squared
		mean_abundance = mean(observed, na.rm=T)
		abundance = ifelse(mean_abundance < .005, .005, mean_abundance)
		q1 = quantile(observed, .25)
		q3 = quantile(observed, .75)
		IQR = q3-q1
		RMSE.iqr = RMSE/IQR
		RMSE.norm = RMSE/abundance
	}, error = function(e) {
		cat("Warning: Some metrics calculation failed, setting to NA. Error:", e$message, "\n")
		RMSE <<- NA; RSQ.1 <<- NA; BIAS <<- NA; MAE <<- NA; RSQ <<- NA
		mean_abundance <<- NA; abundance <<- NA; q1 <<- NA; q3 <<- NA; IQR <<- NA
		RMSE.iqr <<- NA; RMSE.norm <<- NA
	})

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
		residual_variance = var(observed - mean_predicted, na.rm=TRUE),
		predictive_variance = var(mean_predicted, na.rm=TRUE),
		total_PL = var(observed - mean_predicted, na.rm=TRUE) + var(mean_predicted, na.rm=TRUE)
	)

	# Lower limit if RSQ 1:1 is 0
	out_df$RSQ.1 = ifelse(out_df$RSQ.1 < 0, 0, out_df$RSQ.1)

	# Upper limit of RMSE.normalized is 5
	out_df$RMSE.norm = ifelse(out_df$RMSE.norm > 5, 5, out_df$RMSE.norm)

	return(out_df)
}



#'  @title 			Tukey
#'  @description run tukey test
#'
#' @export
tukey <- function(x, y, extra_info = NULL, y.offset = .3){
	new.df <- cbind.data.frame("x" = x, "y" = y)
	abs_max <- max(new.df[,2])
	maxs <- new.df %>%
		group_by(x) %>%
		summarise(tot=max(y)+ y.offset * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>%
		dplyr::select(-y) %>%
		left_join(maxs, by="x")
	if (!is.null(extra_info)){
		Tukey_test <- cbind.data.frame(Tukey_test)
	}
	return(Tukey_test)
}



#'  @title 			get_sin_cos
#'  @description return sin and cosine components from month or day
#' @param input_dates a vector of dates in the format "20210130" or "202101"
#'
#' @return  a vector of values on the interval (0,1)
#' @examples
#' get_sin_cos("20210130")
#' get_sin_cos("202101")
#' @export
#'
get_sin_cos <- function(input_dates) {

	# if input is month, divide by 12
	if (class(input_dates[[1]]) == "character"
			& nchar(input_dates[[1]]) == 6) {
		mo <- lubridate::month(as.Date(paste0(input_dates, "01"), format="%Y%m%d"))
		y_sin = sin((2*pi*mo)/12)
		y_cos = cos((2*pi*mo)/12)

		# if input is day, divide by 365
	} else if (class(input_dates[[1]]) == "character" &
						 nchar(input_dates[[1]]) == 8) {
		doy <- lubridate::yday(as.Date(input_dates, format="%Y%m%d"))
		y_sin = sin((2*pi*doy)/365.25)
		y_cos = cos((2*pi*doy)/365.25)
	}  else {
		message("Inputs must be in the character format '201601' or date format '20160101'")
		return()
	}
	return(list(sin=y_sin, cos=y_cos))
}



#' @title crib_fun
#' stolen from colin's NEFI_microbe repo
#' converts a vector of [0,1] values to (0,1) a la Cribari-Neto & Zeileis 2010
#' @param x a vector of values on the interval [0,1]
#' @param N alternative sample size. This is useful when tranforming a matrix in the dirchlet case, rather than just a vector as in the beta case.
#'
#' @return  a vector of values on the interval (0,1)
#' @export
#'
crib_fun <- function(x,N = NA){
	#default use length of vector.
	if( is.na(N)){
		out <- (x * (length(x) - 1) + 0.5) / length(x)
	}
	#custom- useful when I am in multivariate case.
	if(!is.na(N)){
		out <- (x * (N - 1) + 0.5) / N
	}
	return(out)
}



#' @title getMaxMin
#' @description # get maximum amplitude from sin and cosine components
#'
#' @export
#'
getMaxMin <- function(sin, cos, T = 12, max_only = T) {
	sin <- sin[[1]]
	cos <- cos[[1]]

	if (sin == 0 & cos == 0) {
		return(list("min" = NA, "max" = NA))
	}

	# Closed-form peak of y(t) = sin*sin(2*pi*t/T) + cos*cos(2*pi*t/T) using atan2.
	# atan2 returns the angle in (-pi, pi] for any quadrant including the sin=0 axes,
	# avoiding the edge-case bug the old quadrant-by-quadrant branch had when sin==0.
	# The peak is at the angle whose direction is (cos, sin); the trough is half a period later.
	max <- (atan2(sin, cos) * T / (2 * pi)) %% T
	min <- (max + T / 2) %% T

	if (max_only) {
		return(max)
	} else {
		return(list("min" = min, "max" = max))
	}
}




#' @title invlogit
#' @description # invlogit
#'
#' @export
invlogit = function(x) exp(x)/(1+exp(x))


#' @title sin_cos_to_seasonality
#' @description Convert sin and cos effect sizes to seasonal amplitude parameter
#' @export
sin_cos_to_seasonality <- function(sin, cos){
	if (sin==0 & cos==0|is.na(sin)|is.na(cos)) {return(cbind.data.frame(max=NA,
																								amplitude_orig=NA,
																								amplitude = NA))}
	min_max <- getMaxMin(sin, cos, max_only = F)
	amplitude <- sqrt(sin^2 + cos^2)

	t=seq(0,12,0.1)
	monthly_vals = sin*sin(2*pi*t/12)+cos*cos(2*pi*t/12)
	max_val = max(monthly_vals)
	avg_val <- mean(min_max[[1]], min_max[[2]])
	out <- cbind.data.frame(max=min_max[[1]],
													amplitude_orig=amplitude,
													amplitude = max_val)
	return(out)
}


#' @title predictive_loss
#' @description Calculate predictive loss decomposed into predictive variance and residual variance
#' @export
predictive_loss = function(observed, predicted, predicted_sd){
	npred = length(predicted)
	predictive_variance = predicted_sd^2
	residual_variance = (predicted - observed)^2
	P = sum(predictive_variance, na.rm=T)/npred
	G = sum(residual_variance, na.rm=T)/npred
	total_PL = P+G
	data.frame(total_PL = total_PL, predictive_variance=P, residual_variance=G)
}


#' @title pivot_metrics
#' @description Pivot scoring metrics from wide to long format
#' @export
pivot_metrics = function(df) {
	result <- df %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS, CRPS_truncated, RSQ, RSQ.1,
															 RMSE.norm, residual_variance, predictive_variance, total_PL),
											names_to = "metric", values_to = "score")
	preserve_cols <- c("model_id", "fcast_type", "pretty_group", "model_name", "pretty_name", "rank_name", "taxon", "site_prediction", "mean_crps_sample")
	existing_preserve_cols <- preserve_cols[preserve_cols %in% colnames(df)]
	if(length(existing_preserve_cols) > 0) {
		result <- result %>% left_join(
			df %>% select(all_of(existing_preserve_cols)),
			by = existing_preserve_cols
		)
	}
	return(result)
}


#' @title calc_cv
#' @description Calculate coefficient of variation (as percentage)
#' @export
calc_cv <- function(x) sd(x, na.rm = T) / mean(x, na.rm = T) * 100


#' @title tag_facet
#' @description Add letter tags to ggplot facet panels
#' @export
tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf,
											hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
	gb <- ggplot_build(p)
	lay <- gb$layout$layout
	tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
	p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
								vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

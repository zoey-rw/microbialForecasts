
#' @title 	add_scoring_metrics
#' @description Score models and forecasts
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

	if(sum(is.na(observed )) > 0){stop('Error: NAs in observed vector.' )}
	if(sum(is.na(mean_predicted)) > 0){stop('Error: NAs in predicted vector.')}

	# These CRPS stats require distributions for each forecast
	# and cannot be calculated from median
	out_df1 <- cbind.data.frame(observed,mean_predicted,sd_predicted) %>%
		summarise(CRPS = mean(
			crps_norm(observed, mean_predicted, sd_predicted)),
			CRPS_truncated = mean(
				crps(observed,
						 family = "tnorm",
						 location = mean_predicted,
						 scale = sd_predicted,
						 lower = 0, upper = 1)))

	# The rest of these metrics can use forecast median as the best estimate

	if (use_median==TRUE) mean_predicted = median_predicted

	out_df2 = cbind.data.frame(observed,mean_predicted,sd_predicted) %>%
		summarise(
			RMSE = rmse(actual = observed, predicted = mean_predicted),
			RSQ.1 = 1 - (RMSE^2)/var(observed),
			RSQ.1.colin = rsq_1.1(observed, mean_predicted),
			predictive_loss(observed, mean_predicted, sd_predicted),
			BIAS = bias(actual = observed, predicted = mean_predicted),
			MAE = mae(actual = observed, predicted = mean_predicted),
			MAPE = mape(actual = observed, predicted = mean_predicted),
			RSQ = summary(lm(observed ~ mean_predicted))$r.squared,
			mean_abundance = mean(observed, na.rm=T),
			abundance = ifelse(mean_abundance < .005, .005, mean_abundance),
			q1 = quantile(observed, .25),
			q3 = quantile(observed, .75),
			IQR = q3-q1,
			RMSE.iqr = RMSE/IQR,
			RMSE.norm = RMSE/abundance)

	out_df = cbind.data.frame(out_df1, out_df2)

	# Lower limit if RSQ 1:1 is 0
	out_df$RSQ.1 = ifelse(out_df$RSQ.1 < 0, 0, out_df$RSQ.1)

	# Upper limit of RMSE.normalized is 5
	out_df$RMSE.norm.orig = out_df$RMSE.norm
	out_df$RMSE.norm = ifelse(out_df$RMSE.norm > 5, 5, out_df$RMSE.norm)

	out_df <- out_df %>%
		select(!!type)

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

	print(sin[[1]]);
	print(cos[[1]])

	print(length(sin));
	print(length(cos))

	sin <- sin[[1]]
	cos <- cos[[1]]

	t <- atan(sin/cos) * T/(2*pi)

	if (sin==0 & cos==0) {
		return(list("min" = NA, "max" = NA))
	}

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

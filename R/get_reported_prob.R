#' Returns probability of a collected sample having been reported by the analysis's time zero
#'
#' @param sampling_times numeric vector; sampling times
#' @param historic_sample_time numeric vector with historic times samples were collected (GISAID data base can be useful for obtaining such data for a given location)
#' @param historic_report_time numeric vector with historic times sequenced samples were reported (GISAID data base can be useful for obtaining such data for a given location)
#' @param lengthout numeric; length of time grid which defines granularity of reported probabilities, same length as using in BNPR estimation
#'
#' @importFrom dplyr right_join
#' @importFrom mgcv gam
#'
#' @return vector of probabilities of a sample having been reported by the analysis's time zero; given the sampling time
#' @export
get_reported_prob <- function(
    sampling_times, 
    historic_sample_time,
    historic_report_time,
    lengthout
  ){
  
  lengthout <- lengthout + 1 # To account for time zero
  
  historic_data <- data.frame(
    "sample_time" = historic_sample_time,
    "report_time" = historic_report_time,
    "reported" = historic_report_time >= 0 
  )
  
  logistic_fit <- mgcv::gam(
    formula = reported ~ s(sample_time), 
    data = historic_data, 
    family = "binomial"
  )
  
  time_grid <- seq(0, max(sampling_times), length = lengthout)
  midpoints <- time_grid[-lengthout] + diff(time_grid)/2
  
  reported_prob_int_df <- data.frame(
    "interval" = cut(
      x = midpoints, 
      breaks = time_grid, 
      right = FALSE,
      include.lowest = TRUE
    ),
    "reported_prob" = predict(
      logistic_fit, 
      newdata = list("sample_time" = midpoints),
      type = "response"
    )
  )
  
  sampling_time_int_df <- data.frame(
    "sampling_time" = sampling_times,
    "interval" = cut(
      x = sampling_times, 
      breaks = time_grid, 
      right = FALSE,
      include.lowest = TRUE
    )
  )
  
  joined_df <- dplyr::right_join(
    reported_prob_int_df, 
    sampling_time_int_df,
    by = "interval"
  )
  
  joined_df$reported_prob[joined_df$reported_prob == 0] <- 0.0001
  
  joined_df$reported_prob
}

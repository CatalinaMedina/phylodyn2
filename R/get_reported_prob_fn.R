#' Get reported probability function
#' 
#' Creates function that returns probability of a collected sample having been reported by the analysis's time zero
#'
#' @param historic_sample_time numeric vector with historic times samples were collected (sequencing data bases can be useful for obtaining such data for a given location)
#' @param historic_report_time numeric vector with historic times sequenced samples were reported (sequencing data bases can be useful for obtaining such data for a given location)
#' @param time_grid numeric vector; regular time grid used in BNPR estimation
#' @param return_log_rd_fn logical indicating if function should return log probabilities. This should be FALSE only if the reporting delay probability is being incorporated as a covariate and the log will be calculated in infer_coal_samp()
#' @param min_probability minimum probability of being reported, should never be exactly 0 since it will be logged,
#'
#' @importFrom mgcv gam
#' @importFrom stats stepfun
#'
#' @return function that returns vector of probabilities of a sample having been reported by the analysis's time zero; given the sampling time
#' @export
get_reported_prob_fn <- function(
    historic_sample_time,
    historic_report_time,
    time_grid,
    return_log_rd_fn = TRUE,
    min_probability = 0.0001
  ){
  
  if (min_probability == 0) {
    stop("min_probability cannot be exactly 0")
  }
  
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
  
  midpoints <- time_grid[-1] - diff(time_grid) / 2
  
  reported_prob_int_df <- data.frame(
    "interval" = cut(
      x = midpoints, 
      breaks = time_grid, 
      right = FALSE,
      include.lowest = TRUE
    ),
    midpoints = midpoints,
    "reported_prob" = predict(
      logistic_fit,
      newdata = list("sample_time" = midpoints),
      type = "response"
    )
  )
  
  reported_probs <- c(min_probability, reported_prob_int_df$reported_prob, 1)
  
  if (return_log_rd_fn) reported_probs <- log(reported_probs)
  
  prob_reported_fun <- stats::stepfun(
    x = time_grid, 
    y = reported_probs, 
    right = TRUE
  )
  
  prob_reported_fun
}

#' Get reported probability function
#' 
#' Creates function that returns probability of a collected sample having been reported by the analysis's time zero
#'
#' @param historic_reporting_delays numeric vector with historic reporting delays. This must be on in the same units as time grid with all nonnegative values. (sequencing data bases can be useful for obtaining such data for a given location)
#' @param time_grid numeric vector; regular time grid used in BNPR estimation
#' @param return_log_rd_fn logical indicating if function should return log probabilities. This should be FALSE only if the reporting delay probability is being incorporated as a covariate and the log will be calculated in infer_coal_samp()
#' @param min_probability minimum probability of being reported, should never be exactly 0 since it will be logged.
#'
#' @importFrom mgcv gam
#' @importFrom stats stepfun
#'
#' @return function that returns vector of probabilities of a sample having been reported by the analysis's time zero; given the sampling time
#' @export
get_reported_prob_fn <- function(
    historic_reporting_delays,
    time_grid,
    return_log_rd_fn = TRUE,
    min_probability = 0.00001
){
  
  if (min_probability == 0) {
    stop("min_probability cannot be exactly 0")
  }
  
  if (min(time_grid) != 0) {
    warning("Time grid must begin with 0!")
  }
  
  if (max(time_grid) < max(historic_reporting_delays)) { # Only happens in special case of simulations
    grid <- c(time_grid, max(historic_reporting_delays))
  } else {
    grid <- time_grid
  }
  
  sorted_delays <- sort(historic_reporting_delays)
  
  intervals_matched <- data.frame(
    "reporting_delay" = sorted_delays,
    "time_interval" = cut(
      x = sorted_delays,
      breaks = grid,
      right = FALSE,
      include.lowest = TRUE
    )
  )
  
  # From end interval repeat
  intervals_matched <- intervals_matched[!is.na(intervals_matched$time_interval), ]
 
  num_samples <- length(sorted_delays)
  
  int_tabulated <- data.frame(
    "interval" = names(table(intervals_matched$time_interval)),
    "num_reported_in_int" = as.numeric(table(intervals_matched$time_interval))
  )
  
  int_tabulated$per_reported_in_int <- int_tabulated$num_reported_in_int / length(sorted_delays)
  int_tabulated$prob_reported <- cumsum(int_tabulated$per_reported_in_int)
  int_tabulated$prob_reported <- replace(
    int_tabulated$prob_reported, 
    int_tabulated$prob_reported == 0,
    values = min_probability
  )

  reported_probs <- c(min_probability, int_tabulated$prob_reported, 1)
  
  if (return_log_rd_fn) reported_probs <- log(reported_probs)
  
  prob_reported_fun <- stats::stepfun(
    x = grid,
    y = reported_probs,
    right = FALSE
  )
  
  prob_reported_fun
}

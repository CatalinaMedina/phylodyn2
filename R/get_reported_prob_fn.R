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
#' @importFrom tidyr replace_na
#' @import dplyr
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
    time_grid <- c(time_grid, max(historic_reporting_delays))
  }
  
  sorted_delays <- data.frame(
    "reporting_delay" = historic_reporting_delays,
    "time_interval" = cut(
      x = historic_reporting_delays,
      breaks = time_grid,
      right = FALSE,
      include.lowest = TRUE
    )
  ) %>% 
    arrange(time_interval) %>% 
    filter(!is.na(time_interval))
  
  num_samples <- length(historic_reporting_delays)
  
  reported_prob_cumsum <- sorted_delays %>% 
    group_by(time_interval) %>% 
    summarize(n = n(), percent_reported = n()/num_samples) %>% 
    full_join( # Account for intervals with no observations
      data.frame(
        "time_interval" = cut(
          x = time_grid[-length(time_grid)],
          breaks = time_grid,
          right = FALSE,
          include.lowest = TRUE
        )
      ), 
      by = "time_interval"
    ) %>% 
    arrange(time_interval) %>% 
    mutate(percent_reported = tidyr::replace_na(percent_reported, replace = 0)) %>%
    mutate(n = tidyr::replace_na(n, replace = 0)) %>%
    mutate(prob_reported = cumsum(percent_reported)) %>% 
    mutate(prob_reported = replace(prob_reported, prob_reported == 0, values = min_probability))
  
  reported_probs <- c(min_probability, reported_prob_cumsum$prob_reported, 1)
  
  if (return_log_rd_fn) reported_probs <- log(reported_probs)
  
  prob_reported_fun <- stats::stepfun(
    x = time_grid,
    y = reported_probs,
    right = FALSE
  )
  
  prob_reported_fun
}

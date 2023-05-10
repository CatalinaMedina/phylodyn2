#' samp_stats Collect sampling data for estimation accounting for preferential sampling, 
#'
#' @param grid numeric vector; regular grid across genealogy
#' @param samp_times numeric vector; sampling times 
#' @param n_sampled numeric vector; number of sampling events at each sampling time, length matching the length of samp_times
#' @param trim_end logical; TODO
#' 
#' @import stats
#' @importFrom utils head
#'
#' @return data frame with \describe{
#'   \item{time}{Midpoints of time grid intervals}
#'   \item{count}{Number of sampling times within each grid interval}
#'   \item{E}{Vector of differences between consecutive grid points}
#'   \item{E_log}{Log of E}
#' }
#'
samp_stats <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE){
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid) / 2
  E <- diff(grid)
  
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)
  
  if (is.null(n_sampled)) {
    count <- as.vector(table(bins))
    
  } else {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
    
  }
  
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end) {
    result <- result[stats::complete.cases(result), ]
  }
  
  result
  
}

#' joint_stats Prepares data for INLA for case of preferential sampling model
#'
#' @param coal_data data frame with columns: time, event, E_log
#' @param samp_data data frame with columns: time, count, E_log
#'
#' @return list with \describe{
#'   \item{Y}{matrix with two columns. The first column contains the coalescence events padded with NA's at end. The second column contains the sampling counts, padded by NA's at the beginning.}
#'   \item{beta0}{numeric vector with 0's for length of coal_data$time and 1's for length of samp_data$time; for sampling intensity function}
#'   \item{time}{coalescence times padded with 0's at end}
#'   \item{time2}{sampling times padded with 0's at beginning}
#'   \item{w}{numeric vector; 1's repeated for length of coal_data$time followed by -1's repeated for length of samp_data$time}
#'   \item{E_log}{coalescence E_log followed by sampling E_log}
#' }
#'
joint_stats <- function(coal_data, samp_data){
  n1 <- length(coal_data$time)
  n2 <- length(samp_data$time)
  
  beta0 <- c(rep(0, n1), rep(1, n2))
  
  E_log <- c(coal_data$E_log, samp_data$E_log)
  
  Y <- matrix(
    c(coal_data$event, rep(NA, n2), rep(NA, n1), samp_data$count),
    nrow = n1 + n2, 
    byrow = FALSE
  )
  
  w <- c(rep(1, n1), rep(-1, n2))
  
  time  <- c(coal_data$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), samp_data$time)
  
  list(
    Y = Y, 
    beta0 = beta0, 
    time = time, 
    time2 = time2, 
    w = w, 
    E_log = E_log
  )
}

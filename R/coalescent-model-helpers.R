#' coal_stats Returns data frame on midpoints of regular grid with record of coalescent events and Poisson regression offset
#'
#' @param grid numeric vector; regular grid over entire genealogy
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; number of sampling events at each sampling time, length matching the length of samp_times
#' @param log_zero numeric; value to store for E of log(0)
#'
#' @importFrom stats stepfun
#'
#' @return a data frame with \describe{
#'   \item{time}{Vector of midpoints of regular grid where each value represents that either a grid end point, sampling time, or coalescent time occured in that regular grid interval}
#'   \item{event}{Logical vector recording coalescent events with a 1 and regurlar grid points or sampling times with a 0}
#'   \item{E}{Coalescent factor times width of interval (on combined grid level)}
#'   \item{E_log}{Possion regression offset}
#' }
#'
coal_stats <- function(
    grid, 
    samp_times, 
    coal_times, 
    n_sampled = NULL,
    log_zero = -100
){
  
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid) / 2
  
  if (is.null(n_sampled)) {
    n_sampled <- rep(1, length(samp_times))
  }
  
  args <- gen_INLA_args(
    samp_times = samp_times, 
    n_sampled = n_sampled,
    coal_times = coal_times
  )
  
  coal_factor <- args$coal_factor
  s <- args$s # combined and sorted sampling and colescent times
  event <- args$event # 0 = sampling time, 1 = coalescent time
  
  grid_trimmed <- setdiff(x = grid, y = s) # keep elements in x not in y
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE) # identify interval
  time <- field[time_index] # repeat grid midpoints for # of times sampling, coalescent, and or grid point falls in interval
  
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering] # 0 = grid or sampling time, 1 = coalescent
  
  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid) * Cvec
  
  E_log = log(E)
  E_log[E == 0] = log_zero
  
  data.frame(
    time = time, 
    event = event_out[-1], 
    E = E, 
    E_log = E_log
  )
  
}

#' gen_INLA_args Returns data frame on combined grid with record of coalescent events and coalescent factor
#'
#' @param samp_times numeric vector; sampling times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#' @param coal_times numeric vector; coalescence times
#'
#' @importFrom utils head
#'
#' @return data frame with \describe{
#'   \item{coal_factor}{Number of active lineages choose 2}
#'   \item{s}{Sampling and coalescence times combined and sorted} 
#'   \item{event}{Logical vector with 0 = sampling time and 1 = coalescent time}
#'   \item{lineages}{Numeric vector of number of active lineages}
#' }
#'
gen_INLA_args <- function(samp_times, n_sampled, coal_times){
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")
  
  if (length(intersect(coal_times, samp_times)) > 0)
    warning(
      "Coincident sampling event and coalescent event: results may be unpredictable."
    )
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)
  
  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix] # change in number of active lineages (coalescent event results in -1)
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages * (lineages - 1) / 2 # Number of active lineages choose 2
  
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  
  list(
    coal_factor = coal_factor, 
    s = sorting$x, 
    event = event, 
    lineages = lineages
  )
  
}

#' condense_stats Sums all Poisson regression offsets (product of coalescent and combined interval width) within each regular grid interval
#'
#' @param time numeric vector; time of each event
#' @param event numeric vector; number of events at a given time
#' @param E numeric vector; products of coalescent factor and width of interval (on combined grid level)
#' @param log_zero numeric; value to store for E of log(0)
#' 
#' @importFrom stats aggregate
#'
#' @return data frame with \describe{
#'   \item{time}{Vector of midpoints of regular grid where each value represents that either a grid end point, sampling time, or coalescent time occured in that regular grid interval}
#'   \item{event}{Number of coalescent events within each regular grid interval}
#'   \item{E}{Vector of the sum of products of coalescent factor and width of interval (on combined grid level), for each regular grid interval}
#'   \item{E_log}{Possion regression offset}
#' }
#'
condense_stats <- function(time, event, E, log_zero = -100){
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E
  
  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log
  
  result
}

#' Simulate from inhomogeneous, heterochronous coalescent
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param method which sampling method to use. "tt" invoke time-transformation
#'   method, "thin" invokes thinning method.
#' @param val_upper numeric used by time-transformation method to set a starting
#'   point for its dynamic numerical integration upper bound.
#' @param lower_bound numeric lower limit of \code{traj} function on its
#'   support.  Used only by thinning method.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim(0:2, n_sampled = 3:1, traj = unif_traj, lower_bound = 10)
coalsim <- function(
    samp_times, 
    n_sampled, 
    traj, 
    method = "tt", 
    val_upper = 10, 
    lower_bound = 1, 
    ...
  ){
  
  if (method == "tt") {
    result = coalsim_tt(samp_times, n_sampled, traj, val_upper, ...)

  } else if (method == "thin") {
    result = coalsim_thin(samp_times, n_sampled, traj, lower_bound, ...)
    
  } else {
    stop("Argument method not recognized.")
    
  }
  
  result
}


#' Coalescent time transfermation method
#'
#' @param samp_times numeric vector of sampling times
#' @param n_sampled numeric vector of number of samples collected at each sampling time
#' @param traj function that returns effective population size at time t.
#' @param val_upper numeric used by time-transformation method to set a starting
#'   point for its dynamic numerical integration upper bound. 
#' @param ... additional arguments to be passed to \code{traj} function 
#' 
#' @import stats
#'
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#'
coalsim_tt <- function(samp_times, n_sampled, traj, val_upper = 10, ...){
  if (stats::is.stepfun(traj)) {
    knots <- knots(traj)
    midpts <- c(min(knots) - 1, knots[-1] - diff(knots)/2, max(knots) + 1)
    traj_inv <- stats::stepfun(x = knots, y = 1/traj(midpts))
    
    hazard <- function(t, lins, start, target){ 
      int <- integrate_step_fun(traj_inv, start, start + t)
      .5 * lins * (lins - 1) * int - target
    }
    
    is_stepfun <- TRUE
    
  } else {
    traj_inv <- function(t) 1/traj(t, ...)
    hazard <- function(t, lins, start, target) { 
      int <- stats::integrate(
        traj_inv, 
        start, 
        start + t, 
        stop.on.error = FALSE
      )
      
      .5 * lins * (lins - 1) * int$value - target
      
    }
    
    is_stepfun <- FALSE
    
  }
  
  coal_times <- NULL
  lineages <- NULL
  
  curr <- 1
  active_lineages <- n_sampled[curr]
  time <- samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1) {
    if (active_lineages == 1) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    
    target <- stats::rexp(1)
    
    if (is_stepfun) {
      y <- hazard_uniroot_stepfun(
        traj_inv_stepfun = traj_inv,
        lineages = active_lineages,
        start = time, target = target
      )
      
    } else {
      y <- stats::uniroot(
        hazard, 
        lins=active_lineages, 
        start=time, target=target,
        lower=0, upper=val_upper, extendInt = "upX"
      )$root
      
    }
    
    while(curr < length(samp_times) && time + y >= samp_times[curr+1]) {
      target <- -hazard(
        t = samp_times[curr+1] - time, 
        lins = active_lineages,
        start = time, target = target
      )
      
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      
      if (is_stepfun) {
        y <- hazard_uniroot_stepfun(
          traj_inv_stepfun = traj_inv,
          lineages = active_lineages,
          start = time, target = target
        )
        
      } else {
        y <- stats::uniroot(
          hazard, 
          lins = active_lineages, 
          start = time, target = target,
          lower = 0, upper = val_upper, extendInt = "upX"
        )$root
        
      }
      
    }
    
    time <- time + y
    coal_times <- c(coal_times, time)
    lineages <- c(lineages, active_lineages)
    active_lineages <- active_lineages - 1
    
  }
  
  list(
    coal_times = coal_times, 
    lineages = lineages,
    intercoal_times = c(coal_times[1], diff(coal_times)),
    samp_times = samp_times, 
    n_sampled = n_sampled
  )
}


#' Coalescent thinning method
#'
#' @param samp_times numeric vector of sampling times
#' @param n_sampled numeric vector of number of samples collected at each sampling time
#' @param traj function that returns effective population size at time t.
#' @param lower_bound numeric lower limit of \code{traj} function on its
#'   support.  Used only by thinning method.
#' @param ... additional arguments to be passed to \code{traj} function 
#' 
#' @import stats
#' 
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#'
coalsim_thin <- function(samp_times, n_sampled, traj, lower_bound, ...){
  coal_times <- NULL
  lineages <- NULL
  
  curr <- 1
  active_lineages <- n_sampled[curr]
  time <- samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1) {
    if (active_lineages == 1) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    
    time <- time + stats::rexp(
      1, 
      0.5 * active_lineages * (active_lineages - 1) / lower_bound
    )
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1]) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      
    } else if (stats::runif(1) <= lower_bound / traj(time, ...)) {
      coal_times <- c(coal_times, time)
      lineages <- c(lineages, active_lineages)
      active_lineages <- active_lineages - 1
      
    }
    
  }
  
  list(
    coal_times = coal_times, 
    lineages = lineages,
    intercoal_times = c(coal_times[1], diff(coal_times)),
    samp_times = samp_times, 
    n_sampled = n_sampled
  )
}



#' Hazard uniroot step function
#'
#' @param traj_inv_stepfun inverse of the function that returns effective population size at time t
#' @param lineages number of active lineages
#' @param start time to start at
#' @param target TODO
#' 
#' @importFrom stats knots
#'
#' @return result vector minus start
#'
hazard_uniroot_stepfun <- function(traj_inv_stepfun, lineages, start, target){
  knots <- stats::knots(traj_inv_stepfun)
  lin_factor <- 0.5 * lineages * (lineages - 1)
  t <- start
  
  while (target > 0 && sum(knots > t) > 0) {
    next_knot <- min(knots[knots > t])
    
    taj_inv_val <- traj_inv_stepfun(mean(c(t, next_knot)))
    
    if ((next_knot - t) * lin_factor * taj_inv_val > target) {
      result <- t + target / (lin_factor * taj_inv_val)
      target <- 0
      
    } else {
      target <- target - (next_knot - t) * lin_factor * taj_inv_val
      t <- next_knot
      
    }
    
  }
  
  if (sum(knots > t) < 1) {
    result <- t + target / (lin_factor * traj_inv_stepfun(t + 1))
  }
  
  result - start
}


#' Integrates step function
#'
#' @param fun function
#' @param from start
#' @param to end
#' 
#' @importFrom stats knots
#'
#' @return numeric; integration value
#'
integrate_step_fun <- function(fun, from, to) {
  breaks <- stats::knots(fun)

  mask <- breaks > from & breaks < to
  pts <- c(from, breaks[mask], to)
  diffs <- diff(pts)
  midpts <- pts[-1] - diffs / 2
  
  segments <- diffs * fun(midpts)
  
  sum(segments)
}
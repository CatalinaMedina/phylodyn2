#' Sample preferentially
#' 
#' Sample values from an inhomogeneous Poisson process with intensity
#'   c*f(t)^beta.
#' 
#' @param f function to be sampled from preferentially.
#' @param c numeric constant of proportionality for intensity. Defaults to 1.
#' @param beta numeric power to which f is raised for intensity. Defaults to 1.
#' @param upper numeric upper bound for function between limits defined by lim.
#'   Used in thinning algorithm.
#' @param lim lower and upper sampling bounds. Defaults to c(0,1).
#' @param grid.len integer if upper is undefined, automatically evaluates f at
#'   grid.len points between min(lim) and max(lim) to find upper bound.
#'
#' @importFrom stats rpois
#' @importFrom stats runif
#' 
#' @return numeric vector of sampling times
#'
#' @export
pref_sample <- function(
    f,  
    c = 1, beta = 1, upper = NULL, 
    lim = c(0, 1), grid.len = 1000
  ){
  if (is.null(upper)) {
    grid <- seq(lim[1], lim[2], length.out = grid.len)
    eval <- f(grid)^beta
    upper <- max(eval)
  }
  
  nfull <- stats::rpois(1, c * upper * diff(lim))
  pfull <- sort(stats::runif(nfull, lim[1], lim[2]))
  
  pthin <- NULL
  
  for (p in pfull) {
    if (stats::runif(1) < f(p)^beta / upper){
      pthin <- c(pthin, p)
    }
  }

  pthin
}

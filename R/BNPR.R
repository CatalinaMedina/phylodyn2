#' Bayesian nonparametric phylodynamic reconstruction.
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used?
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\kappa}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
#' @param fns list containing functions of covariates.
#' @param log_fns logical whether or not to to apply a log-transformation to
#'   the output of the functions in \code{fns}.
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the 
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'   
#' @return Phylodynamic reconstruction of effective population size at grid 
#'   points. \code{result} contains the INLA output, \code{data} contains the 
#'   information passed to INLA, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop} contains a vector 
#'   of the posterior median effective population size estimates, 
#'   \code{effpop025} and \code{effpop975} contain the 2.5th and 97.5th 
#'   posterior percentiles, \code{summary} contains a data.frame of the 
#'   estimates, and \code{derivative} (if \code{derivative = TRUE}) contains a
#'   data.frame summarizing the log-derivative.
#' @export
#' 
#' @examples
#' \dontrun{
#' data("NY_flu")
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'  res = BNPR(NY_flu)
#'  plot_BNPR(res)
#' }
#' }
BNPR <- function(
    data, lengthout = 100, pref=FALSE, 
    prec_alpha=0.01, prec_beta=0.01, beta1_prec = 0.001, 
    rd_prob_fn = NULL, fns = NULL, log_fns = TRUE,
    simplify = TRUE, derivative = FALSE, forward = TRUE, link = 1
  ){
  
  if (class(data) == "phylo") {
    phy <- summarize_phylo(data)
    
  } else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data))) {
    phy <- with(
      data, 
      list(
        samp_times = samp_times, 
        coal_times = coal_times,
        n_sampled = n_sampled)
    )
    
  }
  
  result <- infer_coal_samp(
    samp_times = phy$samp_times, coal_times = phy$coal_times,
    n_sampled = phy$n_sampled, 
    rd_prob_fn = rd_prob_fn, fns = fns, lengthout = lengthout,
    prec_alpha = prec_alpha, prec_beta = prec_beta,
    beta1_prec = beta1_prec, use_samp = pref, log_fns = log_fns,
    simplify = simplify, derivative = derivative, link = link
  )
  
  result$samp_times <- phy$samp_times
  result$n_sampled  <- phy$n_sampled
  result$coal_times <- phy$coal_times
  
  result$effpop     <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975  <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025  <- exp(-result$result$summary.random$time$`0.975quant`)
  
  result$summary <- with(
    result$result$summary.random$time,
    data.frame(
      time = ID, 
      mean = exp(-mean),
      sd = sd * exp(-mean),
      quant0.025 = exp(-`0.975quant`),
      quant0.5 = exp(-`0.5quant`),
      quant0.975 = exp(-`0.025quant`)
    )
  )
  
  if (derivative) {
    if (forward) {
      ind <- c(1:(lengthout - 1), (lengthout - 1))
    } else {
      ind <- c(1, 1:(lengthout - 1))
    }
    
    result$derivative <- with(
      result$result$summary.lincomb,
      data.frame(
        time = result$x, 
        mean = -mean[ind], 
        sd = sd[ind],
        quant0.025 = -`0.975quant`[ind],
        quant0.5   = -`0.5quant`[ind],
        quant0.975 = -`0.025quant`[ind]
      )
    )
    
  }
  
  if (pref) {
    result$beta0     <- result$result$summary.fixed["beta0","0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0",]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1     <- result$result$summary.hyperpar[2,"0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2,]
    rownames(result$beta1summ) <- "Beta1"
    
  }
  
  result
  
}


#' @describeIn BNPR Uses preferential sampling model.
#' @export
BNPR_PS <- function(
    data, lengthout = 100, 
    prec_alpha=0.01, prec_beta=0.01, beta1_prec = 0.001, 
    rd_prob_fn = NULL, fns = NULL, log_fns = TRUE,
    simplify = TRUE, derivative = FALSE, forward = TRUE, link = 1
  ){
  
  BNPR(
    data = data, lengthout = lengthout, pref = TRUE,
    prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec, 
    rd_prob_fn = rd_prob_fn, fns = fns, log_fns = log_fns,
    simplify = simplify, derivative = derivative, forward = forward, link = link
  )
  
}

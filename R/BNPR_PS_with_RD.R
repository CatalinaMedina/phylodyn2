#' Bayesian nonparametric phylodynamic reconstruction with reporting delay in preferential sampling model
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param historic_reporting_delays numeric vector with historic reporting delays. This must be on in the same units as time grid with all nonnegative values. (sequencing data bases can be useful for obtaining such data for a given location)
#' @param rd_as_offset logical whether reporting delay adjustment should be implemented as offset (TRUE), or as a covariate function with a steep prior (FALSE)
#' @param prec_alpha numeric; hyperparameter alpha for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param prec_beta numeric; hyperparameter beta for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param beta1_mean numeric; mean of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param beta1_prec numeric; precision of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param time0_offset_from_sim_rd numeric; For BNPR_PS_with_RD only. The time between the first reported sampling time and the true first sampling time, time zero, in the case of simulating reporting delays. Generally should be left as NULL. Time zero is considered to be the first sampling time, but when simulating reporting delays and dropping tips from the tree, time zero shifts to the first reported sampling time. 
#' @param fns function; list of covariate functions for the sampling intensity
#' @param log_fns logical; specifies if the log of the covariate functions, fns, needs to be taken. FALSE indicates that the covariate function already returns log transformed values
#' @param fns_coeff_prior_mean numeric vector; normal prior mean for fns coefficient(s). If non NULL, must match length of fns
#' @param fns_coeff_prior_prec numeric vector; normal prior precision for fns coefficient(s). If non NULL, must match length of fns
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#' @param link link for INLA "regression"
#'   
#' @return Phylodynamic reconstruction of effective population size at grid points:\describe{ 
#'   \item{result}{contains the INLA output}
#'   \item{data}{contains the information passed to INLA}
#'   \item{grid}{contains the grid end points}
#'   \item{x}{contains the grid point centers} 
#'   \item{effpop}{contains a vector of the posterior median effective population size estimates}
#'   \item{effpop025}{The 2.5th posterior percentiles}
#'   \item{effpop975}{The 97.5th posterior percentiles} 
#'   \item{summary}{contains a data.frame of the estimates}
#'   \item{derivative}{(if \code{derivative = TRUE}) contains a data frame summarizing the log-derivative}
#'   \item{rd_prob_fn}{Function which computes probability of a sample having been reported based on it's sampling time. Will be log probabilities unless reporting delay function was specified as covariate and log_fns is FALSE.}
#'  }
#'   
#' @export
#' 
BNPR_PS_with_RD <- function(
    data, 
    historic_reporting_delays,
    rd_as_offset = TRUE, time0_offset_from_sim_rd = NULL,
    lengthout = 100, 
    prec_alpha = 0.01, prec_beta = 0.01, beta1_mean = 0, beta1_prec = 0.001, 
    fns = NULL, log_fns = TRUE, 
    fns_coeff_prior_mean = NULL, fns_coeff_prior_prec = NULL,
    simplify = TRUE, derivative = FALSE, 
    forward = TRUE, link = 1
){
  
  # Get grid for rd_fn
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
  
  grid <- seq(min(phy$samp_times), max(phy$coal_times), length.out = lengthout + 1)
  
  rd_fn_will_be_logged_later <- !rd_as_offset & log_fns
  
  rd_fn <- get_reported_prob_fn(
    historic_reporting_delays = historic_reporting_delays,
    time_grid = grid,
    return_log_rd_fn = !rd_fn_will_be_logged_later
  )
  
  if (is.null(time0_offset_from_sim_rd)) {
    rd_prob_fn <- rd_fn
  } else if (is.numeric(time0_offset_from_sim_rd) & time0_offset_from_sim_rd > 0) {
    rd_prob_fn <- function (x) rd_fn(x + time0_offset_from_sim_rd)
  } else {
    stop("time0_offset_from_sim_rd should be NULL or a positive number.")
  }
  
  if (rd_as_offset) {
    res <- BNPR_PS(
      data, lengthout = lengthout, 
      prec_alpha = prec_alpha, prec_beta = prec_beta, 
      beta1_mean = beta1_mean, beta1_prec = beta1_prec, 
      rd_prob_fn = rd_prob_fn, 
      fns = fns, log_fns = log_fns, 
      fns_coeff_prior_mean = fns_coeff_prior_mean, 
      fns_coeff_prior_prec = fns_coeff_prior_prec,
      simplify = simplify, derivative = derivative, 
      forward = forward, link = link
    )
    
  } else {
    rd_fn_mean <- 1
    rd_fn_prec <- 1000
    
    fns_coeff_prior_mean <- c(list(rd_prob_fn = rd_fn_mean), as.list(fns_coeff_prior_mean))
    fns_coeff_prior_prec <- c(list(rd_prob_fn = rd_fn_mean), as.list(fns_coeff_prior_prec))
    
    fns <- c(list(rd_prob_fn), fns)
    
    # Functions need to be named so priors can be specified later
    if (length(fns) > 1 & length(fns_coeff_prior_mean) > 1) {
      names(fns_coeff_prior_mean) <- c(
        "rd_fn",
        paste0("fn", 1:(length(fns_coeff_prior_mean) - 1))
      )
      
      if(length(fns) != length(fns_coeff_prior_mean)){
        stop("length(fns_coeff_prior_mean) must equal length(fns) if nonNULL")
      }
    }
    
    if (length(fns) > 1 & length(fns_coeff_prior_prec) > 1) {
      names(fns_coeff_prior_prec) <- c(
        "rd_fn",
        paste0("fn", 1:(length(fns_coeff_prior_prec) - 1))
      )
      
      if(length(fns) != length(fns_coeff_prior_prec) - 1){
        stop("length(fns_coeff_prior_prec) must equal length(fns) if nonNULL")
      }
    }
    
    fns_coeff_prior_mean$default <- 0
    fns_coeff_prior_prec$default <- 0.001
    
    res <- BNPR_PS(
      data, lengthout = lengthout, 
      prec_alpha = prec_alpha, prec_beta = prec_beta, 
      beta1_mean = beta1_mean, beta1_prec = beta1_prec,
      rd_prob_fn = NULL,
      fns = fns, log_fns = log_fns, 
      fns_coeff_prior_mean = fns_coeff_prior_mean, 
      fns_coeff_prior_prec = fns_coeff_prior_prec,
      simplify = simplify, derivative = derivative, 
      forward = forward, link = link
    )
    
  }
  
  res$rd_prob_fn <- rd_prob_fn
  
  if (is.null(time0_offset_from_sim_rd)) {
    return(res)
    
  } else {
    adj_BNPR_output <- res
    adj_BNPR_output$grid <- res$grid + time0_offset_from_sim_rd
    adj_BNPR_output$x <- res$x + time0_offset_from_sim_rd
    adj_BNPR_output$samp_times <- res$samp_times + time0_offset_from_sim_rd
    adj_BNPR_output$coal_times <- res$coal_times + time0_offset_from_sim_rd
    return(adj_BNPR_output)
    
  }
}

#' Bayesian nonparametric phylodynamic reconstruction.
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used? If so use BNPR_PS for a better report of estimation
#' @param historic_sample_time numeric vector with historic times samples were collected (GISAID data base can be useful for obtaining such data for a given location)
#' @param historic_report_time numeric vector with historic times sequenced samples were reported (GISAID data base can be useful for obtaining such data for a given location)
#' @param rd_as_offset logical whether reporting delay adjustment should be implememented as offset (TRUE), or as a covariate function with a steep prior (FALSE)
#' @param prec_alpha numeric; hyperparameter alpha for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param prec_beta numeric; hyperparameter beta for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param beta1_mean numeric; mean of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param beta1_prec numeric; precision of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param rd_prob_fn function; a function that takes in a vector of sampling times and returns the probability of a collected sample having been reported
#' @param fns function; list of covariate functions for the sampling intensity
#' @param log_fns Boolean; specifies if the log of the covariate functions, fns, needs to be take. FALSE indicates that the covriate function already returns log transformed values
#' @param fns_coeff_prior_mean value to specify mean for normal prior for fn coefficient
#' @param fns_coeff_prior_prec value to specify precision for normal prior for fn coefficient
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
#'  }
#'   
#' @export
#' 
BNPR <- function(
    data, lengthout = 100, pref = FALSE, 
    prec_alpha = 0.01, prec_beta = 0.01, beta1_mean = 0, beta1_prec = 0.001, 
    rd_prob_fn = NULL, 
    fns = NULL, log_fns = TRUE,
    fns_coeff_prior_mean = NULL, fns_coeff_prior_prec = NULL,
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
    rd_prob_fn = rd_prob_fn, 
    fns = fns, log_fns = log_fns,
    fns_coeff_prior_mean = fns_coeff_prior_mean, fns_coeff_prior_prec = fns_coeff_prior_prec,
    lengthout = lengthout,
    prec_alpha = prec_alpha, prec_beta = prec_beta,
    beta1_mean = beta1_mean, beta1_prec = beta1_prec, 
    use_samp = pref, simplify = simplify, derivative = derivative, link = link
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
    prec_alpha = 0.01, prec_beta = 0.01, beta1_mean = 0, beta1_prec = 0.001, 
    rd_prob_fn = NULL, 
    fns = NULL, log_fns = TRUE, 
    fns_coeff_prior_mean = 0, fns_coeff_prior_prec = 0.01,
    simplify = TRUE, derivative = FALSE, forward = TRUE, link = 1
  ){
  
  BNPR(
    data = data, lengthout = lengthout, pref = TRUE,
    prec_alpha = prec_alpha, prec_beta = prec_beta, 
    beta1_mean = beta1_mean, beta1_prec = beta1_prec, 
    rd_prob_fn = rd_prob_fn, 
    fns = fns, log_fns = log_fns,
    fns_coeff_prior_mean = fns_coeff_prior_mean, fns_coeff_prior_prec = fns_coeff_prior_prec,
    simplify = simplify, derivative = derivative, forward = forward, link = link
  )
  
}

#' @describeIn BNPR Uses preferential sampling model with adjustment for reporting delay
#' @export
BNPR_PS_with_RD <- function(
    data, 
    historic_sample_time, historic_report_time,
    rd_as_offset = TRUE, lengthout = 100,
    prec_alpha = 0.01, prec_beta = 0.01, beta1_mean = 0, beta1_prec = 0.001, 
    fns = NULL, #Don't have the ability to handle rd as func with other covariates
    log_fns = TRUE, 
    fns_coeff_prior_mean = 0, fns_coeff_prior_prec = 0.01,
    simplify = TRUE, derivative = FALSE, forward = TRUE, link = 1  
  ){
  
  get_reported_prob_fn <- function(sampling_times) {
    get_reported_prob(
      sampling_times, 
      historic_sample_time = historic_sample_time,
      historic_report_time = historic_report_time,
      lengthout = lengthout
    )
  }
  
  get_log_reported_prob_fn <- function(sampling_times) {
    log(get_reported_prob_fn(sampling_times))
  }
  
  if (rd_as_offset) {
    res <- BNPR_PS(
      data = data, lengthout = lengthout, 
      prec_alpha = prec_alpha, prec_beta = prec_beta, 
      beta1_mean = beta1_mean, beta1_prec = beta1_prec, 
      rd_prob_fn = get_reported_prob_fn, 
      fns = fns, log_fns = log_fns, 
      fns_coeff_prior_mean = fns_coeff_prior_mean, 
      fns_coeff_prior_prec = fns_coeff_prior_prec,
      simplify = simplify, derivative = derivative, 
      forward = forward, link = link
    )
    
  } else {
    if (!is.null(fns)) {
      stop("We currently cannot handle reporting probability funciton in addition to other sampling intensity formula covariates.")
    }
    
    if (log_fns) {
      fns_list <- list(get_reported_prob_fn)
    } else {
      fns_list <- list(get_log_reported_prob_fn)
    }
    
    res <- BNPR_PS(
      data = data, lengthout = lengthout, 
      prec_alpha = prec_alpha, prec_beta = prec_beta, 
      beta1_mean = beta1_mean, beta1_prec = beta1_prec,
      rd_prob_fn = NULL,
      fns = fns_list,
      log_fns = log_fns, 
      fns_coeff_prior_mean = 1, 
      fns_coeff_prior_prec = 1000,
      simplify = simplify, derivative = derivative, 
      forward = forward, link = link
    )
    
  }
  
  res
}

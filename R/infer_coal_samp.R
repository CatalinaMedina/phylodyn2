#' Infer coalescence sample
#'
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#' @param lengthout numeric; length of time for the estimation
#' @param rd_prob_fn function; a function that takes in a vector of sampling times and returns the probability of a collected sample having been reported
#' @param fns function; list of covariate functions for the sampling intensity
#' @param log_fns logical; specifies if the log of the covariate functions, fns, needs to be take. FALSE indicates that the covriate function already returns log transformed values
#' @param fns_coeff_prior_mean numeric vector; normal prior mean for fns coefficient(s). If non NULL, must match length of fns
#' @param fns_coeff_prior_prec numeric vector; normal prior precision for fns coefficient(s). If non NULL, must match length of fns
#' @param prec_alpha numeric; hyperparameter alpha for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param prec_beta numeric; hyperparameter beta for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param beta1_mean numeric; mean of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param beta1_prec numeric; precision of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param use_samp logical; Indicates if the preferential sampling model should be used
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the log-derivative.
#' @param events_only logical; TODO
#' @param link link for INLA "regression"
#' 
#' @import INLA
#' @import utils
#'
#' @return list with \describe{
#'   \item{result}{INLA call return}
#'   \item{data}{data frame; data fed to the INLA call}
#'   \item{grid}{vector; sequence from minimum sample time to maximum sample time, of length one more than lengthout.}
#'   \item{x}{vector; coal data time}
#'}
#'
infer_coal_samp <- function(
  samp_times, coal_times, n_sampled=NULL, lengthout = 100,
  rd_prob_fn = NULL, 
  fns = NULL, log_fns = TRUE,
  fns_coeff_prior_mean = NULL, fns_coeff_prior_prec = NULL,
  prec_alpha = 0.01, prec_beta = 0.01, 
  beta1_mean = 0, beta1_prec = 0.001, 
  use_samp = FALSE, simplify = FALSE, 
  events_only = FALSE, derivative = FALSE, link = 1 
){
  
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      'INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
      call. = FALSE
    )
    
  }
  
  if (min(coal_times) < min(samp_times)) {
    stop("First coalescent time occurs before first sampling time")
  }
  
  if (max(samp_times) > max(coal_times)) {
    stop("Last sampling time occurs after last coalescent time")
  }
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout + 1)
  
  if (is.null(n_sampled)) {
    n_sampled <- rep(1, length(samp_times))
  }
  
  coal_data <- coal_stats(
    grid = grid, 
    samp_times = samp_times, 
    n_sampled = n_sampled,
    coal_times = coal_times
  )
  
  if (simplify) { 
    coal_data <- with(
      coal_data, 
      condense_stats(time = time, event = event, E = E)
    )
    
  }
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (!use_samp) {
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    
    formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
    
    family <- "poisson"
    
  } else if (use_samp) {
    if (events_only) {
      samp_data <- samp_stats(
        grid = grid, 
        samp_times = samp_times
      )
      
    } else {
      samp_data <- samp_stats(
        grid = grid, 
        samp_times = samp_times, 
        n_sampled = n_sampled
      )
      
    }
    
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
    
    if (is.null(fns)) {
      formula <- Y ~ -1 + beta0 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed = FALSE, param = c(beta1_mean, beta1_prec))
      
    } else {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      
      for (fni in fns) {
        if (log_fns) {
          vals <- cbind(
            vals, 
            c(rep(0, bins), log(fni(samp_data$time)))
          )
          
        } else {
          vals <- cbind(
            vals, 
            c(rep(0, bins), fni(samp_data$time))
          )
          
        }
      }
      
      data$fn <- vals
      
      if ("rd_fn" %in% names(fns_coeff_prior_mean)) {# already in correct format
        if (ncol(vals) == 1) {
          colnames(data$fn) <- c("rd_fn")
        } else {
          colnames(data$fn) <- c("rd_fn" , paste0("fn", 1:(ncol(vals) - 1)))
        }
        
      } else {
        colnames(data$fn) <- paste0("fn", 1:ncol(vals))
        
        if (!is.null(fns_coeff_prior_mean)) {
          names(fns_coeff_prior_mean) <- c(paste0("fn", 1:(length(fns_coeff_prior_mean))))
          fns_coeff_prior_mean <- as.list(fns_coeff_prior_mean)
        }
        
        if (!is.null(fns_coeff_prior_prec)) {
          names(fns_coeff_prior_prec) <- c(paste0("fn", 1:(length(fns_coeff_prior_prec))))
          fns_coeff_prior_prec <- as.list(fns_coeff_prior_prec)
        }
        
        fns_coeff_prior_mean <- c(fns_coeff_prior_mean, list(default = 0))
        fns_coeff_prior_prec <- c(fns_coeff_prior_prec, list(default = 0.001))
      }
      
      formula <- Y ~ -1 + beta0 + fn +
        f(time, model = "rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy = "time", fixed = FALSE, param = c(beta1_mean, beta1_prec))
     
    }
    
    family <- c("poisson", "poisson")
    
  } else {
    stop("Invalid use_samp value, should be logical.")
    
  }
  
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid) / 2
    A <- diag(1 / diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
    
  } else {
    lc_many <- NULL
    
  }
  
  if (is.null(rd_prob_fn)) {
    regression_offset <- data$E_log
    
  } else {
    data$log_rd_prob <- c(
      rep(0, length(coal_data$time)),
      rd_prob_fn(samp_data$time)
    )
    
    regression_offset <- data$E_log + data$log_rd_prob
    
  }
    
  if (is.null(fns)) {
    mod <- INLA::inla(
      formula, 
      family = family, 
      data = data,
      lincomb = lc_many, 
      offset = regression_offset,
      control.predictor = list(compute = TRUE, link = link)
    )
    
  } else {
    mod <- INLA::inla(
      formula, 
      family = family, 
      data = data,
      lincomb = lc_many, 
      offset = regression_offset,
      control.predictor = list(compute = TRUE, link = link),
      control.fixed = list(
        mean = fns_coeff_prior_mean,
        prec = fns_coeff_prior_prec
      )
    )
    
  }
  
  list(
    result = mod, 
    data = data, 
    grid = grid, 
    x = coal_data$time
  )
}

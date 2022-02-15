#' gen_INLA_args Title TODO
#'
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#'
#' @return data frame with \describe{
#'   \item{coal_factor}{TODO}
#'   \item{s}{sampling and coalescence times combined and sorted} 
#'   \item{event}{TODO}
#'   \item{lineages}{TODO}
#' }
#'
gen_INLA_args <- function(samp_times, coal_times, n_sampled)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")
  
  if (length(intersect(coal_times, samp_times)) > 0)
    warning(
      "Coincident sampling event and coalescent event: results may be unpredictable."
    )
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)
  
  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages * (lineages - 1) / 2
  
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  
  list(
    coal_factor = coal_factor, 
    s = sorting$x, 
    event = event, 
    lineages = lineages
  )
  
}


#' Title TODO
#'
#' @param grid numeric vector; TODO
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#' @param log_zero numeric; TODO
#'
#' @return a data frame with \describe{
#'   \item{time}{TODO}
#'   \item{event}{TODO}
#'   \item{E}{TODO}
#'   \item{E_log}{TODO}
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
  s <- args$s
  event <- args$event
  
  grid_trimmed <- setdiff(x = grid, y = s)
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering]
  
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


#' Infer coalescence sample
#'
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#' @param lengthout numeric; length of time for the estimation
#' @param rd_prob_fn function; a function that takes in a vector of sampling times and returns the probability of a collected sample having been reported
#' @param fns function; list of covariate functions for the sampling intensity
#' @param prec_alpha numeric; hyperparameter alpha for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param prec_beta numeric; hyperparameter beta for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param beta1_prec numeric; precision of the normal prior assigned to the coefficient of the log effective population size in the sampling intensity formula
#' @param use_samp Boolean; TODO
#' @param log_fns Boolean; specifies if the log of the covariate functions, fns, needs to be take. FALSE indicates that the covriate function already returns log transformed values
#' @param simplify Boolean; TODO
#' @param events_only Boolean; TODO
#' @param derivative Boolean; TODO
#' @param link link for INLA "regression"
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
    rd_prob_fn = NULL, fns = NULL,
    prec_alpha = 0.01, prec_beta = 0.01, beta1_prec = 0.001, 
    use_samp = FALSE, log_fns = TRUE, simplify = FALSE, 
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
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
      
    } else {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      
      for (fni in fns) {
        if (log_fns) {
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        } else {
          vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
        }
        
      }
      
      data$fn <- vals
      
      formula <- Y ~ -1 + beta0 + fn +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
      
    }
    
    family <- c("poisson", "poisson")
    
  } else {
    stop("Invalid use_samp value, should be boolean.")
    
  }
  
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
    
  } else {
    lc_many <- NULL
    
  }
  
  if (is.null(rd_prob_fn)) {
    mod <- INLA::inla(
      formula, 
      family = family, 
      data = data,
      lincomb = lc_many, 
      offset = data$E_log,
      control.predictor = list(compute = TRUE, link = link)
    )
    
  } else {
    if (data$E_log == length(samp_data$time)) {
      data$rd_prob <- log(rd_prob_fn(samp_data$time))
      
    } else if (data$E_log == 2 * length(samp_data$time)) {
      data$rd_prob <- rep(log(rd_prob_fn(samp_data$time)), 2)
      
    } else {
      stop("E_log length unpredictable")
      
    }
    
    mod <- INLA::inla(
      formula, family = family, data = data,
      lincomb = lc_many, offset = data$E_log + data$rd_prob,
      control.predictor = list(compute = TRUE, link = link)
    )
    
  }
  
  list(
    result = mod, 
    data = data, 
    grid = grid, 
    x = coal_data$time
  )
}


#' Infer coalescence
#'
#' @param samp_times numeric vector; sampling times
#' @param coal_times numeric vector; coalescence times
#' @param n_sampled numeric vector; The number of sampling events at each sampling time, length matching the length of samp_times
#' @param lengthout numeric; length of time for the estimation 
#' @param prec_alpha numeric; hyperparameter alpha for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param prec_beta numeric; hyperparameter beta for the gamma prior of kappa, the precision of the Gaussian random walk prior
#' @param simplify Boolean; TODO
#' @param derivative Boolean; TODO
#'
#' @return list with \describe{
#'   \item{result}{INLA call return}
#'   \item{data}{data frame; data fed to the INLA call}
#'   \item{grid}{vector; sequence from minimum sample time to maximum sample time, of length one more than lengthout.}
#'   \item{x}{vector; coal data time}
#' }
#'
infer_coal <- function(
    samp_times, coal_times, n_sampled = NULL, lengthout = 100,
    prec_alpha = 0.01, prec_beta = 0.01, 
    simplify = FALSE, derivative = FALSE
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
  
  if (max(samp_times) > max(coal_times)){
    stop("Last sampling time occurs after last coalescent time")
  }
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout + 1)
  
  if (is.null(n_sampled)){
    n_sampled <- rep(1, length(samp_times))
  }
  
  coal_data <- coal_stats(
    grid = grid, 
    samp_times = samp_times, 
    n_sampled = n_sampled,
    coal_times = coal_times
  )
  
  if (simplify){
    coal_data <- with(
      coal_data, 
      condense_stats(time = time, event = event, E = E)
    )
    
  }
  
  data <- with(
    coal_data, 
    data.frame(y = event, time = time, E_log = E_log)
  )
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
  
  if (derivative) {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid) / 2
    A <- diag(1 / diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
    
    mod <- INLA::inla(
      formula, 
      family = "poisson", 
      data = data, 
      lincomb = lc_many,
      control.predictor = list(compute = TRUE),
      control.inla = list(lincomb.derived.only = FALSE)
    )
    
  } else {
    mod <- INLA::inla(
      formula, 
      family = "poisson", 
      data = data, 
      offset = data$E_log,
      control.predictor = list(compute=TRUE)
    )
    
  }
  
  list(
    result = mod, 
    data = data, 
    grid = grid, 
    x = coal_data$time
  )
  
}

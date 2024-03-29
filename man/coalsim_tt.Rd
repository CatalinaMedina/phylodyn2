% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coalsim.R
\name{coalsim_tt}
\alias{coalsim_tt}
\title{Coalescent time transfermation method}
\usage{
coalsim_tt(samp_times, n_sampled, traj, val_upper = 10, ...)
}
\arguments{
\item{samp_times}{numeric vector of sampling times}

\item{n_sampled}{numeric vector of number of samples collected at each sampling time}

\item{traj}{function that returns effective population size at time t.}

\item{val_upper}{numeric used by time-transformation method to set a starting
point for its dynamic numerical integration upper bound.}

\item{...}{additional arguments to be passed to \code{traj} function}
}
\value{
A list containing vectors of coalescent times \code{coal_times}, 
  intercoalescent times \code{intercoal_times}, and number of active lineages
  \code{lineages}, as well as passing along \code{samp_times} and
  \code{n_sampled}.
}
\description{
Coalescent time transfermation method
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coalsim.R
\name{coalsim}
\alias{coalsim}
\title{Simulate from inhomogeneous, heterochronous coalescent}
\usage{
coalsim(
  samp_times,
  n_sampled,
  traj,
  method = "tt",
  val_upper = 10,
  lower_bound = 1,
  ...
)
}
\arguments{
\item{samp_times}{numeric vector of sampling times.}

\item{n_sampled}{numeric vector of samples taken per sampling time.}

\item{traj}{function that returns effective population size at time t.}

\item{method}{which sampling method to use. "tt" invoke time-transformation
method, "thin" invokes thinning method.}

\item{val_upper}{numeric used by time-transformation method to set a starting
point for its dynamic numerical integration upper bound.}

\item{lower_bound}{numeric lower limit of \code{traj} function on its
support.  Used only by thinning method.}

\item{...}{additional arguments to be passed to \code{traj} function.}
}
\value{
A list containing vectors of coalescent times \code{coal_times}, 
  intercoalescent times \code{intercoal_times}, and number of active lineages
  \code{lineages}, as well as passing along \code{samp_times} and
  \code{n_sampled}.
}
\description{
Simulate from inhomogeneous, heterochronous coalescent
}
\examples{
coalsim(0:2, n_sampled = 3:1, traj = unif_traj, lower_bound = 10)
}

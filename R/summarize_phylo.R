#' Summarizes phylogeny, connecting sampling times to tips
#'
#' @param phy \code{phylo} object containing a phylogeny
#' @param backwards Boolean; TODO
#' 
#' @importFrom ape dist.nodes
#' @import utils
#'
#' @return list with \describe{
#'   \item{ungrouped_samp_df}{data frame connecting sampling times to corresponding tip labels}
#'   \item{samp_times}{sampling times}
#'   \item{n_sampled}{number sampled at each sampling time}
#'   \item{coal_times}{coalescence times}
#' }
#' @export
#'
summarize_phylo2 <- function(phy, backwards = TRUE){
  if (class(phy) != "phylo"){
    stop("object \"phy\" is not of class \"phylo\"")
  }
  
  n_nodes <- phy$Nnode
  n_tips <- length(phy$tip.label)
  
  root_node <- phy$edge[1,1]
  raw_times <- ape::dist.nodes(phy)[root_node, ]
  
  raw_samp_df <- data.frame(
    "raw_samp_times" = utils::head(raw_times, n_tips), 
    "tip_labels" = phy$tip.label
  )
  
  raw_coal_times <- utils::tail(raw_times, n_nodes)
  
  if (backwards) {
    max_samp_time <- max(raw_samp_df$raw_samp_times)
    raw_coal_times <- max_samp_time - raw_coal_times
    raw_samp_df$raw_samp_times <- max_samp_time - raw_samp_df$raw_samp_times
    
  }
  
  raw_samp_df <- raw_samp_df[order(raw_samp_df$raw_samp_times), ]
  samp_tab <- table(raw_samp_df$raw_samp_times)
  samp_times <- as.numeric(names(samp_tab))
  n_sampled <- as.numeric(samp_tab)
  
  coal_times <- sort(as.numeric(raw_coal_times))
  
  # Assumes sampling times start at t = 0
  list(
    ungrouped_samp_df = raw_samp_df,
    samp_times = samp_times,
    n_sampled  = n_sampled,
    coal_times = coal_times
  )
}


#' branching_sampling_times Title TODO
#'
#' @param tr \code{phylo} object containing a phylogeny
#' 
#' @importFrom ape node.depth.edgelength
#'
#' @return numeric vector TODO
#'
branching_sampling_times <- function(tr){
  ##Updated by Julia Sep 1, 2021. The previous function assumed internal nodes are ordered
  if (class(tr) != "phylo"){
    stop("object \"tr\" is not of class \"phylo\"")
  }
  
  edge.mat <- tr$edge
  n.sample <- tr$Nnode + 1
  
  t.tot <- max(ape::node.depth.edgelength(tr))
  n.t <- t.tot - ape::node.depth.edgelength(tr)
  
  xx <- as.numeric(rep(NA, 2 * n.sample - 1))
  names(xx) <- as.character(c(-(1:(n.sample - 1)), 1:n.sample))
  
  xx[1:(n.sample - 1)] <- sort(
    n.t[(n.sample + 1):length(n.t)], 
    decreasing = TRUE
  )
  xx[n.sample:length(xx)] <- sort(
    n.t[1:n.sample],
    decreasing = TRUE
  )
  
  xx
}


#' heterochronous_gp_stat Title TODO
#'
#' @param phy a \code{phylo} object containing a phylogeny
#' @param tol numeric; tolerance
#'
#' @return list with \describe{
#'   \item{coal_times}{numeric vector; coalescent times}
#'   \item{samp_times}{numeric vector; sampling times}
#'   \item{n_sampled}{numeric vector; number sampled at each sampling time}
#' }
#' 
heterochronous_gp_stat <- function(phy, tol = 0.0){
  #Update Aug 2015 by Julia. Adhoc for simulation with a tolerance parameters
  b.s.times <- branching_sampling_times(phy)
  
  int.ind <- which(as.numeric(names(b.s.times)) < 0)
  tip.ind <- which(as.numeric(names(b.s.times)) > 0)
  
  num.tips <- length(tip.ind)
  
  num.coal.events <- length(int.ind)
  
  sampl.suf.stat <- rep(NA, num.coal.events)
  coal.interval <- rep(NA, num.coal.events)
  coal.lineages <- rep(NA, num.coal.events)
  
  sorted.coal.times <- sort(b.s.times[int.ind])
  names(sorted.coal.times) <- NULL
  sampling.times <- sort((b.s.times[tip.ind]))
  
  for (i in 2:length(sampling.times)) {
    if ((sampling.times[i] - sampling.times[i - 1]) < tol) {
      sampling.times[i] <- sampling.times[i - 1]
    }
    
  }
  
  unique.sampling.times <- unique(sampling.times)
  sampled.lineages <- NULL
  
  for (sample.time in unique.sampling.times) {
    sampled.lineages <- c(sampled.lineages, sum(sampling.times == sample.time))
  }
  
  list(
    coal_times = sorted.coal.times,
    samp_times = unique.sampling.times,
    n_sampled = sampled.lineages
  )
}


#' Summarize a phylogeny
#' 
#' @param phy a \code{phylo} object containing a phylogeny.
#'   
#' @return A list containing vectors of sampling times \code{samp_times}, number 
#'   sampled per sampling time \code{n_sampled}, and coalescent times
#'   \code{coal_times}.
#' @export
#' 
summarize_phylo <- function(phy) {
  hgpstat <- heterochronous_gp_stat(phy)
  
  list(
    samp_times = hgpstat$samp_times,
    n_sampled  = hgpstat$n_sampled,
    coal_times = hgpstat$coal_times
  )
}
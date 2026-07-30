#' Trim phylogeny with reporting delays
#'
#' Function to simulate reporting delays from provided \code{calc_report_prob()} according to sampling times. Useful for subsetting tree to only observed samples and keeping track of original time zero.
#'
#' @param full_phylo a phylo object containing a phylogeny
#' @param calc_report_prob a function which returns probability of a sample having been reported by time zero of analysis given the time at which the sample was collected 
#' @param original_samp_times a numeric vector of the sampling times for the \code{full_tree}
#'
#' @importFrom ape drop.tip
#' @importFrom stats rbinom
#'
#' @returns A list with two elements \code{obs_phylo} contains the tree in phylo format only containing samples that were simulated as observed, \code{time0_offset} a numeric of the most recently observed sampling time. This is useful in simulations because the phylo object will automatically consider time zero to be the most recently observed sampling time. 
#' @export
#'
trim_tree_with_delays <- function(
    full_phylo,
    calc_report_prob,
    original_samp_times
) {
  
  # Match coalescence to sampling times
  tips_matched <- data.frame(
    connect_sample_to_tips(full_phylo)$ungrouped_samp_df,
    "true0_samp_times" = original_samp_times
  )
  # simulate reporting delays for our simulated sampling times using the Bernoulli distribution
  tips_matched$prob_reported <-
    calc_report_prob(tips_matched$true0_samp_times)
  
  tips_matched$reported <- mapply(
    FUN = rbinom, 
    n = 1, 
    p = tips_matched$prob_reported,  
    MoreArgs = list(size = 1)
  )
  
  tips_matched$drop <- ifelse(
    tips_matched$reported == 0, yes = TRUE, no = FALSE)
  
  # trim the tree 
  obs_phylo <- ape::drop.tip(
    full_phylo,
    tip = tips_matched$tip_labels[tips_matched$drop]
  )
  
  # Match observed tree tips to sampling times
  obs_tips_matched <- data.frame(
    connect_sample_to_tips(obs_phylo)$ungrouped_samp_df,
    "true0_samp_times" = tips_matched$true0_samp_times[!tips_matched$drop]
  )
  
  list(
    "obs_phylo" = obs_phylo,
    "time0_offset" = min(obs_tips_matched$true0_samp_times)
  )
}
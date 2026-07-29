#' Title
#'
#' @param full_tree `phylo` object containing a phylogeny
#' @param calc_report_prob 
#' @param original_samp_times 
#'
#' @returns
#' @export
#'
#' @examples
trim_tree_with_delays <- function(full_tree, calc_report_prob, original_samp_times) {
  
  # Match coalescence to sampling times
  tips_matched <- data.frame(
    connect_sample_to_tips(full_tree$newick)$ungrouped_samp_df,
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
  
  tips_matched$drop <- ifelse(tips_matched$reported == 0, yes = TRUE, no = FALSE)
  
  # trim the tree 
  obs_tree <- ape::drop.tip(
    full_tree$newick,
    tip = tips_matched$tip_labels[tips_matched$drop]
  )
  
  # Match observed tree tips to sampling times
  obs_tips_matched <- data.frame(
    connect_sample_to_tips(obs_tree)$ungrouped_samp_df,
    "true0_samp_times" = tips_matched$true0_samp_times[!tips_matched$drop]
  )
  
  list(
    "obs_tree" = obs_tree,
    "time0_offset" = min(obs_tips_matched$true0_samp_times)
  )
}
#' Returns a phylo object from the arguments generated with coalsim
#' 
#' @param args is a list containing vectors of coalescent times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}, etc. This list is the output of coalsim
#' 
#' @import ape
#' 
#' @return A list with two elements \code{newikck} contains the tree in phylo format, \code{labels} a vector with tip labels 
#' @export
#' 
#' @examples
#' constant <- function(x){return (rep(1,length(x)))}
#' simulation1 <- coalsim(samp_times = 0, n_sampled = 10, traj = constant)
#' tree1 <- generate_newick(simulation1)
#' plot(tree1$newick)
generate_newick <- function(args) {
  n <- sum(args$n_sampled)
  
  labels <- paste(
    rep("t", n), seq(1, n, 1), rep("_", n),
    rep(args$samp_times[1], args$n_sampled[1]),
    sep = ""
  )
  
  #we could chose labels at random to coalesce but since the process is exchangeable, we don't care. At least not for now
  tb <- args$n_sampled[1] #Total branches (initial)
  s <- 0 #time for branch lengths
  temp_labels <- labels[1:tb]
  temp_times <- rep(args$samp_times[1], args$n_sampled[1])
  initial.row <- 2
  args2 <- gen_INLA_args(
    samp_times = args$samp_times, 
    n_sampled = args$n_sampled, 
    coal_times = args$coal_times
  )
  
  for (j in 2:length(args2$event)) {
    if (args2$event[j] == 1) {
      s <- args2$s[j]
      ra <- sample(tb, 1) #choose at random one of them, the other is the one to the right so not really random
    
      if (ra < tb) {
        new_label <- paste(
          "(", temp_labels[ra], ":", s - temp_times[ra], ",", 
          temp_labels[ra + 1], ":", s - temp_times[ra + 1], ")",
          sep = ""
        )
        temp_labels[ra] <- new_label
        temp_labels <- temp_labels[-(ra + 1)]
        temp_times[ra] <- s
        temp_times <- temp_times[-(ra + 1)]
        
      } else {
        new_label <- paste(
          "(",temp_labels[ra], ":", s - temp_times[ra], ",",
          temp_labels[1], ":", s - temp_times[1], ")",
          sep = ""
        )
        temp_labels[1] <- new_label
        temp_labels <- temp_labels[-(ra)]
        temp_times[1] <- s
        temp_times <- temp_times[-(ra)]
        
      }
      tb <- tb - 1
      
    } else { #I will be adding samples at 
      s <- args2$s[j]; 
      
      if (args$n_sample[initial.row] == 1) {
        temp_labels <- c(
          temp_labels,
          labels[cumsum(args$n_sampled)[initial.row]]
        )
        initial.row <- initial.row + 1
        tb <- tb + 1
        temp_times <- c(temp_times, s) 
        
      } else {
        end <- cumsum(args$n_sampled)[initial.row]
        ini <- cumsum(args$n_sampled)[initial.row - 1] + 1
        
        for (k in ini:end){
          temp_labels <- c(temp_labels, labels[k])
          tb <- tb + 1
          temp_times <- c(temp_times, s)      
        }
        
        initial.row <- initial.row + 1
        
      }
      
    }
    
  }  
  
  out.tree <- ape::read.tree(text = paste(temp_labels, ";", sep = ""))
  
  list(newick = out.tree, labels = labels)
}
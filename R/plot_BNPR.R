#' Plot a BNPR output
#' 
#' @param BNPR_out output of BNPR or BNPR_PS.
#' @param traj function summarizing the true effective population size 
#'   trajectory.
#' @param xlim numeric x-axis interval.
#' @param ylim numeric y-axis interval.
#' @param nbreaks integer number of bins for sampling heatmap.
#' @param lty line type for estimated trajectory.
#' @param lwd line width for estimated trajectory.
#' @param col color for estimated trajectory.
#' @param main character main plot title.
#' @param log_y boolean whether to plot y coordinates on log scale.
#' @param ylab character y-axis label.
#' @param xlab character x-axis label.
#' @param xmarline numeric if not using default x-axis labels, how far to put 
#'   the labels from the axis.
#' @param axlabs character vector x-axis labels.
#' @param traj_lty,traj_lwd,traj_col line type, line width, and line color for 
#'   the true trajectory.
#' @param newplot boolean whether to create a new plot or superimpose over a 
#'   previously open plot.
#' @param credible_region logical whether to display pointwise credible region.
#' @param heatmaps boolean whether to display sampling and coalescent heatmaps.
#' @param heatmap_labels boolean whether to display labels on heatmaps.
#' @param heatmap_labels_side string which side of plot to display heatmaps.
#' @param heatmap_labels_cex numeric scaling factor for heatmap labels.
#' @param heatmap_width numeric how wide heatmaps should be.
#' @param yscale numeric scaling applied to all effective population
#'   calculations.
#' @param ... additional arguments to be passed onto plot().
#'   
#' @import graphics
#'   
#' @export
plot_BNPR <- function(
    BNPR_out, traj = NULL, xlim = NULL, ylim = NULL, nbreaks = 40,
    lty = 1, lwd = 2, col = "black", main = "", log_y = TRUE,
    ylab = "Effective Population Size",
    xlab = "Time", xmarline = 3, axlabs = NULL,
    traj_lty = 2, traj_lwd = 2, traj_col = col,
    newplot = TRUE, credible_region = TRUE,
    heatmaps = TRUE, heatmap_labels = TRUE,
    heatmap_labels_side = "right", heatmap_labels_cex = 0.7,
    heatmap_width = 7, yscale = 1, ...
  ){
  
  grid <- BNPR_out$grid
  
  if (is.null(xlim)) {
    xlim <- c(max(grid), min(grid))
  }
  
  mask <- BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
  
  t <- BNPR_out$x[mask]
  
  y <- BNPR_out$effpop[mask] * yscale
  yhi <- BNPR_out$effpop975[mask] * yscale
  ylo <- BNPR_out$effpop025[mask] * yscale
  
  if (newplot) {
    if (is.null(ylim)) {
      ymax <- max(yhi)
      ymin <- min(ylo)
      
    } else {
      ymin <- min(ylim)
      ymax <- max(ylim)
      
    }
    
    if (heatmaps) {
      yspan <- ymax / ymin
      yextra <- yspan^(1/10)
      ylim <- c(ymin / (yextra^1.35), ymax)
      
    } else {
      ylim <- c(ymin, ymax)
      
    }
    
    if (is.null(axlabs)) {
      if (log_y) {
        graphics::plot(
          1, 1, type = "n", log = "y",
          xlab = xlab, ylab = ylab, main = main,
          xlim = xlim, ylim = ylim, ...
        )
        
      } else {
        graphics::plot(
          1, 1, type = "n",
          xlab = xlab, ylab = ylab, main = main,
          xlim = xlim, ylim = ylim, ...
        )
        
      }
      
    } else {
      if (log_y) {
        graphics::plot(
          1, 1, type = "n", log = log,
          xlab = "", ylab = ylab, main = main,
          xlim = xlim, ylim = ylim, xaxt = "n", ...
        )
        graphics::axis(1, at = axlabs$x, labels = axlabs$labs, las = 2)
        graphics::mtext(text = xlab, side = 1, line = xmarline)
        
      } else {
        graphics::plot(
          1, 1, type = "n",
          xlab = "", ylab = ylab, main = main,
          xlim = xlim, ylim = ylim, xaxt = "n", ...
        )
        graphics::axis(1, at = axlabs$x, labels = axlabs$labs, las = 2)
        graphics::mtext(text = xlab, side = 1, line = xmarline)
        
      }
      
    }
    
  }
  
  if (credible_region) {
    shade_band(x = t, ylo = ylo, yhi = yhi, col = "lightgray")
  }
  
  if (!is.null(traj)) {
    graphics::lines(t, traj(t), lwd = traj_lwd, lty = traj_lty, col = traj_col)
  }
  
  if (newplot) {
    if (heatmaps) {
      samps <- rep(BNPR_out$samp_times, BNPR_out$n_sampled)
      samps <- samps[samps <= max(xlim) & samps >= min(xlim)]
      
      coals <- BNPR_out$coal_times
      coals <- coals[coals <= max(xlim) & coals >= min(xlim)]
      
      breaks <- seq(min(xlim), max(xlim), length.out=nbreaks)
      h_samp <- graphics::hist(samps, breaks=breaks, plot=FALSE)
      h_coal <- graphics::hist(coals, breaks=breaks, plot=FALSE)
      
      hist2heat(h_samp, y= ymin / yextra^0.5, wd = heatmap_width)
      hist2heat(h_coal, y= ymin / yextra, wd = heatmap_width)
      
      if (heatmap_labels) {
        if (heatmap_labels_side == "left") {
          lab_x <- max(xlim)
          lab_adj <- 0
          
        } else if (heatmap_labels_side == "right") {
          lab_x <- min(xlim)
          lab_adj <- 1
          
        } else {
          warning('heatmap_labels_side not "left" or "right", defaulting to right')
          lab_x <- min(xlim)
          lab_adj <- 1
          
        }
        
        graphics::text(
          x = lab_x, y = ymin/(yextra^0.20), labels = "Sampling events",
          adj = c(lab_adj, 0), cex = heatmap_labels_cex
        )
        graphics::text(
          x = lab_x, y = ymin/(yextra^1.25), labels = "Coalescent events",
          adj = c(lab_adj, 1), cex = heatmap_labels_cex
        )
        
      }
    }
  }
  
  graphics::lines(t, y, lwd = lwd, col = col, lty = lty)
}


#' Plot shade bands
#'
#' @param x x values corresponding to the band bounds
#' @param ylo lower band bound
#' @param yhi upper band bound
#' @param xlim x plotting limits
#' @param col color
#' 
#' @importFrom graphics polygon
#'
shade_band <- function(x, ylo, yhi, xlim = NULL, col = "gray"){
  if (is.null(xlim)) {
    xlim <- c(0, Inf)
  }
  
  mask <- x >= min(xlim) & x <= max(xlim)
  
  x <- x[mask]
  ylo <- ylo[mask]
  yhi <- yhi[mask]
  
  graphics::polygon(c(x, rev(x)), c(yhi, rev(ylo)), col = col, border = NA)
}


#' Plot heatmap of a histogram
#' 
#' @param hist \code{histogram} object to be displayed.
#' @param y numeric y-coordinate to display heatmap.
#' @param wd numeric width of heatmap in y-units.
#' 
#' @importFrom graphics segments
#'
hist2heat <- function(hist, y, wd) {
  breaks <- hist$breaks
  counts <- hist$counts
  upper <- max(counts)
  n <- length(counts)
  
  f <- counts / upper
  cols <- sprintf(
    "#%02x%02x%02x", 
    floor((1 - f) * 255), 
    floor((1 - f) * 255), 
    floor((1 - f) * 255)
  )

  graphics::segments(
    x0 = breaks[1:n], 
    y0 = y, 
    x1 = breaks[2:(n + 1)], 
    y1 = y, 
    lwd = wd, col = cols, lend = 1
  )
}

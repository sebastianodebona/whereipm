#' Plot spatial distribution through time
#' 
#' This function plots the distribution of the population across space over time
#' 
#' @export
#' 
#' @param IPM.model Matrix of population sizes in space and across the time of the simulation. This should be an output from the function \code\{\link{runIPM}}.
#' 
#' @param Hue Character value defining which colour hue to use for the lines. The colours are generated with the \code\{\link[RColorBrewer]{brewer.pal}} function. Default is \code{'Blues'}. To see more colour combinations see \url{http://colorbrewer2.org/}
#' @param breaks How many lines to be drawn, each representing a different time point along the population dynamics. The point will be equally spaced. The defaul is 10.
#' 

draw.pop.dyn <- function(IPM.model, Hue='Blues', breaks=10){
  
  # creating an empty plot to draw lines in
  plot(0, xlim=c(min(IPM.model$y), max(IPM.model$y)), ylim=c(0, max(IPM.model$pop.size)), col='white', ylab='population size', xlab='space')
  
  # creating a vector to select columns from the population projection
  selector <- round(seq(1,IPM.model$n.iterations, length.out = breaks))
  
  # determining the vector of colour values
  colorz <- as.list(grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,Hue)[2:9])(breaks))
  
  # turning columns of population size into items of a list, to feed them into mapply
  popsizes <- split(IPM.model$pop.size[, selector], 
                    rep(1:ncol(IPM.model$pop.size[, selector]), each = nrow(IPM.model$pop.size[, selector])))
  
  # using mapply to draw the lines                                                         
  invisible(mapply(draw.pop.distr, pop = popsizes, col = colorz, MoreArgs = list(x = IPM.model$y,  new.plot = FALSE)))
  
}



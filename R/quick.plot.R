#' Generate quick plots to visualise IPM results
#' 
#' Plots change in population distribution through space, population spread through time and population size through time.
#' 
#' @export
#' 
#' @param IPM.model A fitted IPM object, produced using the function \code{\link{runIPM}}.
#' @param legend Legend to be displayed in the top right corner of the plot, to identify the IPM 
#' 
#' 

quick.plot <- function(IPM.model, legend){
  
  # storing original parameters
  op <- par()
  # setting new parameters
  par(mfrow = c(3,2), mar = c(4,4,2,1), oma = c(2,2,2,2))
  # defining layout matrix
  m <- matrix(data = c(1,1,2,1,1,3), nrow = 3)
  layout(mat = m)
  
  # plotting population through space and time
  draw.pop.dyn(IPM.model, breaks=IPM.model$n.iterations+1, Hue='Greys')
  legend('topright', legend = legend, bty = 'n', cex = 1.5)
  legend('topleft', legend = paste(
    paste("n = ", round(colSums(IPM.model$pop.size)[IPM.model$n.iterations+1]), sep = ""),
    paste("iter = ", IPM.model$n.iterations, sep = ""),
    paste("lambda = ", round(colSums(IPM.model$pop.size)[IPM.model$n.iterations+1]/
                                colSums(IPM.model$pop.size)[IPM.model$n.iterations], 3), sep = ""),
    sep = "\n"
    ), bty = 'n')
  
  # plotting spread through time
  plot(IPM.model$pop.spread + 0.5, col='black', type = 'o', xlab = 'Time', ylab = 'spread')
  # plotting population growth through time
  plot(colSums(IPM.model$pop.size), col='black', type = 'o', xlab = 'Time', ylab = 'population size')
  
  # resetting graphic parameters
  suppressWarnings(par(op))
  
}
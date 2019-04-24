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
  par(mfrow = c(3,2), mar = c(3,4,2,1), oma = c(2,2,0,2))
  # defining layout matrix
  m <- matrix(data = c(1,1,2,4,1,1,3,5), nrow = 4)
  layout(mat = m)
  
  # plotting population through space and time
  draw.pop.dyn(IPM.model, breaks=IPM.model$n.iterations+1, Hue='Greys')
  legend('topright', legend = legend, bty = 'n', cex = 1.5)
  legend('topleft', legend = paste(
    paste("n = ", round(colSums(IPM.model$pop.size)[IPM.model$n.iterations+1]), sep = ""),
    paste("iter = ", IPM.model$n.iterations, sep = ""),
    paste("avg_lambda = ", round(IPM.model$lambda_avg, 3), sep = ""),
    sep = "\n"
    ), bty = 'n')
  
  
  # plotting population growth through time
  plot(colSums(IPM.model$pop.size), col='black', type = 'o', xlab = '', 
       ylab = "", main = 'population size')
  # plotting lambda through time
  plot(IPM.model$lambdas, col = "red", type = "o", xlab = "", 
       ylab = "", main = "lambda")
  # plotting spread through time
  plot(IPM.model$pioneer, col='black', type = 'o', xlab = '', 
       ylab = "", main = 'position of pioneer (m)')
  title(xlab = "time", line = 2)
  # plotting front behaviour
  plot(IPM.model$front, col = 'black', type = 'o', xlab = "", 
       ylab = "", main = "position of front (m)")
  title(xlab = "time", line = 2)
  
  
  # resetting graphic parameters
  suppressWarnings(par(op))
  
}


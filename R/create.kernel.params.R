#' Make a list of Kernel Parameters
#' 
#' Creates an object list to store kernel-defining parameters
#' 
#' @export
#' 
#' @param L lowest value the state variable can take. Default is 0.
#' @param U highest value the state variable can take. Default is 100.
#' @param n.bins number of bins defining the kernel. Determines how much the state variable is discretised. Default is 100
#' 

create.kernel.params <- function(L=0, U=100, n.bins=100){
  
  # deriving missing kernel parameters
  # defining meshpoints (separating bins) and midpoints (central point of each bin)
  b <- L+c(0:n.bins)*(U-L)/n.bins # meshpoints
  y <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
  h <- (U-L)/n.bins
  
  list(L=L,
       U=U,
       n.bins=n.bins,
       b=b,
       y=y,
       h=h)
  
}

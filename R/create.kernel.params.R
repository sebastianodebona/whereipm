#' Make a list of Kernel Parameters
#' 
#' Creates an object list to store kernel-defining parameters
#' 
#' @export
#' 
#' @param L lowest value the state variable can take.
#' @param U highest value the state variable can take.
#' @param n.bins number of bins defining the kernel. Determines how much the state variable is discretised.
#' 

create.kernel.params <- function(L, U, n.bins){
  
  # deriving missing kernel parameters
  # calculating bin size
  h <- (U-L)/n.bins
  
  # defining meshpoints (separating bins) and midpoints (central point of each bin)
  b <- seq(L, U, length.out = n.bins+1) # meshpoints
  y <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
  
  list(L=L,
       U=U,
       n.bins=n.bins,
       b=b,
       y=y,
       h=h)
  
}

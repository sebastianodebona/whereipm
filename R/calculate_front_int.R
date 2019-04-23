#' Calculates front position
#' 
#' The function calculates the point \code{x} in space for which the integral of the population distribution from \{x} to \infty equals \code{1/n}, where n is the total population size.
#' 
#' @param pop_vector Vector containing population size in each bin
#' @param kernel.params Parameter kernels used to characterise the population distribution, obtained from the IPM.
#' @param how_many Number of individuals to be found past x. Default is 1.
#' 
#' @export
#' 

calculate_front_int <- function(pop_vector, kernel.param, how_many = 1){
  
  # calculating cumulative population size from infinity back towards zero.
  index_front <- 
    max(which(rev(cumsum(rev(pop_vector))) >= how_many))
  
  x <- kernel.param$y[index_front]
  return(x)
  
}
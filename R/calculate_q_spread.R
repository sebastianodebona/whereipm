#' Calculates spread
#' 
#' The function calculates how many meters around a midpoint the proportion q of a given population occupies
#' 
#' @export
#' 
#' @param pop_vector Vector of the population
#' @param kernel.params Parameter kernels used to characterise the population distribution, often obtained from the IPM.
#' @param q Quantile to be evaluated. The default \code{q = 0.95} evaluates how much space is occupied by 95% of the population.

calculate_q_spread <- function(pop_vector, kernel.params, q = 0.95){
  
  # transforming pop into proportional pop
  prop_pop <- pop_vector/sum(pop_vector)
  
  # finding how many bins contain q*100% of the population
  # calculating midpoint
  midpoint <- kernel.params$n.bins/2
  # creating sequence from 0 to midpoint (to be used in pop subsetting)
  step_seq <- 0:midpoint
  # calculating what's the proportion of population contained in x steps form center (towards both directions)
  quant_pop <- sapply(step_seq, function(x){sum(prop_pop[(midpoint-x):(midpoint+x)])})
  # finding the number of bins containing q proportion of the pop
  steps <- min(which(quant_pop > q))*2
  # finding how many meters this corresponds to
  spread <- steps*kernel.params$h
   
  return(spread)
}

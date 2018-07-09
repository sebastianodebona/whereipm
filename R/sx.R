#' Survival function
#' 
#' Calculates the average survival probability at a given position x. Survival can be 
#' density dependent. The function takes a vector of local densities at the midpoints of the bins,
#' and parameters to define the slope and intercept of the survival per position curve
#' 
#' @export
#' 
#' @param density  A vector of length equal to the number of bins in the IPM kernel, containing local density at the midpoints
#' @param s.int  The intercept of the logistic regression determining density dependent survival (i.e. survival probability at density 0)
#' @param s.slo  The slope of the logistic regression determining density dependent survival (i.e. increment in survival probability per unit of density)
#' @param s.multiplier  Optional argument to inflate survival probability 
#' @return Returns a vector of average density-dependent survival probability at the midpoints of the bins of the IPM kernel


sx <- function(density, s.int, s.slo, s.multiplier=1) {
  
  if(plogis(s.int+s.slo*density)*s.multiplier<=1){
    plogis(s.int+s.slo*density)*s.multiplier
  } else {
    1
  }

}

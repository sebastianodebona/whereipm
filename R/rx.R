#'Fertility function
#'
#'Calculates the average number of offspring produced at a vector of given positions (in a vector form). Number of offspring can be density dependent. The function takes a vector of local densities at the midpoints of the bins, and parameters that determine slope and intercept of the number of offspring per position
#' @param density  A vector of length equal to the number of bins in the IPM kernel, containing local density at the midpoints
#' @param r.int  The intercept of the logistic regression determining density dependent offspring number (i.e. offspring number at density 0)
#' @param r.slo  The slope of the logistic regression determining density dependent offspring number (i.e. increment in offspring number per unit of density)
#' @param r.multiplier  Optional argument to inflate number of offspring
#' @return Returns a vector of average (density-dependent) number of offspring at the midpoints of the bins of the IPM kernel

# fertility function
rx <- function(density, r.int, r.slo, r.multiplier=1) {

  exp(r.int + r.slo * density)*r.multiplier
  
}
#'Fertility function
#'
#'Calculates the average number of offspring produced at a vector of given positions (in a vector form). Number of offspring can be density dependent. The function takes a vector of local densities at the midpoints of the bins, and parameters that determine slope and intercept of the number of offspring per position
#'
#' @export
#' @param density  A vector of length equal to the number of bins in the IPM kernel, containing local density at the midpoints
#' @param r.int  The intercept of the log regression determining density dependent offspring number (i.e. offspring number at density 0). If \code{distribution = 'ztpoisson'}, this becomes the intercept of the logistic regression determining the probability of producing NO offpsring.
#' @param r.slo  The slope of the log regression determining density dependent offspring number (i.e. increment in offspring number per unit of density). If \code{distribution = 'ztpoisson'}, this becomes the slope of the logistic regression determining the probability of producing NO offpsring.
#' @param f.int The intercept of a zero-truncated log regression determining the number of offspring (i.e. offspring number at density 0).
#' @param f.slo The slope of a zero-truncated log regression determining the number of offspring (i.e. offspring number at density 0).
#' @param r.multiplier  Optional argument to inflate number of offspring.
#' @param f.distribution Argument defining which error distribution was fitted for the number of offspring. The default is a Poisson distribution (\code{'poisson'}). The alternative is a zero truncated Poisson (\code{'ztpoisson'}), modeled as the probability to produce any offspring time the estimated average offspring number if any are produced.
#' @return Returns a vector of average (density-dependent) number of offspring at the midpoints of the bins of the IPM kernel

# fertility function
rx <- function(density, r.int, r.slo, r.multiplier = 1, f.distribution = 'ztpoisson', f.int = NA, f.slo = NA) {

  if(f.distribution=='poisson'){
    
    exp(r.int + r.slo * density)*r.multiplier
    
  } else if (f.distribution=='ztpoisson'){

    (1-plogis((r.int + r.slo * density)))*untruncZTP(f.int + f.slo*density)*r.multiplier
    }
}
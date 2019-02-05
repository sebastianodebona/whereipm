#' Calulate neighbourhood density
#' 
#' Creates a vector of neighbourhood densities (mean density around 100m)
#' 
#' @export
#' 
#' @param local_density vector of local densities in each bins of the space
#' @param midpoints vector of midpoints, of the same length as density, from the kernel parameters
#' @param bin.size size in m of a single bin, from the kernel.parameters
#' 

calculate_nb_density <- function(local_density, midpoints, bin.size){
  
  # using the size of the bin to establish where exactly the 100 m range sits
  # in the vector of positions. I cannot just add exactly +/- 50 because the
  # steps might not be integer. Multiplying an integer number of bins that fit
  # into 50 by the size of the bin will give me the closest value to 50 that I
  # can add or remove to the position x to calculate neighbouring density.
  multiplier <- round(50/bin.size)
  
  # creating vectors of lower and upper boundaries of the neighbourhood
  los <- sapply(midpoints - multiplier*bin.size, function(x)max(x, min(midpoints), na.rm = T))
  ups <- sapply(midpoints + multiplier*bin.size, function(x)min(x, max(midpoints), na.rm = T))
  
  # transforming those into indexes
  lo_index <- match(los, midpoints)
  up_index <- match(ups, midpoints)
  
  # calculating the vector of neighbourhood density
  nb_density <- mapply(function(lo, up, local_density){mean(local_density[lo:up])},
                       lo = lo_index, up = up_index,
                       MoreArgs = list(local_density = local_density))
  
  return(nb_density)
}

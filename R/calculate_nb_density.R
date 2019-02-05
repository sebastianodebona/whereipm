#' Calulate neighbourhood density
#' 
#' Creates a vector of neighbourhood densities (mean density around 100m)
#' 
#' @export
#' 
#' @param local_density vector of local densities in each bins of the space
#' @param midpoints vector of midpoints, of the same length as density, from the kernel parameters
#' 
#' 

calculate_nb_density <- function(local_density, midpoints){
  
  # creating vectors of lower and upper boundaries of the neighbourhood
  los <- sapply(midpoints - 50, function(x)max(x, min(midpoints), na.rm = T))
  ups <- sapply(midpoints + 50, function(x)min(x, max(midpoints), na.rm = T))
  
  # transforming those into indexes
  lo_index <- match(los, midpoints)
  up_index <- match(ups, midpoints)
  
  # calculating the vector of neighbourhood density
  nb_density <- mapply(function(lo, up, local_density){mean(local_density[lo:up])},
                       lo = lo_index, up = up_index,
                       MoreArgs = list(local_density = local_density))
  
  return(nb_density)
}

#' Calulate neighbourhood density
#' 
#' Creates a vector of neighbourhood densities (mean density around 100m)
#' 
#' @export
#' 
#' @param local_density vector of local densities in each bins of the space
#' @param bin.size size in m of a single bin, from the kernel.parameters
#' 

calculate_nb_density <- function(local_density, bin.size){
  
  # using the size of the bin to establish where exactly the 100 m range sits
  # in the vector of positions. I cannot just add exactly +/- 50 because the
  # steps might not be integer. I will check what is the integer number of bins 
  # that fit into 50m, based on the the size of the bin. This will correspond to
  # the number of steps needed to go back/forth 50m.
  steps <- round(50/bin.size)
  
  # creating vectors of indexes for the lower and upper boundaries of the neighbourhood
  lo_index <- sapply(seq_along(local_density) - steps,
                     function(x)max(x, 1, na.rm = T))
  
  up_index <- sapply(seq_along(local_density) + steps, 
                     function(x)min(x, length(local_density), na.rm = T))
  
  # calculating the vector of neighbourhood density
  nb_density <- mapply(function(lo, up, local_density){mean(local_density[lo:up])},
                       lo = lo_index, up = up_index,
                       MoreArgs = list(local_density = local_density))
  
  return(nb_density)
}

#' Calculate mean of a ZTP
#' 
#' Calculate the mean estimated count value of a zero-truncated Poisson defined by lambda
#' 
#' @export
#' @param lambda Estimated lambda (e.g. estimate from a \code{\link[glmmTMB:glmmTMB]{glmmTMB}} model)
#' @details See \url{https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution} for more info
#'  

untruncZTP <- function(lambda){
  
  exp(lambda)/(1-exp(-exp(lambda)))
  
  }

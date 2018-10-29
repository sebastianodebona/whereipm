#' Calculate mean of a ZTP
#' 
#' Calculate the mean estimated count value of a zero-truncated Poisson defined by lambda
#' 
#' @export
#' @param lambda Estimated lambda (e.g. estimate from a \code{\link[glmmTMB:glmmTMB]{glmmTMB}} model)
#' @details See \url{https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution} for more info. Here \lambda needs to be transformed into e^\lambda to solve the log link from the glmmTMB model. The truncated poisson starts behaving funnily if \lambda is < -20, so a threshold sets the result to 1 if \lambda < -20.
#'  

untruncZTP <- function(lambda){
  
  unTr <- exp(lambda)/(1-exp(-exp(lambda)))
  unTr[lambda < -20] = 1
  return(unTr)
  
  }

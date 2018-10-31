#' Return if existing
#' 
#' Returns a value (e.g. parameter estimate) if this exists. Otherwise, it returns zero as defaul or any other preferred value.
#' 
#' @param object An object, typically a call to extract a parameter from a model.
#' @param subs What should be returned if x is \code{object} is empty. Default is 0.
#' 

returnIfExists <- function(object, subs = 0){
  
  ifelse(length(object)==0,
         subs,
         object)
  
}

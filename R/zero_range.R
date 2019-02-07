#' Determine if a vector has range zero.
#' 
#' Checks whether all elements of a vector are the same number. All credits for this function to Hadley Wickham.
#' 
#' @export
#' @param x Numeric vector on which to perform the check
#' @param tol Tolerance (based on machine power)
#' @details All credits to Hadley Wickham. Source: https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector.
#' 

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x, na.rm = T) / mean(x, na.rm = T)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
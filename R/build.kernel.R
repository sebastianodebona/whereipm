#' Build Kernel
#' 
#' Build a kernel given a population size, kernel parameters and vital rate parameters
#' 
#' @export
#' 
#' @param kernel.params list of parameters to build the kernel (must include L, H and n.bins). See output of \code{create.kernel.params} for an example.
#' @param vital.params list of parameters to determine vital rate functions. Follow format of the output of \code{create.vital.params}
#' @param pop.size Vector of initial population size per position. The lenght of the vector must be equal to the number of bins specified in \code{kernel.params}. There is no default (as opposed to \code{\link{test.init.kernel}}), so population size \emph{must} be specified.
#' @param ... additional arguments, such as \code{unidirectional} for movement functions (see \code{\link{gx}})
#' @details As a preliminary step before, run a test with \code{\link{test.init.kernel}} with the same parameters, to check visually if the kernel builds without problems.
#' 

build.kernel <- function(kernel.params, vital.params, pop.size, ...){
  
  # deriving missing kernel parameters
  # defining meshpoints (separating bins) and midpoints (central point of each bin)
  b <- kernel.params$L+c(0:kernel.params$n.bins)*(kernel.params$U-kernel.params$L)/kernel.params$n.bins # meshpoints
  y <- 0.5*(b[1:kernel.params$n.bins]+b[2:(kernel.params$n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
  h <- (kernel.params$U-kernel.params$L)/kernel.params$n.bins

  # determining density
  density <- pop.size/kernel.params$h
  
  # survival kernel
  S <- diag(sx(density=density, 
               s.int=vital.params$s.int, 
               s.slo=vital.params$s.slo, 
               s.multiplier=vital.params$s.multiplier))
  # movement kernel
  G <- outer(kernel.params$y, kernel.params$y, 
             FUN=gx, p.stat=vital.params$p.stat, 
                     s.stat=vital.params$s.stat, 
                     s.mob=vital.params$s.mob,
             ...) # the dots allow to specify directionality (uni- or bi-) of the movement kernel
  # reproduction
  R <- diag(rx(density=density, 
               r.int=vital.params$r.int, 
               r.slo=vital.params$r.slo, 
               r.multiplier=vital.params$r.multiplier))
  # recruited position (natal dispersal)
  N <- outer(kernel.params$y, kernel.params$y, 
             FUN=gx, p.stat=vital.params$p.stat.nat, 
                     s.stat=vital.params$s.stat.nat, 
                     s.mob=vital.params$s.mob.nat,
                     ...) # the dots allow to specify directionality (uni- or bi-) of the movement
  
  # P kernel
  Pmatrix <- kernel.params$h*as.matrix(G%*%S)
  # F kernel
  Fmatrix <- kernel.params$h*as.matrix(N%*%R)
  
  # Whole kernel
  K<- Pmatrix + Fmatrix
  
  return(K)
  
}
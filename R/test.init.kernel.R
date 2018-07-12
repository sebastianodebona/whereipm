#' Generate initial kernel for visual examination
#' 
#' Builds and plots a full starting kernel, and all of its components (S, G, R, I, P and F kernels), using pre-defined parameters.
#' 
#' @export
#' 
#' @param kernel.params list of parameters to build the kernel. Can be an output of `create.kernel.params`.
#' @param vital.params list of parameters to determine vital rate functions. Can be an output of `create.vital.params`
#' @param init.density Optional vector of initial population densities. If provided, the lenght of the vector must be equal to the number of bins specified in `kernel.params`. Default is a population with 10 individuals at the spatial origin.
#' @details This function is to be used as a preliminary step before running an IPM simulation. It allows to evaluate visually the vital rates kernels and how they come together to provide the first full kernel.
#' 


test.init.kernel<-function(kernel.params, vital.params, init.density){
  
  # deriving missing kernel parameters
  # defining meshpoints (separating bins) and midpoints (central point of each bin)
  b <- kernel.params$L+c(0:kernel.params$n.bins)*(kernel.params$U-kernel.params$L)/kernel.params$n.bins # meshpoints
  y <- 0.5*(b[1:kernel.params$n.bins]+b[2:(kernel.params$n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
  h <- (kernel.params$U-kernel.params$L)/kernel.params$n.bins
  
  
  if(missing(init.density)){
    init.density <- c(10, rep(0, kernel.params$n.bins-1))/h
  } # this sets default initial density if not provided

  
  # survival kernel
  S <- diag(sx(density=init.density, 
               s.int=vital.params$s.int, 
               s.slo=vital.params$s.slo, 
               s.multiplier=vital.params$s.multiplier))
  # movement kernel
  G <- outer(y,y, FUN=gx, p.stat=vital.params$p.stat, 
                      s.stat=vital.params$s.stat, 
                      s.mob=vital.params$s.mob)
  # reproduction
  R <- diag(rx(density=init.density, 
               r.int=vital.params$r.int, 
               r.slo=vital.params$r.slo, 
               r.multiplier=vital.params$r.multiplier))
  # recruited position (natal dispersal)
  N <- outer(y,y, FUN=gx, p.stat=vital.params$p.stat.nat, 
             s.stat=vital.params$s.stat.nat, 
             s.mob=vital.params$s.mob.nat)
  
  # P kernel
  Pmatrix <- h*as.matrix(G%*%S)
  
  # F kernel
  Fmatrix <- h*as.matrix(N%*%R)

  # Whole kernel
  K<- Pmatrix + Fmatrix
  
  
  op <- par()
  par(mfrow=c(2,2))
  
  image(t(S), main='survival')
  image(t(G),main='dispersal')
  image(t(R), main='reproduction')
  image(t(N),main='natal disp')
  
  par(ask=T)
  image(t(Pmatrix),main='P kernel')
  abline(0,1,lty=2)
  image(t(Fmatrix), main='F kernel')
  abline(0,1,lty=2)
  image(t(K), main='full kernel') 
  
  suppressWarnings(par(op))
}  
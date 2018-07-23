#' Generate initial kernel for visual examination
#' 
#' Builds and plots a full starting kernel, and all of its components (S, G, R, I, P and F kernels), using pre-defined parameters.
#' 
#' @export
#' 
#' @param kernel.params list of parameters to build the kernel. \emph{Must} be an output of \code{create.kernel.params}.
#' @param vital.params list of parameters to determine vital rate functions. Can be an output of \code{create.vital.params}.
#' @param init.pop.size Optional vector of initial population size per position. If provided, the lenght of the vector must be equal to the number of bins specified in \code{kernel.params}. Default is a population with 10 individuals at the spatial origin 0.
#' @param ... additional parameters, such as \code{unidirectional} for movement functions (see \code{\link{gx}})
#' @details This function is to be used as a preliminary step before running an IPM simulation. It allows to evaluate visually the vital rates kernels and how they come together to provide the first full kernel. It is very similar in its body to \code{\link{build.kernel}}, but the latter doesn't allow automatic plotting.
#' 


test.init.kernel<-function(kernel.params, vital.params, init.pop, ...){
  
  
  if(missing(init.pop)){  # this sets default initial density if not provided
    init.density <- c(10, rep(0, kernel.params$n.bins-1))/kernel.params$h
  } else {
    init.density <- init.pop/kernel.params$h
  }

  
  # survival kernel
  S <- diag(sx(density=init.density, 
               s.int=vital.params$s.int, 
               s.slo=vital.params$s.slo, 
               s.multiplier=vital.params$s.multiplier))
  # movement kernel
  G <- outer(kernel.params$y,kernel.params$y, 
             FUN=gx, p.stat=vital.params$p.stat, 
                      s.stat=vital.params$s.stat, 
                      s.mob=vital.params$s.mob)
  # reproduction
  R <- diag(rx(density=init.density, 
               r.int=vital.params$r.int, 
               r.slo=vital.params$r.slo, 
               r.multiplier=vital.params$r.multiplier))
  # recruited position (natal dispersal)
  N <- outer(kernel.params$y,kernel.params$y, 
             FUN=gx, p.stat=vital.params$p.stat.nat, 
                     s.stat=vital.params$s.stat.nat, 
                     s.mob=vital.params$s.mob.nat)
  
  # P kernel
  Pmatrix <- kernel.params$h*as.matrix(G%*%S)
  
  # F kernel
  Fmatrix <- kernel.params$h*as.matrix(N%*%R)

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
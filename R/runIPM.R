#' Project Population with IPM
#' 
#' The function runs an IPM to project the spread of the population through space in time
#' 
#' @export
#' 
#' @inheritParams build.kernel
#' 
#' @param init.pop Vector of lenght equal to the number of bins defininf the starting population size (\emph{not density!})
#' @param niter Number of iterations, equivalent to time steps, the population should be projected for. Default is 10 (way too small!)
#' @param store.k Logic value describing whether the kernels at each iteration need to be saved and store in a list
#'

runIPM <- function(kernel.params, vital.params, init.pop, niter=10,
                   sk = "full", rk = "full", gk = "full", nk = "full",
                   store.k=FALSE, ...){ 
  
  # allocating space for population size over time and for spread over time
  pop.size <- matrix(NA, nrow=kernel.params$n.bins, ncol=niter + 1)
  pop.size[,1] <- init.pop
  
  pop.spread <- vector('numeric', niter + 1)
  pop.spread[1] <- NA
  
  #  allocating space for kernels, if specified
  if(store.k){
    K.list <- list()
  }
  
  
  for(i in 1:niter){
    
    # re-building the kernel at each time step based on population density, which affects vital rates
    K<-build.kernel(kernel.params=kernel.params,
                    vital.params=vital.params,
                    pop.size=pop.size[,i],
                    sk = sk, rk = rk, gk = gk, nk = nk)
    
    # storing the kernel if asked to 
    if(store.k){
      K.list[[i]]<-K
    }
    
    # determining next population size by multiplying current kernel and current population size
    next.pop.size <- K%*%pop.size[,i]
    pop.size[,i+1] <- next.pop.size
    pop.spread[i+1] <- if(length(which(next.pop.size > 1)) == 0){NA} else {
      kernel.params$y[max(which(next.pop.size > 1))]}
  }
  
  # # claculating speed, accelaration and population size at the introduction site at the end of the simulation
  # speed <- pop.spread[niter+1] - pop.spread[niter]
  # 
  # # running a linear model to estimate acceleration
  # t <- 1:niter
  # acceleration <- summary(lm(pop.spread[-1] ~ t + I(t^2)))$coefficients['I(t^2)', 'Estimate']
  # 
  
  if(store.k){
    
  return(list(pop.size=pop.size,
              pop.spread=pop.spread,
              y=kernel.params$y,
              n.iterations=niter,
              kernel.params = kernel.params,
              vital.params = vital.params,
              # speed = speed,
              # acceleration = acceleration,
              Kernels = K.list))  
    
  } else {
    
    return(list(pop.size = pop.size,
                pop.spread = pop.spread,
                y = kernel.params$y,
                n.iterations = niter,
                kernel.params = kernel.params,
                vital.params = vital.params
                # speed = speed,
                # acceleration = acceleration
                ))
    
  }
  
}
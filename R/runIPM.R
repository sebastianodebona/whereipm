#' Project Population with IPM
#' 
#' The function runs an IPM to project the spread of the population through space in time
#' 
#' @export
#' 
#' @inheritParams build.kernel
#' 
#' @param niter Number of iterations, equivalent to time steps, the population should be projected for. Default is 10 (way too small!)
#'

runIPM <- function(kernel.params, vital.params, init.pop, niter=10, store.k=FALSE, ...){ 
  
  # allocating space for population size over time and for spread over time
  pop.size <- matrix(NA, nrow=kernel.params$n.bins, ncol=niter + 1)
  pop.size[,1] <- init.pop
  
  pop.spread <- vector('numeric', niter + 1)
  
  #  allocating space for kernels, if specified
  if(store.k){
    K.list <- list()
  }
  
  
  for(i in 1:niter){
    
    # re-building the kernel at each time step based on population density, which affects vital rates
    K<-build.kernel(kernel.params=kernel.params,
                    vital.params=vital.params,
                    pop.size=pop.size[,i])
    
    # storing the kernel if asked to 
    if(store.k){
      K.list[[i]]<-K
    }
    
    # determining next population size by multiplying current kernel and current population size
    next.pop.size <- K%*%pop.size[,i]
    pop.size[,i+1]<-next.pop.size
    pop.spread[i+1] <- max(sum(((kernel.params$y-mean(kernel.params$y))^2)*next.pop.size)/(sum(next.pop.size)-1),0)
  }
  
  
  if(store.k){
    
  return(list(pop.size=pop.size,
              pop.spread=pop.spread,
              y=kernel.params$y,
              n.iterations=niter,
              Kernels <- K.list))  
    
  } else {
    
    return(list(pop.size=pop.size,
                pop.spread=pop.spread,
                y=kernel.params$y,
                n.iterations=niter))
    
  }
  
}
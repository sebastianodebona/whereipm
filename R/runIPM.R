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
  
  pop.spread <- front <- vector('numeric', niter + 1)
  pop.spread[1] <- calculate_q_spread(init.pop, kernel.params)
  front[1] <- if(length(which(init.pop > 1*kernel.params$h)) == 0){NA} else {
    kernel.params$y[max(which(init.pop > 1*kernel.params$h))]}
  
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
    pop.spread[i+1] <- calculate_q_spread(next.pop.size, kernel.params, q = 0.95)
    front[i+1] <- if(length(which(next.pop.size > 1*kernel.params$h)) == 0){NA} else {
                          kernel.params$y[max(which(next.pop.size > 1*kernel.params$h))]}
  }
  
  # claculating speed at each time step, acceleration and lambda
  speed <- (front[2:niter] - front[1:(niter-1)])/1

  # acceleration
  acceleration_avg <- mean((speed[2:(niter-1)] - speed[1:(niter - 2)])/1, na.rm = T)
  acceleration_tot <- (speed[max(which(!is.na(speed)))] - speed[1])/max(which(!is.na(speed)))
  # running a linear model to estimate acceleration
  # I will run the model at two conditions:
  # 1) the final speed is not NA, meaning the final or semifinal spread is NA
  # 2) the values for spread are not all equal (degenerate lm)
  # to check for the latter I am using this approach, found here from hadley wickham
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  if(sum(is.na(speed)) != 0){
    acceleration_est <- NA
  } else if (zero_range(front)){
    acceleration_est <- NA
  } else {
    t <- 1:niter
    acceleration_est <- summary(lm(front[-1] ~ t + I(t^2)))$coefficients['I(t^2)', 'Estimate']
  }
  acceleration = list(acceleration_avg,
                      acceleration_tot,
                      acceleration_est)

  # lambdas
  lambdas <- colSums(pop.size)[2:niter]/colSums(pop.size)[1:(niter-1)]
  lambda_avg = mean(lambdas)
  
  if(store.k){
    
  return(list(pop.size=pop.size,
              pop.spread=pop.spread,
              front = front,
              y=kernel.params$y,
              n.iterations=niter,
              kernel.params = kernel.params,
              vital.params = vital.params,
              speed = speed,
              acceleration = acceleration,
              lambdas = lambdas,
              lambda_avg = lambda_avg,
              Kernels = K.list))  
    
  } else {
    
    return(list(pop.size = pop.size,
                pop.spread = pop.spread,
                front = front,
                y = kernel.params$y,
                n.iterations = niter,
                kernel.params = kernel.params,
                vital.params = vital.params,
                speed = speed,
                acceleration = acceleration,
                lambdas = lambdas,
                lambda_avg = lambda_avg
                ))
    
  }
  
}
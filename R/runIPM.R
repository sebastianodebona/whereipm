#' Project Population with IPM
#' 
#' The function runs an IPM to project the spread of the population through space in time
#' 
#' @export
#' 
#' @inheritParams build.kernel
#' 
#' @param init.pop Vector of lenght equal to the number of bins defining the starting population size (\emph{not density!})
#' @param niter Number of iterations, equivalent to time steps, the population should be projected for. Default is 10 (way too small!)
#' @param store.k Logic value describing whether the kernels at each iteration need to be saved and store in a list
#'

runIPM <- function(kernel.params, vital.params, init.pop, niter=10,
                   sk = "full", rk = "full", gk = "full", nk = "full",
                   store.k=FALSE, ...){ 
  
  # allocating space for population size over time and for spread over time
  pop.size <- matrix(NA, nrow=kernel.params$n.bins, ncol=niter + 1)
  pop.size[,1] <- init.pop
  
  # Spread is calculated three ways, 
  # - pop.spread: space containing 95% of the population
  # - front: most advanced position with at least one individual
  # - pioneer: most advanced position past which there is at least one cumulative individual
  
  pop.spread <- front <- pioneer <- vector('numeric', niter + 1)
  pop.spread[1] <- calculate_q_spread(init.pop, kernel.params)
  front[1] <- if(length(which(init.pop > 1*kernel.params$h)) == 0){NA} else {
    kernel.params$y[max(which(init.pop > 1*kernel.params$h))]}
  pioneer[1] <- calculate_front_int(init.pop, kernel.params, how_many = 1)
  
  
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
    
    pioneer[i+1] <- calculate_front_int(next.pop.size, kernel.params, how_many = 1)
  }
  
  # claculating speed at each time step, acceleration and lambda
  speed_front <- (front[2:(niter + 1)] - front[1:(niter)])/1
  speed_pioneer <- (pioneer[2:(niter + 1)] - pioneer[1:(niter)])/1

  # acceleration
  acceleration_tot_front <- (speed_front[max(which(!is.na(speed_front)))] -
                               speed_front[1])/max(which(!is.na(speed_front)))
  acceleration_tot_pioneer <- (speed_pioneer[max(which(!is.na(speed_pioneer)))] - 
                                 speed_pioneer[1])/max(which(!is.na(speed_pioneer)))
  
  # running a linear model to estimate acceleration
  # I will run the model at two conditions:
  # 1) the final speed is not NA, meaning the final or semifinal spread is NA
  # 2) the values for spread are not all equal (degenerate lm)
  # to check for the latter I am using this approach, found here from hadley wickham
  # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
  # front acceleration
  if(sum(is.na(speed_front)) != 0){
    acceleration_est_front <- NA
  } else if (zero_range(front)){
    acceleration_est_front <- NA
  } else {
    t <- 1:niter
    acceleration_est_front <- summary(lm(front[-1] ~ t + I(t^2)))$coefficients['I(t^2)', 'Estimate']
  }
  # pioneer acceleration
  if(sum(is.na(speed_pioneer)) != 0){
    acceleration_est_pioneer <- NA
  } else if (zero_range(pioneer)){
    acceleration_est_pioneer <- NA
  } else {
    t <- 1:niter
    acceleration_est_pioneer <- summary(lm(pioneer[-1] ~ t + I(t^2)))$coefficients['I(t^2)', 'Estimate']
  }

  
  # lambdas
  lambdas <- colSums(pop.size)[2:(niter + 1)]/colSums(pop.size)[1:(niter)]
  lambda_avg = mean(lambdas)
  
  # more outputs
  final_pop_size = sum(pop.size[,niter+1])
  final_max_density = max(pop.size[,niter+1])/kernel.params$h
  final_front_position = front[niter + 1]
  final_pioneer_position = pioneer[niter + 1]
  avg_front_speed = mean(speed_front, na.rm = T)
  avg_pioneer_speed = mean(speed_pioneer, na.rm = T)
  initial_front_speed = speed_front[1]
  final_front_speed = speed_front[length(speed_front)]
  initial_pioneer_speed = speed_pioneer[1]
  final_pioneer_speed = speed_pioneer[length(speed_pioneer)]
  
  if(store.k){
    
  return(list(pop.size=pop.size,
              pop.spread=pop.spread,
              front = front,
              pioneer = pioneer,
              y=kernel.params$y,
              n.iterations=niter,
              kernel.params = kernel.params,
              vital.params = vital.params,
              speed_front = speed_front,
              speed_pioneer = speed_pioneer,
              acceleration_tot_front = acceleration_tot_front,
              acceleration_tot_pioneer = acceleration_tot_pioneer,
              acceleration_est_front = acceleration_est_front,
              acceleration_est_pioneer = acceleration_est_pioneer,
              final_pop_size = final_pop_size,
              final_max_density = final_max_density,
              final_front_position = final_front_position,
              final_pioneer_position = final_pioneer_position,
              avg_front_speed = avg_front_speed,
              avg_pioneer_speed = avg_pioneer_speed,
              initial_front_speed = initial_front_speed,
              final_front_speed = final_front_speed,
              initial_pioneer_speed = initial_pioneer_speed,
              final_pioneer_speed = final_pioneer_speed,
              lambdas = lambdas,
              lambda_avg = lambda_avg,
              Kernels = K.list))  
    
  } else {
    
    return(list(pop.size=pop.size,
                pop.spread=pop.spread,
                front = front,
                pioneer = pioneer,
                y=kernel.params$y,
                n.iterations=niter,
                kernel.params = kernel.params,
                vital.params = vital.params,
                speed_front = speed_front,
                speed_pioneer = speed_pioneer,
                acceleration_tot_front = acceleration_tot_front,
                acceleration_tot_pioneer = acceleration_tot_pioneer,
                acceleration_est_front = acceleration_est_front,
                acceleration_est_pioneer = acceleration_est_pioneer,
                final_pop_size = final_pop_size,
                final_max_density = final_max_density,
                final_front_position = final_front_position,
                final_pioneer_position = final_pioneer_position,
                avg_front_speed = avg_front_speed,
                avg_pioneer_speed = avg_pioneer_speed,
                initial_front_speed = initial_front_speed,
                final_front_speed = final_front_speed,
                initial_pioneer_speed = initial_pioneer_speed,
                final_pioneer_speed = final_pioneer_speed,
                lambdas = lambdas,
                lambda_avg = lambda_avg
                ))
    
  }
  
}
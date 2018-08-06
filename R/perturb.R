#' Perturbation analyses
#' 
#' Runs an IPM based on a given IPM, perturbing one of the parameters defining the Kernel
#' 
#' @export
#' 
#' @param orig.IPM Non-perturbed IPM result. Must contain population size array and population spread through time. Fitted IPM from \code{\link{runIPM}}.
#' @param epsilon Proportional change to apply to the vital parameter of choice. The default is the inverse of the sensitivity of the Machine to the cube (based on Data-driven modelling of Structured Populations, page 219, Box 7.2)
#' @param par Parameter to perturb.
#' @return Produces a proportional change in the accelaration of spread, speed of spread at the end of the simulation, local population size at the origin at the end of the simulation, and lambda at the end of the simulation.
#' 

perturb <- function(orig.IPM, epsilon, par){
  
  # extracting parameters to be identical to those of original IPM
  kernel.params <- orig.IPM$kernel.params
  vital.params <- orig.IPM$vital.params
  niter <- orig.IPM$n.iterations
  init.pop <- orig.IPM$pop.size[,1]
  
  # modifying the chosen parameter
  vital.params[[par]] <- vital.params[[par]] + epsilon*vital.params[[par]]
  
  #running IPM
  perturb.IPM <- runIPM(kernel.params = kernel.params,
                        vital.params = vital.params,
                        init.pop = init.pop,
                        niter = niter)
  
  # comparing final speed
  delta.speed <- perturb.IPM$speed - orig.IPM$speed
  delta.accelaration <- perturb.IPM$acceleration - orig.IPM$acceleration
  delta.lambda <- sum(perturb.IPM$pop.size[,niter+1])/sum(perturb.IPM$pop.size[,niter]) - sum(orig.IPM$pop.size[,niter+1])/sum(orig.IPM$pop.size[,niter])
  delta.pop <- perturb.IPM$pop.size[kernel.params$n.bins,niter+1] - orig.IPM$pop.size[kernel.params$n.bins,niter+1]
  
  return(list(delta.speed = delta.speed,
              delta.accelaration = delta.accelaration,
              delta.lambda = delta.lambda,
              delta.pop = delta.pop))
  
}
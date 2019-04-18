#' Perturbation analyses
#' 
#' Runs an IPM based on a given IPM, perturbing one of the parameters defining the Kernel
#' 
#' @export
#' 
#' @param orig.IPM Non-perturbed IPM result. Must contain population size array and population spread through time. Fitted IPM from \code{\link{runIPM}}.
#' @param epsilon Proportional change to apply to the vital parameter of choice. The default is the inverse of the sensitivity of the Machine to the cube (based on Data-driven modelling of Structured Populations, page 219, Box 7.2)
#' @param par Parameter to perturb.
#' @param outputs Defines which outputs are compared between original and perturbed IPMs.
#' @return Produced the perturbed IPM and a list of proportional changes in the output variables, as resulting from the perturbation.
#' 

perturb <- function(orig.IPM, epsilon, par, outputs){
  
  # extracting parameters to be identical to those of original IPM
  kernel.params <- orig.IPM$kernel.params
  vital.params <- orig.IPM$vital.params
  niter <- orig.IPM$n.iterations
  init.pop <- orig.IPM$pop.size[,1]
  
  # modifying the chosen parameter
  vital.params[[par]] <- vital.params[[par]] + epsilon*vital.params[[par]]
  
  # running IPM
  perturb.IPM <- runIPM(kernel.params = kernel.params,
                        vital.params = vital.params,
                        init.pop = init.pop,
                        niter = niter)
  
  # creating a list to store relative deltas
  deltas <- list()
  # comparing output parameters
  for(i in outputs){
    deltas[[paste("delta", i, sep = ".")]] <- (perturb.IPM[[i]] - orig.IPM[[i]])/abs(orig.IPM[[i]])
  }
  
  return(list(newIPM = perturb.IPM,
              deltas = deltas))
  
}
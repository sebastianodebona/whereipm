#' Build Kernel
#' 
#' Build a kernel given a population size, kernel parameters and vital rate parameters
#' 
#' @export
#' 
#' @param kernel.params list of parameters to build the kernel (must include L, H and n.bins). See output of \code{create.kernel.params} for an example.
#' @param vital.params list of parameters to determine vital rate functions. Follow format of the output of \code{create.vital.params}
#' @param pop.size Vector of initial population size per position. The lenght of the vector must be equal to the number of bins specified in \code{kernel.params}. There is no default (as opposed to \code{\link{test.init.kernel}}), so population size \emph{must} be specified.
#' @param test Logic value defining whether all matrices (S, R, G, N, P and F) need to be saved for testing. Default if \code{FALSE}.
#' @details As a preliminary step before, run a test with \code{\link{test.init.kernel}} with the same parameters, to check visually if the kernel builds without problems.
#' 

build.kernel <- function(kernel.params, vital.params, pop.size, test=F){
  

  # determining density
  density <- pop.size/kernel.params$h
  
  # determining landscape density
  Land.density <- mean(density[density>=1])
  
  # survival kernel
  S <- diag(sx(density=density,
               s.int=vital.params$s.int, 
               s.slo=vital.params$s.slo, 
               s.multiplier=vital.params$s.multiplier))
  # movement kernel
  combs <- merge(data.frame(y=kernel.params$y), cbind(x=kernel.params$y, density.x=density))
  combs$p <- gx_ln(y = combs[['y']], x = combs[['x']],
                   density.x = combs[['density.x']],
                   Land.density = Land.density,
                   p.stat.int = vital.params$p.stat.int.ad,
                   p.stat.slo = vital.params$p.stat.slo.ad,
                   sigma.stat = vital.params$sigma.stat.ad, 
                   mu.mob.int = vital.params$mu.mob.int.ad,
                   mu.mob.slo = vital.params$mu.mob.slo.ad,
                   sigma.mob = vital.params$sigma.mob.ad)
  G <- matrix(NA, kernel.params$n.bins, kernel.params$n.bins)
  G[1:length(unique(combs$x)),1:length(unique(combs$y))] <- combs$p
  
  # reproduction
  R <- diag(rx(density = density, 
               r.int = vital.params$r.int, 
               r.slo = vital.params$r.slo, 
               r.multiplier = vital.params$r.multiplier,
               f.distribution = ifelse(is.null(vital.params$f.distribution),
                                             'poisson',
                                             vital.params$f.distribution),
               f.int = ifelse(is.null(vital.params$f.int), 
                              NA,
                              vital.params$f.int),
               f.slo = ifelse(is.null(vital.params$f.slo), 
                              NA,
                              vital.params$f.slo)))
  # recruited position (natal dispersal)
  nombs <- merge(data.frame(y=kernel.params$y), cbind(x=kernel.params$y, density.x=density))
  nombs$p <- gx_ln(y = combs[['y']], x = combs[['x']],
                   density.x = combs[['density.x']],
                   Land.density = Land.density,
                   p.stat.int = vital.params$p.stat.int.nat,
                   p.stat.slo = vital.params$p.stat.slo.nat,
                   p.stat.sloLand = vital.params$p.stat.sloLand.nat,
                   p.stat.sloInt = vital.params$p.stat.sloInt.nat,
                   sigma.stat = vital.params$sigma.stat.nat, 
                   mu.mob.int = vital.params$mu.mob.int.nat,
                   mu.mob.slo = vital.params$mu.mob.slo.nat,
                   mu.mob.sloLand = vital.params$mu.mob.sloLand.nat,
                   mu.mob.sloInt = vital.params$mu.mob.sloInt.nat,
                   sigma.mob = vital.params$sigma.mob.nat)
  N <- matrix(NA, kernel.params$n.bins, kernel.params$n.bins)
  N[1:length(unique(nombs$x)),1:length(unique(nombs$y))] <- nombs$p
  
  # P kernel
  Pmatrix <- kernel.params$h*as.matrix(G%*%S)
  # F kernel
  Fmatrix <- kernel.params$h*as.matrix(N%*%R)
  
  # Whole kernel
  K<- Pmatrix + Fmatrix
  
  if(!test){
    return(K)
  } else if(test){
    return(list(S = S, G = G, R = R, N = N, Pmatrix = Pmatrix, Fmatrix = Fmatrix, K = K))
  }
  
}
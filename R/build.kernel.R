#' Build Kernel
#' 
#' Build a kernel given a population size, kernel parameters and vital rate parameters
#' 
#' @export
#' 
#' @param kernel.params list of parameters to build the kernel (must include L, H and n.bins). See output of \code{create.kernel.params} for an example.
#' @param vital.params list of parameters to determine vital rate functions. Follow format of the output of \code{create.vital.params}
#' @param pop.size Vector of initial population size per position. The lenght of the vector must be equal to the number of bins specified in \code{kernel.params}. There is no default (as opposed to \code{\link{test.init.kernel}}), so population size \emph{must} be specified.
#' @param sk String expression defining how the S matrix should be built. \code{"1"} for everyone surviving; \code{"coinflip"} for constant 50% chances of individuals surviving; \{"full"} for the full density-dependent model, using parameters in \code{"vital.params"} and \code{\link{sx}}.
#' @param rk String expression defining how the R matrix should be built. \code{"0"} for no reproduction; \code{"small"} for constant but little reproduction, 1 individual recruited every 10; \code{"full"} for the full density-dependent model, using parameters in \code{"vital.params"} and the function \code{\link{rx}}.
#' @param gk String expression defining how the G (adult dispersal) matrix should be built. \code{"still"} for no dispersal; \code{"gaussian"} for redistribution (a.k.a. dispersal) curve being a gaussian with mean 0 and small variance; \code{"lognormal"} for redistribution kernel following two back-to-back lognormal distributions; \code{"rw"} for a density-dependent redistribution curve following a Radinger-Wolter curve (sum of gaussians), defined in \code{\link{gl_rw}}; \code{"full"} for a density-dependent redistribution curve as a combination of a gaussian and two back-to-back lognormal, defined in \code{\link{gx_ln}}.
#' @param nk String expression defining how the N (natal dispersal) matrix should be built. \code{"still"} for no dispersal; \code{"gaussian"} for redistribution (a.k.a. dispersal) curve being a gaussian with mean 0 and small variance; \code{"lognormal"} for redistribution kernel following two back-to-back lognormal distributions; \code{"rw"} for a density-dependent redistribution curve following a Radinger-Wolter curve (sum of gaussians), defined in \code{\link{gl_rw}}; \code{"full"} for a density-dependent redistribution curve as a combination of a gaussian and two back-to-back lognormal, defined in \code{\link{gx_ln}}.
#' @param test Logic value defining whether all matrices (S, R, G, N, P and F) need to be saved for testing. Default if \code{FALSE}.
#' @details As a preliminary step before, run a test with \code{\link{test.init.kernel}} with the same parameters, to check visually if the kernel builds without problems.
#' 

build.kernel <- function(kernel.params, vital.params, pop.size, test=F){
  

  # determining density
  density <- pop.size/kernel.params$h
  
  # determining landscape density
  Land.density <- mean(density[density>=1])
  
  # Survival kernel
  if(sk == "1"){
    # "everyone survives"
    S <- diag(kernel.params$n.bins)
    
  } else if (sk == "coinflip") {
    # half the population survives
    S <- 0.5*diag(kernel.params$n.bins)
    
  } else if (sk == "full"){
    # full, density-dependent dispersal kernel, based on vital parameters
    S <- diag(sx(density = density,
                 s.int = vital.params$s.int, 
                 s.slo = vital.params$s.slo,
                 s.slo2 = vital.params$s.slo2,
                 s.multiplier = vital.params$s.multiplier))
    }
    
    
  # Movement kernel
  if(gk == "still"){
    # no dispersal
    G <- diag(kernel.params$n.bins)
    
  } else if (gk == "gaussian"){
    # redistribution following a gaussian distribution
    G <- mapply(gx_intgr_normal, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  sd = 3 
                ))
    
  } else if (gk == "lognormal") {
    # redistribution following a gaussian distribution
    G <- mapply(gx_intgr_ln, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  meanlog = 2.5,
                  sdlog = 0.6 
                ))
    
  } else if (gk == "rw"){
    # redistribution following a Radinger-Wolter curve, superimposition of two gaussians with different sd, describing a resident component and a mobile component of the population. Here density-indepenent, with made up values.
    G <- mapply(gx_intgr_rw, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  sd.stat = 2,
                  sd.mob = 10,
                  p.leave = 0.3
                ))
    
  } else if (gk == "full"){
    # density-dependent redistribution curves made of a compoud gaussian + mirrored lognormal, built from vital paramters
    G <- mapply(gx_ln, x = kernel.params$y, density.x = density,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  Land.density = Land.density,
                  p.stat.int = vital.params$p.stat.int.ad,
                  p.stat.slo = vital.params$p.stat.slo.ad,
                  sigma.stat = vital.params$sigma.stat.ad, 
                  mu.mob.int = vital.params$mu.mob.int.ad,
                  mu.mob.slo = vital.params$mu.mob.slo.ad,
                  sigma.mob = vital.params$sigma.mob.ad
                ))
  }
  
  
  

  # Recruitment kernel
  if(rk == "0"){
    # no recruitment
    R <- 0*diag(kernel.params$n.bins)
    
  } else if (rk == "small"){
    # 1 recruited every 1 individuals
    R <- 0.1*diag(kernel.params$n.bins)
    
  } else if (rk == "full"){
    # full density-dependent recruitment based on vital parameters
    R <- diag(rx(density = density, 
                 r.int = vital.params$r.int, 
                 r.slo = vital.params$r.slo, 
                 r.multiplier = vital.params$r.multiplier,
                 f.distribution = vital.params$f.distribution,
                 f.int = vital.params$f.int,
                 f.slo = vital.params$f.slo))
    }
 
  
  # Recruited position kernel (natal dispersal)
  if(nk == "still"){
    # no dispersal
    N <- diag(kernel.params$n.bins)
    
  } else if (nk == "gaussian"){
    # redistribution following a gaussian distribution
    N <- mapply(gx_intgr_normal, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  sd = 2
                ))
    
  } else if (nk == "lognormal") {
    # redistribution following a gaussian distribution
    N <- mapply(gx_intgr_ln, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  meanlog = 1.5,
                  sdlog = 0.5 
                ))
    
  } else if (nk == "rw"){
    # redistribution following a Radinger-Wolter curve, superimposition of two gaussians with different sd, describing a resident component and a mobile component of the population. Here density-indepenent, with made up values.
    N <- mapply(gx_intgr_rw, x = kernel.params$y,
                MoreArgs = list(
                  meshpoints = kernel.params$b,
                  sd.stat = 1,
                  sd.mob = 5,
                  p.leave = 0.6
                ))
    
  } else if (nk == "full"){
    # density-dependent redistribution curves made of a compoud gaussian + mirrored lognormal, built from vital paramters
    N <- mapply(gx_ln, x = kernel.params$y, density.x = density,
                MoreArgs = list(meshpoints = kernel.params$b,
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
                                sigma.mob = vital.params$sigma.mob.nat
                ))
  }
  
  

  # P kernel
  Pmatrix <- as.matrix(G%*%S)
  # F kernel
  Fmatrix <- as.matrix(N%*%R)
  # WARNING, not multiplying by h, since G and N store already probabilities
  
  # Whole kernel
  K <- Pmatrix + Fmatrix
  
  if(!test){
    return(K)
  } else if(test){
    return(list(S = S, G = G, R = R, N = N, Pmatrix = Pmatrix, Fmatrix = Fmatrix, K = K))
  }
  
}
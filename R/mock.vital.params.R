#' Make a list of default parameters for Vital Rates
#' 
#' Creates an object list to store all vital rates parameters, using default values. Useful to test if a specific kernel works and runs.
#' 
#' @export
#' 
#' @param s.int,s.slo Intercept and slope of the linear regression between survival (logit transformed) and density. The intercepts determines the average survival probability at density = 0; the slope defines how strongly density affects survival probability. Defaults are 0 (50% chances of survival) and 0 (no effect of density) respectively.
#' @param r.int,r.slo Intercept and slope of linear regression between number of offpsring produced (log transformed) and density. The intercept defines how many offspring on average an individual produces at density = 0; the slope defines how strongly density affects the number of offspring. Defaults are 2 offspring per individual and no effect of density.
#' @param s.multiplier,r.multiplier Optional arguments to inflate survival and reproduction probability.
#' @param p.stat,p.stat.nat Proportion of adults (and juveniles) representing the static component of the population. 1-p.stat represents the proportion of mobile (i.e. dispersing) individuals. Deafault is 0.90
#' @param s.stat,s.stat.nat Standard deviation of the gaussian curve centred on zero, describing the movement of the static component of adults (and juveniles). Default is 1.
#' @param s.mob,s.mob.nat Standard deviation of the gaussian curve centred on zero, describing the movement of the mobile component (i.e. dispersing individuals) among adults (and juveniles). Default is 1.
#' @param noise Logical value that allows to add random noise to the otherwise fixed values
#' 

mock.vital.params <- function(s.int = 0,
                              s.slo = -0.01,
                              s.multiplier = 1,
                              r.int = log(2),
                              r.slo = -.01,
                              r.multiplier = 1,
                              f.distribution = 'poisson',
                              p.stat.int.ad = 1.5,
                              p.stat.slo.ad = 0,
                              sigma.stat.ad = 0.25,
                              mu.mob.int.ad = 0.15,
                              mu.mob.slo.ad = 0,
                              sigma.mob.ad = 1,
                              p.stat.int.nat = 1,
                              p.stat.slo.nat = 0,
                              sigma.stat.nat = 0.25,
                              mu.mob.int.nat = 0.2,
                              mu.mob.slo.nat = 0.01,
                              sigma.mob.nat = 1,
                              noise = FALSE){
  
  
  if(!noise){
  list(s.int = s.int,
       s.slo = s.slo,
       s.multiplier = s.multiplier,
       r.int = r.int,
       r.slo = r.slo,
       r.multiplier = r.multiplier,
       f.distribution = f.distribution,
       r.multiplier = r.multiplier,
       p.stat.int.ad = p.stat.int.ad,
       p.stat.slo.ad = p.stat.slo.ad,
       sigma.stat.ad = sigma.stat.ad,
       mu.mob.int.ad = mu.mob.int.ad,
       mu.mob.slo.ad = mu.mob.slo.ad,
       sigma.mob.ad = sigma.mob.ad,
       p.stat.int.nat = p.stat.int.nat,
       p.stat.slo.nat = p.stat.slo.nat,
       sigma.stat.nat = sigma.stat.nat,
       mu.mob.int.nat = mu.mob.int.nat,
       mu.mob.slo.nat = mu.mob.slo.nat,
       sigma.mob.nat = sigma.mob.nat)
    
    
  } else {
    list(s.int = s.int + rnorm(1, 0, 0.1),
         s.slo = s.slo + rnorm(1, 0, abs(s.slo)*0.1),
         s.multiplier = s.multiplier,
         r.int = r.int + rnorm(1, 0, 0.1),
         r.slo = r.slo + rnorm(1, 0, abs(r.slo)*0.1),
         r.multiplier = r.multiplier,
         f.distribution = f.distribution,
         r.multiplier = r.multiplier,
         p.stat.int.ad = p.stat.int.ad + rnorm(1, 0, 0.1),
         p.stat.slo.ad = p.stat.slo.ad,
         sigma.stat.ad = sigma.stat.ad,
         mu.mob.int.ad = mu.mob.int.ad,
         mu.mob.slo.ad = mu.mob.slo.ad,
         sigma.mob.ad = sigma.mob.ad,
         p.stat.int.nat = p.stat.int.nat + rnorm(1, 0, 0.05),
         p.stat.slo.nat = p.stat.slo.nat,
         sigma.stat.nat = sigma.stat.nat,
         mu.mob.int.nat = mu.mob.int.nat,
         mu.mob.slo.nat = mu.mob.slo.nat,
         sigma.mob.nat = sigma.mob.nat)
  }
  
}
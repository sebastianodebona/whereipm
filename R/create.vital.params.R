#' Make a list of parameters for Vital Rates
#' 
#' Creates an object list to store all vital rates parameters
#' 
#' @export
#' 
#' @param s.int,s.slo Intercept and slope of the linear regression between survival (logit transformed) and density. The intercepts determines the average survival probability at density = 0; the slope defines how strongly density affects survival probability. Defaults are 0 (50% chances of survival) and 0 (no effect of density) respectively.
#' @param r.int,r.slo Intercept and slope of linear regression between number of offpsring produced (log transformed) and density. The intercept defines how many offspring on average an individual produces at density = 0; the slope defines how strongly density affects the number of offspring. Defaults are 2 offspring per individual and no effect of density.
#' @param s.multiplier,r.multiplier Optional arguments to inflate survival and reproduction probability.
#' @param p.stat,p.stat.nat Proportion of adults (and juveniles) representing the static component of the population. 1-p.stat represents the proportion of mobile (i.e. dispersing) individuals. Deafault is 0.90
#' @param s.stat,s.stat.nat Standard deviation of the gaussian curve centred on zero, describing the movement of the static component of adults (and juveniles). Default is 1.
#' @param s.mob,s.mob.nat Standard deviation of the gaussian curve centred on zero, describing the movement of the mobile component (i.e. dispersing individuals) among adults (and juveniles). Default is 1.
#' 

create.vital.params <- function(s.int=0,
                                s.slo=0,
                                r.int=log(2),
                                r.slo=0,
                                s.multiplier=1,
                                r.multiplier=1,
                                p.stat=0.9,
                                p.stat.nat=0.9,
                                s.stat=1,
                                s.stat.nat=1,
                                s.mob=1,
                                s.mob.nat=1){
  
  
  list(s.int=s.int,
       s.slo=s.slo,
       r.int=r.int,
       r.slo=r.slo,
       s.multiplier=s.multiplier,
       r.multiplier=r.multiplier,
       p.stat=p.stat,
       p.stat.nat=p.stat.nat,
       s.stat=s.stat,
       s.stat.nat=s.stat.nat,
       s.mob=s.mob,
       s.mob.nat=s.mob.nat)
  
}
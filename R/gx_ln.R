#' Movement function using Lognormal kernel
#'
#' Calculates the probability to move from x to y, from a symmetric lognormal probability density function (i.e. a dispersal kernel curve). Dispersal can be density dependent.
#' 
#' @export
#' 
#' @param x Current position (single value).
#' @param y Vector of possible positions at the next time step, towards which the function predicts the probability to disperse to.
#' @param density.x Density for the position x.
#' @param Land.density Density at the landscape level (i.e. average local landscape among occupied bins)
#' @param p.stat.int Intercept for the logistic regression defining the proportion of static individuals (not dispersing, but just performing short distance movements) in the population based on local density (i.e. proportion of dispersers at density 0). If \code{p.stat.slo} = 1, no individuals disperse for long distances when density is close to 0.
#' @param p.stat.slo Slope for the logistic regression describing how the proportion of static individuals changes with density. If \code{p.stat.slo} is large in absolute value, the strength of density dependent dispersal probability is high.
#' @param p.stat.sloLand Effect of landscape density over the proportion of static individuals. If \code{p.stat.sloLand} is large, landscape density has a strong effect on the probability to disperse.
#' @param p.stat.sloInt Magnitude of the interaction between local and Landscape density in determining probability to disperse. If positive, means higher landscape density increase th positive effect of local density. If negative, means higher landscape density lessen the effect of local density.
#' @param sigma.stat Standard deviation of the gaussian curve describing small movements of the static component of the population. This is small, to grant some small movement of individuals. By default this is 0.25, so that ~95% of the population will stay within half unit of space in one time step. The mean of this Gaussian is 0.
#' @param mu.mob.int Intercept for the mean of the Gaussian curve which defines the lognormal distribution describing the movement of the mobile (dispersing) component of the population. If \code{l} is the distance moved, the distribution of the movement distances following a lognormal distribution is \code{log(l) ~ N(mu, sigma)}. \code{mu.mob} is the mean of that Gaussian distribution \code{N}.
#' @param mu.mob.slo Slope for the mean of the Gaussian distribution, defining how that changes with density.
#' @param mu.mob.sloLand Effect of landscape density on the mean of the Gaussian distribution.
#' @param mu.mob.sloInt Interaction between local and landscape density affecting the mean of the Gaussian ditribution.
#' @param sigma.mob Standard deviation of the Gaussian curve which defines the lognormal sistribution describing the movement of the mobile (dispersing) component of the population.
#' @return The function returns a value for the probability of an individual being in y at t1 if it was in x at t0.


gx_ln <- function(y, x, density.x, Land.density,
                  p.stat.int, p.stat.slo, 
                  p.stat.sloLand = 0, p.stat.sloInt = 0, sigma.stat = 0.25, 
                  mu.mob.int, mu.mob.slo,
                  mu.mob.sloLand = 0, mu.mob.sloInt = 0, sigma.mob) {
  
  stopifnot(length(x) == 1)
  stopifnot(length(density.x)==1)
  # this produces a density dependent movement probability.
  # NOTE: the redistribution kernel of individuals originating form position x is only dependent on density at that position and on Landscape density. Also, the redistribution kernel is symmetrical, so what is affected is the distance dispersed (regardless of the direction).
  
  # The redistribution kernel of individuals through time is a composed distribution of three other: 
  #1) a central gaussian distribution with mean 0 (so that it's centered on the starting point) and a small variance.
  #2-3) two lognormal distributions, one to the right and one to the left of the central gaussian distribution.
  
  # Here I define all needed parameters for the redistribution kernel
  # probability to leave 
  #   (defining the proportion of Gaussian over Lognormal distributed values)
  p.leave = plogis(p.stat.int + p.stat.slo*density.x +
                      p.stat.sloLand*Land.density + 
                      p.stat.sloInt*density.x*Land.density)
  # mean of the LogNormal distribution
  meanlog = mu.mob.int + mu.mob.slo*density.x +
    mu.mob.sloLand*Land.density + mu.mob.sloInt*density.x*Land.density
  # sd of the LogNormal distribution
  sdlog = sigma.mob
  # mean of the Normal distribution
  meann = 0
  # sd of the Normal distribution
  sdn = sigma.stat
  

  # defining function for normal distribution
  Normal <- function(i){
    (1/(sqrt(2*pi*(sdn^2))))*exp(-((i-meann)^2/(2*(sdn)^2)))
  }
  
  # defininf function for lognormal distribution
  logNormal <- function(i){
    (1/(i*sdlog*(sqrt(2*pi))))*exp(-((log(i)-meanlog)^2/(2*(sdlog^2))))
  }
  
  D <- function(i){
    y <- ifelse(i==0,
                (1-p.leave)*Normal(i),
                (1-p.leave)*Normal(i) + p.leave*0.5*logNormal(abs(i)))
    # p.leave is multiply by 0.5 since the logNormal are two and mirrored respective to 0
    return(y)
  }
  
  # Making sure D Integrates to 1 (with some tolerance error)
  stopifnot(abs(integrate(D, -Inf, Inf)$value-1) < 1e-3)
  
  # calculating redistribution from x to each value of y
  d <- D(y-x)
  
  return(d)
}

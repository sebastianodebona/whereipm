#' Movement function using Lognormal kernel
#'
#' Calculates the probability to move from x to y, from a symmetric lognormal probability density function (i.e. a dispersal kernel curve). Dispersal can be density dependent.
#' 
#' @export
#' 
#' @param x Current position.
#' @param y Position at the next time step towards which the function predicts the probability to disperse to.
#' @param density.x Densities for the position x.
#' @param p.stat.int Intercept for the logistic regression defining the proportion of static individuals (not dispersing, but just performing short distance movements) in the population based on local density (i.e. proportion of dispersers at density 0). If \code{p.stat.slo} = 1, no individuals disperse for long distances when density is close to 0.
#' @param p.stat.slo Slope for the logistic regression describing how the proportion of static individuals changes with density. If \code{p.stat.slo} is large in absolute value, the strength of density dependent dispersal probability is high.
#' @param sigma.stat Standard deviation of the gaussian curve describing small movements of the static component of the population. This is small, to grant some small movement of individuals. By default this is 0.25, so that ~95% of the population will stay within half unit of space in one time step. The mean of this Gaussian is 0.
#' @param mu.mob.int Intercept for the mean of the Gaussian curve which defines the lognormal distribution describing the movement of the mobile (dispersing) component of the population. If \code{l} is the distance moved, the distribution of the movement distances following a lognormal distribution is \code{log(l) ~ N(mu, sigma)}. \code{mu.mob} is the mean of that Gaussian distribution \code{N}.
#' @param mu.mob.slo Slope for the mean of the Gaussian distribution, defining how that changes with density.
#' @param sigma.mob Standard deviation of the Gaussian curve which defines the lognormal sistribution describing the movement of the mobile (dispersing) component of the population.
#' @return The function returns a value for the probability of an individual being in y at t1 if it was in x at t0.


gx_ln <- function(y, x, density.x, p.stat.int, p.stat.slo,  sigma.stat = 0.25, mu.mob.int, mu.mob.slo, sigma.mob) {
  # this produces for now a density dependent movement probability.
  # The dispersal Kernel follows an overlay of a gaussian distribution and two symmetrical lognormal distributions, centered on 0.
  
  p.leave <- plogis(p.stat.int + p.stat.slo*density.x) # defining what is the p to leave
  mean.lognorm <- mu.mob.int + mu.mob.slo*density.x # and the mean of the lognormal describing movement
  
  d <- (1- p.leave)*dnorm(y, mean = x, sd = sigma.stat) + # gaussian of static individuals
    #defines te probability of a value occurring in y if the gaussian is centered in x. This is equivalent to defining the probability of having a value = |y-x| from a Gaussian centered at 1.
    + ((p.leave)*dlnorm(abs(y-x), meanlog = mean.lognorm, sdlog = sigma.mob))/2 
  # lognormal part. This part of the probability is halved, because if moved, there is half chance to move left as there is to move right. It defines what is the probability of obtaining a movement of |x-y| from a lognormal distribution build from a Gaussian with mean given mean and sdev.
  
  return(d)
}
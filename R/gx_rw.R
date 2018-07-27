#' Movement function more Radinger & Wolter
#'
#' Calculates the probability to move from x to y, from a probability density function (i.e. a dispersal kernel curve) defined as an overalp of two gaussian distributions. Dispersal can be density dependent (not yet implemented).
#' 
#' @export
#' 
#' @param x current position
#' @param y position at the next time step towards which the function predicts the probability to disperse
#' @param p.stat Proportion of static individuals (not dispersing, or short-distance dispersers) in the population. If p.stat = 1, no individuals disperse for long distances and the kernel becomes a gaussian curve with mean 0 and var = s.stat^2
#' @param s.stat Standard deviation of the gaussian curve describing the movement of the static component of the population
#' @param s.mob Standard deviation of the gaussian curve describing the movement of the mobile (dispersing) component of the population
#' @param unidirectional Logic parameter describing whether movement is unidirectional (never moving from x1 to x2 if x2 < x1), or birirectional. _Bidirectional_ is default.
#' @return The function returns a value for the probability of an individual being in y at t1 if it was in x at t0


gx_rw <- function(y, x, p.stat, s.stat, s.mob, unidirectional=FALSE) {
  # this produces for now a density INDEPENDENT movement probability.
  # The dispersal Kernel follows Radinger & Wolter 2014, Fish and Fisheries
  
  if(!unidirectional){
    
    p.stat*(1/(sqrt(2*pi*s.stat^2)))*exp(-((y-x)^2/(2*s.stat^2))) +
      (1-p.stat)*(1/(sqrt(2*pi*s.mob^2)))*exp(-((y-x)^2/(2*s.mob^2)))
    
  } else if(unidirectional){
      
    2*(p.stat*(1/(sqrt(2*pi*s.stat^2)))*exp(-((y-x)^2/(2*s.stat^2))) +
         (1-p.stat)*(1/(sqrt(2*pi*s.mob^2)))*exp(-((y-x)^2/(2*s.mob^2))))
    
    } else {print('Failed')}
  
}
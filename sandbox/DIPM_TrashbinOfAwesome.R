# TRASHBIN of AWESOMENESS....To COME
# Here I am dumping a piece of code where I am trying to introduce density dependence in dispersal decision. This part was removed from the DIPM_Simulation.R script. I cannot quite figure out what I was trying to do, but I was onto something so I want to remove it, but not discard it.
# Figure this out at some point!


# breeding dispersal function
mox <- function(y, x, density) {
  # for the movement function I have a density INdependent probability distribution of movements, that has to be centered each time on the current position
  # the function calculates the probability of moving to y, if the distribution kernel is centred on x
  
  # WARNING!!!
  # I am here implementing a half dispersal kernel, with individuals moving only in one direction. To do so, I force y to be always > x, and double the halfed probability to force the integral to be one.
  # MOREOVER, I want to have a density dependent proportion of static vs mobile individuals. In order to do this, I need to add a modifier to the proportion (ps$disp.p). This cannot be added straightaway. Being p a proportion, I need to transform it with a logit transformation, to create a linear relationship with density first, and then backtransform all. Altoghether, p should become plogis(logit(p)+slope*density). Logit(p) would be the intercept, i.e. the proportion of stationary individuals at density 0).
  ifelse(y<x,
         0,
         2*(plogis(log(ps$disp.p[k]/(1-ps$disp.p[k]))+disp.p.slo*density)*(1/(sqrt(2*pi*ps$disp.stat[k]^2)))*exp(-((y-x)^2/(2*ps$disp.stat[k]^2))) +
              (plogis(log((1-ps$disp.p[k])/ps$disp.p[k])+disp.p.slo*density))*(1/(sqrt(2*pi*ps$disp.mob[k]^2)))*exp(-((y-x)^2/(2*ps$disp.mob[k]^2))))
  )
}
### * Description

# Toy script to play with the idea of stochastic (agent-based) model

### * Setup

library(tidyverse)
library(whereipm)
library(profvis)
library(microbenchmark)

devtools::load_all() # To get the source code in profiling results

### * Functions

### ** plotHistory

#' Plot a history of positions
#'
#' @param xHistory A list with x positions for each time step
#' @param multFactor Numeric, factor used to adjust height for visualization purposes
#'
#' @return The numerical integration for the last step; should be close to one.

plotHistory <- function(xHistory, xMin = -50, xMax = 50, nGrid = 1024, bandwith = 0.1,
                        yMin = 1, yMax = NULL, multFactor = 0.1, reverse = TRUE,
                        col = adjustcolor("lightgrey", alpha.f = 0.8),
                        border = adjustcolor("darkgrey", alpha.f = 0.9)) {
    if (is.null(yMax)) { yMax <- length(xHistory) }
    plot(0, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax),
         las = 1, xlab = "Space", ylab = "Time")
    sequence <- seq_along(xHistory)
    if (reverse) { sequence <- rev(sequence) }
    for (i in sequence) {
        nInds <- length(xHistory[[i]])
        if (nInds == 1 && is.na(xHistory[[i]])) {
            lines(c(xMin, xMax), c(i, i), col = border)
        } else {
            d <- density(xHistory[[i]], bw = bandwidth, n = nGrid, from = xMin, to = xMax)
            polygon(d$x, d$y * multFactor * nInds + i, col = col, border = border)
        }
    }
    # Check
    sum(d$y * (xMax - xMin)/nGrid) # Should be close to one
}

### ** runIPM.vanilla

#' This is a simplified version of the original runIPM() function by Seba.
#' 
#' The function runs an IPM to project the spread of the population through space in time
#' 
#' @export
#' 
#' @inheritParams build.kernel
#' 
#' @param init.pop Vector of lenght equal to the number of bins defininf the starting population size (\emph{not density!})
#' @param niter Number of iterations, equivalent to time steps, the population should be projected for. Default is 10 (way too small!)
#' @param store.k Logic value describing whether the kernels at each iteration need to be saved and store in a list
#'

runIPM.vanilla <- function(kernel.params, vital.params, init.pop, niter=10,
                           sk = "full", rk = "full", gk = "full", nk = "full"){ 
    
    # allocating space for population size over time and for spread over time
    pop.size <- matrix(NA, nrow=kernel.params$n.bins, ncol=niter + 1)
    pop.size[,1] <- init.pop
    
    pop.spread <- vector('numeric', niter + 1)
    pop.spread[1] <- NA
    
    for(i in 1:niter){
        
        # re-building the kernel at each time step based on population density, which affects vital rates
        K<-build.kernel(kernel.params=kernel.params,
                        vital.params=vital.params,
                        pop.size=pop.size[,i],
                        sk = sk, rk = rk, gk = gk, nk = nk)
        
        # determining next population size by multiplying current kernel and current population size
        next.pop.size <- K%*%pop.size[,i]
        pop.size[,i+1] <- next.pop.size
        pop.spread[i+1] <- if(length(which(next.pop.size > 1)) == 0){NA} else {
                                                                           kernel.params$y[max(which(next.pop.size > 1))]}
    }
    
    # calculating speed, accelaration and population size at the introduction site at the end of the simulation
    speed <- pop.spread[niter+1] - pop.spread[niter]

    
    # running a linear model to estimate acceleration
    # I will run the model at two conditions:
    # 1) the final speed is not NA, meaning the final or semifinal spread is NA
    # 2) the values for spread are not all equal (degenerate lm)
    # to check for the latter I am using this approach, found here from hadley wickham
    # https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
    
    if(is.na(speed)){
        acceleration <- NA
    } else if (zero_range(pop.spread)){
        acceleration <- NA
    } else {
        t <- 1:niter
        acceleration <- summary(lm(pop.spread[-1] ~ t + I(t^2)))$coefficients['I(t^2)', 'Estimate']
    }

    # Return
    return(list(pop.size = pop.size,
                pop.spread = pop.spread,
                y = kernel.params$y,
                n.iterations = niter,
                kernel.params = kernel.params,
                vital.params = vital.params
                # speed = speed,
                # acceleration = acceleration
                ))
    
}

### * Main

### ** Prepare data for runs (taken from vignettes/DIPM-030_IPM_Phases.Rmd)

# Using data from the vignettes (so the vignettes should have been run before to generate the data files)
data <- ipmdata
vitalRateModels <- readRDS(file.path(here::here(), "results", "vitalRateModels.rds"))

# kernel parameters
kp <- create.kernel.params(L=-300, U=300, n.bins = 800)
# initial population size
pop0 <- c(rep(0, kp$n.bins/2), rep(10,1), rep(0, kp$n.bins/2-1))

Models <- vitalRateModels %>%
    transmute(
        survival = map(survival, broom::tidy),
        recruit = map(recruit, broom.mixed::tidy),
        adDispersal = map(adDispersal, broom::tidy),
        adDistance = map(adDistance, broom::tidy),
        natDispersal = map(natDispersal, broom::tidy),
        natDistance = map(natDistance, broom::tidy))

# Storing values into vital params list, divided into early, mid and late
# given local density is scaled (but not centered) the estimates for the slope have to be transformed back to the non-scaled value.

# Using a general loop to do so for the three phases:

for(i in 1:dim(vitalRateModels)[1]){
 

  temp.vp <- list(# Survival
                 s.int = Models$survival[[i]]%>%
                   filter(term == "(Intercept)") %>%
                   pull(estimate),
                 s.slo = (Models$survival[[i]]%>%
                   filter(term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 s.slo2 = returnIfExists(Models$survival[[i]]%>%
                   filter(term == "I(s_loc_density^2)") %>%
                   pull(estimate)),
                 s.multiplier = 1,
                 
                # Recruitment
                 r.int = Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "zi",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 r.slo = (Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "zi",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 r.multiplier = 1,
                 f.int = Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "cond",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 f.slo = (Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "cond",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 f.distribution = 'ztpoisson',
                 
                # Adult dispersal
                 p.stat.int.ad = Models$adDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 p.stat.slo.ad = (Models$adDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 sigma.stat.ad = 0.25, 
                 mu.mob.int.ad = Models$adDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 mu.mob.slo.ad = (Models$adDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 sigma.mob.ad = Models$adDistance[[i]]%>%
                   filter(effect == "ran_pars") %>%
                   pull(estimate) %>%
                   sum %>%
                   sqrt,
                 
                # Natal dispersal
                 p.stat.int.nat = Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 p.stat.slo.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 p.stat.sloLand.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_land_density.mid") %>%
                   pull(estimate))/scalingAttr$land_density$scale,
                 p.stat.sloInt.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid:s_land_density.mid") %>%
                   pull(estimate))/(scalingAttr$loc_density$scale*scalingAttr$land_density$scale),
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 mu.mob.slo.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 mu.mob.sloLand.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_land_density.mid") %>%
                   pull(estimate))/scalingAttr$land_density$scale,
                 mu.mob.sloInt.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid:s_land_density.mid") %>%
                   pull(estimate))/(scalingAttr$loc_density$scale*scalingAttr$land_density$scale),
                 sigma.mob.nat = Models$natDistance[[i]]%>%
                   filter(effect == "ran_pars") %>%
                   pull(estimate) %>%
                   sum %>%
                   sqrt)
 
  # replacing empty elements with zeros
  temp.vp <- lapply(temp.vp, function(x) ifelse(length(x)==0,0,x))
  
  # renaming list
  assign(paste("vp", vitalRateModels$phase[i], sep = "."), temp.vp)
}

### ** Profiling

# Profiling runIPM.vanilla
prof <- profvis({
    IPM.early <- runIPM.vanilla(kp, vp.early, init.pop = pop0, niter=12)    
})
# To visualize the results, type: print(prof)
# It looks like a bit more than half the run time is spent in calls to plnorm() and pnorm(), so there is not much room for improvement here.

### ** Stochastic approach

# What we are tracking is x, a vector containing individual positions.

# We could imagine also keeping track of e.g. ages and such, but for now it will do.

### *** Just a random walk

# Nobody dies, individuals just move around randomly.

# Initialize x
nStart <- 30
x <- rnorm(n = nStart, mean = 0, sd = 2)
t <- 1 # Time indices start at 1, not 0, so that we can use them for list indexing to store history.
xHistory <- list() # That's where we are going to store the x vectors at different time points.
xHistory[[t]] <- x

# Let's just assume that individuals move randomly for now
nIters <- 30
for (t in 2:nIters) {
    indMove <- sapply(xHistory[[t-1]], function(pos) rnorm(n = 1, mean = 0, sd = 1))
    xHistory[[t]] <- xHistory[[t-1]] + indMove
}

plotHistory(xHistory, multFactor = 1/30)

### *** Random survival and random walk

# Initialize x
nStart <- 30
x <- rnorm(nStart, 0, 2)
xHistory <- list(x)

# Let's assume that individuals can (1) survive with prob pSurv and then (2) move randomly
nIters <- 30
pSurv <- 0.9
for (t in 2:nIters) {
    indSurvival <- rbinom(n = length(xHistory[[t-1]]), size = 1, prob = pSurv)
    currentX <- xHistory[[t-1]][indSurvival > 0]
    if (length(currentX) > 0) {
        indMove <- sapply(currentX, function(pos) rnorm(n = 1, mean = 0, sd = 1))
        xHistory[[t]] <- currentX + indMove
    } else {
        xHistory[[t]] <- NA
    }
}

plotHistory(xHistory, multFactor = 1/10)

### *** Random reproduction and random walk

# Initialize x
nStart <- 30
x <- rnorm(nStart, 0, 2)
xHistory <- list(x)

# Let's assume that individuals can (1) reproduce with prob pRepr and then (2) move randomly (including babies born at the same turn)
nIters <- 30
pRepr <- 0.1
for (t in 2:nIters) {
    indRepr <- rbinom(n = length(xHistory[[t-1]]), size = 1, prob = pRepr)
    currentX <- c(xHistory[[t-1]], xHistory[[t-1]][indRepr > 0])
    if (length(currentX) > 0) {
        indMove <- sapply(currentX, function(pos) rnorm(n = 1, mean = 0, sd = 1))
        xHistory[[t]] <- currentX + indMove
    } else {
        xHistory[[t]] <- NA
    }
}

plotHistory(xHistory, multFactor = 1/30)

### *** Random death, reproduction and random walk

# Initialize x
nStart <- 30
x <- rnorm(nStart, 0, 2)
xHistory <- list(x)
# Let's assume that individuals can (1) reproduce and die with prob pRepr and pSurv (you can die and reproduce at the same turn, but you cannot be born and die at the same turn), and then (2) move randomly (including babies born at the same turn)
nIters <- 50
pRepr <- 0.1
pSurv <- 1/(1 + pRepr) # This doesn't seem to balance the population exactly
pSurv <- 0.905 # This seems like a sweet spot
for (t in 2:nIters) {
    indRepr <- rbinom(n = length(xHistory[[t-1]]), size = 1, prob = pRepr)
    indSurv <- rbinom(n = length(xHistory[[t-1]]), size = 1, prob = pSurv)
    currentX <- c(xHistory[[t-1]][indSurv > 0], xHistory[[t-1]][indRepr > 0])
    if (length(currentX) > 0) {
        indMove <- sapply(currentX, function(pos) rnorm(n = 1, mean = 0, sd = 1))
        xHistory[[t]] <- currentX + indMove
    } else {
        xHistory[[t]] <- NA
    }
}

plotHistory(xHistory, multFactor = 1/30)

### *** Introducing density dependence for random walk

# Explore how density() handles rectangular kernels
x <- runif(5, min = 0, max = 20)
z <- density(x, width = 2, kernel = "rectangular") # So width is defining how wide is the "spreading" of each point
plot(z$x, z$y * length(x), type = "l")
rug(x, lwd = 3)

# IPM
nStart <- 60
x <- rnorm(nStart, 0, 2)
xHistory <- list(x)
# How far an individual moves is density-dependent
nIters <- 30
kernelWidth <- 5 # Each individual has an effect spread on 5 meters, i.e. contributes 0.2 density along 5 meters
for (t in 2:nIters) {
    # Evaluate local densities using a rectangular kernel
    z <- density(xHistory[[t-1]], width = kernelWidth, kernel = "rectangular",
                 n = 1024)
    xDensities <- approx(z$x, z$y, xout = xHistory[[t-1]])$y * length(xHistory[[t-1]])
    indMoves <- rnorm(n = length(xHistory[[t-1]]), mean = 0, sd = xDensities * 1/2)
    xHistory[[t]] <- xHistory[[t-1]] + indMoves
}

plotHistory(xHistory)

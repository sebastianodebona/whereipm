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

plotHistory <- function(xHistory, xMin = -50, xMax = 50, nGrid = 1024, bandwidth = 0.1,
                        yMin = 1, yMax = NULL, multFactor = 1, reverse = TRUE,
                        col = adjustcolor("lightgrey", alpha.f = 0.8),
                        border = adjustcolor("darkgrey", alpha.f = 0.9),
                        stepShift = 0) {
    if (is.null(yMax)) { yMax <- length(xHistory) }
    plot(0, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax),
         las = 1, xlab = "Space", ylab = "Time")
    sequence <- seq_along(xHistory)
    if (reverse) { sequence <- rev(sequence) }
    for (i in sequence) {
        nInds <- length(xHistory[[i]])
        if (nInds == 0) {
            lines(c(xMin, xMax) + stepShift * i, c(i, i), col = border)
        } else {
            d <- density(xHistory[[i]], bw = bandwidth, n = nGrid, from = xMin, to = xMax)
            polygon(d$x + stepShift * i, d$y * 0.1 * multFactor * nInds + i,
                    col = col, border = border)
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
        xHistory[[t]] <- as.numeric()
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
        xHistory[[t]] <- as.numeric()
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

#' Basic wrapper to run a stochastic model
#'
#' This function accepts four optional functions:
#'
#' - survFun (default: everyone survives): function giving the survival status of an individual at a given turn. The function should be of the form f(localDensities) and return a vector of 0/1 of the same length as the input.
#'
#' - reprFun (default: no one reproduces): function giving the reproduction status of an individual at a given turn, independent of whether it survives or not until the end of the turn. The function should be of the form f(localDensities) and return a vector of 0/1 of the same length as the input.
#'
#' - offspringFun (default: one offspring per individual): function giving the number of offspring produced per individual. In the current implementation of runStochasticModel(), this function is called only on individuals with 1 in the output of reprFun. The function should be of the form f(localDensities) and return a vector of (non-strictly) positive integers of the same length as the input.
#'
#' - adultDispFun (default: nobody moves): function giving the disperser status of an adult individual at a given turn. Dispersal decisions are taken at the beginning of the turn, i.e. using the initial densities before anybody dies. The function should be of the form f(localDensities) and return a vector of 0/1 of the same length as the input.
#'
#' - adultDistFun (default: Gaussian kernel following a Z distribution): function giving the distance moved for a disperser adult. The function should be of the form f(localDensities) and return a vector of numeric values of the same length as the input. Those values will be randomly multiplied by -1 or 1 by runStochasticModel(), so it is ok to return only positive values (i.e. the returned kernel is automatically mirrored by runStochasticModel()).
#'
#' @param x0 Numeric vector of individual starting positions
#' @param presenceWidth Numeric value, width of the segment over which the presence of an individual will increase density for density-dependent effects. For example, if set to 5, an individual will increase density by 0.2 over a 5-unit wide segment centered on itself.
#' @param nIters Integer, number of iterations to run. The starting condition x0 is returned as the first element of the simulation, so the total length of a simulation is nIters+1.
#' @param maxInds Integer, maximum number of individuals allowed per time step. The function will raise an error if this number+1 is reached.
#' @return A list of numeric vectors. Each element contains the positions of individuals alive at a given iteration.
#'
#' @examples
#' # Boring: no death, no birth, no move
#' x0 <- runif(n = 50, min = -20, max = 20)
#' sim <- runStochasticModel(x0, 5, nIters = 48)
#' plotHistory(sim)
#'
#' # Survival: constant pSurv
#' survFun <- function(localDensities) {
#'   pSurv <- 0.9
#'   return(rbinom(length(localDensities), size = 1, prob = pSurv))
#' }
#' sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun)
#' plotHistory(sim)
#'
#' # Survival: density-dependent (higher mortality at high density)
#' survFun <- function(localDensities) {
#'   pSurv <- 1/(1 + exp(-8 + 3*localDensities))
#'   return(rbinom(length(localDensities), size = 1, prob = pSurv))
#' }
#' x0 <- rnorm(n = 50, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun)
#' plotHistory(sim)
#'
#' # Survival as above plus everybody has a constant chance of reproducing
#' reprFun <- function(localDensities) {
#'   pRepr <- 0.2
#'   return(rbinom(length(localDensities), size = 1, prob = pRepr))
#' }
#' x0 <- rnorm(n = 50, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun, reprFun = reprFun)
#' plotHistory(sim, multFactor = 0.5, stepShift = 1/2)
#'
#' # Survival and reproduction as above, but instead of default (1 offspring per parent) nOffspring depends on density (less offspring when density increases)
#' offspringFun <- function(localDensities) {
#'   lambdas <- exp(-0.5 - 0.5*localDensities)
#'   return(rpois(length(localDensities), lambda = lambdas))
#' }
#' x0 <- rnorm(n = 50, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun, reprFun = reprFun,
#'          offspringFun = offspringFun)
#' plotHistory(sim, multFactor = 0.5, stepShift = 1/2)
#'
#' # No death, no birth, random walk independent of density
#' adultDispFun <- function(localDensities) {
#'   pDisp <- 0.5
#'   return(rbinom(length(localDensities), size = 1, prob = pDisp))
#' }
#' x0 <- rnorm(n = 5, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, adultDispFun = adultDispFun)
#' plotHistory(sim, multFactor = 0.5)
#'
#' # Same as above, but dispersal probability is higher when density increases
#' adultDispFun <- function(localDensities) {
#'   pDisp <- 1/(1 + exp(4 - 2*localDensities))
#'   return(rbinom(length(localDensities), size = 1, prob = pDisp))
#' }
#' x0 <- rnorm(n = 50, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, adultDispFun = adultDispFun)
#' plotHistory(sim, multFactor = 0.5)
#'
#' # Same as above, but dispersal distance is higher when density increases
#' adultDistFun <- function(localDensities) {
#'   return(rnorm(length(localDensities), mean = 0, sd = localDensities))
#' }
#' x0 <- rnorm(n = 50, mean = 0, sd = 5)
#' sim <- runStochasticModel(x0, 5, nIters = 48, adultDispFun = adultDispFun,
#'          adultDistFun = adultDistFun)
#' plotHistory(sim, multFactor = 0.5)
#'
#' # Putting everything together
#' x0 <- seq(-5, 5, length.out = 30) # 30 individuals evenly spaced over 10 meters
#' sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun, reprFun = reprFun,
#'          offspringFun = offspringFun, adultDispFun = adultDispFun,
#'          adultDistFun = adultDistFun)
#' plotHistory(sim, multFactor = 0.5)
#'
#' # Running many simulations
#' nSims <- 100
#' mySims <- lapply(seq_len(nSims), function(i) {
#'   sim <- runStochasticModel(x0, 5, nIters = 48, survFun = survFun, reprFun = reprFun,
#'            offspringFun = offspringFun, adultDispFun = adultDispFun,
#'            adultDistFun = adultDistFun)
#' })
#' # Calculate average history
#' avgHist <- averageHistory(mySims)
#' plotAvgHistory(avgHist)
#' 

runStochasticModel <- function(x0, presenceWidth, nIters, maxInds = 1e3, 
                               survFun = NULL,
                               reprFun = NULL, offspringFun = NULL,
                               adultDispFun = NULL, adultDistFun = NULL) {
    N_GRID <- 1024 # Used for density estimation; using a power of 2 is nice for the FFT.
    # Parse vital functions
    if (is.null(survFun)) {
        survFun <- function(localDensities) {
            # Everybody survives
            return(rep(1, length(localDensities)))
        }
    }
    if (is.null(reprFun)) {
        reprFun <- function(localDensities) {
            # Nobody reproduces
            return(rep(0, length(localDensities)))
        }
    }
    if (is.null(offspringFun)) {
        offspringFun <- function(localDensities) {
            # One offspring per individual
            return(rep(1, length(localDensities)))
        }
    }
    if (is.null(adultDispFun)) {
        adultDispFun <- function(localDensities) {
            # Nobody moves
            return(rep(0, length(localDensities)))
        }
    }
    if (is.null(adultDistFun)) {
        adultDistFun <- function(localDensities) {
            # Gaussain kernel using Z distribution
            return(rnorm(length(localDensities), mean = 0, sd = 1))
        }
    }
    # Run the model
    xHistory <- list(x0)
    for (i in seq_len(nIters)) {
        prevX <- xHistory[[i]]
        if (length(prevX) > maxInds) {
            stop("Maximum number of individuals reached (", length(prevX),
                 " reached, max. allowed is ", maxInds, ").")
        }
        if (length(prevX) > 0) {
            # Calculate local densities at the beginning of the turn
            # (a) Calculate density profile
            z <- density(prevX, width = presenceWidth, kernel = "rectangular", n = 1024)
            # (b) Interpolate the profile at each individual location
            prevLocalDensities <- approx(z$x, z$y, xout = prevX)$y * length(prevX)
            # Stochastic survival
            survAtX <- survFun(prevLocalDensities)
            curX <- prevX[survAtX == 1]
            # Stochastic reproduction (everybody alive at the beginning of the turn gets a chance)
            reprAtX <- reprFun(prevLocalDensities)
            # Stochastic number of offspring for reproducing individuals
            nJuveniles <- offspringFun(prevLocalDensities[reprAtX == 1])
            if (any(is.na(nJuveniles))) {
                stop("\"NA\" number of juveniles returned from offspringFun().\n",
                     "Are you sure your offspringFun() can handle the density values reached?\n",
                     "Here is a summary of the densities reached:\n",
                     paste(capture.output(summary(prevLocalDensities)), collapse = "\n"))
            }
            if ((sum(nJuveniles)+length(curX)) > maxInds) {
                stop("Maximum number of individuals reached (", (sum(nJuveniles)+length(curX)),
                     " reached, max. allowed is ", maxInds, ").")

            }
            newRecruits <- rep(prevX[reprAtX == 1], times = nJuveniles)
            # Stochastic adult dispersal
            ## Here we assume the dispersal decision was made at the beginning of the turn, i.e. the local densities used are the ones at the beginning of the turn
            disperserAtX <- adultDispFun(prevLocalDensities[survAtX == 1])
            directionMoved <- sample(c(-1, 1), size = sum(disperserAtX > 0), replace = TRUE)
            distanceMoved <- adultDistFun(prevLocalDensities[survAtX == 1][disperserAtX == 1])
            curX[disperserAtX == 1] <- curX[disperserAtX == 1] + distanceMoved
            # Assemble the population
            curX <- c(curX, newRecruits)
        } else {
            # Nobody was alive at the previous iteration
            curX <- as.numeric()
        }
        # Store history
        xHistory[[i+1]] <- curX
    }
    # Return
    return(xHistory)
}

#' Average histories from several simulations
#'
#' @param sims A list containing outputs from runStochasticModel(). All outputs should have the same number of iterations.
#' @param xMin,xMax,bandwidth Passed to density()
#'
#' @return A list containing one density profile for each time step
#'

averageHistory <- function(sims, xMin = -50, xMax = 50, bandwidth = 0.2) {
    N_GRID <- 1024
    nTimesteps <- length(sims[[1]])
    densitiesPerTimestep <- lapply(seq_len(nTimesteps), function(i) {
        densities <- lapply(1:length(sims), function(j) {
            density(sims[[j]][[i]], bw = bandwidth, n = N_GRID, from = xMin, to = xMax)
        })
        x <- densities[[1]]$x
        y <- do.call(rbind, lapply(densities, "[[", "y"))
        y <- colSums(y) / length(sims)
        return(list(x = x, y = y))
    })
    return(densitiesPerTimestep)
}

#' Plot average history
#'
#' @param avgHistory Output from averageHistory()
#'

plotAvgHistory <- function(avgHistory, multFactor = 1, reverse = TRUE,
                           col = adjustcolor("lightgrey", alpha.f = 0.8),
                           border = adjustcolor("darkgrey", alpha.f = 0.9),
                           stepShift = 0) {
    xMin <- min(avgHistory[[1]]$x)
    xMax <- max(avgHistory[[1]]$x)
    yMin <- 0
    yMax <- length(avgHistory) + 1
    plot(0, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax),
         las = 1, xlab = "Space", ylab = "Time")
    sequence <- seq_along(avgHistory)
    if (reverse) { sequence <- rev(sequence) }
    for (i in sequence) {
        polygon(avgHistory[[i]]$x + stepShift * i,
                avgHistory[[i]]$y * 5 * multFactor + i,
                col = col, border = border)
    }
}


plotHistory <- function(xHistory, xMin = -50, xMax = 50, nGrid = 1024, bandwidth = 0.1,
                        yMin = 1, yMax = NULL, multFactor = 1, reverse = TRUE,
                        col = adjustcolor("lightgrey", alpha.f = 0.8),
                        border = adjustcolor("darkgrey", alpha.f = 0.9),
                        stepShift = 0) {
    if (is.null(yMax)) { yMax <- length(xHistory) }
    plot(0, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax),
         las = 1, xlab = "Space", ylab = "Time")
    sequence <- seq_along(xHistory)
    if (reverse) { sequence <- rev(sequence) }
    for (i in sequence) {
        nInds <- length(xHistory[[i]])
        if (nInds == 0) {
            lines(c(xMin, xMax) + stepShift * i, c(i, i), col = border)
        } else {
            d <- density(xHistory[[i]], bw = bandwidth, n = nGrid, from = xMin, to = xMax)
            polygon(d$x + stepShift * i, d$y * 0.1 * multFactor * nInds + i,
                    col = col, border = border)
        }
    }
    # Check
    sum(d$y * (xMax - xMin)/nGrid) # Should be close to one
}

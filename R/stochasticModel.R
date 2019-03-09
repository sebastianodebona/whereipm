### * Description

# Functions related to stochastic model running and plotting

### * runStochasticModel()

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
#' @export

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

### * averageHistory()

#' Average histories from several simulations
#'
#' @param sims A list containing outputs from runStochasticModel(). All outputs should have the same number of iterations.
#' @param xMin,xMax,bandwidth Passed to density()
#'
#' @return A list containing one density profile for each time step
#'
#' @export

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

### * plotHistory()

#' Plot a history of positions from a single simulation
#'
#' @param xHistory Output from runStochasticModel()
#' @param multFactor Numeric, factor used to adjust height for visualization purposes
#'
#' @return The numerical integration for the last step; should be close to one.
#'
#' @export

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

### * plotAvgHistory()

#' Plot average history
#'
#' @param avgHistory Output from averageHistory()
#'
#' @export

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

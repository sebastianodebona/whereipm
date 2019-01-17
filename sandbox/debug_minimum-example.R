### * Description

# Script to implement a minimal example reproducing some issues with the IPM.

### * Simplified functions

### ** runIPM

runIPM <- function(kernel.params, init.pop, niter) {
    # allocating space for population size over time and for spread over time
    pop.size <- matrix(NA, nrow=kernel.params$n.bins, ncol=niter + 1)
    pop.size[,1] <- init.pop
    pop.spread <- vector('numeric', niter + 1)
    pop.spread[1] <- NA
    # Loop
    for(i in 1:niter){
        # re-building the kernel at each time step based on population density,
        # which affects vital rates
        K <- build.kernel(kernel.params=kernel.params,
                          pop.size=pop.size[,i])
        # determining next population size by multiplying current kernel and
        # current population size
        next.pop.size <- K%*%pop.size[,i]
        pop.size[,i+1] <- next.pop.size
    }
    # Return
    return(pop.size)
}

### ** makeInitPop

#' Create an initial population vector
#'
#' Given a kp (with U, L and n.bins), a number of individuals to introduce, and
#' a width of the introduction area, this function calculates how many bins are
#' covered by the introduction area, rounds this up if needed, and distributes
#' evenly the introduced individuals on it. This makes it easier to change the
#' number of bins of a kp, but to keep the same spatial distribution of
#' introduced individuals.
#'
#' @param kp Output from create.kernel.params()
#' @param nIntro Total number of individuals to introduce
#' @param widthIntro Width of the introduction section, in the same unit as L
#'     and U in kp
#'
#' @return An initial pop vector of length kp$n.bins, such that the sum of it
#'     is equal to nIntro
#'
#' @examples
#' makeInitPop(create.kernel.params(-10, 10, 20), nIntro = 10, widthIntro = 2)
#' makeInitPop(create.kernel.params(-10, 10, 40), nIntro = 10, widthIntro = 2)
#' # The number of bins over which the individuals are distributed is not always
#' # matching the width what is asked for because there is some rounding
#' # happening:
#' makeInitPop(create.kernel.params(-10, 10, 10), nIntro = 10, widthIntro = 2)

makeInitPop <- function(kp, nIntro, widthIntro) {
    widthKp <- kp$U - kp$L
    stopifnot(widthIntro <= widthKp)
    borderWidth <- (widthKp-widthIntro)/2
    nBinPerBorder <- floor(borderWidth/kp$h)
    indicesIntroduction <- (nBinPerBorder+1):(kp$n.bins-nBinPerBorder)
    stopifnot(length(indicesIntroduction)>0)
    introPerCell <- nIntro/length(indicesIntroduction)
    init.pop <- rep(0, kp$n.bins)
    init.pop[indicesIntroduction] <- introPerCell
    stopifnot(abs(sum(init.pop) - nIntro) < 2e-16)
    return(init.pop)
}

### ** normalizeLeftStochasticMatrix

#' Normalize a left stochastic matrix (i.e. each column sums up to 1)
#'
#' A transition matrix normalized this way will keep the number of individuals
#' constant before and after matrix multiplication with the population column
#' vector. If you want to allow for individuals to "fall over the edge", you
#' should normalize in a different way.
#'
#' @param x Square matrix with non-negative real coefficients
#'
#' @return Normalized matrix (coefficients in each column sum up to one)

normalizeLeftStochasticMatrix <- function(x) {
    x.cols.sums <- apply(x, 2, sum)
    x.norm <- t(t(x) / x.cols.sums)
    stopifnot(max(abs(apply(x.norm, 2, sum)-1)) < 2e-16)
    return(x.norm)
}

### ** gaussianKernel

#' Build a gaussian kernel
#'
#' Return a non-normalized Gaussian distribution, in a vector of size 2*n.bins,
#' centered on the middle of the vector, with an sd of 10% the number of bins.
#'
#' examples
#' # Test for independence between n.bins and kernel shape
#' plot(gaussianKernel(n.bins = 10), type = "l")
#' plot(gaussianKernel(n.bins = 50), type = "l")
#' plot(gaussianKernel(n.bins = 100), type = "l")

gaussianKernel <- function(n.bins) {
    bin.width <- 1 / n.bins
    sd <- 0.1 # in relative terms
    meshpoints <- seq(-1, 1, length.out = 2*n.bins + 1)
    midpoints <- meshpoints[1:(length(meshpoints)-1)] + diff(meshpoints)/2
    k <- dnorm(midpoints, sd = sd)
    return(k)
}

### ** build.kernel

build.kernel <- function(kernel.params, pop.size){
    n.bins <- kernel.params$n.bins
    # determining density
    density <- pop.size/kernel.params$h
    # survival
    S <- diag(n.bins)
    # adult dispersal
    gaussianTemplate <- gaussianKernel(n.bins)
    G.cols <- lapply(1:n.bins, function(i) {
        gaussianTemplate[(n.bins+2-i):(2*n.bins+1-i)]
    })
    G <- do.call(cbind, G.cols) %>%
        normalizeLeftStochasticMatrix()
    # reproduction
    R <- 0*diag(n.bins)
    # natal dispersal
    N <- diag(n.bins) %>%
        normalizeLeftStochasticMatrix()
    # P kernel
    Pmatrix <- as.matrix(G%*%S)
    # F kernel
    Fmatrix <- as.matrix(N%*%R)
    # Whole kernel
    K <- Pmatrix + Fmatrix
 }

### * Run

library(whereipm)
library(tidyverse)
library(txtplot)

# kernel parameters
kps <- list("kp200" = create.kernel.params(L=-200, U=200, n.bins = 200),
           "kp400" = create.kernel.params(L=-200, U=200, n.bins = 400),
           "kp800" = create.kernel.params(L=-200, U=200, n.bins = 800))
# IPMs
IPMs <- lapply(kps, function(x) {
    runIPM(kernel.params = x,
           # I build init pop size on the fly
           init.pop = makeInitPop(kp = x, nIntro = 10, widthIntro = 4),
           niter = 12)
})
# Pop size versus time (txtplot and regular plot)
par(mfrow = c(2, 2))
ylim <- c(0, 12)
for (i in seq_along(IPMs)) {
    IPM <- IPMs[[i]]
    time <- seq_len(ncol(IPM))
    totalPop <- apply(IPM, 2, sum)
    txtplot(time, totalPop, ylim = ylim)
    plot(time, totalPop, main = names(kps)[i], las = 1,
         type = "o", pch = 19, ylim = ylim)
}

# Checking that the density profiles are the same, whatever the number of bins
# (note that here the G normalization ensures that nobody is falling over the
# edge, so the population always remains at 10).
densities <- lapply(seq_along(IPMs), function(i) {
    return(IPMs[[i]]/kps[[i]]$h)
})
par(mfrow = c(2, 2))
ylim <- c(0, 0.1)
cols <- viridis::viridis(12)
for (i in seq_along(IPMs)) {
    density <- densities[[i]]
    plot(0, type = "n", xlab = "", ylab = "Density",
         xlim = extendrange(1:nrow(density)), ylim = ylim,
         main = names(kps)[i])
    for (t in 1:12) {
        lines(seq_len(nrow(density)), density[, t], col = cols[t])
    }
}

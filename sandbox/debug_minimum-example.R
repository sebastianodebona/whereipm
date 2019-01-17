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
    stopifnot(all(apply(x.norm, 2, sum) == 1))
    return(x.norm)
}

### ** build.kernel

build.kernel <- function(kernel.params, pop.size){
    n.bins <- kernel.params$n.bins
    # determining density
    density <- pop.size/kernel.params$h
    # survival
    S <- diag(n.bins)
    # adult dispersal
    G <- diag(n.bins) %>%
        normalizeLeftStochasticMatrix()
    # reproduction
    R <- 0 * diag(n.bins)
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
           init.pop = c(rep(0, x$n.bins/2), rep(10,1), rep(0, x$n.bins/2-1)),
           niter = 12)
})
# Pop size versus time
for (IPM in IPMs) {
    time <- seq_len(ncol(IPM))
    totalPop <- apply(IPM, 2, sum)
    txtplot(time, totalPop)
}

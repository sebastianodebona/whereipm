### * Description

# This file defines functions to calculate probabilities of dispersal along a
# discretized unidimensional space, from a starting position.

### * gx_intr_normal

#' Calculate the probabilities of dispersal over a unidimensional grid (normal distribution)
#'
#' Given a unidimensional grid and a starting coordinate, calculate the
#' probability of dispersal from this cell over all the cells of the
#' grid. Probabilities are calculated rigorously, and no further adjustment
#' using e.g. cell width is necessary.
#'
#' @param meshpoints Vector of meshpoints defining the coordinates of the
#'     borders of the grid cells (should be ordered from smaller to larger)
#' @param x A single numerical value defining the starting coordinate for
#'     dispersal (same coordinate system as for "meshpoints")
#' @param sd Standard deviation of the normal distribution describing the
#'     dispersal around the starting coordinate (same unit as for meshpoints
#'     and x)
#'
#' @return A vector of length (length(meshpoints)-1), giving the probability of
#'     dispersing to each grid cell defined by the meshpoints
#'
#' @examples
#' par(mfrow = c(2, 3))
#' ylim <- c(0, 0.5)
#' # Coarse grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 20)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_normal(meshpoints, x = -200, sd = 10)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, sd = 10")
#' probs <- gx_intgr_normal(meshpoints, x = -50, sd = 10)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "20 bins, sd = 10")
#' probs <- gx_intgr_normal(meshpoints, x = 20, sd = 50)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "20 bins, sd = 50")
#' # Finer grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 50)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_normal(meshpoints, x = -200, sd = 10)
#' stickplot(midpoints, probs, ylim = ylim, main = "50 bins, sd = 10")
#' probs <- gx_intgr_normal(meshpoints, x = -50, sd = 10)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "50 bins, sd = 10")
#' probs <- gx_intgr_normal(meshpoints, x = 20, sd = 50)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "50 bins, sd = 50")
#'
#' @export

gx_intgr_normal <- function(x, meshpoints, sd) {
    # Calculate CDF at each meshpoint
    cdfAtMesh <- pnorm(meshpoints - x, mean = 0, sd = sd)
    # Calculate the probability of dispersing between each pair of meshpoints
    probs <- diff(cdfAtMesh)
    stopifnot(sum(probs) <= 1)
    # Return
    return(probs)
}

### * gx_intr_rw

#' Calculate the probabilities of dispersal over a unidimensional grid (Radinger-Wolter distribution)
#'
#' Similar to \code{\link{gx_intr_normal}}, but uses a Radinger-Wolter
#' distribution instead of a normal distribution.
#'
#' @param meshpoints Vector of meshpoints defining the coordinates of the
#'     borders of the grid cells (should be ordered from smaller to larger)
#' @param x A single numerical value defining the starting coordinate for
#'     dispersal (same coordinate system as for "meshpoints")
#' @param sd1 Standard deviation of the first normal distribution describing
#'     the dispersal around the starting coordinate (same unit as for
#'     meshpoints and x)
#' @param sd2 Standard deviation of the second normal distribution describing
#'     the dispersal around the starting coordinate (same unit as for
#'     meshpoints and x)
#' @param p Probability of dispersing according to the first normal
#'     distribution (sd=sd1). Conversely, the probability of dispersing
#'     according to the second normal distribution will be (1-p).
#'
#' @return A vector of length (length(meshpoints)-1), giving the probability of
#'     dispersing to each grid cell defined by the meshpoints
#'
#' @examples
#' par(mfrow = c(2, 3))
#' # Coarse grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 20)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_rw(meshpoints, x = -200, sd1 = 10, sd2 = 50, p = 0.8)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, p = 0.8")
#' probs <- gx_intgr_rw(meshpoints, x = -50, sd1 = 10, sd2 = 50, p = 0.5)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, p = 0.5")
#' probs <- gx_intgr_rw(meshpoints, x = 20, sd1 = 10, sd2 = 50, p = 0.2)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, p = 0.2")
#' # Finer grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 50)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_rw(meshpoints, x = -200, sd1 = 10, sd2 = 50, p = 0.8)
#' stickplot(midpoints, probs, ylim = ylim, main = "50 bins, p = 0.8")
#' probs <- gx_intgr_rw(meshpoints, x = -50, sd1 = 10, sd2 = 50, p = 0.5)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "50 bins, p = 0.5")
#' probs <- gx_intgr_rw(meshpoints, x = 20, sd1 = 10, sd2 = 50, p = 0.2)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "50 bins, p = 0.2")
#' 
#' @export

gx_intgr_rw <- function(meshpoints, x, sd.stat, sd.mob, p.leave) {
    # Calculate CDF at each meshpoint
    cdfAtMesh <- ((1-p.leave) * pnorm(meshpoints - x, mean = 0, sd = sd1) +
                  (p.leave) * pnorm(meshpoints - x, mean = 0, sd = sd2))
    # Calculate the probability of dispersing between each pair of meshpoints
    probs <- diff(cdfAtMesh)
    stopifnot(sum(probs) <= 1)
    # Return
    return(probs)
}

### * gx_intgr_ln

#' Calculate the probabilities of dispersal over a unidimensional grid (log-normal distribution)
#'
#' Similar to \code{\link{gx_intr_normal}}, but uses a (mirrored) log-normal
#' distribution instead of a normal distribution.
#'
#' @param meshpoints Vector of meshpoints defining the coordinates of the
#'     borders of the grid cells (should be ordered from smaller to larger)
#' @param x A single numerical value defining the starting coordinate for
#'     dispersal (same coordinate system as for "meshpoints")
#' @param meanlog Mean (on the log scale) of the distribution describing the
#'     dispersal on each side of the starting coordinate (same unit as for
#'     meshpoints and x)
#' @param sdlog Standard deviation (on the log scale) of the distribution
#'     describing the dispersal on each side of the starting coordinate (same
#'     unit as for meshpoints and x)
#' 
#' @return A vector of length (length(meshpoints)-1), giving the probability of
#'     dispersing to each grid cell defined by the meshpoints
#'
#' @examples
#' # How to calculate the CDF for a log-normal distribution mirrored on both
#' # sides of x=0?
#' x <- seq(-5, 5, length.out = 256)
#' plot(x, dlnorm(x) + dlnorm(-x), type = "l")
#' # wrong
#' plot(x, 1/2 * (plnorm(x) + plnorm(-x)), type = "l")
#' # right
#' plot(x, 1/2 * (plnorm(x) + plnorm(-x, lower.tail = FALSE)), type = "l")
#'
#' # Test gx_intgr_ln()
#' par(mfrow = c(2, 3))
#' # Coarse grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 20)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_ln(meshpoints, x = -200, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, meanlog = 3.5, sdlog = 0.6")
#' probs <- gx_intgr_ln(meshpoints, x = -50, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "20 bins")
#' probs <- gx_intgr_ln(meshpoints, x = 20, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "20 bins")
#' # Finer grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 50)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_ln(meshpoints, x = -200, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, main = "50 bins, meanlog = 3.5, sdlog = 0.6")
#' probs <- gx_intgr_ln(meshpoints, x = -50, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "50 bins")
#' probs <- gx_intgr_ln(meshpoints, x = 20, meanlog = 3.5, sdlog = 0.6)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "50 bins")
#' 
#' @export

gx_intgr_ln <- function(meshpoints, x, meanlog, sdlog) {
    # Calculate CDF at each meshpoint
    cdfAtMesh <- 1/2*(plnorm(meshpoints - x, meanlog = meanlog, sdlog = sdlog) +
                      plnorm((x - meshpoints), meanlog = meanlog, sdlog = sdlog,
                             lower.tail = FALSE))
    # Calculate the probability of dispersing between each pair of meshpoints
    probs <- diff(cdfAtMesh)
    stopifnot(sum(probs) <= 1)
    # Return
    return(probs)
}

### * gx_intgr_nln

#' Calculate the probabilities of dispersal over a unidimensional grid (log-normal distribution)
#'
#' Similar to \code{\link{gx_intr_normal}}, but uses a mixture of a normal and
#' a (mirrored) log-normal distribution instead of a simple normal
#' distribution.
#'
#' @param meshpoints Vector of meshpoints defining the coordinates of the
#'     borders of the grid cells (should be ordered from smaller to larger)
#' @param x A single numerical value defining the starting coordinate for
#'     dispersal (same coordinate system as for "meshpoints")
#' @param sd Standard deviation of the normal distribution (same unit as for
#'     meshpoints and x)
#' @param meanlog Mean (on the log scale) of the log-normal distribution (same
#'     unit as for meshpoints and x)
#' @param sdlog Standard deviation (on the log scale) of the log-normal
#'     distribution (same unit as for meshpoints and x)
#' @param p Probability of dispersing according to the normal distribution.
#'     Conversely, the probability of dispersing according to the (mirrored)
#'     log-normal distribution will be (1-p).
#' 
#' @return A vector of length (length(meshpoints)-1), giving the probability of
#'     dispersing to each grid cell defined by the meshpoints
#'
#' @examples
#' par(mfrow = c(2, 3))
#' # Coarse grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 20)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_nln(meshpoints, x = -200, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.9)
#' stickplot(midpoints, probs, ylim = ylim, main = "20 bins, p = 0.9")
#' probs <- gx_intgr_nln(meshpoints, x = -50, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.5)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "20 bins, p = 0.5")
#' probs <- gx_intgr_nln(meshpoints, x = 20, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.2)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "20 bins, p = 0.2")
#' # Finer grid
#' kp <- create.kernel.params(L = -200, U = 200, n.bins = 50)
#' meshpoints <- kp$b
#' midpoints <- kp$y
#' probs <- gx_intgr_nln(meshpoints, x = -200, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.9)
#' stickplot(midpoints, probs, ylim = ylim, main = "50 bins, p = 0.9")
#' probs <- gx_intgr_nln(meshpoints, x = -50, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.5)
#' stickplot(midpoints, probs, ylim = ylim, col = "blue", main = "50 bins, p = 0.5")
#' probs <- gx_intgr_nln(meshpoints, x = 20, sd = 5, meanlog = 3.5, sdlog = 0.6, p = 0.2)
#' stickplot(midpoints, probs, ylim = ylim, col = "magenta", main = "50 bins, p = 0.2")
#' 
#' @export

gx_intgr_nln <- function(meshpoints, x, sd, meanlog, sdlog, p.leave) {
    # Calculate CDF at each meshpoint
    cdfAtMesh <- ((1-p.leave) * pnorm(meshpoints - x, mean = 0, sd = sd) +
                  p.leave * 1/2*(plnorm(meshpoints - x, meanlog = meanlog, sdlog = sdlog) +
                               plnorm((x - meshpoints), meanlog = meanlog, sdlog = sdlog,
                                      lower.tail = FALSE)))
    # Calculate the probability of dispersing between each pair of meshpoints
    probs <- diff(cdfAtMesh)
    stopifnot(sum(probs) <= 1)
    # Return
    return(probs)
}

### * Helper function: stickplot()

#' Produce a stick plot
#'
#' @param x X coordinates
#' @param y Y quantities
#'
#' @examples
#' par(mfrow = c(1, 2))
#'
#' x <- seq(-4, 4, length.out = 11)
#' y <- dnorm(x)
#' stickplot(x, y, las = 1)
#'
#' x <- seq(-4, 4, length.out = 31)
#' y <- dnorm(x)
#' stickplot(x, y, las = 1)
#'
#' @export

stickplot <- function(x, y, ylim = NULL, ...) {
    if (is.null(ylim)) {
        ylim <- extendrange(y)
        if (max(ylim) < 0) { ylim[2] <- 0 }
        if (min(ylim) > 0) { ylim[1] <- 0 }
    }
    plot(x, y, type = "n", ylim = ylim, ...)
    for (i in seq_along(x)) {
        lines(rep(x[i], 2), c(0, y[i]), ...)
    }
}

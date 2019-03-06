#' Create an initial population vector
#'
#' Given a kp (with U, L and n.bins), a number of individuals to introduce, and
#' a width of the introduction area, this function calculates how many bins are
#' covered by the introduction area, rounds this up if needed, and distributes
#' evenly the introduced individuals on it. This makes it easier to change the
#' number of bins of a kp, but to keep the same spatial distribution of
#' introduced individuals.
#'
#' @export
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
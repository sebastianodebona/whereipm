#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Modelling population dynamics through space using IPMs
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Defining vital rate parameters

vp <- create.vital.params(s.int = 0,
                    s.slo = -0.05,
                    r.int = log(2),
                    r.slo = -0.01,
                    p.stat = 0.8,
                    s.stat = 1,
                    s.mob = 2,
                    p.stat.nat = .5,
                    s.stat.nat = 0.6,
                    s.mob.nat = 4)


# Defining kernel parameters

kp <- create.kernel.params(L=0, U=100, n.bins=200)


# Testing the initial kernel

#test.init.kernel(kernel.params = kp, vital.params = vp)


# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 30, rep(0, kp$n.bins/2))


# Running the IPM

IPM <- runIPM(kp, vp, init.pop = pop0, niter=20, unidirectional=F)

plotly::plot_ly(x=1:ncol(IPM$pop.size), y=rev(kp$y), z=IPM$pop.size, type="contour")

draw.pop.dyn(IPM.model = IPM, breaks=20)




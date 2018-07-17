#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Modelling population dynamics through space using IPMs
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Defining vital rate parameters

vp <- create.vital.params(s.int = 1.6,
                    s.slo = -0.05,
                    r.int = -1,
                    r.slo = -0.01,
                    p.stat = 0.8,
                    s.stat = 1,
                    s.mob = 15,
                    p.stat.nat = .5,
                    s.stat.nat = 0.6,
                    s.mob.nat = 9.5)


# Defining kernel parameters

kp <- create.kernel.params(L=0, U=100, n.bins=100)


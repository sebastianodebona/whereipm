#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Modelling population dynamics through space using IPMs
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(whereipm)

# Defining vital rate parameters

vp <- mock.vital.params()


# Defining kernel parameters
kp <- create.kernel.params(L = -50, U = 50, n.bins = 100)
# check why weird stuff happens when U/n.bins is changed


# Testing the initial kernel

#test.init.kernel(kernel.params = kp, vital.params = vp)

K <- build.kernel(kp, vp, pop.size = pop0, test=T)

# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 30, rep(0, kp$n.bins/2))


# Running the IPM

IPM <- runIPM(kp, vp, init.pop = pop0, niter=25, unidirectional=F)
# 3d plot of population size through time and space
plotly::plot_ly(x=1:ncol(IPM$pop.size), y=rev(kp$y), z=IPM$pop.size, type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM, breaks=25)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Now with REAL DATA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# loading data
load('data/data.Rda')
head(data, 3)
str(data)

# only individuals that did produce offspring have a value for recruited. I need to manually add zeros to all females (remember this is a subset) that have no information == 0 rectruitment
data$recruited <- ifelse(is.na(data$recruited), 0, data$recruited)

# subsetting for sex, habitat and time (here I am building a model for females in pools, only in the early phase)
data.early <- subset(data, Habitat=='P' & sex=='F' & time>=1 & time<=12)
data.mid <- subset(data, Habitat=='P' & sex=='F' & time>=26 & time>=37)
data.late <- subset(data, Habitat=='P' & sex=='F' & time>=51 & time<=62)

source('DIPM_RunModels.R')

# parameters for Kernel
# kernel parameters
kp <- create.kernel.params(L=-500, U=500, n.bins = 1000)
# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 10, rep(0, kp$n.bins/2))

#### EARLY IPM ####
# Storing values into vital params list, divided into early and late
# given local density is scaled (but not centered) the estimates for the slope have to be transformed back to the non-scaled value.
vp.early <- list(s.int = summary(e.survival)$coefficients['(Intercept)','Estimate'],
                 s.slo = 0,
                 s.multiplier = 1,
                 
                 r.int = summary(e.recruit)$coefficients$zi['(Intercept)','Estimate'],
                 r.slo = 0,
                 r.multiplier = 1,
                 f.int = summary(e.recruit)$coefficients$cond['(Intercept)','Estimate'],
                 f.slo = summary(e.recruit)$coefficients$cond['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 distribution = 'ztpoisson',
                 
                 p.stat.int.ad = summary(e.dispersal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.ad = 0,
                 sigma.stat.ad = 0.25, 
                 mu.mob.int.ad = summary(e.dispersal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.ad = 0,
                 sigma.mob.ad = sqrt(sum(as.data.frame(lme4::VarCorr(e.dispersal.dist))[,'sdcor'])),
                 
                 p.stat.int.nat = summary(e.natal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.nat = 0,
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = summary(e.natal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.nat = 0,
                 sigma.mob.nat = sqrt(sum(as.data.frame(lme4::VarCorr(e.natal.dist))[,'sdcor']))
)

# Running the IPM for the early stage
IPM.early <- runIPM(kp, vp.early, init.pop = pop0, niter=12)
plotly::plot_ly(x=1:ncol(IPM.early$pop.size), y=rev(kp$y), z=log(IPM.early$pop.size+1), type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM.early, breaks=15, Hue='Reds')
# population explodes after 15 months: we probably need to set a limit to population growth rate manually, since survival is not density dependent.
####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#### MID IPM ####
vp.mid <- list(s.int = summary(m.survival)$coefficients['(Intercept)','Estimate'],
                s.slo = summary(m.survival)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                s.multiplier = 1,
                
                r.int = summary(m.recruit)$coefficients$zi['(Intercept)','Estimate'],
                r.slo = summary(m.recruit)$coefficients$zi['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                r.multiplier = 1,
                f.int = summary(m.recruit)$coefficients$cond['(Intercept)','Estimate'],
                f.slo = summary(m.recruit)$coefficients$cond['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                f.distribution = 'ztpoisson',
                
                p.stat.int.ad = summary(m.dispersal.p)$coefficients['(Intercept)','Estimate'],
                p.stat.slo.ad = summary(m.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                sigma.stat.ad = 0.25, 
                mu.mob.int.ad = summary(m.dispersal.dist)$coefficients['(Intercept)','Estimate'],
                mu.mob.slo.ad = summary(m.dispersal.dist)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                sigma.mob.ad = sqrt(sum(as.data.frame(lme4::VarCorr(m.dispersal.dist))[,'sdcor'])),
                
                p.stat.int.nat = summary(m.natal.p)$coefficients['(Intercept)','Estimate'],
                p.stat.slo.nat = summary(m.natal.p)$coefficients['s_loc_density.mid','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                sigma.stat.nat = 0.25, 
                mu.mob.int.nat = summary(m.natal.dist)$coefficients['(Intercept)','Estimate'],
                mu.mob.slo.nat = summary(m.natal.dist)$coefficients['s_loc_density.mid','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                sigma.mob.nat = sqrt(sum(as.data.frame(lme4::VarCorr(m.natal.dist))[,'sdcor']))
)


# Running the IPM for the early stage
IPM.mid <- runIPM(kp, vp.mid, init.pop = pop0, niter = 36)
plotly::plot_ly(x=1:ncol(IPM.mid$pop.size), y=rev(kp$y), z=IPM.mid$pop.size, type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM.mid, breaks=48, Hue='Blues')
plot(IPM.mid$pop.spread + 0.5, col='darkblue', type='o')
plot(colSums(IPM.mid$pop.size), col='darkblue')
plot(IPM.mid$pop.size[which(IPM.late$y==0.5),], col='darkblue')
####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#### LATE IPM ####
vp.late <- list(s.int = summary(l.survival)$coefficients['(Intercept)','Estimate'],
                 s.slo = summary(l.survival)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 s.multiplier = 1,
                 
                 r.int = summary(l.recruit)$coefficients$zi['(Intercept)','Estimate'],
                 r.slo = 0,
                 r.multiplier = 1,
                 f.int = summary(l.recruit)$coefficients$cond['(Intercept)','Estimate'],
                 f.slo = 0,
                 f.distribution = 'ztpoisson',
                 
                 p.stat.int.ad = summary(l.dispersal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.ad = summary(l.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 sigma.stat.ad = 0.25, 
                 mu.mob.int.ad = summary(l.dispersal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.ad = 0,
                 sigma.mob.ad = sqrt(sum(as.data.frame(lme4::VarCorr(l.dispersal.dist))[,'sdcor'])),
                 
                 p.stat.int.nat = summary(l.natal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.nat = summary(l.natal.p)$coefficients['(Intercept)','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = summary(l.natal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.nat = 0,
                 sigma.mob.nat = sqrt(sum(as.data.frame(lme4::VarCorr(l.natal.dist))[,'sdcor']))
)

# increasing the space will also increase the maximum population size. This is because at the moment, individuals that overcome the boundaries are lost forever, whereby they could and should move back into the frame we are looking at or produce offspring that would do so.
# a smart way around it is to change the scale of movement distances (mu of the lognormal distributions), which can be lowered to a tenth or hundredth of its current value. This way we'll be projecting movements in tens or hundreds of meters instead of actual meters.


# Running the IPM for the LATE stage
IPM.late <- runIPM(kp, vp.late, init.pop = pop0, niter=36)
plotly::plot_ly(y=kp$y, x=1:ncol(IPM.late$pop.size), z=IPM.late$pop.size, type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM.late, breaks=48, Hue='Greens')
plot(IPM.late$pop.spread + 0.5, col='darkgreen')
plot(colSums(IPM.late$pop.size), col='darkgreen')
plot(IPM.late$pop.size[which(IPM.late$y==0.5),], col='darkgreen')



# creating array to store results of perturbation analyses
perturbation <- array(NA, dim = c(4,14,3), dimnames = list(c('speed', 'acceleration', 'lambda', 'pop'), c('s.int', 's.slo',  'r.int',  'r.slo', 'f.int',  'f.slo',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'mu.mob.slo.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat', 'mu.mob.slo.nat'), c('IPM.early', 'IPM.mid', 'IPM.late')))

phases <- list(IPM.early, IPM.mid, IPM.late)
params <- list(early = c('s.int',  'r.int', 'f.int', 'f.slo', 'p.stat.int.ad',  'mu.mob.int.ad', 'p.stat.int.nat', 'mu.mob.int.nat'),
               mid = c('s.int', 's.slo',  'r.int',  'r.slo', 'f.int',  'f.slo',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'mu.mob.slo.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat', 'mu.mob.slo.nat'), 
               late = c('s.int', 's.slo',  'r.int', 'f.int',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat'))

start <- Sys.time()
for(I in 1:length(phases)){
  for(j in 1:length(params[[I]])){
    temp <- perturb(orig.IPM = phases[[I]], epsilon = 0.01, par = params[[I]][j])
    perturbation[, params[[I]][[j]], I] <- unlist(temp)
  }
}
Sys.time() - start


save(perturbation, file = 'sensitivity.rds')




oopt <- ani.options(interval = 0.5, nmax = 48)
for (i in 1:ani.options("nmax")) {
  if(i==1){
    plot(0, xlim=range(IPM.mid$y), ylim=c(0, max(IPM.mid$pop.size)), col='white')
  }
  draw.pop.distr(x=IPM.mid$y, pop=IPM.mid$pop.size[,i], col='darkgreen', new.plot = F)
  ani.pause()
}
ani.options(oopt)

saveHTML({ani.options(interval = 0.05, nmax = 50)
  brownian.motion(pch = 21, cex = 5, col = "red", bg = "yellow")},
  img.name = "bm_plot", ani.height = 300, ani.width = 550,
    title = "Demonstration of the Brownian Motion",
  description = c("Random walk on the 2D plane: for each point",
                  "(x, y), x = x + rnorm(1) and y = y + rnorm(1)."))

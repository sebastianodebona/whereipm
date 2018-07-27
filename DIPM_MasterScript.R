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
kp <- create.kernel.params(L = 0, U = 100, n.bins = 200)
# check why weird stuff happens when U/n.bins is changed


# Testing the initial kernel

#test.init.kernel(kernel.params = kp, vital.params = vp)


# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 30, rep(0, kp$n.bins/2))


# Running the IPM

IPM <- runIPM(kp, vp, init.pop = pop0, niter=20, unidirectional=F)
# 3d plot of population size through time and space
plotly::plot_ly(x=1:ncol(IPM$pop.size), y=rev(kp$y), z=IPM$pop.size, type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM, breaks=20)



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
data.early <- subset(data, Habitat=='P' & sex=='F' & time<35 & time>=12)
data.late <- subset(data, Habitat=='P' & sex=='F' & time>=61-23)

# running models to extract vital rate parameters for Early and Late phase

#### Survival ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.survival <- lme4::glmer(survived ~ 1 +
                                   + (1|FishID) + (1|reach) + (1|sampling),
                        data=data.early,
                        family=binomial,
                        lme4::glmerControl(optimizer=c('bobyqa')))
summary(e.survival)
# in early phase survival is not affectd by density
l.survival <- lme4::glmer(survived ~ s_loc_density +
                            + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                          data=data.late,
                          family=binomial,
                          lme4::glmerControl(optimizer=c('bobyqa')))
summary(l.survival)
# in late phase survival is positively affected by density

#### Recruitment ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.recruit <- glmmTMB::glmmTMB(recruited ~ 1 +
                              + (1|FishID) + (1|reach) + (1|sampling),
                            zi = ~ 1 +
                              + (1|FishID) + (1|reach)  + (1|sampling),
                            family=list(family="truncated_poisson", link="log"),
                            data=rm.NAs(data.early, 'recruited', fixEff=c('s_loc_density')))
# in the early phase, density affects the number of offspring but not the probability to produce any.
l.recruit <- glmmTMB::glmmTMB(recruited ~ 1 +
                                + (1|FishID) + (1|reach) + (1|sampling),
                              zi = ~ s_loc_density +
                                + (1|FishID) + (1|reach)  + (s_loc_density|sampling),
                              family=list(family="truncated_poisson", link="log"),
                              data=rm.NAs(data.late, 'recruited', fixEff=c('s_loc_density')))
# the late phase is not affected by density
summary(e.recruit)
summary(l.recruit)

#### Adult dispersal probability ####+++++++++++++++++++++++++++++++++++++++++++++++++
e.dispersal.p <- lme4::glmer(p_moving ~ s_loc_density +
                                 + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                               family = 'binomial',
                               data = data.early,
                               lme4::glmerControl(optimizer=c('bobyqa')))
# small negative effect of density on dispersal (high density, low dispersal) in early phase
l.dispersal.p <- lme4::glmer(p_moving ~ s_loc_density +
                               + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                             family = 'binomial',
                             data = data.late,
                             lme4::glmerControl(optimizer=c('bobyqa')))
# strong positive effect of density on dispersal in the late phase
summary(e.dispersal.p)
summary(l.dispersal.p)

#### Adult dispersal distance ####+++++++++++++++++++++++++++++++++++++++++++++++++++++++
e.dispersal.dist <- lme4::lmer(log(dist_moved) ~ s_loc_density +
                                        + (1|FishID) + (1|reach) + (s_loc_density|sampling),
                               data = data.early)
car::Anova(e.dispersal.dist, type='III')
# early phase distance does not depend on density
l.dispersal.dist <- lme4::lmer(log(dist_moved) ~ 1 +
                                 + (1|FishID) + (1|reach) + (1|sampling),
                               data = data.late)
car::Anova(l.dispersal.dist, type='III')
# neither late phase dispersal distance does
summary(e.dispersal.dist)
summary(l.dispersal.dist)

#### Natal dispersal probability ####++++++++++++++++++++++++++++++++++++++++++++++++
e.natal.p <- lme4::glmer(p_natal ~ 1 +
                                   + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                                 family = 'binomial',
                                 data = data.early,
                                 lme4::glmerControl(optimizer=c('bobyqa')))
# no effect of local density early on
l.natal.p <- lme4::glmer(p_natal ~ 1 +
                           + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                         family = 'binomial',
                         data = data.late,
                         lme4::glmerControl(optimizer=c('bobyqa')))
# no effect late either
summary(e.natal.p)
summary(l.natal.p)

#### Natal dispersal distance ####++++++++++++++++++++++++++++++++++++++++++++++++
e.natal.dist <- lme4::lmer(log(dist_natal) ~ 1 + 
                                     + (1|Maternity) + (1|reach.birth) + (1|sampling.birth),
                                   data = data.early)
car::Anova(e.natal.dist, type = 'III')
# no effect of density on distance in early phase
l.natal.dist <- lme4::lmer(log(dist_natal) ~ s_loc_density.mid  + 
                             + (1|Maternity) + (1|reach.birth) + (s_loc_density.mid|sampling.birth),
                           data = data.late)
car::Anova(l.natal.dist, type = 'III')
# neither in late
summary(e.natal.dist)
summary(l.natal.dist)

# Storing values into vital params list, divided into early and late
# given local density is scaled (but not centered) the estimates for the slope have to be transformed back to the non-scaled value.
vp.early <- list(s.int = summary(e.survival)$coefficients['(Intercept)','Estimate'],
                 s.slo = 0,
                 s.multiplier = 1,
                 
                 r.int = summary(e.recruit)$coefficients$zi['(Intercept)','Estimate'],
                 r.slo = 0,
                 r.multiplier = 1,
                 f.int = summary(e.recruit)$coefficients$cond['(Intercept)','Estimate'],
                 f.slo = 0,
                 distribution = 'ztpoisson',
                 
                 p.stat.int.ad = summary(e.dispersal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.ad = summary(e.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 sigma.stat.ad = 0.25, 
                 mu.mob.int.ad = summary(e.dispersal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.ad = summary(e.dispersal.dist)$coefficients['(Intercept)','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 sigma.mob.ad = sqrt(sum(as.data.frame(lme4::VarCorr(e.dispersal.dist))[,'sdcor'])),
                 
                 p.stat.int.nat = summary(e.natal.p)$coefficients['(Intercept)','Estimate'],
                 p.stat.slo.nat = 0,
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = summary(e.natal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.nat = 0,
                 sigma.mob.nat = sqrt(sum(as.data.frame(lme4::VarCorr(e.natal.dist))[,'sdcor']))
)

kp <- create.kernel.params()


# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 10, rep(0, kp$n.bins/2))

# Running the IPM for the early stage
IPM <- runIPM(kp, vp.early, init.pop = pop0, niter=15)
plotly::plot_ly(x=1:ncol(IPM$pop.size), y=rev(kp$y), z=IPM$pop.size, type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM, breaks=15)

# population explodes after 15 months: we probably need to set a limit to population growth rate manually, since survival is not density dependent.
vp.late <- list(s.int = summary(l.survival)$coefficients['(Intercept)','Estimate'],
                 s.slo = summary(l.survival)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
                 s.multiplier = 1,
                 
                 r.int = summary(l.recruit)$coefficients$zi['(Intercept)','Estimate'],
                 r.slo = summary(l.recruit)$coefficients$zi['s_loc_density','Estimate'],
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
                 p.stat.slo.nat = 0,
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = summary(l.natal.dist)$coefficients['(Intercept)','Estimate'],
                 mu.mob.slo.nat = summary(l.natal.dist)$coefficients['s_loc_density.mid','Estimate'],
                 sigma.mob.nat = sqrt(sum(as.data.frame(lme4::VarCorr(l.natal.dist))[,'sdcor']))
)

kp <- create.kernel.params(L=-150, U=150, n.bins = 300)
# increasing the space will also increase the maximum population size. This is because at the moment, individuals that overcome the boundaries are lost forever, whereby they could and should move back into the frame we are looking at or produce offspring that would do so.
# a smart way around it is to change the scale of movement distances (mu of the lognormal distributions), which can be lowered to a tenth or hundredth of its current value. This way we'll be projecting movements in tens or hundreds of meters instead of actual meters.

# Defining the initial population size
pop0 <- c(rep(0, kp$n.bins/2-1), 30, rep(0, kp$n.bins/2))

# Running the IPM for the LATE stage
IPM <- runIPM(kp, vp.late, init.pop = pop0, niter=120)
plotly::plot_ly(x=kp$y, y=1:ncol(IPM$pop.size), z=t(IPM$pop.size), type="contour")
# another way to do it
draw.pop.dyn(IPM.model = IPM, breaks=50)

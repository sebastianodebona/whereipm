#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# PART 2: PARAMETERS
# this code uses the prepped data to fit models in order to get needed parameters to describe vital rates
#
#___________________________________________________________________________________


setwd('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Dispersal/Chapter I - Dispersal IPM')

source('DIPM_DataPrepping.R')


########################################
#
# modelling regressions of vital rates
#
########################################
# 1) SURVIVAL
# 2) REPRODUCTION
# 3) ADULT DISPERSAL
# 4) NATAL DISPERSAL


# loading needed packages
library(lme4)
#library(VGAM)
#library(mgcv)
#library(glmmADMB)


# 1) SURVIVAL regression
# survival is here density dependent (local density)
# location (specific reach), sampling (time variable), section (spatial variable) and FishID are included as random effects to the model
s0<-glmer(surv~1 + (1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=binomial)
s1<-glmer(surv~spot.dens+(1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=binomial)
#s.null<-glm(surv~spot.dens, data=d, family=binomial)
#AIC(s0,s1)
#summary(s1)
surmod<-s1



# 2) REPRODUCTION regression
# Here modelled as a Poisson. See below for alternatives (2-way split and zero truncated poisson)
# offspring corresponds to the number of offspring recruited (lots of zeros)
# the offspring size is considered density dependent (both local and global density)
# location (specific reach), sampling (time variable), section (spatial variable) and FishID are included as random effects to the model
off0 <- glmer(offspring~1+(1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=poisson)
off1 <- glmer(offspring~spot.dens+(1|sampling) + (1|section) + (1|location)+ (1|FishID), data=d, family=poisson)
#AIC(off0,off1)
#summary(off1)
repmod<-off1



# 3) ADULT DISPERSAL
#I will use the function in Radinger and Wolter 2013 (package fishmove)
library(fishmove)
# Radinger and Wolter describe a leptokurtic dispersal kernel as the combination of two normal distributions with same mean but different variance, as if the population was characterized by a proportion of static individuals and one of moving individuals. This creates the 'fat-tails' typical of dispersal kernels
# The function 'fishmove.estimates' I am about to use takes as an imput a vector of (absolute value) distances and return the three parameters describing a leptokurtic dispersal kernel:
# sigma.mov=variance of distance moved by proportion of moving fish in the pop
# sigma.stat=variance of distance moved by proportion of static fish in the pop
krnl.parameters <- fishmove.estimate(abs(subset(d,!is.na(dist.moved))$dist.moved+0.01))
# function to determine the dispersal kernel
# this function calculates the probability fx of a certain dispersal movement. This functio is symmetrical, and centered on mu (in this specific case, 0). Changing mu I can change where it is centered. I can also interpret it otherwise: if I think of mu as my position now, and x as my position next, the function calculates the probability of the movement x-mu occurring. This is why I am using this function also down below, when contructing the movement function mox.
DK<-function(mu,x, param.estimates){
  # extracting parameters
  sigma_stat<-as.numeric(param.estimates[['estimate']]['sigma_stat'])
  sigma_mob<-as.numeric(param.estimates[['estimate']]['sigma_mob'])
  p<-as.numeric(param.estimates[['estimate']]['p'])
  fx = p*(1/(sqrt(2*pi*sigma_stat^2)))*exp(-((x-mu)^2/(2*sigma_stat^2))) +
    (1-p)*(1/(sqrt(2*pi*sigma_mob^2)))*exp(-((x-mu)^2/(2*sigma_mob^2)))
  return(fx)
}

movmod <- krnl.parameters



# 4) NATAL DISPERSAL
# (see MOVEMENT for rational over splitting)
# Again, doing this in two ways, with Radinger and Wolter function and with splitting
krnl.parameters.natal <- fishmove.estimate(abs(subset(d,!is.na(off.dist))$off.dist+0.001))
natmod<-krnl.parameters.natal








#$$$$$$$$$$$$$$$$$$
# PLOTTING
# creating discrete categories form density for plotting purposes
op<-par()
library(plyr)
d$dens_cat<-cut(d$spot.dens, seq(0,max(d$spot.dens, na.rm=T), by=2), right=FALSE, labels=1:25)
toplot<-ddply(d,.(dens_cat),summarize,mean.s=mean(surv, na.rm=T),mean.rep=mean(repr, na.rm=T),mean.rec=mean(recruited,na.rm=T),mean.off=mean(offspring,na.rm=T),p.mov=mean(moved, na.rm=T),dir.mov=mean(direction, na.rm=T), mean.mov=mean(dist.moved, na.rm=T))
toplot$dens_cat<-as.numeric(toplot$dens_cat)*2

par(mfrow=c(2,2))
plot(toplot$mean.s~toplot$dens_cat, ylim=c(0,1), xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='survival prob', main='Survival')
#abline(a=plogis(coef(summary(s1))['(Intercept)','Estimate']), b=plogis(coef(summary(s1))['spot.dens','Estimate']))
xvalues<-seq(0,50,2)
yvalues<-exp(coef(summary(s1))['(Intercept)','Estimate']+coef(summary(s1))['spot.dens','Estimate']*xvalues)/(1+exp(coef(summary(s1))['(Intercept)','Estimate']+coef(summary(s1))['spot.dens','Estimate']*xvalues))
lines(xvalues,yvalues, col='red')


plot(toplot$mean.off~toplot$dens_cat, ylim=c(0,0.5),xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='mean offsrping n', main='Recruitment')
yvalues<-exp(coef(summary(off1))['(Intercept)','Estimate']+coef(summary(off1))['spot.dens','Estimate']*xvalues)
lines(xvalues,yvalues, col='red')


x<-seq(-50,50,by=0.5)
ys<-DK(param.estimates=krnl.parameters,mu=0,x=x)
hist(d$dist.moved,breaks=200,xlim=c(-60,60),ylim=c(0,7000),main='Breeding Dispersal',xlab='space')
lines(x,ys*nrow(subset(d,!is.na(dist.moved))),type='l',xlim=c(-60,60),col='red',ylim=c(0,7000))


ys.natal<-DK(param.estimates=krnl.parameters.natal,mu=0,x=x)
hist(d$off.dist,breaks=200,xlim=c(-60,60),ylim=c(0,500), main='Natal Dispersal',xlab='space')
lines(x,ys.natal*nrow(subset(d,!is.na(off.dist))),type='l',xlim=c(-60,60),ylim=c(0,500),col='red',xlab='space',ylab='frequency')


par(op)


# ALTERNATIVE FOR REPRODUCTION: 2-way split of offspring size
# REPRODUCTION PROBABILITY regression
# rp0<-glmer(repr~1+(1|sampling) + (1|section) + (1|habitat), data=d, family=binomial)
# rp1<-glmer(repr~spot.dens+(1|sampling) + (1|section) + (1|habitat), data=d, family=binomial)
# AIC(rp0,rp1)
# summary(rp1)
# repmod<-rp0
# 
# plot(toplot$mean.rep~toplot$dens_cat, ylim=c(0,1),xlim=c(0,50), pch=19, xlab='density (inds/m)', ylab='reproduction prob', main='Reproduction')
# #lines(xvalues,predict(rp1, list(spot.dens=xvalues), type='response'), col='red')
# 
# # RECRUITMENT
# #fitting a zero truncated poisson distribution
# rc0<-vglm(recruited~1, data=d, family=pospoisson())
# rc1<-vglm(recruited~spot.dens, data=d, family=pospoisson())
# summary(rc1)
# 
# recmod <- rc1
# 


# ALTERNATIVE FOR AD DISPERSAL
# #******* UNCOMMENT TO TILL NEXT ASTERISCS TO RUN ALTERNATIVE ANALYSES ******
# # Alternatively
# I will consider movement as distance moved, and break it into three
# 1) probability of moving
# 2) distance moved (in absolute value)
# 3) direction moved

# # probability of moving
# m0 <- glmer(moved~1+(1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=binomial)
# m1 <- glmer(moved~spot.dens+(1|sampling) + (1|section) + (1|location) +(1|FishID), data=d, family=binomial)
# AIC(m0,m1)
# summary(m1)
# 
# plot(toplot$p.mov~toplot$dens_cat, ylim=c(0,1), xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='prob of moving', main='Movement (prob)')
# yvalues<-exp(coef(summary(m1))['(Intercept)','Estimate']+coef(summary(m1))['spot.dens','Estimate']*xvalues)/(1+exp(coef(summary(m1))['(Intercept)','Estimate']+coef(summary(m1))['spot.dens','Estimate']*xvalues))
# lines(xvalues,yvalues, col='red')
# 
# # distance moved
# #!!! This has to be changed to accommodate the positive only dispersal. see what distributions dispersal is modelled with in published data and apply that
# # not much information. some on geometric or exponential distributions (Starrfeldt and Kokko 2010)
# # other big problem I had before. When modelling the distance moved I have to EXCLUDE individuals who did not move.
# hist(subset(d,dist.moved!=0)$dist.moved,breaks=100)
# abline(v=0,col='red')
# dm0 <- lmer(I(log(abs(dist.moved)))~1+(1|sampling)+(1|section) + (1|location)+ (1|FishID), data=subset(d, dist.moved!=0))
# dm1 <- lmer(I(log(abs(dist.moved)))~spot.dens+(1|sampling) + (1|section) + (1|location)+(1|FishID), data=subset(d, dist.moved!=0))
# AIC(dm0,dm1)
# summary(dm0)
# summary(dm1)
# #plot to see the goodness of fit
# x <- seq(0,20,by=0.1)
# hx <- coef(summary(dm0))[1,1]*exp(-coef(summary(dm0))[1,1]*x)
# lines(x,hx*length(subset(d,dist.moved!=0)$dist.moved), col='red')
# 
# plot(abs(d$dist.moved)~d$spot.dens, xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='meters moved', main='Movement (m)')
# #lines(xvalues,predict(dm1, list(spot.dens=xvalues), type='response'), col='red')
# 
# # direction moved
# dir0 <-glmer(direction~1+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# dir1 <- glmer(direction~spot.dens+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# AIC(dir0,dir1)
# summary(dir1)
# # it seems like the higher the density, the higer the chance to move upstream
# plot(toplot$dir.mov~toplot$dens_cat, ylim=c(0,1), xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='prob of moving downstream', main='Movement (direction)')
# yvalues<-exp(coef(summary(dir1))['(Intercept)','Estimate']+coef(summary(dir1))['spot.dens','Estimate']*xvalues)/(1+exp(coef(summary(dir1))['(Intercept)','Estimate']+coef(summary(dir1))['spot.dens','Estimate']*xvalues))
# lines(xvalues,yvalues, col='red')
# 
# #***************************************************************************




# ALTERNATIVE FOR NATAL DISPERSAL
# #******* UNCOMMENT TO TILL NEXT ASTERISCS TO RUN ALTERNATIVE ANALYSES ******
# # ALTERNATIVELY
# # probability of offspring moving from natal patch
# n0 <- glmer(off.moved~1+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# n1 <- glmer(off.moved~spot.dens+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# AIC(n0,n1)
# summary(n1)
# # distance moved by offspring
# dn0 <- lmer(I(abs(off.moved))~1+(1|sampling) + (1|section) + (1|location), data=d)
# dn1 <- lmer(I(abs(off.moved))~spot.dens+(1|sampling) + (1|section) + (1|location), data=d)
# AIC(dn0,dn1)
# summary(dn1)
# # direction moved
# ndir0 <-glmer(off.direction~1+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# ndir1 <- glmer(off.direction~spot.dens+(1|sampling) + (1|section) + (1|location), data=d, family=binomial)
# AIC(dir0,dir1)
# summary(dir1)
# #****************************************************************************



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# this code constructs an IPM where position in the stream is the state variable. Position is then linked to vital rates such as survival, reproduction, movement (i.e. growth) and natal dispersal (i.e. inheritance) through the effects of local density on these traits.
# this code makes use of a dataset previosly constructed

rm(list=ls())

setwd('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Dispersal/Chapter I - Dispersal IPM')

# reading dataset 
load('DIPMdata.Rda')
head(DIPMdata)

d<-DIPMdata
# removing some typos in distance.next
d<-d[-which(d$position.next>110),]
# also, adding NA when recruited is 0, to model truncated poisson
d$offspring <- d$recruited
d$recruited<-ifelse(d$recruited==0,NA,d$recruited)
head(d)


# I am here creating a variable called SECTION to account for heterogeneity in the stream (to be included as a random effect in the model)
# how long the sections will be is determined by the following parameter, to be adjusted
how.fine<-5
# the following command creates sections of the stream cut each "how.fine" meters
d$section<-cut(d$position, seq(0,max(d$position, na.rm=T), by=how.fine), right=FALSE)


# since the reach labels change every month (different interns or actual morphometric change in the stream), the position is showing small changes without movement
# I need to remove small changes due to variations in stream morphology
# I calculate the distance moved, and remove the negligible small scale movements likely due to changes in how the distance is calculated
# dist.moved is therefore defined as the actualy distance moved if this is not smaller than a certain value, either absolute or relative to the reach length
neglect <- 0
# as an alternative, the threashold could be a fourth of the reach legth, meaning that I will ignore a 2 meters movement in a 8 meters reach, but consider a 1 meter movement in a 3 meter reach
#neglect <- d$reach.l/4

d$dist.moved<-ifelse(abs(d$position.next-d$position)>neglect,
                     d$position.next-d$position,
                     0)
# moved corresponds to a binomial where 0 is stay and 1 is moved
d$moved<-as.numeric(d$dist.moved!=0)
# finally, direction corresponds to 0 whether the movement was upstream and 1 whether it was downstream
d$direction <- sign(d$dist.moved)
d$direction <- ifelse(d$direction==0,NA,d$direction)
d$direction <- ifelse(d$direction==-1,0,d$direction)


#should do the same for natal dispersal
#distance the offsrping moved from the mother in the first month (calculated as halfdistance of the total after 2 months)
d$off.dist<-ifelse(abs(d$halfpos_off-d$position)>neglect,
                     d$halfpos_off-d$position,
                     0)
# moved or not
d$off.moved<-as.numeric(d$off.dist!=0)

# direction of natal dispersal
d$off.direction <- sign(d$off.dist)
d$off.direction <- ifelse(d$off.direction==0,NA,d$off.direction)
d$off.direction <- ifelse(d$off.direction==-1,0,d$off.direction)








########################################
#
# modelling regressions of vital rates
#
########################################

library(lme4)
#library(VGAM)
#library(mgcv)
library(glmmADMB)


# SURVIVAL regression
s0<-glmer(surv~1 + (1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=binomial)
s1<-glmer(surv~spot.dens+(1|sampling) + (1|section) + (1|location) + (1|FishID), data=d, family=binomial)
#s.null<-glm(surv~spot.dens, data=d, family=binomial)
AIC(s0,s1)
summary(s1)
surmod<-s1

# creating discrete categories form density for plotting purposes
library(plyr)
d$dens_cat<-cut(d$spot.dens, seq(0,max(d$spot.dens, na.rm=T), by=2), right=FALSE, labels=1:25)
toplot<-ddply(d,.(dens_cat),summarize,mean.s=mean(surv, na.rm=T),mean.rep=mean(repr, na.rm=T),mean.rec=mean(recruited,na.rm=T),mean.off=mean(offspring,na.rm=T),p.mov=mean(moved, na.rm=T),dir.mov=mean(direction, na.rm=T), mean.mov=mean(dist.moved, na.rm=T))
toplot$dens_cat<-as.numeric(toplot$dens_cat)*2

par(mfrow=c(2,3))
plot(toplot$mean.s~toplot$dens_cat, ylim=c(0,1), xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='survival prob', main='Survival')
#abline(a=plogis(coef(summary(s1))['(Intercept)','Estimate']), b=plogis(coef(summary(s1))['spot.dens','Estimate']))
xvalues<-seq(0,50,2)
yvalues<-exp(coef(summary(s1))['(Intercept)','Estimate']+coef(summary(s1))['spot.dens','Estimate']*xvalues)/(1+exp(coef(summary(s1))['(Intercept)','Estimate']+coef(summary(s1))['spot.dens','Estimate']*xvalues))
lines(xvalues,yvalues, col='red')







# trying for a moment to simplify things and model number of offsrping as a simple poisson or neg bin to avoid having to deal with complications of truncated distributions and multiple random effects
# COMMENTED OUT
# # REPRODUCTION PROBABILITY regression
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
plot(toplot$mean.off~toplot$dens_cat, ylim=c(0,0.5),xlim=c(0,50), pch=19, xlab='density(inds/m)', ylab='mean offsrping n', main='Recruitment')
yvalues<-exp(coef(summary(off1))['(Intercept)','Estimate']+coef(summary(off1))['spot.dens','Estimate']*xvalues)
lines(xvalues,yvalues, col='red')

#using a new variable to collapses recruited and repr (i.e. number of offsrping) -> d$offsrping
off0 <- glmer(offspring~1+(1|sampling) + (1|section) + (1|location), data=d, family=poisson)
off1 <- glmer(offspring~spot.dens+(1|sampling) + (1|section) + (1|location), data=d, family=poisson)
AIC(off0,off1)
summary(off1)
repmod<-off1








# MOVEMENT
#I will use the function in Radinger and Wolter 2013 (package fishmove)
# RIGHT HERE ->

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

x<-seq(-50,50,by=0.5)
ys<-DK(param.estimates=krnl.parameters,mu=0,x=x)

hist(d$dist.moved,breaks=200,xlim=c(-60,60),ylim=c(0,7000))
plot(x,ys*nrow(subset(d,!is.na(dist.moved))),type='l',xlim=c(-60,60),col='red',ylim=c(0,7000),xlab='space')
abline(h=0, lty=2)

movmod <- krnl.parameters






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







# NATAL DISPERSAL (see MOVEMENT for rational over splitting)
# Again, doing this in two ways, with Radinger and Wolter function and with splitting
hist(d$off.dist,breaks=100,xlim=c(-60,60),ylim=c(0,500))
krnl.parameters.natal <- fishmove.estimate(abs(subset(d,!is.na(off.dist))$off.dist+0.001))
ys.natal<-DK(param.estimates=krnl.parameters.natal,mu=0,x=x)
plot(x,ys.natal,type='l',xlim=c(-60,60),col='red',xlab='space',ylab='frequency')

natmod<-krnl.parameters.natal






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












#######################################
#
# Vital rates functions
#
######################################

# survival function
sx <- function(density) { 
  # x is the vector of midpoint values in the bins
  # density is a vector of the SAME length, with local density for the midpoints xs
  # I cannot use the function predict, since I have random effects in the models and those should be specified. Although, I can manually calculate the probability of survival per location (mediated by density) given the vector of spot densities
  #plogis(predict(surmod, data.frame(spot.dens=density)))
  plogis(coef(summary(surmod))['(Intercept)','Estimate']+coef(summary(surmod))['spot.dens','Estimate']*density)
}

# fertility function
fx <- function(density) {
  #mean number of offspring produced:
  # I cannot use the function predict, since I have random effects in the models and those should be specified. Although, I can manually calculate the probability of survival per location (mediated by density) given the vector of spot densities
  exp(coef(summary(repmod))['(Intercept)','Estimate']+coef(summary(repmod))['spot.dens','Estimate']*density)*3

  #!!! given the high number of orphans, this still has to be corrected for underestimated reproductive output
}


# movement function
mox <- function(y, x) {
  # for the movement function I have a density INdependent probability distribution of movements, that has to be centered each time on the current position
  DK(mu=x,x=y,param.estimates=movmod)
}



# natal dispersal function (new recruits position)
nax <- function(y, x) {
  DK(mu=x,x=y,param.estimates=natmod)
}



#######################################
#
# Kernel assembling
#
######################################

# number of bins
n.bins=200

# min and max value of the state variable
# I know it can't be less than 0 or higher than, lets say 110
L<-0
U<-110

# minsize here is negative because size is on the log scale in the dataset
b = L+c(0:n.bins)*(U-L)/n.bins # meshpoints
y = 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # midponts
h = y[2]-y[1] # bin width

# starting density vector
# I have 20 individuals in the fist 3 meters, and that is my starting population. To calculate density, I have to calculate how many ind per bin length, and divide that number BY bin length (because my density is individual/meter)
starting.ind.per.bin<-19/(3/h)
pop.size0<-c(rep(0,(n.bins-6)/2),10,10,10,10,10,10,rep(0,(n.bins-6)/2))
density0<-pop.size0/h


# first of all, let have a look at the kernels for the initial density
# survival kernel
  S <- diag(sx(density0))
  image(t(S))
# movement kernel
  G <- outer(y, y, mox)
  image(t(G))
  
  
  
  
# P kernel
  Pmatrix <- h*as.matrix(G%*%S)
  image(t(Pmatrix))
  abline(0,1,lty=2)
  
  
# reproduction
  R <- diag(fx(density0))
  image(t(R))
# recruited position (natal dispersal)
  N <- outer(y, y, nax)
  image(t(N))
  
# F kernel
  Fmatrix <- h*as.matrix(N%*%R)
  image(t(Fmatrix))
  abline(0,1,lty=2)

# Whole kernel
  K<- Pmatrix + Fmatrix
  image(t(K)) 
  
  
  
  
# how long to run the simulation for
  nsims<-60
# creating objects to fill  
  pop.size.tosave<-array(NA,dim=c(n.bins,nsims+1))
# filling with first initial densities
  pop.size.tosave[,1]<-pop.size0
  
  K.list<-list()
  

# Function for calculating Kernel
K.fun <- function(y,densities){

  S <- diag(sx(densities)) #surv kernel
  G <- outer(y, y, mox) # movement kernel
  # P kernel
  Pmatrix <- h*as.matrix(G%*%S)
  
  R <- diag(fx(densities)) # reproduction
  N <- outer(y, y, nax) # recruited position (natal dispersal)
  # F kernel
  Fmatrix <- h*as.matrix(N%*%R)

  # Whole kernel
  K<- Pmatrix + Fmatrix
  
  return(K)
}

# Loop for simulation
for(i in 1:nsims){
  densities<-pop.size.tosave[,i]/h
  K.fun(y=y,densities=densities)
  K.list[[i]]<-K
  new.pop.size<-K%*%pop.size.tosave[,i]
  pop.size.tosave[,i+1]<-new.pop.size
}

# plotting total pop size through time
par(mfrow=c(1,1))
plot(colSums(pop.size.tosave), type='l', ylab='Total population size', xlab='months')


# plotting population spread through time
par(mfrow=c(2,3), oma=c(3,4,3,3),mar=c(3,3,1,2))
for(m in c(2,3,4,6,12,24)){
plot(pop.size.tosave[,m]~seq(1:n.bins), type='l', ylim=c(0,6))
legend('top', legend=paste('month= ',m, sep=''), bty='n')
abline(h=0,lty=2)
}
mtext('position (m/2)',1, outer=T)
mtext('local population size', 2, outer=T)


#library(fields)
par(mfrow=c(2,3), oma=c(3,4,3,3),mar=c(3,3,1,2))
for(m in c(2,3,4,6,12,24)){
  image(t(K.list[[m]]))
  par(new=T)
  contour(t(K.list[[m]]), axes=F)
  legend('top', legend=paste('month= ',m, sep=''), bty='n', text.col='white')
}



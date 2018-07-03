#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# PART 4: PARAMETER EXTRACTION and SIMULATION
# this code sources two other parts (parameters and function building), build the IPM and runs simulations
#
#___________________________________________________________________________________

rm(list=ls())


#Start<-Sys.time()
#source('DIPM_Parameters.R') # subsourcing DIPM_DataPrepping
#Sys.time()-Start

#!!!!!! WARNING! Some models do not converge


#######################################
#
# ALL PARAMETERS & KERNEL SPECS
#
######################################

# BASIC PARAMETERS


# Survival
surv.int <- 1.6 #intercept of linear regression between survival (logit transformed) and density 
surv.slo <- -0.05 #slope of linear regression between survival (logit transformed) and density 
surv.multiplier <- 1 # this value is multiplied to the mean survival probability (survival buffering)
# Reproduction
repr.int <- -1.00 #intercept of linear regression (log transformed) and density 
repr.slo <- -0.01 #slope of linear regression (log transformed) and density 
repr.multiplier <- 1 # this value is multiplied to the mean number of offspring produced
# Breeding Dispersal
disp.stat <- 1.00 #variance of stationary pop fraction
disp.mob <- 15.00 #variance of mobile pop fraction
disp.p <- 0.8 # proportion of stationary adults at density 0
disp.p.slo <- 0 # density dependence of proportion of adult dispersers
# Natal Dispersal
nat.stat <- 0.6 #variance of stationary pop fraction
nat.mob <- 9.50 #variance of mobile pop fraction
nat.p <- 0.5 # proportion of stationary juveniles at density 0
nat.p.slo <- 0 # density dependence of proportion of natal dispersers


#$$$$$$$$$$$
#
# MODIFIED PARAMETERS
#
#$$$$$$$$$$$

#here
# changed.param <-repr.multiplier<-4
changed.paramII <-surv.slo<-seq(-0.2,0.2, length.out = 3)
changed.paramIII <- disp.p <- seq(0,1, length.out = 3)


ps<-expand.grid(surv.int=surv.int,surv.slo=surv.slo,surv.multiplier=surv.multiplier,
                         repr.int=repr.int,repr.slo=repr.slo,repr.multiplier=repr.multiplier,
                         disp.stat=disp.stat,disp.mob=disp.mob,disp.p=disp.p,
                         nat.stat=nat.stat,nat.mob=nat.mob,nat.p=nat.p)


# General Kernel specifications
# number of bins (how many the intervals to approximate continuity are)
n.bins <- 100
# min and max value of the state variable
# I know it can't be less than 0 or higher than, lets say 100
L <- 0
U <- 100
# defining meshpoints (separating bins) and midpoints (central point of each bin)
b <- L+c(0:n.bins)*(U-L)/n.bins # meshpoints
y <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
h <- y[2]-y[1] # bin width (in meters)

# SEEDING POPULATION
# starting density vector
pop.size0<-c(30,rep(0,(n.bins-1)))
density0<-pop.size0/h


#######################################
#
# Specs for Simluation
#
#######################################


# how long to run the simulation for
nsims<-10 # months

# OBJECTS TO SAVE
# 1) Population size
pop.size.tosave<-array(NA,dim=c(n.bins,nsims+1,nrow(ps))) # creating array to fill wit population size
pop.size.tosave[,1,]<-pop.size0 # filling with first initial densities

# 2) Kernels
K.list<-list()

# 3) Population Spread
pop.spread<-array(NA,dim=c(nsims+1,nrow(ps)))
pop.spread[1,]<-max(sum(((y-mean(y))^2)*pop.size0)/(sum(pop.size0)-1),0)



##############################################
#
# SIMULATION OVER PARAMETER SPACE
#
##############################################

Start<-Sys.time()
for(k in 1:nrow(ps)){
  source('DIPM_Simulation.R')
}
Sys.time()-Start


# PLOTTING
colourz<-colorRampPalette(c("red","yellow","springgreen","royalblue"))

# # IF ONE PARAMETER WAS CHANGED
# # plotting spread through time
# matplot(pop.spread, type='l',col=colourz(nrow(ps)))
# # plotting total pop size through time
# matplot(colSums(pop.size.tosave), type='l', ylab='Total population size', xlab='months', col=colourz(nrow(ps)))
# # plotting final pop spread aross parameter space
# plot(1:nrow(ps),pop.spread[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
# axis(1, at=1:nrow(ps), labels=changed.param)
# # plotting final pop size aross parameter space
# plot(1:nrow(ps),colSums(pop.size.tosave)[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
# axis(1, at=1:nrow(ps), labels=changed.param)



# IF TWO PARAMETERS CHANGED
# creating matrix for countour plot to plot pop spread according to two params
param.mat<- matrix(NA, nrow=length(changed.paramII), ncol=length(changed.paramIII))
for(i in 1:length(changed.paramII)){
  for(j in 1:length(changed.paramIII)){
    param.mat[i,j]<-pop.spread[nsims+1,which(ps$surv.slo==changed.paramII[i] & ps$disp.p==changed.paramIII[j])]
  }
}
filled.contour(x=changed.paramIII,y=changed.paramII,z=t(param.mat), xlab='proportion of static (adults)',ylab='strength of dd for survival', main="Population spread")#,levels=c(0,2000,4000,6000,10000,15000,20000,25000,50000))
legend('topleft',legend='repr mult = 4',cex=0.8, box.lwd = 0, bg='transparent')

# # creating matrix for countour plot to plot pop SIZE according to two params
# pop.mat<-matrix(NA, nrow=length(changed.paramII), ncol=length(changed.paramIII))
# for(i in 1:length(changed.paramII)){
#   for(j in 1:length(changed.paramIII)){
#     pop.mat[i,j]<-colSums(pop.size.tosave)[nsims+1,which(ps$surv.slo==changed.paramII[i] & ps$disp.p==changed.paramIII[j])]
#   }
# }
# filled.contour(x=changed.paramIII,y=changed.paramII,z=t(pop.mat), xlab='proportion of static (adults)',ylab='strength of dd for survival', main="Final population size")

# #the following bit of code fixes one of the two parameters and plots results when the other parameter is changed
# sub.ps<-which(ps$repr.multiplier==2.3)
# matplot(pop.spread[,sub.ps], type='l',col=colourz(length(sub.ps)),main='Population spread through time',ylab=expression(sigma^2),xlab='months',xlim=c(0,80))
# legend('right',legend=changed.paramII,col=colourz(length(sub.ps)), lty=1, cex=0.8, title='survival slope')
# legend('topleft',legend='repr.mult=2.3',cex=0.8)
# # plotting total pop size through time
# matplot(colSums(pop.size.tosave[,,sub.ps]), type='l', main='Total population size', xlab='months',ylab='n individuals', col=colourz(length(sub.ps)),xlim=c(0,80))
# legend('right',legend=changed.param,col=colourz(length(sub.ps)), lty=1, cex=0.8,title='survival slope')
# legend('top',legend='repr.mult=2.3',cex=0.8)


# # plotting final pop spread aross parameter space
# plot(1:nrow(ps),pop.spread[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
# axis(1, at=1:nrow(ps), labels=changed.param)
# # plotting final pop size aross parameter space
# plot(1:nrow(ps),colSums(pop.size.tosave)[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
# axis(1, at=1:nrow(ps), labels=changed.param)
# 







# # plotting population spatial spread through time
# par(mfrow=c(2,3), oma=c(3,4,3,3),mar=c(3,3,1,2))
# for(m in c(2,10,25,50,57,100)){
# plot(pop.size.tosave[,m,1]~seq(1:n.bins), type='l')
# legend('top', legend=paste('year = ',m, sep=''), bty='n')
# abline(h=0,lty=2)
# }
# mtext('position (m/2)',1, outer=T)
# mtext('local population size', 2, outer=T)



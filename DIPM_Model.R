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

setwd('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Dispersal/Chapter I - Dispersal IPM')

Start<-Sys.time()
source('DIPM_Parameters.R') # subsourcing DIPM_DataPrepping
Sys.time()-Start

#!!!!!! WARNING! Some models do not converge


#######################################
#
# ALL PARAMETERS & KERNEL SPECS
#
######################################

# BASIC PARAMETERS

# Survival
surv.int <- coef(summary(surmod))['(Intercept)','Estimate'] #intercept of linear regression
surv.slo <- coef(summary(surmod))['spot.dens','Estimate'] #slope of linear regression
surv.multiplier <- 1 # this value is multiplied to the mean survival probability
# Reproduction
repr.int <- coef(summary(repmod))['(Intercept)','Estimate'] #intercept of linear regression
repr.slo <- coef(summary(repmod))['spot.dens','Estimate'] #slope of linear regression
repr.multiplier <- 1 # this value is multiplied to the mean number of offspring produced
# Breeding Dispersal
disp.stat <- as.numeric(movmod[['estimate']]['sigma_stat']) #variance of stationalry pop fraction
disp.mob <- as.numeric(movmod[['estimate']]['sigma_mob']) #variance of mobile pop fraction
disp.p <- as.numeric(movmod[['estimate']]['p']) # proportion of stationary individuals
# Natal Dispersal
nat.stat <- as.numeric(natmod[['estimate']]['sigma_stat']) #variance of stationalry pop fraction
nat.mob <- as.numeric(natmod[['estimate']]['sigma_mob']) #variance of mobile pop fraction
nat.p <- as.numeric(natmod[['estimate']]['p']) # proportion of stationary individuals

#$$$$$$$$$$$
#
# MODIFIED PARAMETERS
#
#$$$$$$$$$$$

#here
changed.param <-repr.multiplier<-seq(0.5,4,by=0.25)


param.space<-expand.grid(surv.int=surv.int,surv.slo=surv.slo,surv.multiplier=surv.multiplier,
                         repr.int=repr.int,repr.slo=repr.slo,repr.multiplier=repr.multiplier,
                         disp.stat=disp.stat,disp.mob=disp.mob,disp.p=disp.p,
                         nat.stat=nat.stat,nat.mob=nat.mob,nat.p=nat.p)


# General Kernel specifications
# number of bins (how many the intervals to approximate continuity are)
n.bins <- 200
# min and max value of the state variable
# I know it can't be less than 0 or higher than, lets say 110, since it's the position in the stream
L <- 0
U <- 110
# defining meshpoints (separating bins) and midpoints (central point of each bin)
b <- L+c(0:n.bins)*(U-L)/n.bins # meshpoints
y <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # midponts (basically, mean between consecutive meshpoints)
h <- y[2]-y[1] # bin width (in meters)
# starting density vector
# I have 30 individuals in the mid ~3 meters (6 bins), and that is my starting population. To calculate density, I have to calculate how many ind per bin length, and divide that number BY bin length (because my density is individual/meter)
starting.ind.per.bin<-19/(3/h)
pop.size0<-c(rep(0,(n.bins-6)/2),5,5,5,5,5,5,rep(0,(n.bins-6)/2))
density0<-pop.size0/h


#######################################
#
# Specs for Simluation
#
#######################################


# how long to run the simulation for
nsims<-60 # months

# OBJECTS TO SAVE
# 1) Population size
pop.size.tosave<-array(NA,dim=c(n.bins,nsims+1,nrow(param.space))) # creating array to fill wit population size
pop.size.tosave[,1,]<-pop.size0 # filling with first initial densities

# 2) Kernels
K.list<-list()

# 3) Population Spread
pop.spread<-array(NA,dim=c(nsims+1,nrow(param.space)))
pop.spread[1,]<-max(sum(((y-mean(y))^2)*pop.size0)/(sum(pop.size0)-1),0)



##############################################
#
# SIMULATION OVER PARAMETER SPACE
#
##############################################

for(k in 1:nrow(param.space)){
  source('DIPM_Simulation.R')
}

# PLOTTING
colourz<-colorRampPalette(c("red","yellow","springgreen","royalblue"))

# plotting spread through time
matplot(pop.spread, type='l',col=colourz(nrow(param.space)))

# plotting total pop size through time
matplot(colSums(pop.size.tosave), type='l', ylab='Total population size', xlab='months', col=colourz(nrow(param.space)))

# plotting final pop spread aross parameter space
plot(1:nrow(param.space),pop.spread[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
axis(1, at=1:nrow(param.space), labels=changed.param)
# plotting final pop size aross parameter space
plot(1:nrow(param.space),colSums(pop.size.tosave)[nsims+1,], type='l', xaxt='n',xlab='parameter value',ylab='final population spread')
axis(1, at=1:nrow(param.space), labels=changed.param)




# plotting population spatial spread through time
par(mfrow=c(2,3), oma=c(3,4,3,3),mar=c(3,3,1,2))
for(m in c(2,3,4,6,12,24)){
plot(pop.size.tosave[,m]~seq(1:n.bins), type='l', ylim=c(0,max(pop.size.tosave, na.rm=T)+1))
legend('top', legend=paste('month= ',m, sep=''), bty='n')
abline(h=0,lty=2)
}
mtext('position (m/2)',1, outer=T)
mtext('local population size', 2, outer=T)



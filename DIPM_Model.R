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



#######################################
#
# ALL PARAMETERS & KERNEL SPECS
#
######################################


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
# VITAL RATES FUNCTIONS
#
######################################

# survival function
sx <- function(density) {
  # mean probability of survival
  # density is a vector of the SAME length of the bins vector, with local density for the midpoints
  # I can manually calculate the probability of survival per location (mediated by density) given the vector of spot densities
  plogis(surv.int+surv.slo*density)*surv.multiplier
}

# fertility function
fx <- function(density) {
  # mean number of offspring produced:
  # I can manually calculate the probability of survival per location (mediated by density) given the vector of spot densities
  exp(repr.int+repr.slo*density)*repr.multiplier
  #!!! given the high number of orphans, this still has to be corrected for underestimated reproductive output. This can be done by changing the repr.multiplier value
}

# breeding dispersal function
mox <- function(y, x) {
  # for the movement function I have a density INdependent probability distribution of movements, that has to be centered each time on the current position
  # the function calculates the probability of moving to y, if the distribution kernel is centred on x-mu (so on x if mu is 0, as by default)
  disp.p*(1/(sqrt(2*pi*disp.stat^2)))*exp(-((y-x)^2/(2*disp.stat^2))) +
    (1-disp.p)*(1/(sqrt(2*pi*disp.mob^2)))*exp(-((y-x)^2/(2*disp.mob^2)))
}


# natal dispersal function (new recruits position)
nax <- function(y, x) {
  # the function calculates the probability of an foosrping moving to y, if distribution kernel is centred on x-mu (so on x if mu is 0, as by default), and if x-mu is the position of the mother
  nat.p*(1/(sqrt(2*pi*nat.stat^2)))*exp(-((y-x)^2/(2*nat.stat^2))) +
    (1-nat.p)*(1/(sqrt(2*pi*nat.mob^2)))*exp(-((y-x)^2/(2*nat.mob^2)))
}



#######################################
#
# KERNEL ASSEMBLING
#
######################################
# Constructing starting kernel to check everything works (initial density)
check.init.kernels<-function(){
par(mfrow=c(2,2))
# survival kernel
  S <- diag(sx(density0))
  image(t(S), main='survival')
# movement kernel
  G <- outer(y,y, mox)
  image(t(G),main='dispersal')
# reproduction
  R <- diag(fx(density0))
  image(t(R), main='reproduction')
# recruited position (natal dispersal)
  N <- outer(y, y, nax)
  image(t(N),main='natal disp')
  
# P kernel
  Pmatrix <- h*as.matrix(G%*%S)
  image(t(Pmatrix),main='P kernel')
  abline(0,1,lty=2)
# F kernel
  Fmatrix <- h*as.matrix(N%*%R)
  image(t(Fmatrix), main='F kernel')
  abline(0,1,lty=2)

# Whole kernel
  K<- Pmatrix + Fmatrix
  image(t(K), main='full kernel') 
}  
#check.init.kernels()
##############################  

# Function for creating Kernel within simulation
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





#######################################
#
# SIMULATION
#
######################################

# how long to run the simulation for
  nsims<-60 # months

# OBJECTS TO SAVE
# 1) Population size
  pop.size.tosave<-array(NA,dim=c(n.bins,nsims+1)) # creating array to fill wit population size
  pop.size.tosave[,1]<-pop.size0 # filling with first initial densities

# 2) Kernels
  K.list<-list()
  
# 3) Population Spread
  pop.spread<-rep(NA,nsims+1)
  pop.spread[1]<-sum(((y-mean(y))^2)*pop.size0)
  
  
  

#Start <- Sys.time()
# Loop for simulation
for(i in 1:nsims){
  densities<-pop.size.tosave[,i]/h
  K.fun(y=y,densities=densities)
  K.list[[i]]<-K
  new.pop.size<-K%*%pop.size.tosave[,i]
  pop.size.tosave[,i+1]<-new.pop.size
  pop.spread[i]<-
}
#Sys.time() - Start


# plotting total pop size through time
par(mfrow=c(1,1))
plot(colSums(pop.size.tosave), type='l', ylab='Total population size', xlab='months')


# plotting population spread through time
par(mfrow=c(2,3), oma=c(3,4,3,3),mar=c(3,3,1,2))
for(m in c(2,3,4,6,12,24)){
plot(pop.size.tosave[,m]~seq(1:n.bins), type='l', ylim=c(0,max(pop.size.tosave, na.rm=T)+1))
legend('top', legend=paste('month= ',m, sep=''), bty='n')
abline(h=0,lty=2)
}
mtext('position (m/2)',1, outer=T)
mtext('local population size', 2, outer=T)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# PART 3: (dependent on PART 4)
# KERNEL CONSTRUCTION FOR SIMULATION
#
#___________________________________________________________________________________



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
  plogis(param.space$surv.int[k]+param.space$surv.slo[k]*density)*param.space$surv.multiplier[k]
}

# fertility function
fx <- function(density) {
  # mean number of offspring produced:
  # I can manually calculate the probability of survival per location (mediated by density) given the vector of spot densities
  exp(param.space$repr.int[k]+param.space$repr.slo[k]*density)*param.space$repr.multiplier[k]
  #!!! given the high number of orphans, this still has to be corrected for underestimated reproductive output. This can be done by changing the repr.multiplier value
}

# breeding dispersal function
mox <- function(y, x) {
  # for the movement function I have a density INdependent probability distribution of movements, that has to be centered each time on the current position
  # the function calculates the probability of moving to y, if the distribution kernel is centred on x-mu (so on x if mu is 0, as by default)
  param.space$disp.p[k]*(1/(sqrt(2*pi*param.space$disp.stat[k]^2)))*exp(-((y-x)^2/(2*param.space$disp.stat[k]^2))) +
    (1-param.space$disp.p[k])*(1/(sqrt(2*pi*param.space$disp.mob[k]^2)))*exp(-((y-x)^2/(2*param.space$disp.mob[k]^2)))
}


# natal dispersal function (new recruits position)
nax <- function(y, x) {
  # the function calculates the probability of an foosrping moving to y, if distribution kernel is centred on x-mu (so on x if mu is 0, as by default), and if x-mu is the position of the mother
  param.space$nat.p[k]*(1/(sqrt(2*pi*param.space$nat.stat[k]^2)))*exp(-((y-x)^2/(2*param.space$nat.stat[k]^2))) +
    (1-param.space$nat.p[k])*(1/(sqrt(2*pi*param.space$nat.mob[k]^2)))*exp(-((y-x)^2/(2*param.space$nat.mob[k]^2)))
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

#Start <- Sys.time()
# Loop for simulation
for(i in 1:nsims){
  densities<-pop.size.tosave[,i,k]/h
  K<-K.fun(y=y,densities=densities)
  #K.list[[i]]<-K
  new.pop.size<-K%*%pop.size.tosave[,i,k]
  pop.size.tosave[,i+1,k]<-new.pop.size
  pop.spread[i+1,k]<- sum(((y-mean(y))^2)*new.pop.size)/(sum(new.pop.size)-1)
}
#Sys.time() - Start

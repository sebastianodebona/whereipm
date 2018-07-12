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



# Now specified as package functions! sx, fx and gx



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
  S <- diag(sx(s.int=surv.int, s.slo=surv.slo, density=densities)) #surv kernel
  G <- outer(y, y, gx, p.stat=disp.p, s.stat=disp.stat, s.mob=disp.mob, unidirectional=T) # movement kernel
  # P kernel
  Pmatrix <- h*as.matrix(G%*%S)
  
  R <- diag(fx(f.int = repr.int, f.slo=repr.slo, density=densities)) # reproduction
  N <- outer(y, y, gx, p.stat=nat.p, s.stat=nat.stat, s.mob=nat.slo, unidirectional=T) # recruited position (natal dispersal)
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
  pop.spread[i+1,k]<- max(sum(((y-mean(y))^2)*new.pop.size)/(sum(new.pop.size)-1),0)
}
#Sys.time() - Start

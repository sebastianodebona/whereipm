#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# PART 1: DATA PREPPING
# this code modifies the dataset created in PART 0 to get it ready for analyses
#
#___________________________________________________________________________________


setwd('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Dispersal/Chapter I - Dispersal IPM')

# reading dataset 
load('DIPMdata.Rda')
head(DIPMdata)
d<-DIPMdata
# removing some typos in distance.next
d<-d[-which(d$position.next>110),]
# also, adding NA when recruited is 0, to model truncated poisson
# to avoid losing actual number of recruited even when it is 0, I am creating a duplicate called offspring
d$offspring <- d$recruited
d$recruited<-ifelse(d$recruited==0,NA,d$recruited)
head(d)


# I am here creating a variable called SECTION to account for heterogeneity in the stream (to be included as a random effect in the model)
# how long the sections will be is determined by the following parameter, to be adjusted
how.fine<-5
# the following command creates sections of the stream cut each "how.fine" meters
d$section<-cut(d$position, seq(0,max(d$position, na.rm=T), by=how.fine), right=FALSE)


# since the reach labels change every month (different interns or actual morphometric change in the stream), the position may be showing small changes without actual movement. I want to introduce a way to remove small changes due to variations in stream morphology
# I calculate the distance moved, and remove the negligible small scale movements likely due to changes in how the stream is sampled
# dist.moved is therefore defined as the actualy distance moved if this is not smaller than a certain value, either absolute or relative to the reach length
neglect <- 0
# as an alternative, the threashold could be a fourth of the reach legth, meaning that I will ignore a 2 meters movement in a 8 meters reach, but consider a 1 meter movement in a 3 meter reach
#neglect <- d$reach.l/4
d$dist.moved<-ifelse(abs(d$position.next-d$position)>neglect,
                     d$position.next-d$position,
                     0)
# the variable moved corresponds to a binomial where 0 is stay and 1 is moved
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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Density-dependent dispersal  during guppies colonisation of low-predation habitats
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# PART 0: DATA
# this code organises the dataset to be used in the model construction of a Dispersal IPM
#
#___________________________________________________________________________________

# clearing env
rm(list=ls(all=TRUE))

# setting working directory and reading datafile
setwd('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Dispersal/Chapter I - Dispersal IPM')
pedigree<-read.csv2('Pedigree.csv',header=T)
# this file contains pedigree information by individual, and uses as a unit 'form' from the field (both marking forms and recapture forms). Useful information derives from marking forms only, when the information on paternity/maternity is first input.  I need to get rid of recap forms
# a quick and dirty way is to only include lines which current form is equal to the cohort (meaning that fish was a new mark) 
pedigree<-pedigree[which(pedigree$Form==pedigree$Cohort),]


# I can extract information in the locality using FishID. The second and third character stand for stream
pedigree$locality<-factor(substr(pedigree$FishID,2,3))
# sex is the first letter
pedigree$sex<-factor(substr(pedigree$FishID,1,1))
# and I want only females here since I am constructing a one sex IPM
pedigree<-subset(pedigree, sex=='F')

# not all individuals have maternity assigned. I can check what percentage. This will later be used to correct for offspring number
table(pedigree$Maternity=='')
p.orphans<-mean(pedigree$Maternity=='') #this represents the probability of being 'orphan' (no maternity assigned)


# I need to match the time specification of the MR dataset, which is by 'sampling'. Given that sampling starts from March 2008 and proceeds monthly, and that the form specification is the form YYMM, I can sort the Form here by numerical order, and assign increasing integers as sampling. 
sampling <- sort(unique(pedigree$Form))
# April 2008 is for some reason missing, so I need to add it manually to avoid the messing up of the sampling order
sampling <- c(803,804,sort(unique(pedigree$Form))[-1])
# going from March 2008 to April 2013 the number of sampling here should be 10+12*4+4=62
# --> checking this
length(sampling)==10+12*4+4

# now the position on the sorted vector corresponds to the sampling month. I can use match to create a sampling column
pedigree$sampling <- match(pedigree$Form, sampling)
# I can do the same for cohort, first and last capture
pedigree$cohort_samp <- match(pedigree$Cohort, sampling)
pedigree$last_samp <- match(pedigree$Last_Capture, sampling)
pedigree$first_samp <- match(pedigree$First_capture, sampling)


# Now read MR datafile
mrdata <- read.csv("/Users/sebastianodebona/Dropbox/GMR/Mark-recapture Data/2016-01d/Files For Use/GuppyData.csv", header=TRUE)
# subsetting for the only stream with pedigree information, and for females only, to make the IPMization easier
mrdata<-subset(mrdata,stream=='LL' & sex=='F')
str(mrdata)


# I want to use MR data only up to when pedigree was calculated (or I will artificially have a bunch of females with 0 recruitment since their daughters haven't been genotyped yet)
#!!!!!!!!!!!!!!!!CHECK THIS STEP
# it is strictly < max(sampling)-1 (and not <=) because I want the pedigree to overextend two months so that I can inlude captured daugthers two months after they were produced
#!!!!!!!!!!!!!!!!
mrdata <- subset(mrdata, sampling<max(pedigree$sampling, na.rm=T)-1)
length(table(mrdata$sampling))

# To make things easier later, I want to assign location at first capture and first size (weight) to all individuals in pedigree
pedigree<-merge(pedigree,mrdata[,c('FishID','sampling','distance','weight')], by=c(as.character('FishID'),'sampling'), all.x=T)
# I am renaming the variable 'distance' as 'postion'
mrdata$position<-mrdata$distance

# RECRUITMENT
# I now have to find a way to assign daughtes produced monthly to each female
mrdata$recruited<-rep(NA, nrow(mrdata))
#!!!!!!!!!!!!!!!!CHECK THIS STEP
# I created a loop that goes through the dataset and assigns to each row (individual) the number of individuals that have that fish as a mother, and were first captured 3 months after the mum was. This is to assign recruitment to the next step (time t+1), and accounts for indivuals being captured only 2 months after being born
#!!!!!!!!!!!!!!!!CHECK THIS STEP
for(i in 1:nrow(mrdata)){
  mum <- as.character(mrdata$FishID[i])
  when <- mrdata$sampling[i]+2
  mrdata$recruited[i] <- nrow(subset(pedigree, as.character(Maternity)==mum & cohort_samp==when))
}


# REPRODUCTION (binary)
# assign reproduced or not (whether recruit was non-zero)
mrdata$repr <- as.numeric(mrdata$recruited!=0)


# OFFSPRING POSITION (to model inheritance)
# I now have to assign the mean position of the offspring to the mother
# this again goes throught the dataseta and assigns mean offspring position (when many), or NA (when o offspring)
mrdata$position_off<-rep(NA, nrow(mrdata))
# mean offspring position
for(i in 1:nrow(mrdata)){
  mum <- as.character(mrdata$FishID[i])
  when <- mrdata$sampling[i]+3
  mrdata$position_off[i] <- ifelse(
    nrow(subset(pedigree, as.character(Maternity)==mum & cohort_samp==when))!=0,
    mean(subset(pedigree, as.character(Maternity)==mum & cohort_samp==when)[,'distance'], na.rm=T),
    NA)
}

# since the first position of the offspring is calculated two months after the offspring is born, I will calculate the half distance as (mother position - offsrping position two months later)/2
mrdata$halfpos_off<-rep(NA, nrow(mrdata))
# mean offspring position
for(i in 1:nrow(mrdata)){
  mrdata$halfpos_off[i] <- mrdata$position[i]+(mrdata$position_off[i]-mrdata$position[i])/2
}



# MOVEMENT
# next time mother's position
mrdata$position.next<-rep(NA,nrow(mrdata))
for(i in 1:nrow(mrdata)){
  mum <- as.character(mrdata$FishID[i])
  nextt <- mrdata$sampling[i]+1
  mrdata$position.next[i] <- ifelse(
    length(which(mrdata$FishID==mum & mrdata$sampling==nextt))==0,
    NA,
    mrdata$position[which(mrdata$FishID==mum & mrdata$sampling==nextt)])
}
table(is.na(mrdata$position.next))
#not very many fish have position next....(1/5)


# SURVIVAL
# I can add a column of individual survival to time t+1. I can make a very rough approximation, and given the very high recpature rate, say that if I never again capture an invidual, it means it's dead
mrdata$surv<-rep(NA,nrow(mrdata))
for(i in 1:nrow(mrdata)){
  mrdata$surv[i]<-ifelse(
    nrow(subset(mrdata, FishID==FishID[i] & sampling>sampling[i]))==0,
    0,
    1)
}


# POP SIZE at each timestep
mrdata$pop_size<-rep(NA, nrow(mrdata))
for(i in 1:nrow(mrdata)){
  mrdata$pop_size[i]<-nrow(subset(mrdata,sampling==sampling[i]))
}
toplot<-as.data.frame(cbind(sampling=mrdata$sampling,pop=mrdata$pop_size)[order(mrdata$sampling),])
plot(toplot$pop~toplot$sampling, type='l')



# POP SIZE at each time point, per location
mrdata$spot.pop<-rep(NA, nrow(mrdata))
for(i in 1:nrow(mrdata)){
  mrdata$spot.pop[i]<-nrow(subset(mrdata,sampling==sampling[i] & position==position[i])) #pop size here&now
}
range(mrdata$spot.pop)



# before the next step, some weird stuff to remove
mrdata<-mrdata[-c(which((mrdata$downstream-mrdata$upstream)>20), which((mrdata$downstream-mrdata$upstream)<=0),which(is.na(mrdata$downstream-mrdata$upstream))),]


# REACH LENGHT useful to calculate densities, but also to establish real movement
mrdata$reach.l <- mrdata$downstream-mrdata$upstream
  
# POP DENSITY at each time point, per location
mrdata$spot.dens<-mrdata$spot.pop/mrdata$reach.l #pop size here&now/reach length


# HABITAT TYPE
mrdata$habitat<-droplevels(mrdata$Habitat)


DIPMdata<-mrdata[c('FishID','sampling','surv','position','position.next','repr','recruited','position_off','halfpos_off','spot.pop','spot.dens','pop_size','habitat','reach.l','location','upstream','downstream')]
head(DIPMdata)
save(DIPMdata, file='DIPMdata.Rda')


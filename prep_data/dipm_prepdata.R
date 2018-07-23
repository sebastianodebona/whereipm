#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Preparing a polished dataset from raw data
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# getting the pedigree datafile
pedigree <- read.csv2('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Raw Data/Pedigree.csv',header=T)
str(pedigree)

# this file is structured according to marking and recapture forms from the field. The information is redundant.
# To solve this, I will only consider the marking form, reporting individuals at their first capture, and containing the informations on parents.
pedigree <- pedigree[which(pedigree$Form==pedigree$Cohort),]

range(pedigree$Form) #checking how far in time this goes
# I need to match the time specification of the MR dataset, which is by 'sampling'. Given that sampling starts from March 2008 and proceeds monthly, and that the form specification is the form YYMM, I can sort the Form here by numerical order, and assign increasing integers as sampling. 
sampling <- sort(unique(pedigree$Form))

# April 2008 is missing from pedigree file (no new borns big enough after only one month from intoduction), so I need to add it manually to avoid the messing up of the sampling order
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


# reading the GMR datafile
gmr <- read.csv("/Users/sebastianodebona/Documents/Research/Projects/Guppies/GuppyData.csv", header=TRUE)


#### SURVIVAL ###############################################
# This step will come in handy later
# I can add a column of individual survival to time t+1. I can make a very rough approximation, and given the very high recapture rate, say that if I never again capture an invidual, it means it's dead
# defining when each individual was last captured (and never after that)
# extracting information on the last capture for each individual
last_capture <- plyr::ddply(gmr, c('FishID'), plyr::summarise, last_capture = max(sampling, na.rm=T))
gmr <- merge(gmr, last_capture, by='FishID')
# checking whether the individual is captured again after every capture instance. if not, considered dead.
gmr$survived <- ifelse(gmr$sampling==gmr$last_capture, 0, 1)
# checking if the already existing survival variable checks with the one I just created
table(gmr$survival==gmr$survived)
# not always...


# I will analyse reproduction (and the effect of local vs global density on it). therefore, I need to focus only on the piece of mark-recapture data for which pedigree analyses are available (only LL, only up to April 2013)
data <- subset(gmr, stream=='LL' & sampling<=max(pedigree$sampling))
str(data)
# checking again survival == survived
table(data$survival==data$survived)
# at least it is working for the time I am interested in


# Below I will calculate the variables I'm interested in: survival (done before), reproduction and dispersal, but first...

#### MATERNITY AND PATERNITY #####################################
# Adding a column for each individuals father and mother
pedigree.Table <- subset(pedigree, select=c('FishID', 'Paternity', 'Maternity'))
# filling empty rows with NA
pedigree.Table$Maternity <- ifelse(pedigree.Table$Maternity=='', NA, as.character(pedigree.Table$Maternity))
pedigree.Table$Paternity <- ifelse(pedigree.Table$Paternity=='', NA, as.character(pedigree.Table$Paternity))
data<-merge(data, pedigree.Table, by='FishID')


#### REPRODUCTION ###############################################
# trying to derive how many offspring are produced by each individual in a month
recruited <- plyr::ddply(data, c('Maternity', 'sampling'), nrow)
recruited <- recruited[recruited$Maternity!='',]
# renaming V1
colnames(recruited)[ncol(recruited)] <- 'recruited'
# merging this information back into the main dataset
data <- merge(data, recruited, by.x = c('FishID', 'sampling'), by.y = c('Maternity', 'sampling'), all.x = T)
table(data$recruited)


#### LOCAL POP SIZE #############################################
# at each time point, per location 
local.pop.size <- plyr::ddply(data, c('sampling', 'location'), nrow)
colnames(local.pop.size)[ncol(local.pop.size)] <- 'local.pop'
# merging into main data
data <- merge(data, local.pop.size, by=c('sampling', 'location'), all.x=T)

# REACH LENGHT useful to calculate densities, but also to establish real movement
table(data$downstream-data$upstream<=0)
# there are a few negative lengths, and only few zeros. I will assign NA to those
data$reach_l <- ifelse(data$downstream-data$upstream>0,data$downstream-data$upstream,NA)
# I will also remove clearly wrong values, above the total length of the stream
data$reach_l <- ifelse(data$reach_l<110, data$reach_l, NA)

# LOCAL POPULATION DENSITY, at each time point (and each location)
data$loc_density <- data$local.pop/data$reach_l

# VARIANCE IN LOCAL POP DENSITY
var_density <- tapply(data$loc_density, data$sampling, var, na.rm=T)
var_density <- cbind(sampling=1:62, density_var=var_density)
# merging
data <- merge(data, var_density, by = c('sampling'), all.x = T)


#### MOVEMENT #################################################
# next time position
start <- Sys.time()
data$distance.next<-rep(NA,nrow(data)) #first, creating an empty column
for(i in 1:nrow(data)){
  data$distance.next[i]<-ifelse(
    length(which(data$FishID==data$FishID[i] & data$sampling==data$sampling[i]+1))==0,
    NA,
    data$distance[which(data$FishID==data$FishID[i] & data$sampling==data$sampling[i]+1)]
  )
}
table(is.na(data$distance.next))
#not very many fish have position next....(1/4)
data$movement <- data$distance.next - data$distance
Sys.time() - start

#### NATAL DISPERSAL #############################################
# in order to estimate natal dispersal I'll take the offspring perspective. I need to find the first capture per each individual. And then check the mother's position two months before first capture (likely the month of birth for this fish). I may also remove all fish above a threshold length, to exclude "old" first captured individuals (which might have be born many months before)
# the first capture corresponds to the lowest sampling value in the mark recapture file
first_capture <- plyr::ddply(data, c('FishID'), plyr::summarize ,first_capture = min(sampling, na.rm=T))
# what I created is a vector of first captures per each individual. I can now use 'match' to associate a first capture value to all fish
data <- merge(data, first_capture, by = 'FishID', all.x = T)

# Now I can check where each individual's mother was two months before the first capture. I will only add the mother's position once, in the row corresponding to the first capture.
# One more step is needed to be able to compare adult and natal dispersal. Adult dispersal is calculated in a projective manner, finding the position at the next time step, but from the 'now' point of view (habitat, local and global density NOW). I then need to find those variables from the inset of natal dispersal (i.e. from the data belonging to the mother two months before first capture)

# Extracting matching New recuit ID and Maternity information
birth.info <- plyr::ddply(subset(data, sampling==first_capture), c('Maternity', 'sampling', 'FishID'), nrow)
# making sure there is not redundant information, and if there is, removing it
table(birth.info$V1)
birth.info <- subset(birth.info, V1==1, select=c('Maternity', 'sampling', 'FishID'))
# Calculating extimated time of birth (two months before first capture)
birth.info$birth <- as.numeric(birth.info$sampling -2)
# and time mid (between estimated birth and first capture)
birth.info$mid <- as.numeric(birth.info$sampling - 1)
# Extracting information on position, location, habitat, reach length and local density of pool in which mother was captured. 
birth.info <- merge(birth.info, subset(data, select = c('FishID', 'sampling', 'distance', 'Habitat', 'location', 'reach_l',  'loc_density')), by.x = c('Maternity', 'birth'), by.y = c('FishID', 'sampling'))
# Renaming
colnames(birth.info)[which(colnames(birth.info)%in% c('birth', 'distance', 'Habitat', 'location', 'reach_l', 'loc_density'))] <- c('sampling.birth', 'distance.birth', 'habitat.birth', 'location.birth', 'reach_l.birth', 'loc_density.birth')
# Now extracting same information for mid month
birth.info <- merge(birth.info, subset(data, select = c('FishID', 'sampling', 'distance', 'Habitat', 'location', 'reach_l',  'loc_density')), by.x = c('Maternity', 'mid'), by.y = c('FishID', 'sampling'))
colnames(birth.info)[which(colnames(birth.info)%in% c('mid', 'distance', 'Habitat', 'location', 'reach_l', 'loc_density'))] <- c('sampling.mid', 'distance.mid', 'habitat.mid', 'location.mid', 'reach_l.mid', 'loc_density.mid')
# now merging back into main dataset
data <- merge(data, birth.info, by = c('FishID', 'sampling'), all.x = T)
# for some reasons the number of rows increased :/

# for someone to have a position at birth the have to have defined maternity and for the mother to be caught two months before the first capture of the offspring
table(!is.na(data$distance.birth))
# even less have a position at birth
data$natal_movement <- data$distance-data$distance.birth

# Purging some unrealistic distance specification before saving dataset
table(data$distance>110)
data<-subset(data, distance<=110)

#### REACH ################################################
# As an alternative to section, I have manually assigned reach level random effect, grouping up reaches that changed slightly through time and chunks of the stream which appear to be self-consistent
reach_key<-read.table('/Users/sebastianodebona/Documents/Research/Projects/Guppies/Raw Data/reach_key.txt', header = T)
data$reach <- as.factor(reach_key$newID[match(data$location, reach_key$oldID)])
data$reach.birth <- as.factor(reach_key$newID[match(data$location.birth, reach_key$oldID)])


##### TIME ##############################################
data$time<-as.numeric(data$sampling)-1 # so that introduction becomes t0
data$time.birth<-as.numeric(data$sampling.birth)-1
data$sampling<-as.factor(data$sampling)
data$sampling.birth<-as.factor(data$sampling.birth)


##### SCALING ##############################################
# to improve model convergence, I will take some steps following 'https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html'
# 1) SCALING/CENTERING local and global density in the model
# 2) INCREASE ITERATIONS (later)
# I will scale the variables but keep the original in the data frame
  
# what needs to be scaled is not the local density vector in the data.frame, since values of local density are repeated once per each individual present in each pool and experiencing the same density. I need to scale it based on a vector of unique values of local densities found at each sampling and location
loc <- plyr::ddply(data, c('sampling', 'location', 'loc_density'), nrow)
loc$s_loc_density <- scale(loc$loc_density)
data <- merge(data, loc[,c("sampling", "location", "s_loc_density")], by = c('sampling', 'location'))
# now transferring attributes (scaling and centering value) from loc$s_loc_density to data$s_loc_density, so they can be retrieved later
# attr(x, 'string') <- value assignes a value as an attribute to the object (e.g. vector) x, to be called 'string':
# see: https://stat.ethz.ch/R-manual/R-devel/library/base/html/attr.html
# the function sticky used here, from the package sticky, allows vector to retain the attributes even after subsetting, appending etc
# see: https://cran.r-project.org/web/packages/sticky/vignettes/introduction.html
data$s_loc_density <- sticky::sticky(data$s_loc_density)
attr(data$s_loc_density, "scaled:center") <- attr(loc$s_loc_density, "scaled:center")
attr(data$s_loc_density, "scaled:scale") <- attr(loc$s_loc_density, "scaled:scale")
   
# time will be scaled but not centered, to keep the intercept at 0 in the models
data$s_time<-scale(data$time, center=F)
data$s_time.birth<-scale(data$time.birth, center=F)


##### HABITAT ##############################################
# removing weird habitat types (probably typos)
data$habitat<-ifelse(data$Habitat%in%c('S','B','N','P','R'), as.character(data$Habitat), NA)
data$habitat<-ifelse(data$habitat=='R', 'R', 'P') # simplifying and calling all non-riffle habitats pools
data$habitat.birth<-ifelse(data$habitat.birth%in%c('S','B','N','P','R'), as.character(data$habitat.birth), NA)
data$habitat.birth<-ifelse(data$habitat.birth=='R', 'R', 'P')


##### REPRODUCTION ##############################################
# I will model offspring number as a binomial + a zero truncated count distribution
# fist considering the probability of producing any offspring
data$reproduced<-ifelse(data$recruited>0,1,data$recruited) #defining reproduced variable (binary outcome)
data$offspring<-ifelse(data$recruited>0,data$recruited,NA)


##### MOVEMENT ##############################################
# I will here model dispersal as a  probability of moving or not (binomial), combined with the direction moved and the absolute distance moved
# an individual cannot have moved more than the total length of the stream (110m), up or downstream
data$movement<-ifelse(abs(data$movement)>110,NA,data$movement)
# given month to month discrepancies in drawing the tape-meter in the field, most fish appeared to have moved very short distances (they actually did not move, but the centerpoint of the reach they live in wobbled), inflating the probability of moving.
# I am here stating here that the movement has to be at least as much as half the reach length.
# the probability of moving or not this is 0 if an individual has a recorded position and position.next, and has not moved more than the cutoff (apparent movement), and 1 if it has moved
data$p_moving<-ifelse(abs(data$movement)>data$reach_l/2,1,0)
# then, absolute distance moved
data$dist_moved<-ifelse(data$p_moving==1,   # if it moved
                        abs(data$movement), # absolute distance
                        NA)                 # NA otherwise
#hist(data$dist_moved, breaks=200)
# finally, direction moved
data$dir_moved<-ifelse(data$p_moving==1, # nested ifelse, to account for distance>cutoff and 
                       ifelse(data$movement>0,1,0), # movement positive
                       NA)


##### NATAL MOVEMENT ##############################################
# same purging as for movement
data$natal_movement<-ifelse(abs(data$natal_movement)>110,NA,data$natal_movement)
# defining the probability of moving form birth
data$p_natal<-ifelse(abs(data$natal_movement)>data$reach_l.birth/2,1,0)
# then, absolute distance moved from birthplace
data$dist_natal<-ifelse(data$p_natal==1,
                        abs(data$natal_movement),
                        NA)
data$dir_natal<-ifelse(data$p_natal==1, # nested ifelse, to account for distance>cutoff and
                       ifelse(data$natal_movement>0,1,0), # movement positive
                       NA)


save(data, file = 'data/data.Rda')

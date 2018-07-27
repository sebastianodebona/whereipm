#' Clear NAs form dataset
#' 
#' Removes rows containing NAs from a dataset before the data is fit to a model
#' 
#' @export
#' 
#' @param MyData Dataset (\code{data.frame}) to be cleansed
#' @param response Character string defining the response variable in the model
#' @param fixEff,RanEff Character vectors containing the fixed and random effects to be used in the model. By default: \code{fixEff=c('loc_density', 'SL')} and code\{ranEff=c('sampling', 'reach', 'FishID')}
#' @param whichFactors Character vector defining which variables have to be converted into factors. By default, this is equal to the random effects.

rm.NAs <- function(myData, 
                   response, 
                   fixEff=c('loc_density', 'SL'), 
                   ranEff=c('sampling', 'reach', 'FishID'), 
                   whichFactors=ranEff){

  myData <- myData[complete.cases(myData[,c(response, fixEff, ranEff)]),]
  for(i in 1:length(whichFactors)){
    myData[,whichFactors[i]]<-as.factor(myData[,whichFactors[i]])
  }
  return(myData)
}
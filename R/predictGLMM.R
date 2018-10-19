#' Predict GLMM
#' 
#' Calculate predicted values from fitted glmm
#' 
#' @export
#' 
#' @param x Vector of values for independent variable to be used to predict response variable.
#' @param Coefs Table of coefficients (extracted as e.g. \code{summary(model)$coefficients}) from a fitted GLMM.
#' @param scalingValue If x contains natural values but the model is fitted with scaled values, this parameter un-scales the coefficients
#' @param centeringValue If x contains natural values but the model is fitted with centered values, this parameter un-centers the coefficients.
#' @param covariateString A string identifying the covariate on which the response variable depends (used by \code{grepl}). Default is local density \code{"loc_density"}
#' @param link Link function used by fitted model: \code{"logit"}, \code{"log"}, \{"truncated-pois"} \code{"n"}.
#' 

predictGLMM <- function(x, Coefs, 
                        scalingValue = 1, 
                        centeringValue = 0,
                        covariateString = "loc_density",
                        link){
  
  # extracting model intercept
  alpha = Coefs['(Intercept)','Estimate']
    
  # establishing whether the covariate is in the model (as fixed effect) or not
  if(sum(grepl(covariateString, rownames(Coefs))) == 0){
    
    beta = 0
  
  } else {
    
    # taking the first estimate which checks with the pattern covariateString (main effect)
    # this averages out possible other covariates
    beta = Coefs[grep(covariateString, rownames(Coefs))[1],'Estimate'] 
    
  }
  
  # now evaluating link type
  if(link == "logit"){
    
    y = plogis(alpha + beta*((x-centeringValue)/scalingValue))
    
  } else if (link == "log"){
    
    y = exp(alpha + beta*((x-centeringValue)/scalingValue))
    
  } else if (link == "n"){
    
    y = alpha + beta*((x-centeringValue)/scalingValue)

  } else if (link == "truncated-pois"){
    
    y = untruncZTP(alpha + beta*((x-centeringValue)/scalingValue))
    
  }

  return(y)
}
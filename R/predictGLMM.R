#' Predict GLMM
#' 
#' Calculate predicted values from fitted glmm
#' 
#' @export
#' 
#' @param x Vector of values for independent variable to be used to predict response variable.
#' @param model Fitted GLMM.
#' @param scalingValue If x contains natural values but the model is fitted with scaled values, this parameter un-scales the coefficients
#' @param centeringValue If x contains natural values but the model is fitted with centered values, this parameter un-centers the coefficients.
#' @param covariateString A string identifying the covariate on which the response variable depends (used by \code{grepl}). Default is local density \code{"loc_density"}
#' @param link Link function used by fitted model: \code{"logit"}, \code{"log"}, \code{"n"}.
#' 

predictGLMM <- function(x, model, 
                        scalingValue = 1, 
                        centeringValue = 0,
                        covariateString = "loc_density",
                        link){
  
  # extracting model intercept
  alpha = summary(model)$coefficients['(Intercept)','Estimate']
    
  # establishing whether the covariate is in the model (as fixed effect) or not
  if(sum(grepl(covariateString, rownames(summary(model)$coefficients))) = 0){
    
    beta = 0
  
  } else {
    
    # taking the first estimate which checks with the pattern covariateString (main effect)
    # this averages out possible other covariates
    beta = summary(model)$coefficients[grep(covariateString, 
                                            rownames(summary(model)$coefficients)))[1],'Estimate'] 
    
  }
  
  # now evaluating link type
  if(link == "logit"){
    
    y = plogis(alpha + beta*x)
    
  } else if (link == "log"){
    
    y = exp(alpha + beta*x)
    
  } else if (link == "n")
    
    y = alpha + beta*x
  
}
---
title: "Perturbation analyses - Phases"
author: "Seba De Bona"
date: "2019-05-18"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: lumen
---

## Rationale and Setup

Here I am performing perturbation analyses to estimate the relative effect of different parameters (and therefore biological processes, such as survival, reproduction, dispersal and their density-dependence) on the spread of the population.

To perform a perturbation analyses on a non deterministic model (which is the case of our IPMs, given density-dependence is included) the parameters are adjusted by a small value, $\epsilon$, one at the time. The model is then run again, with all other paramters unperturbed, to obtain some output variables (e.g. long term population growth rate or speed of spread, in our case). A comparison between these output variables and those originated by the unperturbed model give a measure of the effect of the perturbed variable on the overall population dynamics. For more detail, see chapter XX of Ellner, Reese, Childs 2016. Data-driven modelling of structured populations.

Below, I load packages, model outputs and stored vital parameters.


```r
library(whereipm)
library(tidyverse)
library(furrr)

# loading model outputs (also containing vital params)
phasesIPMs <- readRDS(file.path(here::here(), "results", "phaseModelsIPMs.rds"))
```

And here I allocate space to store the results of the perturbation analyses, by creating an array.


```r
# here I create a list of IPM for different phases
phases <- phasesIPMs$IPMs
names(phases) <- unlist(phasesIPMs$phase)

# here I create a vector of paramters to be perturbed
params <- c('s.int', 's.slo',
            'r.int',  'r.slo', 
            'f.int',  'f.slo',
            'p.stat.int.ad','p.stat.slo.ad', 
            'mu.mob.int.ad', 'sigma.mob.ad',
            'p.stat.int.nat','p.stat.slo.nat', 'p.stat.sloLand.nat', 'p.stat.sloInt.nat',
            'mu.mob.int.nat', 'sigma.mob.nat')

# and a list of output variables
# here I'll have a bunch of variables to evaluate the effect of perturbation on
outputs <- c('final_front_position',
            'final_pioneer_position',
            'initial_front_speed',
            'initial_pioneer_speed',
            'final_front_speed',
            'final_pioneer_speed',
            'avg_front_speed',
            'avg_pioneer_speed',
            'acceleration_tot_front',
            'acceleration_tot_pioneer',
            'acceleration_est_front',
            'acceleration_est_pioneer',
            'lambda_avg',
            'final_pop_size',
            'final_max_density')
# these variables should be in the same order the `perturb` function will provide as output.

# the perturbation array will have output variables as rows,
# perturbed parameters as columns and will store results from different 
# IPM as a third dimension of the array.

perturbationResults <- array(NA, dim = c(length(outputs),
                                         length(params),
                                         length(phases)),
                             dimnames = list(outputs,
                                             params,
                                             names(phases)))
str(perturbationResults)
```

```
##  logi [1:15, 1:16, 1:3] NA NA NA NA NA NA ...
##  - attr(*, "dimnames")=List of 3
##   ..$ : chr [1:15] "final_front_position" "final_pioneer_position" "initial_front_speed" "initial_pioneer_speed" ...
##   ..$ : chr [1:16] "s.int" "s.slo" "r.int" "r.slo" ...
##   ..$ : chr [1:3] "early" "mid" "late"
```

```r
#this is how the array looks like per each phase
perturbationResults[,,1]
```

```
##                          s.int s.slo r.int r.slo f.int f.slo p.stat.int.ad
## final_front_position        NA    NA    NA    NA    NA    NA            NA
## final_pioneer_position      NA    NA    NA    NA    NA    NA            NA
## initial_front_speed         NA    NA    NA    NA    NA    NA            NA
## initial_pioneer_speed       NA    NA    NA    NA    NA    NA            NA
## final_front_speed           NA    NA    NA    NA    NA    NA            NA
## final_pioneer_speed         NA    NA    NA    NA    NA    NA            NA
## avg_front_speed             NA    NA    NA    NA    NA    NA            NA
## avg_pioneer_speed           NA    NA    NA    NA    NA    NA            NA
## acceleration_tot_front      NA    NA    NA    NA    NA    NA            NA
## acceleration_tot_pioneer    NA    NA    NA    NA    NA    NA            NA
## acceleration_est_front      NA    NA    NA    NA    NA    NA            NA
## acceleration_est_pioneer    NA    NA    NA    NA    NA    NA            NA
## lambda_avg                  NA    NA    NA    NA    NA    NA            NA
## final_pop_size              NA    NA    NA    NA    NA    NA            NA
## final_max_density           NA    NA    NA    NA    NA    NA            NA
##                          p.stat.slo.ad mu.mob.int.ad sigma.mob.ad
## final_front_position                NA            NA           NA
## final_pioneer_position              NA            NA           NA
## initial_front_speed                 NA            NA           NA
## initial_pioneer_speed               NA            NA           NA
## final_front_speed                   NA            NA           NA
## final_pioneer_speed                 NA            NA           NA
## avg_front_speed                     NA            NA           NA
## avg_pioneer_speed                   NA            NA           NA
## acceleration_tot_front              NA            NA           NA
## acceleration_tot_pioneer            NA            NA           NA
## acceleration_est_front              NA            NA           NA
## acceleration_est_pioneer            NA            NA           NA
## lambda_avg                          NA            NA           NA
## final_pop_size                      NA            NA           NA
## final_max_density                   NA            NA           NA
##                          p.stat.int.nat p.stat.slo.nat p.stat.sloLand.nat
## final_front_position                 NA             NA                 NA
## final_pioneer_position               NA             NA                 NA
## initial_front_speed                  NA             NA                 NA
## initial_pioneer_speed                NA             NA                 NA
## final_front_speed                    NA             NA                 NA
## final_pioneer_speed                  NA             NA                 NA
## avg_front_speed                      NA             NA                 NA
## avg_pioneer_speed                    NA             NA                 NA
## acceleration_tot_front               NA             NA                 NA
## acceleration_tot_pioneer             NA             NA                 NA
## acceleration_est_front               NA             NA                 NA
## acceleration_est_pioneer             NA             NA                 NA
## lambda_avg                           NA             NA                 NA
## final_pop_size                       NA             NA                 NA
## final_max_density                    NA             NA                 NA
##                          p.stat.sloInt.nat mu.mob.int.nat sigma.mob.nat
## final_front_position                    NA             NA            NA
## final_pioneer_position                  NA             NA            NA
## initial_front_speed                     NA             NA            NA
## initial_pioneer_speed                   NA             NA            NA
## final_front_speed                       NA             NA            NA
## final_pioneer_speed                     NA             NA            NA
## avg_front_speed                         NA             NA            NA
## avg_pioneer_speed                       NA             NA            NA
## acceleration_tot_front                  NA             NA            NA
## acceleration_tot_pioneer                NA             NA            NA
## acceleration_est_front                  NA             NA            NA
## acceleration_est_pioneer                NA             NA            NA
## lambda_avg                              NA             NA            NA
## final_pop_size                          NA             NA            NA
## final_max_density                       NA             NA            NA
```

I will save each of the perturbed IPM in case they are needed later on. So here I will expand on the `phasesIPM` tibble, adding one column per perturbed parameter.


```r
perturbed_phasesIPMs <- phasesIPMs
for(i in params){
  perturbed_phasesIPMs[[i]] <- as.list(rep(NA, length(phases)))
}
```

Finally, I am picking some colours for later plotting.


```r
Cols <- wesanderson::wes_palette("Darjeeling1", 3, type = "discrete")
```



## Perturbations

I am now ready to run the perturbation analyses. For this, I will use the custom made function `perturb`. The following analyses can be time-consuming. Here I also estimate how long it would take to run them.


```r
# to have an idea of how long it will take, we can add up the time it takes
# to run each IPM for a number of time equal to the number of parameters.

phasesIPMs %>%
  select(timetorun) %>%
  unlist() %>%
  sum() %>%
  prod(length(params)) %>%
  round() %>%
  paste(attr(phasesIPMs$timetorun[[1]], "units"),
        "estimated total run time for Perturbation analyses", sep = " ")
```

```
## [1] "954 mins estimated total run time for Perturbation analyses"
```

```r
# It is worth adding some runtime for the loops themselves, 
# so this estimate is on the low range.
```

Being warned of the possibly long runtime, I can move on to the simulations.
The setup allows to run the process on multiple cores, to shorten the runtime.


```r
# creating a tibble of indices
start <- Sys.time()
run <- expand.grid(I = 1:length(phases),
                   j = 1:length(params))

plan(multisession)
run <- run %>%
  mutate(temp = furrr::future_map2(I, j, 
                            function(phase.index, param.index){
                              perturb(orig.IPM = phases[[phase.index]],
                                      epsilon = 0.01,
                                      par = params[param.index],
                                      outputs = outputs)
                            }))
```

```
## Warning in if (tilde_msg != "no issues") {: the condition has length > 1
## and only the first element will be used

## Warning in if (tilde_msg != "no issues") {: the condition has length > 1
## and only the first element will be used

## Warning in if (tilde_msg != "no issues") {: the condition has length > 1
## and only the first element will be used

## Warning in if (tilde_msg != "no issues") {: the condition has length > 1
## and only the first element will be used
```

```r
invisible(pmap(list(run$I, run$j, run$temp), 
     function(phase.index, param.index, object){
       perturbationResults[, params[param.index], phase.index] <<- unlist(object$delta)
       perturbed_phasesIPMs[[params[param.index]]][[phase.index]] <<- object$newIPM
     }))


# start <- Sys.time()
# for(I in 1:length(phases)){
#   for(j in 1:length(params)){
#     temp <- perturb(orig.IPM = phases[[I]],
#                     epsilon = 0.01,
#                     par = params[j],
#                     outputs = outputs)
#     perturbationResults[, params[[j]], I] <- unlist(temp$deltas)
#     perturbed_phasesIPMs[I , params[j]] <- list(list(temp$newIPM))
#   }
# }

# how long it actually took
Sys.time() - start
```

```
## Time difference of 9.000225 hours
```

Note that in the `perturbationResults` array we have stored *relative* changes in the output parameters.


```r
# obtaining the original values for output variables:
phasesIPMs <- phasesIPMs %>%
  mutate(
    orig.output =   map(IPMs, function(x){
      outs <- lapply(outputs, function(i){x[[i]]})
      names(outs) = outputs
      return(outs)
      }
      ))
```


```r
# extracting original parameters and reading them for plotting
orig_outputs <- phasesIPMs %>%
  select(orig.output) %>%
  pull() %>%
  unlist() %>%
  matrix(ncol = length(phases), nrow = length(outputs)) %>%
  t() %>%
  as.data.frame()
colnames(orig_outputs) <- outputs
rownames(orig_outputs) <- names(phases)

# checking
orig_outputs
```

```
##       final_front_position final_pioneer_position initial_front_speed
## early                 62.5                  337.1                 0.4
## mid                    4.3                  104.9                 0.4
## late                   4.5                   99.3                 1.0
##       initial_pioneer_speed final_front_speed final_pioneer_speed
## early                  38.6              10.0                13.2
## mid                    17.2              -0.4                 4.8
## late                   26.0              -0.8                 4.2
##       avg_front_speed avg_pioneer_speed acceleration_tot_front
## early       5.0000000         27.883333             0.80000000
## mid         0.1500000          8.533333            -0.06666667
## late        0.1666667          8.066667            -0.15000000
##       acceleration_tot_pioneer acceleration_est_front
## early                -2.116667             0.41338661
## mid                  -1.033333            -0.06833167
## late                 -1.816667            -0.09250749
##       acceleration_est_pioneer lambda_avg final_pop_size final_max_density
## early               -1.3074925   1.204359      465.56973         13.823161
## mid                 -0.3815185   1.000849       50.51148          1.468948
## late                -0.3244256   1.001839       51.11450          1.217789
```

```r
# looping the plot
par(mfrow = c(round(length(outputs)/3), 3))
for(i in 1:dim(orig_outputs)[2]){
  plot(orig_outputs[[i]], type='l', main = colnames(orig_outputs)[i], xlim = c(0.7, 3.3),
     xaxt = 'n', ylab = '', xlab = '', las = 1)
  axis(1, at = 1:3, labels = c('Initial', 'Middle', 'End'), las=1)
  points(orig_outputs[[i]], col = Cols, pch=19, cex = 1.5)
  text(1:3 + 0.1, orig_outputs[[i]], 
     labels = round(as.numeric(orig_outputs[[i]]), 2), cex = 1)
}
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

Now, I'll plot the results of the perturbation analyses. This will be represented as an horizontal barplot so that the relative effects could be comparable among the three relpicates.


```r
# For ease of plotting, handling the array t make it a single data.frame with outputs for columns and parameters for rows

pResults <- aperm(perturbationResults, perm = c(2,1,3))
pResults <- data.frame(apply(pResults, 2, c),
                       expand.grid(dimnames(perturbationResults)[c(2,3)]))
colnames(pResults)[(ncol(pResults)-1):ncol(pResults)] <- c("params", "phase")

# creating a vector of better names for parameteres
paramnames <- c('survival (int)', 'survival (slo)',
            'repr p (int)',  'repr p (slo)',
            'offspr n (int)',  'offspr n (slo)',
            'ad dispersal (int)','ad dispersal (slo)',
            'ad distance (mean)', 'ad distance (sd)',
            'nat dispersal (int)','nat dispersal (slo)',
            'nat dispersal (landsc)', 'nat dispersal (loc x land)',
            'nat distance (mean)', 'nat distance (sd)')

# checking the length of parameters and parameter names is the same
# checking also identity
stopifnot(length(params) == length(paramnames))
cbind(params, paramnames)
```

```
##       params               paramnames                  
##  [1,] "s.int"              "survival (int)"            
##  [2,] "s.slo"              "survival (slo)"            
##  [3,] "r.int"              "repr p (int)"              
##  [4,] "r.slo"              "repr p (slo)"              
##  [5,] "f.int"              "offspr n (int)"            
##  [6,] "f.slo"              "offspr n (slo)"            
##  [7,] "p.stat.int.ad"      "ad dispersal (int)"        
##  [8,] "p.stat.slo.ad"      "ad dispersal (slo)"        
##  [9,] "mu.mob.int.ad"      "ad distance (mean)"        
## [10,] "sigma.mob.ad"       "ad distance (sd)"          
## [11,] "p.stat.int.nat"     "nat dispersal (int)"       
## [12,] "p.stat.slo.nat"     "nat dispersal (slo)"       
## [13,] "p.stat.sloLand.nat" "nat dispersal (landsc)"    
## [14,] "p.stat.sloInt.nat"  "nat dispersal (loc x land)"
## [15,] "mu.mob.int.nat"     "nat distance (mean)"       
## [16,] "sigma.mob.nat"      "nat distance (sd)"
```

```r
## plotting all outputs
for(i in outputs){
  p <- ggplot(data = pResults, aes_string(x = "params", y = i, fill = "phase")) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge()) +
  scale_fill_manual(values = Cols) +
    coord_flip()
  theme_minimal()
print(p + theme(axis.text.x = element_text(angle = 0)) +
  scale_x_discrete(labels = paramnames) +
  ggtitle(as.character(i)))
}
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-2.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-3.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-4.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-5.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-6.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-7.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-8.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-9.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-10.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-11.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-12.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-13.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-14.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-15.png)

## Saving perturbation results


```r
saveRDS(perturbationResults, file = file.path(here::here(), "results",
                                              "perturbation_result_phase_models.rds"))
saveRDS(perturbed_phasesIPMs, file = file.path(here::here(), "results",
                                               "perturbed_IPMs_phase_models.rds"))
```

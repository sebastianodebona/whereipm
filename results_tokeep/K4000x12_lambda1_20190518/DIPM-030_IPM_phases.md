---
title: "Real data IPM - Phases"
author: "Seba De Bona"
date: "2019-05-18"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: lumen
---

Here, we run an IPM using real data. Vital rates models were fit in the previous vignette, and refer to three distinct phases ("early", "mid" and "late") during the colonisation process.

## Setup

Here we load the data.


```r
library(whereipm)
library(tidyverse)

# loading the data
data <- ipmdata
```

And the models previously run and observation of population parameters.


```r
vitalRateModels <- readRDS(file.path(here::here(), "results", "vitalRateModels.rds"))
observation_summary <- readRDS(file.path(here::here(), "results", "observation_summary.rds"))
```

## Parameters and initial Population

Here we define kernel parameters and initial population size.


```r
# kernel parameters
kp <- create.kernel.params(L=-400, U=400, n.bins = 4000)

# # initial population size for all phases
# pop0.early <- makeInitPop(kp, 
#                           observation_summary%>%
#                             filter(phases == "early", pool == "top") %>%
#                             pull(avg_n) %>%
#                             round(),
#                           observation_summary%>%
#                             filter(phases == "early", pool == "top") %>%
#                             pull(avg_l) %>%
#                             round())
# 
# pop0.mid <- makeInitPop(kp, 
#                           observation_summary%>%
#                             filter(phases == "mid", pool == "top") %>%
#                             pull(avg_n) %>%
#                             round(),
#                           observation_summary%>%
#                             filter(phases == "mid", pool == "top") %>%
#                             pull(avg_l) %>%
#                             round())
# 
# pop0.late <- makeInitPop(kp, 
#                           observation_summary%>%
#                             filter(phases == "late", pool == "top") %>%
#                             pull(avg_n) %>%
#                             round(),
#                           observation_summary%>%
#                             filter(phases == "late", pool == "top") %>%
#                             pull(avg_l) %>%
#                             round())

pop0.early = pop0.mid = pop0.late <- makeInitPop(kp, 50, 5)
```

And we extract the vital rates parameters from the models fitted in the [Vital Rate Models: phases](DIPM-020_Vital_rate_models.html) vignette.


```r
# First I need to extract the info in the model to make it easier to handle
Models <- vitalRateModels %>%
  transmute(
    survival = map(survival, broom::tidy),
    recruit = map(recruit, broom.mixed::tidy),
    adDispersal = map(adDispersal, broom::tidy),
    adDistance = map(adDistance, broom::tidy),
    natDispersal = map(natDispersal, broom::tidy),
    natDistance = map(natDistance, broom::tidy))

# Storing values into vital params list, divided into early, mid and late
# given local density is scaled (but not centered) the estimates for the slope have to be transformed back to the non-scaled value.

# The percentage to which the recruitment needs to be buffed is calculated 
# in the 'Buffing Recruitment' vignette. Here I'll add the phase, specific proportions
buff_recruit <- c(2.86, 1.90, 2.55)

# Using a general loop to do so for the three phases:

for(i in 1:dim(vitalRateModels)[1]){
 

  temp.vp <- list(# Survival
                 s.int = Models$survival[[i]]%>%
                   filter(term == "(Intercept)") %>%
                   pull(estimate),
                 s.slo = (Models$survival[[i]]%>%
                   filter(term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 s.slo2 = returnIfExists(Models$survival[[i]]%>%
                   filter(term == "I(s_loc_density^2)") %>%
                   pull(estimate)),
                 s.multiplier = 1,
                 
                # Recruitment
                 r.int = Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "zi",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 r.slo = (Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "zi",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 r.multiplier = buff_recruit[i],
                 f.int = Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "cond",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 f.slo = (Models$recruit[[i]]%>%
                   filter(effect == "fixed",
                          component == "cond",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 f.distribution = 'ztpoisson',
                 
                # Adult dispersal
                 p.stat.int.ad = Models$adDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 p.stat.slo.ad = (Models$adDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 sigma.stat.ad = 0.25, 
                 mu.mob.int.ad = Models$adDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 mu.mob.slo.ad = (Models$adDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 sigma.mob.ad = Models$adDistance[[i]]%>%
                   filter(effect == "ran_pars") %>%
                   pull(estimate) %>%
                   sum %>%
                   sqrt,
                 
                # Natal dispersal
                 p.stat.int.nat = Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 p.stat.slo.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 p.stat.sloLand.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_land_density.mid") %>%
                   pull(estimate))/scalingAttr$land_density$scale,
                 p.stat.sloInt.nat = (Models$natDispersal[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid:s_land_density.mid") %>%
                   pull(estimate))/(scalingAttr$loc_density$scale*scalingAttr$land_density$scale),
                 sigma.stat.nat = 0.25, 
                 mu.mob.int.nat = Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "(Intercept)") %>%
                   pull(estimate),
                 mu.mob.slo.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid") %>%
                   pull(estimate))/scalingAttr$loc_density$scale,
                 mu.mob.sloLand.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_land_density.mid") %>%
                   pull(estimate))/scalingAttr$land_density$scale,
                 mu.mob.sloInt.nat = (Models$natDistance[[i]]%>%
                   filter(effect == "fixed",
                          term == "s_loc_density.mid:s_land_density.mid") %>%
                   pull(estimate))/(scalingAttr$loc_density$scale*scalingAttr$land_density$scale),
                 sigma.mob.nat = Models$natDistance[[i]]%>%
                   filter(effect == "ran_pars") %>%
                   pull(estimate) %>%
                   sum %>%
                   sqrt)
 
  # replacing empty elements with zeros
  temp.vp <- lapply(temp.vp, function(x) ifelse(length(x)==0,0,x))
  
  # renaming list
  assign(paste("vp", vitalRateModels$phase[i], sep = "."), temp.vp)
}

# visualizing the three lists side by side
tab <- tibble("parameters" = names(vp.early),
                  "early phase" = as.numeric(unlist(vp.early)),
                  "mid phase" = as.numeric(unlist(vp.mid)),
                  "late phase" = as.numeric(unlist(vp.late)))
```

```
## Warning in eval_tidy(xs[[i]], unique_output): NAs introduced by coercion

## Warning in eval_tidy(xs[[i]], unique_output): NAs introduced by coercion

## Warning in eval_tidy(xs[[i]], unique_output): NAs introduced by coercion
```

```r
knitr::kable(tab, digits = 3, col.names = c("Parameters", "Early", "Mid", "Late"), caption = "Comparison of vital rates among phases")
```



|Parameters         |  Early|    Mid|   Late|
|:------------------|------:|------:|------:|
|s.int              |  2.088|  1.999|  1.406|
|s.slo              |  0.029|  0.006|  0.022|
|s.slo2             |  0.000|  0.000|  0.000|
|s.multiplier       |  1.000|  1.000|  1.000|
|r.int              |  2.493|  2.798|  2.591|
|r.slo              | -0.012|  0.002|  0.020|
|r.multiplier       |  2.860|  1.900|  2.550|
|f.int              | -0.357| -1.653| -1.447|
|f.slo              | -0.040| -0.027| -0.044|
|f.distribution     |     NA|     NA|     NA|
|p.stat.int.ad      | -0.126| -0.873| -0.944|
|p.stat.slo.ad      | -0.198| -0.012|  0.135|
|sigma.stat.ad      |  0.250|  0.250|  0.250|
|mu.mob.int.ad      |  1.857|  1.699|  1.598|
|mu.mob.slo.ad      |  0.000|  0.000|  0.000|
|sigma.mob.ad       |  1.381|  1.309|  1.304|
|p.stat.int.nat     | -1.227| -8.703| -3.992|
|p.stat.slo.nat     |  0.204|  0.421| -0.191|
|p.stat.sloLand.nat |  0.530|  1.433|  1.097|
|p.stat.sloInt.nat  | -0.080| -0.068|  0.024|
|sigma.stat.nat     |  0.250|  0.250|  0.250|
|mu.mob.int.nat     |  2.299|  2.055|  2.284|
|mu.mob.slo.nat     |  0.000|  0.000|  0.000|
|mu.mob.sloLand.nat |  0.000|  0.000|  0.000|
|mu.mob.sloInt.nat  |  0.000|  0.000|  0.000|
|sigma.mob.nat      |  1.418|  1.280|  1.078|

## Running the IPM model

Now that all parameters are ready, I can run the Integral Projection models. First, some diagnostic checks using the tools in the package.

### Visual checks


```r
test.init.kernel(kp, vp.early, pop0.early)
```
![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

### Running the IPM

Below, I run the 3 IPM simulations for the early, mid and late phase. The simulations are rather time consuming, so I will store them once run, in the results folder, as part of the tibble containing vital rate models.


```r
run.time = 12

# Early phase ipm
start <- Sys.time()
IPM.early <- runIPM(kp, vp.early, init.pop = pop0.early, niter=run.time)
early.t <- Sys.time() - start
print(early.t)
```

```
## Time difference of 19.98482 mins
```

```r
# Mid phase ipm
start <- Sys.time()
IPM.mid <- runIPM(kp, vp.mid, init.pop = pop0.mid, niter=run.time)
mid.t <- Sys.time() - start
print(mid.t)
```

```
## Time difference of 19.85135 mins
```

```r
# Late phase ipm
start <- Sys.time()
IPM.late <- runIPM(kp, vp.late, init.pop = pop0.late, niter=run.time)
late.t <- Sys.time() - start
print(late.t)
```

```
## Time difference of 19.77368 mins
```

```r
par(mfrow=c(1,3))
plotly::plot_ly(x=1:ncol(IPM.early$pop.size), y=rev(kp$y), z=IPM.early$pop.size, type="contour")
```

```
## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

```
## Warning in normalizePath(f2): path[1]="./webshot34562c3f8681.png": No such
## file or directory
```

```
## Warning in file(con, "rb"): cannot open file './webshot34562c3f8681.png':
## No such file or directory
```

```
## Error in file(con, "rb"): cannot open the connection
```

```r
plotly::plot_ly(x=1:ncol(IPM.mid$pop.size), y=rev(kp$y), z=IPM.mid$pop.size, type="contour")
```

```
## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

```
## Warning in normalizePath(f2): path[1]="./webshot345676e132cf.png": No such
## file or directory
```

```
## Warning in file(con, "rb"): cannot open file './webshot345676e132cf.png':
## No such file or directory
```

```
## Error in file(con, "rb"): cannot open the connection
```

```r
plotly::plot_ly(x=1:ncol(IPM.late$pop.size), y=rev(kp$y), z=IPM.late$pop.size, type="contour")
```

```
## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

```
## Warning in normalizePath(f2): path[1]="./webshot345614159895.png": No such
## file or directory
```

```
## Warning in file(con, "rb"): cannot open file './webshot345614159895.png':
## No such file or directory
```

```
## Error in file(con, "rb"): cannot open the connection
```

```r
# another way to do it
quick.plot(IPM.early, 'Early')
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
quick.plot(IPM.mid, "Mid")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)

```r
quick.plot(IPM.late, 'Late')
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-3.png)


### Storing results


```r
# creating a tibble to store IPM results
phaseModelsIPMs <- tibble(phase = list("early", "mid", "late"),
                          vital.params = list(vp.early, vp.mid, vp.late),
                          IPMs = list(IPM.early, IPM.mid, IPM.late),
                          timetorun = list(early.t, mid.t, late.t))

# storing the updated tibble
saveRDS(phaseModelsIPMs, file = file.path(here::here(), "results", "phaseModelsIPMs.rds"))
```

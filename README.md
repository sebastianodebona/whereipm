[![pipeline status](https://gitlab.com/sebadebona/dispersal-ipm/badges/master/pipeline.svg)](https://gitlab.com/sebadebona/dispersal-ipm/commits/pkgdev)

# Description

This project explores the relationship between density dependence and population spread during an invasion using an Integral Projection Model approach.
The individual based dataset collected in Trinidad during the Guppy CMR project is used to parameterise the model, which predicts how the population will diffuse in space given natal and adult dispersal, survival and reproduction. The code included, altogether, constructs the IPM and evaluates the effect of different values across the parameter space on population spread and the speed of invasion.
The project constitutes part of Sebastiano De Bona's PhD thesis.

You can access a more detailed description of the project, a tutorial, and some result figres here:

- [Where I (p)m?](https://sebadebona.gitlab.io/dispersal-ipm/)

## How to clone this project locally

If you have access to the GitLab repository, you can clone this project
remotely by opening a terminal and typing:

```
git clone git@gitlab.com:sebadebona/dispersal-ipm.git
```

## How to prepare your R for the analyses

A preparation step is required before the analyses, and it entails installing a list of required packages. Copy-paste the following code in R or RStudio and execute:

```
install.packages(c('tidyverse', 'plyr', 'RColorBrewer', 'lme4', 'glmmTMB', 'gridExtra', 'car', 'sticky', 'rtf', 'furrr', 'distr', 'plotly', 'sjPlot'))
```


# How to reproduce this analysis

Please, see the [tutorial](https://sebadebona.gitlab.io/dispersal-ipm/) to know how to reproduce the analyses and run a spatial IPM.

If you want to know more, please contact @sebadebona.


# Data

The package comes with a datafile attached (`ipmdata.Rda`). The script `prep-ipmdata.R` uses raw data files to create the Rda file.

# References

For more information, see published paper
on [LoST group's page](http://ecoevo.pagecloud.com/).

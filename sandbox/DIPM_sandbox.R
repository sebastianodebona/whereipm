x <- seq(0, 80, length.out=100)
y <- dnorm(x)
plot(y~x, type='l')

# figuring out lognormal density function

sigma=1.334262
mu=1.8463
y <- (1/x)*(1/(sigma*sqrt(2*pi)))*exp(-((log(x)-mu)^2)/(2*sigma^2))

y1 <- dlnorm(x, meanlog = mu, sdlog = sigma)
lines(y1~x, col='red')

hist(data.early$dist_moved, breaks=20, ylim=c(0,300))
n <- nrow(data.early)
lines(y1*n~x, col='red')


# testing gx_ln function
x <- rep(0, 200)
y <- seq(-99, 100)

p <- gx_ln(y = combs[['y']], x = combs[['x']],
      density.x = combs[['density.x']],
      p.stat.int = summary(e.dispersal.p)$coefficients['(Intercept)','Estimate'],
      p.stat.slo = summary(e.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'),
      sigma.stat = 0.25, 
      mu.mob.int = summary(e.dispersal.dist)$coefficients['(Intercept)','Estimate'],
      mu.mob.slo = 0,
      sigma.mob = sqrt(sum(as.data.frame(lme4::VarCorr(e.dispersal.dist))[,'sdcor'])))

plot(p~y, type='l')
lines(p~y, col='green')

x <- seq(1, 100)
y <- seq(1, 100)
densities <- rep(0, 100)

combs <- merge(data.frame(y=x), cbind(x=x, density.x=densities))



mat <- matrix(NA, length(y), length(y))
mat[x,y] <- p
image(t(log(mat)))

# Creating animation for population spread
# first, creating series of frames
pdf(file = "Images/toanimate.pdf", width = 4.5, height = 5)
for (i in 1:37){
  plot(0, xlim=range(IPM.early$y), ylim=c(0, 60), col='white', xlab = '', ylab = '', yaxt = 'n')
  axis(2, at = c(0, 20, 40, 60), labels = c(0, 20, 40, 60), las = 2)
  mtext('n', side = 2, las = 2, line = 3)
  mtext('x', side = 1, line = 2)
  I <- ifelse(i>13, 13, i)
  for (j in 1:I){
    draw.pop.distr(x=IPM.early$y, pop=IPM.early$pop.size[,j], col='#AE1A12', new.plot = F)
  }
}
dev.off()



frames <- magick::image_read_pdf('Images/toanimate.pdf', density = 400)
animation <- magick::image_animate(frames, fps = 5)
image_write(animation, "Images/Early.gif")




# some random graphs for plotting
x <- 0:100

png(filename = "Images/Surv", width = 6, height = 7.5, units = "cm", res = 300)
survival <- function(x){plogis(1.5 - 0.02*x)}
curve(survival(x), ylim = c(0,1), xlim = c(0,100), ylab = '', xlab = '', yaxt = 'n', xaxt='n')
dat <- data.frame(y = rbinom(length(x), 1, prob = survival(x) + rnorm(length(x), 0, 0.03)), x = x)
points(y~x, pch = 16, cex = 0.8, data = dat[sample(x, 50),])
axis(2, at = c(0,.5,1), labels = c('0.0','0.5','1.0'), las = 2)
axis(1, at = c(0, 50, 100), labels = c(0, 50, 100))
dev.off()


png(filename = "Images/Fert", width = 6, height = 7.5, units = "cm", res = 300)
offspring <- function(x){exp(0.6 - 0.02*x)}
curve(offspring(x), ylim = c(0,4), xlim = c(0,100), ylab = '', xlab = '', yaxt = 'n', xaxt='n')
dat <- data.frame(y = rpois(length(x), lambda = offspring(x) + rnorm(length(x), 0, 0.04)), x = x)
points(y~x, pch = 16, cex = 0.8, data = dat[sample(x, 50),])
axis(2, at = c(0,2,4), labels = c('0','2','4'), las = 2)
axis(1, at = c(0, 50, 100), labels = c(0, 50, 100))
dev.off()


x <- seq(0, 80, length.out=100)
y <- dnorm(x)
plot(y~x, type='l')

# figuring out lognormal density function

sigma=1
mu=1.5

png(filename = "Images/NDisp", width = 6, height = 7.5, units = "cm", res = 300)
p <- 0.2
x <- seq(-50, 50, by=1)
y1 <- ifelse(x>=0, p*dnorm(x, mean = 0, sd = 0.15) + ((1-p)/2)*dlnorm(x, meanlog = mu, sdlog = sigma), p*dnorm(x, mean = 0, sd = 0.25) + ((1-p)/2)*dlnorm(-x, meanlog = mu, sdlog = sigma))
plot(y1~x, col='black', xlim = c(-55, 55), ylim = c(0,1), type='l', xaxt = 'n', yaxt = 'n', xlab = '', ylab='')
axis(2, at = c(0,.5,1), labels = c('0','.5','1'), las = 2)
axis(1, at = c(-50, 0, 50), labels = c(-50, 0, 50))
dev.off()

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


# drawing vital rates and dispersal through time
# DISPESAL PROBABILITY
dens <- seq(0, 50, length.out = 100)

p.e <- function(x){plogis(summary(e.dispersal.p)$coefficients['(Intercept)','Estimate'] + 0*x)}
p.m <- function(x){plogis(summary(m.dispersal.p)$coefficients['(Intercept)','Estimate'] + (summary(m.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'))*x)}
p.l <- function(x){plogis(summary(l.dispersal.p)$coefficients['(Intercept)','Estimate'] + (summary(l.dispersal.p)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'))*x)}

png(filename = "Images/ADispModel", width = 9, height = 10, units = "cm", res = 300)
plot(p.e(dens) ~ dens, ylim = c(0,1), type = 'l', col = '#AE1A12', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = 'density (n/m)', ylab = '')
lines(p.m(dens) ~ dens, col = '#114085', lwd = 2)
lines(p.l(dens) ~ dens, col = '#258646', lwd = 2)
axis(1, at=c(0,10,20,30,40, 50), labels = c(0,10,20,30,40, 50))
axis(2, at=c(0,.25,.5,.75,1), las = 2)
mtext('p', side = 2, las=2, line=3)
dev.off()


# ADULT DISPERSAL DISTANCE
dens <- seq(0, 50, length.out = 100)
mu.e <- function(x){exp(summary(e.dispersal.dist)$coefficients['(Intercept)','Estimate'] + 0*x)}
mu.m <- function(x){exp(summary(m.dispersal.dist)$coefficients['(Intercept)','Estimate'] + (summary(m.dispersal.dist)$coefficients['s_loc_density','Estimate']/attr(data$s_loc_density, 'scaled:scale'))*x)}
mu.l <- function(x){exp(summary(l.dispersal.dist)$coefficients['(Intercept)','Estimate'] + 0*x)}

png(filename = "Images/ADistModel", width = 9, height = 10, units = "cm", res = 300)
plot(mu.e(dens) ~ dens, type = 'l', ylim = c(3,9), col = '#AE1A12', lwd = 2, xaxt = 'n', yaxt = 'n', xlab = 'density (n/m)', ylab = '')
lines(mu.m(dens) ~ dens, col = '#114085', lwd = 2)
lines(mu.l(dens) ~ dens, col = '#258646', lwd = 2)
axis(1, at=c(0,10,20,30,40, 50), labels = c(0,10,20,30,40, 50))
axis(2, at=c(2, 4, 6, 8, 10), las = 2)
mtext('m', side = 2, las=2, line=2)
dev.off()

# Drawing estimated distribution kernels given parameter estimates
k.e <- function(i){
  y <- ifelse(x>=0,
              (1-p.e(i))*dnorm(x, mean = 0, sd = 0.15) + (p.e(i)/2)*dlnorm(x, meanlog = log(mu.e(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(e.dispersal.dist))[,'sdcor']))), 
              (1-p.e(i))*dnorm(x, mean = 0, sd = 0.25) + (p.e(i)/2)*dlnorm(-x, meanlog = log(mu.e(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(e.dispersal.dist))[,'sdcor']))))
  return(y)
}

k.m <- function(i){
  y <- ifelse(x>=0,
              (1-p.m(i))*dnorm(x, mean = 0, sd = 0.15) + (p.m(i)/2)*dlnorm(x, meanlog = log(mu.m(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(m.dispersal.dist))[,'sdcor']))), 
              (1-p.m(i))*dnorm(x, mean = 0, sd = 0.25) + (p.m(i)/2)*dlnorm(-x, meanlog = log(mu.m(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(m.dispersal.dist))[,'sdcor']))))
  return(y)
}

k.l <- function(i){
  y <- ifelse(x>=0,
              (1-p.l(i))*dnorm(x, mean = 0, sd = 0.15) + (p.l(i)/2)*dlnorm(x, meanlog = log(mu.l(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(l.dispersal.dist))[,'sdcor']))), 
              (1-p.l(i))*dnorm(x, mean = 0, sd = 0.25) + (p.l(i)/2)*dlnorm(-x, meanlog = log(mu.l(i)), sdlog = sqrt(sum(as.data.frame(lme4::VarCorr(l.dispersal.dist))[,'sdcor']))))
  return (y)
  }

space <- -20:20

# EARLY kernels through different densities
png(filename = "Images/KernelEarly", width = 9, height = 9, units = "cm", res = 500)
par(bg=NA)
plot(k.e(5) ~ space, type='l', col = '#AE1A12', ylim = c(-1, 2), xlim = c(-40,60), bty = 'n', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
segments(x0=0, y0=0, y1=2.1, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.e(25) - 0.33) ~ I(space + 10), col = '#AE1A12')
segments(x0=10, y0=-.33, y1=1.7, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.e(50) - 0.66) ~ I(space + 20), col = '#AE1A12')
segments(x0=20, y0=-.66, y1=1.4, lty = 2, col = 'grey88')
dev.off()
# MID KERNELS
png(filename = "Images/KernelMid", width = 9, height = 9, units = "cm", res = 500)
par(bg=NA)
plot(k.m(5) ~ space, type='l', col = '#114085', ylim = c(-1, 2), xlim = c(-40,60), bty = 'n', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
segments(x0=0, y0=0, y1=2.1, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.m(25) - 0.33) ~ I(space + 10), col = '#114085')
segments(x0=10, y0=-.33, y1=1.7, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.m(50) - 0.66) ~ I(space + 20), col = '#114085')
segments(x0=20, y0=-.66, y1=1.4, lty = 2, col = 'grey88')
dev.off()
# LATE KERNELS
png(filename = "Images/KernelLate", width = 9, height = 9, units = "cm", res = 500)
par(bg=NA)
plot(k.l(5) ~ space, type='l', col = '#258646', ylim = c(-1, 2), xlim = c(-40,60), bty = 'n', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
segments(x0=0, y0=0, y1=2.1, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.l(25) - 0.33) ~ I(space + 10), col = '#258646')
segments(x0=10, y0=-.33, y1=1.7, lty = 2, col = 'grey88')
par(new = T)
lines(I(k.l(50) - 0.66) ~ I(space + 20), col = '#258646')
segments(x0=20, y0=-.66, y1=1.4, lty = 2, col = 'grey88')
dev.off()




# creating array to store results of perturbation analyses
perturbation <- array(NA, dim = c(4,14,3), dimnames = list(c('speed', 'acceleration', 'lambda', 'pop'), c('s.int', 's.slo',  'r.int',  'r.slo', 'f.int',  'f.slo',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'mu.mob.slo.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat', 'mu.mob.slo.nat'), c('IPM.early', 'IPM.mid', 'IPM.late')))

phases <- list(IPM.early, IPM.mid, IPM.late)
params <- list(early = c('s.int',  'r.int', 'f.int', 'f.slo', 'p.stat.int.ad',  'mu.mob.int.ad', 'p.stat.int.nat', 'mu.mob.int.nat'),
               mid = c('s.int', 's.slo',  'r.int',  'r.slo', 'f.int',  'f.slo',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'mu.mob.slo.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat', 'mu.mob.slo.nat'), 
               late = c('s.int', 's.slo',  'r.int', 'f.int',  'p.stat.int.ad','p.stat.slo.ad',  'mu.mob.int.ad', 'p.stat.int.nat','p.stat.slo.nat', 'mu.mob.int.nat'))

# running loop to run perturbation analyses and store the results
start <- Sys.time()
for(I in 1:length(phases)){
  for(j in 1:length(params[[I]])){
    temp <- perturb(orig.IPM = phases[[I]], epsilon = 0.01, par = params[[I]][j])
    perturbation[, params[[I]][[j]], I] <- unlist(temp)
  }
}
Sys.time() - start

# turning absolute values for perturbation analyses into proportional ones
orig.early <- c(IPM.early$speed,
                IPM.early$acceleration,
                sum(IPM.early$pop.size[,IPM.early$n.iterations+1])/sum(IPM.early$pop.size[,IPM.early$n.iterations]),
                IPM.early$pop.size[IPM.early$kernel.params$n.bins/2,IPM.early$n.iterations+1])
orig.mid <- c(IPM.mid$speed,
              IPM.mid$acceleration,
              sum(IPM.mid$pop.size[,IPM.mid$n.iterations+1])/sum(IPM.mid$pop.size[,IPM.mid$n.iterations]),
              IPM.mid$pop.size[IPM.mid$kernel.params$n.bins/2,IPM.mid$n.iterations+1])
orig.late <- c(IPM.late$speed,
               IPM.late$acceleration,
               sum(IPM.late$pop.size[,IPM.late$n.iterations+1])/sum(IPM.late$pop.size[,IPM.late$n.iterations]),
               IPM.late$pop.size[IPM.late$kernel.params$n.bins/2,IPM.late$n.iterations+1])

perturbation[,,1] <- perturbation[,,1]/orig.early
perturbation[,,2] <- perturbation[,,2]/orig.mid
perturbation[,,3] <- perturbation[,,3]/orig.late
# save(perturbation, file = 'sensitivity.rds')

load('sensitivity.rds')

outputs <- as.data.frame(rbind(orig.mid, orig.late))
colnames(outputs) <- c('speed', 'acc', 'lambda', 'pop')

#png(filename = "Images/Lambda.png", width = 8, height = 9, units = "cm", res = 300)
#par(mfrow = c(2,2), mar = c(3,4,3,2), oma = c(0,2,1,1))
plot(outputs$lambda, type='l', main='Lambda', ylim = c(0.5, 1.5), xaxt = 'n', ylab = '', xlim = c(0.7, 2.3), xlab = '', yaxt = 'n')
axis(2, at = c(0.6, 0.8, 1, 1.2, 1.4), labels = c('0.6', '0.8', '1.0', '1.2', '1.4'), las=2)
axis(1, at = c(1,2), labels = c('Mid', 'Late'), las=1)
points(outputs$lambda, col = c('#114085', '#258646'), pch=19, cex = 1.5)
text(1:2+0.12, outputs$lambda + 0.1, labels = round(as.numeric(outputs$lambda), 2), cex = 1)
mtext(expression(lambda), 2, line = 2.5)
#dev.off()

#png(filename = "Images/PopCore.png", width = 8, height = 9, units = "cm", res = 300)
plot(outputs$pop, type='l', main='Core population size', xaxt = 'n', ylab = '', xlim = c(0.7, 2.3), xlab = '', ylim = c(0, 80), yaxt = 'n')
points(outputs$pop, col = c('#114085', '#258646'), pch=19, cex = 1.5)
axis(2, at = c(0, 20, 40, 60, 80), labels = c(0, 20, 40, 60, 80), las=2)
text(1:2+0.12, outputs$pop + 5, labels = round(as.numeric(outputs$pop), 0), cex = 1)
axis(1, at = c(1,2), labels = c('Mid', 'Late'), las=1)
mtext('N', 2, line = 2.5)
#dev.off()

#png(filename = "Images/Speed.png", width = 8, height = 9, units = "cm", res = 300)
plot(outputs$speed, type='l', main='Speed at the front', ylim = c(-1,20), xaxt = 'n', xlim = c(.7, 2.3), xlab = '', ylab = '', yaxt = 'n')
points(outputs$speed, col = c('#114085', '#258646'), pch=19, cex = 1.5)
axis(2, at = c(0, 5, 10, 15, 20), labels = c(0, 5, 10, 15, 20), las=2)
axis(1, at = c(1,2), labels = c('Mid', 'Late'), las=1)
text(1:2-0.1, outputs$speed - 2, labels = round(as.numeric(outputs$speed), 2), cex = 1)
mtext('m/month', 2, line = 2.5)
#dev.off()

#png(filename = "Images/Acc.png", width = 8, height = 9, units = "cm", res = 300)
plot(outputs$acc, type='l', ylim = c(0,0.4), main='Acceleration', xaxt = 'n', xlim = c(.7, 2.3), xlab = '', ylab = '', yaxt = 'n')
axis(2, at = c(0, .1, .2, .3, .4), labels = c(0, .1, .2, .3, .4), las=2)
axis(1, at = c(1,2), labels = c('Mid', 'Late'), las=1)
points(outputs$acc, col = c('#114085', '#258646'), pch=19, cex = 1.5)
text(1:3+0.1, outputs$acc + 0.05, labels = round(as.numeric(outputs$acc), 2), cex = 1)
mtext(expression(m/month^2), 2, line = 2.5)
#dev.off()

# plot perturbations
#png(filename = "Images/Speed_el.png", width = 22, height = 12, units = "cm", res = 1000)
par(mfrow = c(1,2), oma = c(2,3,2,2), mar = c(3,3,3,3))
col <- c('red', '#114085', '#258646')
ph <- c('Early', 'Mid', 'Late')
for(i in 2:3){
  pattern <- vector('numeric', length = length(params[[i]]))
  pattern[grep('slo', params[[i]])] <- 2000
  pattern[grep('int', params[[i]])] <- 50
  barplot(rev(perturbation[1,,i]), horiz = T, density = rev(pattern), col = col[i], yaxt = 'n',  main = ph[i], space = rev(c(0,.8,0,.8,0,.8,0,.8,0,.8,0,.8,0,0)), width = 0.1, xlim = c(0.0, 0.1))
}
#dev.off()


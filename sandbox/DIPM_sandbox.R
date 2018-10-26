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





# testing gx_ln function
x <- rep(0, 200)
y <- seq(-99, 100)

p <- gx_ln(y = combs[['y']], x = combs[['x']],
      density.x = combs[['density.x']],
      p.stat.int = summary(e.dispersal.p)$coefficients['(Intercept)','Estimate'],
      p.stat.slo = summary(e.dispersal.p)$coefficients['s_loc_density','Estimate']/scalingAttr$loc_density$scale,
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





##########################
# Using package distr to create custom functions

y <- seq(-500, 500, by = 0.001)
x = 4

N <- distr::Norm(mean = 0, sd = 0.25)
L <- distr::Lnorm(meanlog = 1.83, sd = 1.33)
D <- (- L*0.3) + L*0.3 + N*0.4

sum(distr::d(D)(y-4)*0.001) # important to multiply by h!!! (here 0.001)


distr::plot(D, panel.first = grid())

z <- seq(-500, 500, by = 0.0005)



### Drawing a function by hand (from scratch)
mean = 0
sd = 1
Normal <- function(x){
  (1/(sqrt(2*pi*(sd^2))))*exp(-((x-mean)^2/(2*(sd)^2)))
}

meanlog = 1.8
sdlog = 1
logNormal <- function(x){
  (1/(x*sdlog*(sqrt(2*pi))))*exp(-((log(x)-meanlog)^2/(2*(sdlog^2))))
}

integrate(Normal, -Inf, Inf)

x <- seq(-100, 100, by =0.001)

plot(D(x)~x, type = "l", col = 'red')

D <- function(x){
  y <- ifelse(x==0,
         0.6*Normal(x),
         0.6*Normal(x) + 0.2*logNormal(abs(x)))
  return(y)
}

abs(integrate(D, -Inf, Inf)$value-1) < 1e-4

D(0)

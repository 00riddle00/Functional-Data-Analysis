# vim: set ft=r tw=88 nu ai et ts=2 sw=2:

rm(list=ls())
library(fda)

### Precipitation data

logprecav = CanadianWeather$dailyAv[, , 'log10precip']
yearRng  = c(0,365)
daybasis = create.fourier.basis(yearRng, 365)
Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

lambda     = 1e6
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fd = smooth.basis(day.5, logprecav, fdParobj)$fd

nharm = 4
pcalist = pca.fd(logprec.fd, nharm, centerfns = TRUE)
plot(pcalist)
plot(pcalist$harmonics)

plotscores(pcalist, loc = 5)


# pca graphs explained
c <- 2
mn <- pcalist$meanfd
phi <- pcalist$harmonics[1]
lambda <- pcalist$values[1]

f1 <- mn - c*sqrt(lambda)*phi
f2 <- mn + c*sqrt(lambda)*phi

opar <- par(mfrow=c(1,2))
plot(pcalist, harm = 1)

plot(mn, ylim = c(-1,1))
lines(f1, col = 2)
lines(f2, col = 3)
par(opar)

# PCA restore the original curves
fd.pca1.list <- list()
fd.pca2.list <- list()
fd.pca3.list <- list()
fd.pca4.list <- list()

for(i in 1:5) {
  fd.pca1.list[[i]] <- mean.fd(logprec.fd) +
    pcalist$scores[i,1]*pcalist$harmonics[1]

  fd.pca2.list[[i]] <- mean.fd(logprec.fd) +
    pcalist$scores[i,1]*pcalist$harmonics[1] +
    pcalist$scores[i,2]*pcalist$harmonics[2]

  fd.pca3.list[[i]]<- mean.fd(logprec.fd) +
    pcalist$scores[i,1]*pcalist$harmonics[1] +
    pcalist$scores[i,2]*pcalist$harmonics[2] +
    pcalist$scores[i,3]*pcalist$harmonics[3]

  fd.pca4.list[[i]]<- mean.fd(logprec.fd) +
    pcalist$scores[i,1]*pcalist$harmonics[1] +
    pcalist$scores[i,2]*pcalist$harmonics[2] +
    pcalist$scores[i,3]*pcalist$harmonics[3] +
    pcalist$scores[i,4]*pcalist$harmonics[4]
}

opar <- par(mfrow=c(2,2), ask = TRUE)
for(i in 1:5) {
  plot(fd.pca1.list[[i]], ylim=c(-1, 1), ylab = "1 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.pca2.list[[i]], ylim=c(-1, 1), ylab = "2 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.pca3.list[[i]], ylim=c(-1, 1), ylab = "3 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.pca4.list[[i]], ylim=c(-1, 1), ylab = "4 PC")
  lines(logprec.fd[i], col = 2)
}
par(opar)



#### Rotation
varmx <- varmx.pca.fd(pcalist)
plot(varmx)

plot(varmx$harmonics)

plotscores(varmx, loc = 5)


# PCA restore the original curves
fd.vrm1.list <- list()
fd.vrm2.list <- list()
fd.vrm3.list <- list()
fd.vrm4.list <- list()

for(i in 1:5) {
  fd.vrm1.list[[i]] <- mean.fd(logprec.fd) +
    varmx$scores[i,1]*varmx$harmonics[1]

  fd.vrm2.list[[i]] <- mean.fd(logprec.fd) +
    varmx$scores[i,1]*varmx$harmonics[1] +
    varmx$scores[i,2]*varmx$harmonics[2]

  fd.vrm3.list[[i]]<- mean.fd(logprec.fd) +
    varmx$scores[i,1]*varmx$harmonics[1] +
    varmx$scores[i,2]*varmx$harmonics[2] +
    varmx$scores[i,3]*varmx$harmonics[3]

  fd.vrm4.list[[i]]<- mean.fd(logprec.fd) +
    varmx$scores[i,1]*varmx$harmonics[1] +
    varmx$scores[i,2]*varmx$harmonics[2] +
    varmx$scores[i,3]*varmx$harmonics[3] +
    varmx$scores[i,4]*varmx$harmonics[4]
}

opar <- par(mfrow=c(2,2), ask = TRUE)
for(i in 1:5) {
  plot(fd.vrm1.list[[i]], ylim=c(-1, 1), ylab = "1 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.vrm2.list[[i]], ylim=c(-1, 1), ylab = "2 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.vrm3.list[[i]], ylim=c(-1, 1), ylab = "3 PC")
  lines(logprec.fd[i], col = 2)

  plot(fd.vrm4.list[[i]], ylim=c(-1, 1), ylab = "4 PC")
  lines(logprec.fd[i], col = 2)
}
par(opar)


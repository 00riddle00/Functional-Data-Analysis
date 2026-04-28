# vim: set ft=r tw=88 nu ai et ts=2 sw=2:

rm(list=ls())
library(fda)
library(fda.usc)

###################################################
### depth
###################################################
tt=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
plot(fdataobj)

l <- c(0 ,2 ^ seq(-2, 9, len = 30))
nb <- seq(7, 31, by = 2)
fdataobj2 <- optim.basis(fdataobj, lambda = l, numbasis = nb)
plot(fdataobj2$fdata.est)

# Fraiman-Muniz Depth
out.FM=depth.FM(fdataobj2$fdata.est,trim=0.1,draw=TRUE)
#Modal Depth
out.mode=depth.mode(fdataobj2$fdata.est,trim=0.1,draw=TRUE)
#Random projection
out.RP=depth.RP(fdataobj2$fdata.est,trim=0.1,draw=TRUE)

# boxplot
logprecav = CanadianWeather$dailyAv[, , 'log10precip']

#  set up a saturated basis: as many basis functions as observations

yearRng  = c(0,365)
daybasis = create.fourier.basis(yearRng, 365)

#  define the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

#  smooth data with lambda that minimizes GCV

lambda     = 1e6
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fd = smooth.basis(day.5, logprecav, fdParobj)$fd
plot(logprec.fd)

boxplot(logprec.fd)
boxplot(fdata2fd(fdataobj2$fdata.est))



#########################################
#https://cran.r-project.org/web/packages/fdaoutlier/fdaoutlier.pdf
library(fdaoutlier)

lgp <- eval.fd(tt, logprec.fd)
head(lgp)

tmp <- eval.fd(tt, fdata2fd(fdataobj2$fdata.est))
head(tmp)

#needs transpose
bd <- band_depth(dt = t(lgp))
names(bd) <- colnames(lgp)
bd
plot(bd, type="l")

bd2 <- band_depth(dt = t(tmp))
names(bd2) <- colnames(tmp)
bd2
plot(bd2, type="l")

mbd <- modified_band_depth(t(lgp))
names(mbd) <- colnames(lgp)
mbd
plot(mbd, type="l")

mbd2 <- modified_band_depth(t(tmp))
names(mbd2) <- colnames(tmp)
mbd2
plot(mbd2, type="l")


fbplot_obj <- functional_boxplot(t(lgp), depth_method = "bd")
fbplot_obj$outliers

fbplot_obj2 <- functional_boxplot(t(tmp), depth_method = "bd")
fbplot_obj2$outliers

fbplot_obj <- functional_boxplot(t(lgp), depth_method = "mbd")
fbplot_obj$outliers

fbplot_obj2 <- functional_boxplot(t(tmp), depth_method = "mbd")
fbplot_obj2$outliers

m <- muod(t(lgp), cut_method = c("boxplot"))
m$outliers

m2 <- muod(t(tmp), cut_method = c("boxplot"))
m2$outliers

fbplot(lgp, method="BD2")
fbplot(lgp, method="MBD")

fbplot(tmp, method="BD2")
fbplot(tmp, method="MBD")

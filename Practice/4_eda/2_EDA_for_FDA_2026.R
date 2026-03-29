###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
###
### ch. 6.  Descriptions of Functional Data
###


#  load the fda package

library(fda)

##
## Section 6.1 Some Functional Descriptive Statistics
##

#   ----------  Statistics for the log precipitation data  ---------------

#  using 'logprec.fd' computed in fdarm-ch05.R as follows:

#  organize data to have winter in the center of the plot

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

#  elementary pointwise mean and standard deviation

meanlogprec   = mean.fd(logprec.fd)
stddevlogprec = std.fd(logprec.fd)

lines(meanlogprec, lwd=4, lty=2, col=2)
lines(stddevlogprec, lwd=4, lty=2, col=4)

lines(meanlogprec-stddevlogprec, lwd=4, lty=2, col=6)
lines(meanlogprec+stddevlogprec, lwd=4, lty=2, col=6)


# Section 6.1.1 The Bivariate Covariance Function v(s; t)

logprecvar.bifd = var.fd(logprec.fd)

weektime        = seq(0,365,length=53)
logprecvar_mat  = eval.bifd(weektime, weektime,
                            logprecvar.bifd)

# Figure 6.1

persp(weektime, weektime, logprecvar_mat,
      theta=-45, phi=25, r=3, expand = 0.5,
      ticktype='detailed',
      xlab="Day",
      ylab="Day",
      zlab="variance(log10 precip)")

contour(weektime, weektime, logprecvar_mat,
        xlab="Day",
        ylab="Day")

# Figure 6.2

day5time = seq(0,365,5)
logprec.varmat = eval.bifd(day5time, day5time,
                           logprecvar.bifd)
contour(day5time, day5time, logprec.varmat,
        col=terrain.colors(12),
        xlab="Day",
        ylab="Day", lwd=2,
        labcex=1)

library(fields)
image.plot(day5time, day5time, logprec.varmat)
contour(day5time, day5time, logprec.varmat,
        col= "white", add = TRUE)


###################################################
### Descriptive measures for functional data.
###################################################
library(fda.usc)
data(poblenou)
nox <- poblenou$nox
working <- poblenou$nox[poblenou$df$day.festive == 0 &
                          as.integer(poblenou$df$day.week) < 6]
nonworking <- poblenou$nox[poblenou$df$day.festive == 1 |
                             as.integer(poblenou$df$day.week) > 5]
# Centrality measures (working)

opra <- par( mfrow=c(2, 2) )
plot(func.mean(working), ylim = c(10, 170),
     main = "Centrality measures in working days")
legend(x = 11, y = 170, cex = 1, box.col = "white", lty = 1:5,
       col = c(1:5), legend = c("mean","trim.mode","trim.RP",
                                "median.mode","median.RP"))
lines(func.trim.mode(working, trim = 0.15), col = 2, lty = 2)
lines(func.trim.RP(working, trim = 0.15), col = 3, lty = 3)
lines(func.med.mode(working, trim = 0.15), col = 4, lty = 4)
lines(func.med.RP(working, trim = 0.15), col = 5, lty = 5)

# Centrality measures (non-working)
plot(func.mean(nonworking), ylim = c(10,170),
     main = "Centrality measures in non-working days")
legend(x = 11, y = 170, cex = 1, box.col = "white",lty = 1:5,
       col = c(1:5), legend = c("mean","trim.mode","trim.RP",
                                "median.mode","median.RP"))
lines(func.trim.mode(nonworking, trim = 0.15),col = 2, lty = 2)
lines(func.trim.RP(nonworking, trim = 0.15),col = 3, lty = 3)
lines(func.med.mode(nonworking, trim = 0.15),col = 4, lty = 4)
lines(func.med.RP(nonworking, trim = 0.15),col = 5, lty = 5)

# Measures of dispersion   (working)
plot(func.var(working),
     main = "Dispersion measures in working days", ylim = c(100 ,5500))
legend(x = 11, y = 5300,cex = 1, box.col = "white", lty = 1:3, col = 1:3,
       legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(working,trim = 0.15), col = 2, lty = 2)
lines(func.trimvar.RP(working,trim = 0.15), col = 3, lty = 3)

# Measures of dispersion   (non-working)
plot(func.var(nonworking),
     main = "Dispersion measures in non-working days", ylim = c(100, 5500))
legend(x = 11, y = 5300, cex = 1, box.col = "white", lty = 1:3, col = 1:3,
       legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(nonworking, trim = 0.15), col = 2, lty = 2)
lines(func.trimvar.RP(nonworking, trim = 0.15), col = 3, lty = 3)
par(opra)

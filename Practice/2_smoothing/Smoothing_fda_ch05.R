### Source : https://rdrr.io/rforge/fda/src/inst/scripts/fdarm-ch05.R

###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to
#  otherwise illustrate how R can be used in the analysis of functional
#  data.
#  We do not claim to reproduce the results in the book exactly by these
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves
#    -- we may have changed our minds about how these analyses should be
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

# Some parts of this script corrected by Jurgita Markeviciute (2021)

###
### ch. 5.  Smoothing: Computing Curves from Noisy Data
###

#  load the fda package

library(fda)

##
## Section 5.2.  Data Smoothing with Roughness Penalties
##

# section 5.2.2 The Roughness Penalty Matrix R

#  ---------------  Smoothing the Canadian weather data  ------------------

#  define a Fourier basis for daily temperature data

yearRng   = c(0,365)
nbasis    = 65
tempbasis = create.fourier.basis(yearRng,nbasis)

#  define the harmonic acceleration operator

#  When the coefficients of the linear differential operator
#  are constant, the Lfd object can be set up more simply
#  by using function vec2Lfd as follows:

#  The first argument is a vector of coefficients for the
#  operator, and the second argument is the range over which
#  the operator is defined.

harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,0), yearRng)

#Lx(t) = b_0(t) x(t) + b_1(t)Dx(t) + ... + b_{m-1}(t) D^{m-1} x(t) + D^mx(t)
constb <- create.constant.basis(yearRng)
bwlist <- vector("list", 3)
for (j in 1:3) bwlist[[j]] <- fd(c(0,(2*pi/365)^2,0)[j], constb)
bwlist
lfdobj <- Lfd(3, bwlist)

#  compute the penalty matrix R
# R = integral phi(t) phi^T(t) dt
# Estimator = (Phi^T Phi + lanbda R)^{-1} Phi^T y

Rmat = eval.penalty(tempbasis, harmaccelLfd)

# section 5.2.4 Defining Smoothing by Functional Parameter Objects

##
## 5.3.  Case Study: The Log Precipitation Data
##

#  organize data to have winter in the center of the plot

dayOfYearShifted = c(182:365, 1:181)

dim(CanadianWeather$dailyAv)

logprecav = CanadianWeather$dailyAv[dayOfYearShifted, , 'log10precip']

#  set up a saturated basis: as many basis functions as observations

nbasis   = 365
daybasis = create.fourier.basis(yearRng, nbasis)

#  set up the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

#  step through values of log(lambda)

loglam        = seq(4,9,0.25)
nlam          = length(loglam)
dfsave        = rep(NA,nlam)
names(dfsave) = loglam
gcvsave       = dfsave
for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda        = 10^loglam[ilam]
  fdParobj      = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist    = smooth.basis(day.5, logprecav,
                               fdParobj)
  dfsave[ilam]  = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

# Figure 5.2.

plot(loglam, gcvsave, type='b', lwd=2, ylab='GCV Criterion',
     xlab=expression(log[10](lambda)) )

#  smooth data with minimizing value of lambda

lambda      = 1e6
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprec.fit$fd
fdnames     = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  plot the functional data object

plot(logprec.fd, lwd=2)

fulfd <- smooth.basis(day.5, logprecav, daybasis)

opar <- par(mfrow=c(2,1))
plot(logprec.fd, lwd=2)
plot(fulfd$fd, lwd=2)
par(opar)


# Compare with actual values

plotfit.fd(logprecav, day.5, logprec.fd)

# plotfit.fd:  Pauses between plots
# *** --->>> input required (e.g., click on the plot)
#            to advance to the next plot

##
## Section 5.4 Positive, Monotone, Density
##             and Other Constrained Functions
##

#   ----------------  Positive smoothing of precipitation  ----------------

lambda      = 1e3
WfdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
VanPrec     = CanadianWeather$dailyAv[
  dayOfYearShifted, 'Vancouver', 'Precipitation.mm']
VanPrecPos  = smooth.pos(day.5, VanPrec, WfdParobj)
Wfd         = VanPrecPos$Wfdobj
Wfd$fdnames = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")

precfit = exp(eval.fd(day.5, Wfd))

plot(day.5, VanPrec, type="p", cex=1.2,
     xlab="Day (July 1 to June 30)",
     ylab="Millimeters",
     main="Vancouver's Precipitation")
lines(day.5, precfit,lwd=2)

#  ------------  5.4.2.1  Monotone smoothing  of the tibia data  -----------

#  set up the data for analysis

day    = infantGrowth[, 'day']
tib    = infantGrowth[, 'tibiaLength']
n      = length(tib)

#  a basis for monotone smoothing

nbasis = 42
Wbasis   = create.bspline.basis(c(1,n), nbasis)

#  the fdPar object for smoothing

Wfd0     = fd(matrix(0,nbasis,1), Wbasis)
WfdPar   = fdPar(Wfd0, 2, 1e-4)

#  smooth the data

result   = smooth.monotone(day, tib, WfdPar)
Wfd      = result$Wfd
beta     = result$beta

#  compute fit and derivatives of fit

dayfine  = seq(1,n,len=151)
tibhat   = beta[1]+beta[2]*eval.monfd(dayfine ,Wfd)
Dtibhat  =        beta[2]*eval.monfd(dayfine, Wfd, 1)
D2tibhat =        beta[2]*eval.monfd(dayfine, Wfd, 2)

#  plot height

op = par(mfrow=c(3,1), mar=c(5,5,3,2), lwd=2)

plot(day, tib, type = "p", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Length (mm)")
lines(dayfine, tibhat, lwd=2)

#  plot velocity

plot(dayfine, Dtibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Velocity (mm/day)")

#  plot acceleration

plot(dayfine, D2tibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Acceleration (mm/day/day)")
lines(c(1,n),c(0,0),lty=2)

par(op)

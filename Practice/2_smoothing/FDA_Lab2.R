# =======================================================
# Graphics / Debugging modes
# =======================================================

# Choose ONE mode by uncommenting it. !! Add dev.off() at the end of the script
# for all modes except MODE 1 !!

# ------------------------------------------------
# MODE 1 — Interactive (default, RStudio-friendly)
# ------------------------------------------------
# Do nothing. Plots will appear in RStudio Plots pane.

# ------------------------------------------------
# MODE 2 — External window (clean, fast, no RStudio quirks)
# ------------------------------------------------
# graphics.off()  # Reset all devices (recommended before opening a new one)
# x11()           # For Linux
# # quartz()      # For macOS
# # windows()     # For Windows

# ------------------------------------------------
# MODE 3 — Non-interactive (save plots to file)
# ------------------------------------------------
# graphics.off()
# png("out.png")  # Redirect all plots, overwriting, to a file (no GUI needed)

# ------------------------------------------------
# MODE 4 — Disable plotting (fastest debugging)
# ------------------------------------------------
# graphics.off()
# pdf(NULL)

# ------------------------------------------------
# Common settings (apply in any mode if needed)
# ------------------------------------------------
# par(ask = FALSE)                     # Disable "Hit <Return>" prompts
# options(device.ask.default = FALSE)  # Global fallback
# grDevices::devAskNewPage(FALSE)      # Force-disable paging on current device

# NOTE:
# -----
# Some functions (e.g. in `fda`) override paging with `ask = TRUE`.
# If "Hit <Return>" prompts still appear, pass: ask = FALSE (e.g.
# plotfit.fd(..., ask = FALSE))

####################################################################
# Smoothing Lab: Smoothing Tasks
####################################################################

library(fda)

####################################################################
# TASK 1: Growth data
####################################################################
# (1.) basis systems for growth data
#  ------------------- Smoothing the growth data ---------------

#   define the range of the ages and set up a fine mesh of ages

ageRng  = c(1,18)
age     = growth$age
agefine = seq(1,18,len=501)

#   set up order 6 spline basis with 12 basis functions for
#   fitting the growth data so as to estimate acceleration

nbasis = 12
norder = 6
heightbasis12 = create.bspline.basis(ageRng, nbasis, norder)
plot(heightbasis12)

#   fit the data by least squares

basismat   = eval.basis(age, heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(basismat, heightmat, intercept = FALSE)$coef

fdnames      = vector('list', 3)
fdnames[[1]] = "Age (years)"
fdnames[[2]] = "Child"
fdnames[[3]] = "Height (cm)"

hgtfd <- fd(heightcoef, heightbasis12, fdnames)
plot(hgtfd)
str(hgtfd)
names(hgtfd)

#   fit the data using the function smooth_basis, which does the same thing.

heightList = smooth.basis(age, heightmat, heightbasis12)
heightfd   = heightList$fd
plot(heightfd)
str(heightfd)
names(heightfd)

height.df = heightList$df
height.df

height.gcv = heightList$gcv
height.gcv

heightbasismat = eval.basis(age, heightbasis12)
y2cMap         = solve(crossprod(heightbasismat)) %*% t(heightbasismat)

coefest <- y2cMap %*% heightmat
dcoef <- coefest - heightcoef
dcoef

# Hat matrix H
# Phi (Phi^T Phi)^{-1} Phi^T
# with penalty:  Phi  (Phi^T Phi + lambda R)^{-1} Phi^T
hatH = heightbasismat %*%
  solve(crossprod(heightbasismat)) %*% t(heightbasismat)

# Differential operators

# The first thing we need to do is to define a linear differential operator
# (Lfd) object. This is a list of functional data objects defining the terms
# on the right hand side of
#
#   D^m x + b_{m-1}(t) D^{m-1} x + ... + b_0(t) x

?Lfd

# The simplest Lfds just penalize the derivative x, and we can define
# them by

D2lfd = int2Lfd(2)
D2lfd

# This says that b_0 = b_1 = 0, which we can see by looking internally at
# D2lfd we see

names(D2lfd)

plot(D2lfd$bwtlist[[2]])

# Instead of evaluating a derivative, we can always evaluate an Lfd
# applied to a function

plot(hgtfd, Lfdobj=D2lfd)
plot(deriv.fd(hgtfd,2))

# Now we need to do some smoothing.
# In order to set up the smoothing operator, we need to define an
# fdPar (functional parameter) object. This is just a list of various
# quantities that avoids the need for many argument values

D2fdPar = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e4)

# The functional parameter object holds a basis, a Lfd and a value for lambda.
# We can now make a new call to

syfd = smooth.basis(age, heightmat, D2fdPar)

plotfit.fd(heightmat,age,syfd$fd)
plot(syfd) # to strait line

# Now we can examine how this changes with lambda. It is most useful to
# vary lambda on the logarithmic scale. Note that large values of
# lambda mean more smoothing

# We'll also keep track of gcv, df and sse

gcv = rep(0,21)
df = rep(0,21)
sse = rep(0,21)

for (i in 1:21) {
  lambda = 10^(i-10)
  tD2fdPar = fdPar(heightbasis12,Lfdobj = int2Lfd(2),lambda=lambda)

  tyfd = smooth.basis(age, heightmat, tD2fdPar)

  gcv[i] = sum(tyfd$gcv)
  df[i] = tyfd$df
  sse[i] = tyfd$SSE
}

# And we'll plot some results

plot(-10:0,df[1:11],type='l',xlab='log lambda',ylab='df',cex.lab=1.5)
plot(-10:0,sse[1:11],type='l',xlab='log lambda',ylab='sse',cex.lab=1.5)
plot(-10:0,gcv[1:11],type='l',xlab='log lambda',ylab='gcv',cex.lab=1.5)

# Let's look at some influence functions

y2cMap = syfd$y2cMap

bvals = eval.basis(age,heightbasis12)

hatMat = bvals%*%y2cMap

matplot(age,hatMat,type='l',xlab='time',ylab='influence',
        cex.lab=1.5,cex.axis=1.5)

# We'll do the same thing with a very small lambda

D2fdPar1 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e-6)
syfd1 = smooth.basis(age,heightmat,D2fdPar1)
plotfit.fd(heightmat,age,syfd1$fd)

y2cMap1 = syfd1$y2cMap
hatMat1 = bvals %*% y2cMap1

matplot(age,hatMat1,type='l',xlab='time',ylab='influence',
        cex.lab=1.5,cex.axis=1.5)

# Or very large lambda

D2fdPar2 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e10)
syfd2 = smooth.basis(age,heightmat,D2fdPar2)
plotfit.fd(heightmat,age,syfd2$fd)

y2cMap2 = syfd2$y2cMap
hatMat2 = bvals%*%y2cMap2

matplot(age,hatMat2,type='l',xlab='time',ylab='influence',
        cex.lab=1.5,cex.axis=1.5)

# Or min GCV

D2fdPar3 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=0.001)
syfd3 = smooth.basis(age,heightmat,D2fdPar3)
plotfit.fd(heightmat,age,syfd3$fd)

y2cMap3 = syfd3$y2cMap
hatMat3 = bvals%*%y2cMap3

matplot(age,hatMat3,type='l',xlab='time',ylab='influence',
        cex.lab=1.5,cex.axis=1.5)

# ------- Smoothing the growth data with a roughness penalty -----------

#   set up a basis for the growth data
#   with knots at ages of height measurement

norder      = 6
nbasis      = length(age) + norder - 2
heightbasis = create.bspline.basis(ageRng, nbasis, norder, age)

heightList = smooth.basis(age, heightmat, heightbasis12)
heightfd   = heightList$fd

# define a functional parameter object for smoothing

heightLfd    = 4
heightlambda = 0.01
heightfdPar  = fdPar(heightbasis, heightLfd, heightlambda)

# smooth the data

heightfdSmooth = smooth.basis(age, heightmat, heightfdPar)
heightfdpen       = heightfdSmooth$fd

opar <- par(mfrow=c(3,2))
plot(heightfd)
plot(heightfdpen)
plot(heightfd, Lfdobj = 1)
plot(heightfdpen, Lfdobj = 1)
plot(heightfd, Lfdobj = 2)
plot(heightfdpen, Lfdobj = 2)
par(opar)

# section 5.2.5 Choosing Smoothing Parameter lambda

loglam         = seq(-6, 0, 0.25)
Gcvsave        = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave         = Gcvsave
for (i in 1:length(loglam)) {
  hgtfdPari  = fdPar(heightbasis, Lfdobj=4, 10^loglam[i])
  hgtSm.i    = smooth.basis(age, heightmat, hgtfdPari)
  Gcvsave[i] = sum(hgtSm.i$gcv)
  Dfsave[i]  = hgtSm.i$df
}

# Figure 5.1.

par(mfrow = c(1, 1))
plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )
abline(h=min(Gcvsave), col=2, lty=2)
abline(v=loglam[which.min(Gcvsave)], col=2, lty=2)
abline(v=loglam[which.min(Gcvsave)+1], col=4, lty=2)
abline(v=loglam[which.min(Gcvsave)-2], col=4, lty=2)

####################################################################
# TASK 2: Lfd operator
####################################################################

####################################################
## The Linear Differential Operator or Lfd Class ###
####################################################

# Define a B-spline basis over one year (0–365 days).
# This is needed because the Lfd construction requires a basis object,
# and specifically its range (thawbasis$rangeval).
# The choice of nbasis/norder here is arbitrary for testing purposes.
thawbasis = create.bspline.basis(c(0, 365), nbasis = 25, norder = 4)

omega           = 2*pi/365
thawconst.basis = create.constant.basis(thawbasis$rangeval)

betalist      = vector("list", 3)
betalist[[1]] = fd(0, thawconst.basis)
betalist[[2]] = fd(omega^2, thawconst.basis)
betalist[[3]] = fd(0, thawconst.basis)
harmaccelLfd  = Lfd(3, betalist)

accelLfd = int2Lfd(2)

harmaccelLfd.thaw = vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
all.equal(harmaccelLfd[-1], harmaccelLfd.thaw[-1])

class(accelLfd)
class(harmaccelLfd)

# Create a functional data object (fd) from synthetic temperature data.
# A sinusoidal signal is generated with yearly periodicity (matching omega),
# evaluated at day.5 (midpoints of days of the year).
# smooth.basis() converts discrete observations into basis coefficients.
temp.fd = smooth.basis(
  day.5,
  sin(2 * pi * day.5 / 365),
  thawbasis
)$fd

Ltempmat   = eval.fd(day.5, temp.fd, harmaccelLfd)

D2tempfd = deriv.fd(temp.fd, 2)
Ltempfd  = deriv.fd(temp.fd, harmaccelLfd)

####################################################################
# TASK 3: Monotone smoothing the Berkeley female data
####################################################################

##
## Compute the monotone smoothing of the Berkeley female growth data.
##

#   set up ages of measurement and an age mesh

age     = growth$age
nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length = nfine)

#   the data

hgtf   = growth$hgtf
ncasef = dim(hgtf)[2]

#   an order 6 bspline basis with knots at ages of measurement

norder = 6
nbasis = nage + norder - 2
wbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#   define the roughness penalty for function W

Lfdobj    = 3          #   penalize curvature of acceleration
lambda    = 10^(-0.5)  #   smoothing parameter
cvecf     = matrix(0, nbasis, ncasef)
Wfd0      = fd(cvecf, wbasis)
growfdPar = fdPar(Wfd0, Lfdobj, lambda)

#   monotone smoothing

growthMon = smooth.monotone(age, hgtf, growfdPar)

# (wait for an iterative fit to each of 54 girls)

Wfd       = growthMon$Wfd
betaf     = growthMon$beta
hgtfhatfd = growthMon$yhatfd

#   Set up functional data objects for the acceleration curves
#   and their mean.   Suffix UN means "unregistered".

accelfdUN     = deriv.fd(hgtfhatfd, 2)
accelmeanfdUN = mean.fd(accelfdUN)

#   plot unregistered curves

par(ask = FALSE)
plot(accelfdUN, xlim=ageRng, ylim=c(-4,3), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Acceleration (cm/yr/yr)")

#dev.off()
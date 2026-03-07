###################################################
### tecator dataset
###################################################
library("fda.usc")
data("tecator")
names(tecator)
absorp <- tecator$absorp.fdata
Fat20 <- ifelse(tecator$y$Fat < 20, 0, 1) * 2 + 2
absorp$names$main <- ""

opar <- par(mfrow = c(1 , 2))
# Figure 1 (left panel)
plot(absorp, col = Fat20)
absorp.d1 <- fdata.deriv(absorp, nderiv = 1)
# Figure 1 (right panel)
plot(absorp.d1, col = Fat20)
par(opar)

###################################################
### convert fdata class to fd class
###################################################
class( absorp.fd <- fdata2fd( absorp, type.basis = "fourier", nbasis = 15) )
class( absorp.fdata <- fdata(absorp.fd) )

###################################################
### phoneme data and smoothing
###################################################
data(phoneme)
plot(phoneme$test)

# parametric smoothing
learn <- phoneme$learn
l <- c(0 ,2 ^ seq(-2, 9, len = 30))
nb <- seq(7, 31, by = 2)

# by default "bspline" is used
out0 <- optim.basis(learn, lambda = l, numbasis = nb)
out0$gcv.opt
out0$numbasis.opt

opar <- par(mfrow=c(1,2))
plot(out0$fdataobj)
plot(l, out0$gcv["27",], type = "l", 
     xlab = "Smoothing parameter", ylab = "GCV",
     main = "GCV for 27 basis functions")
par(opar)

# all GCV curves
matplot(nb, out0$gcv, type = "l")

# Non-parametric smoothing
# Nadaraya-Watson kernel estimator
# Ker.norm
out1 <- optim.np(learn, type.S = S.NW, par.CV = list(criteria = "GCV"))

out1$gcv.opt # optimal GCV
out1$h.opt # optimal window h

# Cosine kernel, window grid h taken from normal kernel
# estimated earlier for comparison
out1a <- optim.np(learn, h = out1$h,
                  Ker = Ker.cos, 
                  type.S = S.NW, 
                  par.CV = list(criteria = "GCV"),
                  correl = FALSE,
                  verbose = TRUE)

# Epanechnikov kernel, window grid h taken from normal kernel
# estimated earlier for comparison
out1b <- optim.np(learn, h = out1$h,
                  Ker = Ker.epa, 
                  type.S = S.NW, 
                  par.CV = list(criteria = "GCV"),
                  correl = FALSE,   
                  verbose = TRUE)

# Triweight kernel, window grid h taken from normal kernel
# estimated earlier for comparison
out1c <- optim.np(learn, h = out1$h,
                  Ker = Ker.tri, 
                  type.S = S.NW, 
                  par.CV = list(criteria = "GCV"),
                  correl = FALSE,   
                  verbose = TRUE)

# Quartic kernel, window grid h taken from normal kernel
# estimated earlier for comparison
out1d <- optim.np(learn, h = out1$h,
                  Ker = Ker.quar, 
                  type.S = S.NW, 
                  par.CV = list(criteria = "GCV"),
                  correl = FALSE,   
                  verbose = TRUE)

# Uniform kernel, window grid h taken from normal kernel
# estimated earlier for comparison
out1e <- optim.np(learn, h = out1$h,
                  Ker = Ker.unif, 
                  type.S = S.NW, 
                  par.CV = list(criteria = "GCV"),
                  correl = FALSE,   
                  verbose = TRUE)

# Local Linear Smoothing with bandwidth parameter h.
# Normal kernel
out2 <- optim.np(learn, type.S = S.LLR, par.CV = list(criteria = "GCV"))

out2$gcv.opt
out2$h.opt

###################################################
### plot GCV criteria
###################################################
opar <- par(mfrow = c(1,2))
contour(nb, l, out0$gcv, ylab = "Lambda", xlab = "Number of basis", 
        main = "GCV criteria by optim.basis()")
plot(out1$h, out1$gcv, type = "l",, ylim = c(1.5, 9.5),
     main = "GCV criteria  by optim.np() ", 
     xlab = "Bandwidth (h) values",ylab = "GCV criteria", col = 1, lwd = 2)
legend(x = 3, y = 9, legend = c("Ker.norm-S.NW", "Ker.norm-S.LLR",
                                "Ker.cos-S.NW", "Ker.epa-S.NW",
                                "Ker.tri-S.NW", "Ker.quar-S.NW",
                                "Ker.unif-S.NW"),
       box.col = "white", lwd = c(2, 2, 2), col = 1:7,cex = 0.75)
lines(out2$h,out2$gcv, col = 2, lwd = 2)
lines(out1a$h,out1a$gcv, col = 3, lwd = 2)
lines(out1b$h,out1b$gcv, col = 4, lwd = 2)
lines(out1c$h,out1c$gcv, col = 5, lwd = 2)
lines(out1d$h,out1d$gcv, col = 6, lwd = 2)
lines(out1d$h,out1d$gcv, col = 7, lwd = 2)
par(opar)

###################################################
### smoothing a fdata curve
###################################################
library(RColorBrewer)

cols <- brewer.pal(9, "Set1")
cols

ind <- 11
nam <- expression( paste("Phoneme curve"[11]) )
par(mar = c(5, 4, 4, 10))
plot(learn[ind, ], main = nam, lty = 2, lwd = 2, col = 8)
legend(x = 160, y = 20,, legend = c("Curve","Bspline basis",
                                  "Ker.norm-S.NW", "Ker.norm-S.LLR",
                                    "Ker.cos-S.NW", "Ker.epa-S.NW",
                                    "Ker.tri-S.NW", "Ker.quar-S.NW",
                                    "Ker.unif-S.NW"),
       lty = c(2, rep(1,8)), lwd = 2, col = c(8, cols), box.col = "white",
       xpd = TRUE)
lines(out0$fdata.est[ind, ], col = cols[1], lty = 1, lwd = 2)
lines(out1$fdata.est[ind, ], col = cols[2], lty = 1, lwd = 2)
lines(out2$fdata.est[ind, ], col = cols[3], lty = 1, lwd = 2)
lines(out1a$fdata.est[ind, ], col = cols[4], lty = 1, lwd = 2)
lines(out1b$fdata.est[ind, ], col = cols[5], lty = 1, lwd = 2)
lines(out1c$fdata.est[ind, ], col = cols[6], lty = 1, lwd = 2)
lines(out1d$fdata.est[ind, ], col = cols[7], lty = 1, lwd = 2)
lines(out1e$fdata.est[ind, ], col = cols[8], lty = 1, lwd = 2)


###################################################
### Calculation of the smoothing parameter (h) 
### for a functional data using nonparametric 
### kernel estimation.
###################################################
hgrid <- h.default(learn, prob = c(0.025, 0.5), len = 100)
hgrid

hgrid <- h.default(learn, prob = c(0.025, 0.25), len = 100)
hgrid

hgrid <- h.default(learn, prob = c(0.025, 0.25), len = 100,
                   Ker = Ker.cos)
hgrid

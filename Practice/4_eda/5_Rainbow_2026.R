###################################
# load the package
library("rainbow")

# plot.type = "function", curves are plotted by time
# the most recent curve is shown in purple
# the distant past cure is shown in red
par(mfrow=c(1,2))
plot(Australiafertility, plot.type = "functions", plotlegend = TRUE)
plot(ElNino_ERSST_region_1and2, plot.type = "functions", plotlegend = TRUE)

# plot.type="depth", curves are plotted by halfspace depth
# depth is distance between median and each curve
# median curve shown in black line is the center
plot(ElNino_ERSST_region_1and2, plot.type="depth",plotlegend=TRUE)
plot(Australiafertility, plot.type="depth",plotlegend=TRUE)

# plot.type = "bivariate", the bivariate principal component scores are displayed
# type = "bag" requests the bagplot
# It uses Tukey's halfspace location depths
# Central point Tukey median, an inner region (the "bag")
# and an outer region, beyond which outliers are shown as individual points.
# The  bag  is  defined  as  the  smallest  depth  region  containing
# at  least  50%  of  the  total  number of  observations.
# The  outer  region  (or  "fence")  of  the  bagplot  is  the  convex  hull  of  the  points contained  within  the  region  obtained
# by  inflating  the  bag  (relative  to  the  Tukey  median)  by  a factor rho.
# rho = 2.58, as that will allow the fence to contain
# 99% of the observations when the projected bivariate scores follow
# standard normal distributions.
fboxplot(ElNino_ERSST_region_1and2, plot.type="bivariate",
         type="bag", projmethod = "PCAproj",
         ylim=c(-10,20), xlim=c(-15,20))
fboxplot(Australiafertility, plot.type="bivariate",
         type="bag", projmethod = "PCAproj",
         ylim=c(-200,200), xlim=c(-400,300))

# plot.type = "functional", the bivariate pc scores are matched to corresponding curves
# The functional bagplot is a mapping of the bagplot of the first two robust principal componentscores to the functional curves.
# The functional bagplot displays the median curve (the curve with the greatest depth), and the inner and outer regions.
# The inner region is defined as the region bounded  by  all  curves  corresponding  to  points  in  the  bivariate  bag.
# Thus,  50%  of  curves  are in  the  inner  region.   The  outer  region  is  similarly  defined  as  the  region
# bounded  by  all  curves corresponding to points within the bivariate fence region.
fboxplot(ElNino_ERSST_region_1and2,
         plot.type = "functional", type = "bag", projmethod = "PCAproj")
fboxplot(Australiafertility, plot.type = "functional",
         type = "bag", projmethod = "PCAproj")


# type = "hdr" requests the HDR boxplot
# alpha requests the coverage probability of inner
# and outer HDR regions, customarily c(0.05,0.5)
# The functional HDR boxplot displays themodal  curve
# (the  curve  with  the  highest  density),  and  the  inner  and
# outer  regions.
# The  inner region  is  defined  as  the  region  bounded  by
# all  curves  corresponding  to  points  inside  the  50% bivariate HDR.
# Thus, 50% of curves are in the inner region.
# The outer region is similarly definedas the region bounded
# by all curves corresponding
# to the points within the outer bivariate HDR.
# As  with  any  outlier  detection  method,
# including  bagplots  and  HDR  boxplots,  the  coverageprobability of the outer region needs to be pre-specified.
fboxplot(ElNino_ERSST_region_1and2, plot.type="bivariate",
         type="hdr", alpha=c(0.07,0.5),
         projmethod="PCAproj", ylim=c(-10,20), xlim=c(-15,20))
fboxplot(ElNino_ERSST_region_1and2, plot.type="bivariate",
         type="hdr", alpha=c(0.05,0.5),
         projmethod="PCAproj", ylim=c(-10,20), xlim=c(-15,20))
fboxplot(Australiafertility, plot.type="bivariate",
         type="hdr", alpha=c(0.07,0.5),
         projmethod="PCAproj", ylim=c(-200,200), xlim=c(-400,300))
fboxplot(Australiafertility, plot.type="bivariate",
         type="hdr", alpha=c(0.05,0.5),
         projmethod="PCAproj", ylim=c(-200,200), xlim=c(-400,300))

fboxplot(ElNino_ERSST_region_1and2, plot.type = "functional",
         type = "hdr", alpha = c(0.07,0.5),
         projmethod="PCAproj")
fboxplot(ElNino_ERSST_region_1and2, plot.type = "functional",
         type = "hdr", alpha = c(0.05,0.5),
         projmethod="PCAproj")
fboxplot(Australiafertility, plot.type = "functional",
         type = "hdr", alpha = c(0.07,0.5),
         projmethod="PCAproj")
fboxplot(Australiafertility, plot.type = "functional",
         type = "hdr", alpha = c(0.05,0.5),
         projmethod="PCAproj")

# Functional outliers
foutliers(ElNino_ERSST_region_1and2, method = "robMah")
foutliers(ElNino_ERSST_region_1and2, method = "lrt")
foutliers(ElNino_ERSST_region_1and2, method = "depth.trim")
foutliers(ElNino_ERSST_region_1and2, method = "depth.pond")
foutliers(ElNino_ERSST_region_1and2, method = "HUoutliers")

foutliers(Australiafertility, method = "robMah")
foutliers(Australiafertility, method = "lrt")
foutliers(Australiafertility, method = "depth.trim")
foutliers(Australiafertility, method = "depth.pond")
foutliers(Australiafertility, method = "HUoutliers")


######################################################################
# temperature data
library(fda)
library(fda.usc)

tt=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
plot(fdataobj)

l <- c(0 ,2 ^ seq(-2, 9, len = 30))
nb <- seq(7, 31, by = 2)
fdataobj2 <- optim.basis(fdataobj, lambda = l, numbasis = nb)
plot(fdataobj2$fdata.est)
matplot(t(fdataobj2$fdata.est$data), type="l")

wth1 <- fds(x = tt, y = CanadianWeather$dailyAv[,,1],
           xname = "Day", yname = "Temperature")
wth1

dt <- as.matrix(t(fdataobj2$fdata.est$data))
colnames(dt) <- colnames(CanadianWeather$dailyAv[,,1])
rownames(dt) <- rownames(CanadianWeather$dailyAv[,,1])

wth2 <- fds(x = tt, y = dt,
            xname = "Day", yname = "Temperature")
wth2

op <- par(mfrow=c(1,2))
plot(wth1, plot.type = "functions", plotlegend = TRUE)
plot(wth2, plot.type = "functions", plotlegend = TRUE)
par(op)

plot(wth2, plot.type = "functions", plotlegend = TRUE)
plot(wth2, plot.type = "depth", plotlegend=TRUE)

fboxplot(wth2, plot.type="bivariate",
         type="bag", projmethod = "PCAproj",
         ylim=c(-200,200), xlim=c(-400,300))
fboxplot(wth2, plot.type = "functional",
         type = "bag", projmethod = "PCAproj")

fboxplot(wth2, plot.type="bivariate",
         type="hdr", alpha=c(0.07,0.5),
         projmethod="PCAproj",
         ylim=c(-100,300), xlim=c(-300,500))

fboxplot(wth2, plot.type="bivariate",
         type="hdr", alpha=c(0.05,0.5),
         projmethod="PCAproj",
         ylim=c(-100,300), xlim=c(-300,500))

fboxplot(wth2, plot.type = "functional",
         type = "hdr", alpha = c(0.07,0.5),
         projmethod="PCAproj")
fboxplot(wth2, plot.type = "functional",
         type = "hdr", alpha = c(0.05,0.5),
         projmethod="PCAproj")


foutliers(wth2, method = "robMah")
foutliers(wth2, method = "lrt")
foutliers(wth2, method = "depth.trim")
foutliers(wth2, method = "depth.pond")
foutliers(wth2, method = "HUoutliers")


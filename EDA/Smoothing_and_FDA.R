library(readxl)

# Create output directory if it doesn't exist
out_dir <- file.path("./EDA/outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the data
dat <- read_excel("./EDA/Flanker_stimulus_FC1_channel.xlsx")

time_vec <- dat$time
Y_mat    <- as.matrix(dat[, -1])

pdf(file.path(out_dir, "00_all_curves.pdf"), width = 10, height = 6)
matplot(time_vec, Y_mat, type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.4),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "Flanker Task — Channel FC1 — All 39 Subjects")
abline(v = 0, lty = 2, col = "red", lwd = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topleft",
       legend = c("Individual subjects", "Stimulus onset"),
       col = c("steelblue", "red"),
       lty = c(1, 2), lwd = c(1, 2), cex = 0.8)
dev.off()


################################################################################
# Step 1: Smoothing
# Flanker Task EEG — Channel FC1 — 39 Subjects
################################################################################

library(fda)
library(readxl)

# ==============================================================================
# Load data
# ==============================================================================

dat      <- read_excel("./EDA/Flanker_stimulus_FC1_channel.xlsx")

time_vec <- dat$time
Y_mat    <- as.matrix(dat[, -1])

cat("Time points:", length(time_vec), "\n")
cat("Time range: ", range(time_vec), "\n")
cat("Subjects:   ", ncol(Y_mat), "\n")

# ==============================================================================
# Set up B-spline basis
# ==============================================================================

time_range <- range(time_vec)
norder     <- 4                          # cubic B-splines
nbasis     <- 20
basis_obj  <- create.bspline.basis(time_range, nbasis = nbasis, norder = norder)

cat("Basis functions:", nbasis, "\n")

# ==============================================================================
# GCV to choose lambda
# (from Smoothing_fda_ch05.R lines 110-123)
# ==============================================================================

Lfdobj  <- 2
loglam  <- seq(-5, 10, by = 0.5)
nlam    <- length(loglam)
dfsave  <- rep(NA, nlam)
gcvsave <- rep(NA, nlam)
names(dfsave)  <- loglam
names(gcvsave) <- loglam

for (ilam in 1:nlam) {
  cat(paste("log10 lambda =", loglam[ilam], "\n"))
  lambda     <- 10^loglam[ilam]
  fdParobj   <- fdPar(basis_obj, Lfdobj, lambda)
  smoothlist <- smooth.basis(time_vec, Y_mat, fdParobj)
  dfsave[ilam]  <- smoothlist$df
  gcvsave[ilam] <- sum(smoothlist$gcv)
}

best_idx    <- which.min(gcvsave)
best_lambda <- 10^loglam[best_idx]

cat("\n=== GCV RESULT ===\n")
cat("Optimal log10(lambda):", loglam[best_idx], "\n")
cat("Optimal lambda:       ", best_lambda, "\n")
cat("Optimal df:           ", dfsave[best_idx], "\n")
cat("==================\n")

# Plot GCV curve
pdf(file.path(out_dir, "01_gcv_lambda.pdf"), width = 10, height = 6)
plot(loglam, gcvsave, type = "b", lwd = 2,
     ylab = "GCV Criterion",
     xlab = expression(log[10](lambda)),
     main = "GCV for Smoothing Parameter Selection")
abline(v = loglam[best_idx], lty = 2, col = "red", lwd = 2)
text(loglam[best_idx], max(gcvsave) * 0.9,
     labels = paste0("optimal: ", loglam[best_idx]),
     col = "red", pos = 4)
dev.off()
cat("Saved: 01_gcv_lambda.pdf\n")

# ==============================================================================
# Smooth with optimal lambda
# ==============================================================================

fdParobj_opt  <- fdPar(basis_obj, Lfdobj, best_lambda)
smooth_result <- smooth.basis(time_vec, Y_mat, fdParobj_opt)
fd_smooth     <- smooth_result$fd

# ==============================================================================
# Plot: before vs after smoothing
# ==============================================================================

pdf(file.path(out_dir, "02_before_after_smoothing.pdf"), width = 12, height = 6)
opar <- par(mfrow = c(1, 2))

matplot(time_vec, Y_mat, type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.4),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "Before Smoothing (raw)")
abline(v = 0, lty = 2, col = "red", lwd = 2)

plot(fd_smooth, col = adjustcolor("steelblue", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = paste0("After Smoothing (lambda = ", format(best_lambda, digits = 2), ")"))
abline(v = 0, lty = 2, col = "red", lwd = 2)

par(opar)
dev.off()
cat("Saved: 02_before_after_smoothing.pdf\n")

# ==============================================================================
# Plot: one subject close-up (raw vs smoothed)
# ==============================================================================

pdf(file.path(out_dir, "03_single_subject_fit.pdf"), width = 10, height = 6)
sub_idx <- 2  # SUB2
y_raw   <- Y_mat[, sub_idx]
y_fit   <- eval.fd(time_vec, fd_smooth[sub_idx])

plot(time_vec, y_raw, type = "l", col = "grey60",
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = paste0("Subject ", sub_idx, ": Raw vs Smoothed"))
lines(time_vec, y_fit, col = "steelblue", lwd = 2)
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Raw", "Smoothed"),
       col = c("grey60", "steelblue"),
       lty = 1, lwd = c(1, 2), cex = 0.8)
dev.off()
cat("Saved: 03_single_subject_fit.pdf\n")

# ==============================================================================
# Save the fd object for later use in EDA
# ==============================================================================

saveRDS(fd_smooth, file.path(out_dir, "fd_smooth.rds"))
cat("\nSmoothed fd object saved to:", file.path(out_dir, "fd_smooth.rds"), "\n")
cat("Use fd_smooth <- readRDS(file.path(out_dir, 'fd_smooth.rds')) in the EDA script.\n")
















library(fda)
library(fda.usc)
library(fdaoutlier)
library(rainbow)
library(fields)

# Load the smoothed fd object from Step 1
fd_smooth <- readRDS(file.path(out_dir, "fd_smooth.rds"))

# Evaluate on a fine grid (for fdaoutlier and rainbow functions)
t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 501)
Y_eval <- eval.fd(t_fine, fd_smooth)



# ==============================================================================
# PART 1: MEAN, SD, COVARIANCE
# (from 2_EDA_for_FDA2026.R lines 45-52)
# ==============================================================================

meanfd <- mean.fd(fd_smooth)
sdfd   <- std.fd(fd_smooth)

pdf(file.path(out_dir, "04_mean_sd.pdf"), width = 10, height = 6)
plot(fd_smooth, col = "gray", lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = "Mean and Standard Deviation")
lines(meanfd, lwd = 4, lty = 2, col = 2)
lines(sdfd,   lwd = 4, lty = 2, col = 4)
lines(meanfd - sdfd, lwd = 4, lty = 2, col = 6)
lines(meanfd + sdfd, lwd = 4, lty = 2, col = 6)
abline(v = 0, lty = 2, col = "black")
legend("topright",
       legend = c("Curves", "Mean", "SD", "Mean +/- SD"),
       col = c("gray", 2, 4, 6), lty = c(1, 2, 2, 2),
       lwd = c(1, 4, 4, 4), cex = 0.7)
dev.off()

# --- Covariance ---
# (from 2_EDA_for_FDA2026.R lines 57-90)

covbifd <- var.fd(fd_smooth)
t_grid  <- seq(fd_smooth$basis$rangeval[1],
               fd_smooth$basis$rangeval[2], length.out = 101)
cov_mat <- eval.bifd(t_grid, t_grid, covbifd)

pdf(file.path(out_dir, "05_covariance_persp.pdf"), width = 10, height = 8)
persp(t_grid, t_grid, cov_mat,
      theta = -45, phi = 25, r = 3, expand = 0.5,
      ticktype = "detailed",
      xlab = "Time s", ylab = "Time t", zlab = "Cov(s,t)",
      main = "Covariance Surface")
dev.off()

pdf(file.path(out_dir, "06_covariance_contour.pdf"), width = 10, height = 8)
contour(t_grid, t_grid, cov_mat,
        xlab = "Time s", ylab = "Time t",
        main = "Covariance Contour", lwd = 2)
dev.off()

pdf(file.path(out_dir, "07_covariance_image.pdf"), width = 10, height = 8)
image.plot(t_grid, t_grid, cov_mat,
           xlab = "Time s", ylab = "Time t",
           main = "Covariance Surface")
contour(t_grid, t_grid, cov_mat, col = "white", add = TRUE)
dev.off()


# ==============================================================================
# PART 2: CENTRALITY AND DISPERSION MEASURES
# (from 2_EDA_for_FDA2026.R lines 96-142)
# ==============================================================================

fdataobj <- fdata(t(Y_eval), argvals = t_fine)

pdf(file.path(out_dir, "08_centrality_dispersion.pdf"), width = 14, height = 7)
opar <- par(mfrow = c(1, 2))

plot(func.mean(fdataobj), ylim = range(Y_eval),
main = "Centrality Measures",
xlab = "Time (s)", ylab = "Amplitude (uV)")
legend("topright", cex = 0.7, box.col = "white", lty = 1:5,
col = 1:5,
legend = c("mean", "trim.mode", "trim.RP",
"median.mode", "median.RP"))
lines(func.trim.mode(fdataobj, trim = 0.15), col = 2, lty = 2)
lines(func.trim.RP(fdataobj, trim = 0.15),   col = 3, lty = 3)
lines(func.med.mode(fdataobj, trim = 0.15),  col = 4, lty = 4)
lines(func.med.RP(fdataobj, trim = 0.15),    col = 5, lty = 5)

plot(func.var(fdataobj),
main = "Dispersion Measures",
xlab = "Time (s)", ylab = "Variance")
legend("topright", cex = 0.7, box.col = "white", lty = 1:3, col = 1:3,
legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(fdataobj, trim = 0.15), col = 2, lty = 2)
lines(func.trimvar.RP(fdataobj, trim = 0.15),   col = 3, lty = 3)

par(opar)
dev.off()


# ==============================================================================
# PART 3: FUNCTIONAL DEPTH
# (from 4_Boxplots_and_outliers2026.R lines 17-22)
# ==============================================================================

pdf(file.path(out_dir, "09_depth_FM.pdf"), width = 10, height = 8)
out.FM <- depth.FM(fdataobj, trim = 0.1, draw = TRUE)
dev.off()

pdf(file.path(out_dir, "10_depth_mode.pdf"), width = 10, height = 8)
out.mode <- depth.mode(fdataobj, trim = 0.1, draw = TRUE)
dev.off()

pdf(file.path(out_dir, "11_depth_RP.pdf"), width = 10, height = 8)
out.RP <- depth.RP(fdataobj, trim = 0.1, draw = TRUE)
dev.off()


# ==============================================================================
# PART 4: PCA
# (from 3_PCA_for_FDA2026.R lines 16-21)
# ==============================================================================

nharm   <- 4
pcalist <- pca.fd(fd_smooth, nharm, centerfns = TRUE)

cat("\nVariance proportions:\n")
print(pcalist$varprop)
cat("Cumulative:", cumsum(pcalist$varprop), "\n")

pdf(file.path(out_dir, "12_pca_plot.pdf"), width = 10, height = 8)
plot(pcalist)
dev.off()

pdf(file.path(out_dir, "13_pca_harmonics.pdf"), width = 10, height = 6)
plot(pcalist$harmonics)
dev.off()

pdf(file.path(out_dir, "14_pca_scores.pdf"), width = 8, height = 6)
#plotscores(pcalist, loc = 5)
dev.off()

# --- Perturbation plots ---
# (from 3_PCA_for_FDA2026.R lines 24-39)

c <- 2
mn <- pcalist$meanfd

pdf(file.path(out_dir, "15_pca_perturbation.pdf"), width = 12, height = 8)
opar <- par(mfrow = c(2, 2))
for (k in 1:4) {
  phi    <- pcalist$harmonics[k]
  lambda <- pcalist$values[k]
  f1 <- mn - c * sqrt(lambda) * phi
  f2 <- mn + c * sqrt(lambda) * phi

  plot(mn, ylim = range(eval.fd(t_fine, f1), eval.fd(t_fine, f2)),
       lwd = 2, xlab = "Time (s)", ylab = "Amplitude (uV)",
       main = paste0("PC", k, " (", round(pcalist$varprop[k] * 100, 1), "%)"))
  lines(f1, col = 2)
  lines(f2, col = 3)
  abline(v = 0, lty = 2, col = "grey50")
}
par(opar)
dev.off()

# --- PCA reconstruction of first 5 subjects ---
# (from 3_PCA_for_FDA2026.R lines 42-80)

pdf(file.path(out_dir, "16_pca_reconstruction.pdf"), width = 12, height = 10)
opar <- par(mfrow = c(2, 2))
for (i in 1:5) {
  fd.pca1 <- mean.fd(fd_smooth) +
    pcalist$scores[i, 1] * pcalist$harmonics[1]

  fd.pca2 <- mean.fd(fd_smooth) +
    pcalist$scores[i, 1] * pcalist$harmonics[1] +
    pcalist$scores[i, 2] * pcalist$harmonics[2]

  fd.pca3 <- mean.fd(fd_smooth) +
    pcalist$scores[i, 1] * pcalist$harmonics[1] +
    pcalist$scores[i, 2] * pcalist$harmonics[2] +
    pcalist$scores[i, 3] * pcalist$harmonics[3]

  fd.pca4 <- mean.fd(fd_smooth) +
    pcalist$scores[i, 1] * pcalist$harmonics[1] +
    pcalist$scores[i, 2] * pcalist$harmonics[2] +
    pcalist$scores[i, 3] * pcalist$harmonics[3] +
    pcalist$scores[i, 4] * pcalist$harmonics[4]

  yrng <- range(eval.fd(t_fine, fd_smooth[i]))

  plot(fd.pca1, ylim = yrng, ylab = "1 PC",
       main = paste0("Subject ", i))
  lines(fd_smooth[i], col = 2)

  plot(fd.pca2, ylim = yrng, ylab = "2 PCs")
  lines(fd_smooth[i], col = 2)

  plot(fd.pca3, ylim = yrng, ylab = "3 PCs")
  lines(fd_smooth[i], col = 2)

  plot(fd.pca4, ylim = yrng, ylab = "4 PCs")
  lines(fd_smooth[i], col = 2)
}
par(opar)
dev.off()


# ==============================================================================
# PART 5: VARIMAX ROTATION
# (from 3_PCA_for_FDA2026.R lines 86-91)
# ==============================================================================

varmx <- varmx.pca.fd(pcalist)

cat("\nVARIMAX variance proportions:\n")
print(varmx$varprop)

pdf(file.path(out_dir, "17_varimax_plot.pdf"), width = 10, height = 8)
plot(varmx)
dev.off()

pdf(file.path(out_dir, "18_varimax_harmonics.pdf"), width = 10, height = 6)
plot(varmx$harmonics)
dev.off()

pdf(file.path(out_dir, "19_varimax_scores.pdf"), width = 8, height = 6)
#plotscores(varmx, loc = 5)
dev.off()


# ==============================================================================
# PART 6: FUNCTIONAL BOXPLOTS
# (from 4_Boxplots_and_outliers2026.R lines 44, 99-103)
# ==============================================================================

pdf(file.path(out_dir, "20_boxplot_fd.pdf"), width = 10, height = 6)
boxplot(fd_smooth)
dev.off()

pdf(file.path(out_dir, "21_fbplot_MBD.pdf"), width = 10, height = 6)
fbplot(Y_eval, method = "MBD",
       xlab = "Time (s)", ylab = "Amplitude (uV)",
       main = "Functional Boxplot (MBD)")
dev.off()

pdf(file.path(out_dir, "22_fbplot_BD2.pdf"), width = 10, height = 6)
fbplot(Y_eval, method = "BD2",
       xlab = "Time (s)", ylab = "Amplitude (uV)",
       main = "Functional Boxplot (BD2)")
dev.off()


# ==============================================================================
# PART 7: BAND DEPTH AND MODIFIED BAND DEPTH
# (from 4_Boxplots_and_outliers2026.R lines 60-88)
# ==============================================================================

bd <- band_depth(dt = t(Y_eval))
names(bd) <- colnames(fd_smooth$coefs)
cat("\nBand Depth:\n")
print(bd)

mbd <- modified_band_depth(t(Y_eval))
names(mbd) <- colnames(fd_smooth$coefs)
cat("\nModified Band Depth:\n")
print(mbd)

pdf(file.path(out_dir, "23_band_depths.pdf"), width = 12, height = 6)
opar <- par(mfrow = c(1, 2))
plot(bd,  type = "l", main = "Band Depth",         xlab = "Subject", ylab = "BD")
plot(mbd, type = "l", main = "Modified Band Depth", xlab = "Subject", ylab = "MBD")
par(opar)
dev.off()

pdf(file.path(out_dir, "24_fdaoutlier_fbplot_bd.pdf"), width = 10, height = 6)
fbplot_bd <- functional_boxplot(t(Y_eval), depth_method = "bd")
dev.off()
cat("\nOutliers (BD):", fbplot_bd$outliers, "\n")

pdf(file.path(out_dir, "25_fdaoutlier_fbplot_mbd.pdf"), width = 10, height = 6)
fbplot_mbd <- functional_boxplot(t(Y_eval), depth_method = "mbd")
dev.off()
cat("Outliers (MBD):", fbplot_mbd$outliers, "\n")


# ==============================================================================
# PART 8: MUOD OUTLIER DETECTION
# (from 4_Boxplots_and_outliers2026.R lines 93-94)
# ==============================================================================

m <- muod(t(Y_eval), cut_method = "boxplot")
names(m)
str(m)
m <- muod(t(Y_eval), cut_method = "boxplot")
cat("\nMUOD outliers:\n")
print(m$outliers)

pdf(file.path(out_dir, "26_muod.pdf"), width = 12, height = 4)
opar <- par(mfrow = c(1, 3))
plot(m$indices$shape,     type = "h", main = "Shape Index",     xlab = "Subject", ylab = "IS")
plot(m$indices$magnitude, type = "h", main = "Magnitude Index", xlab = "Subject", ylab = "IM")
plot(m$indices$amplitude, type = "h", main = "Amplitude Index", xlab = "Subject", ylab = "IA")
par(opar)
dev.off()


# ==============================================================================
# PART 9: RAINBOW PLOTS
# (from 5_Rainbow2026.R lines 118-136)
# ==============================================================================

fds_obj <- fds(x = t_fine, y = Y_eval,
               xname = "Time (s)", yname = "Amplitude (uV)")

pdf(file.path(out_dir, "27_rainbow_functions.pdf"), width = 10, height = 6)
plot(fds_obj, plot.type = "functions", plotlegend = TRUE)
dev.off()

pdf(file.path(out_dir, "28_rainbow_depth.pdf"), width = 10, height = 6)
plot(fds_obj, plot.type = "depth", plotlegend = TRUE)
dev.off()


# ==============================================================================
# PART 10: FUNCTIONAL BAGPLOT
# (from 5_Rainbow2026.R lines 138-142)
# ==============================================================================

pdf(file.path(out_dir, "29_bagplot_bivariate.pdf"), width = 8, height = 8)
fboxplot(fds_obj, plot.type = "bivariate",
         type = "bag", projmethod = "PCAproj")
dev.off()

pdf(file.path(out_dir, "30_bagplot_functional.pdf"), width = 10, height = 6)
fboxplot(fds_obj, plot.type = "functional",
         type = "bag", projmethod = "PCAproj")
dev.off()


# ==============================================================================
# PART 11: HDR BOXPLOT
# (from 5_Rainbow2026.R lines 63-87)
# ==============================================================================

pdf(file.path(out_dir, "31_hdr_bivariate_007.pdf"), width = 8, height = 8)
fboxplot(fds_obj, plot.type = "bivariate",
         type = "hdr", alpha = c(0.07, 0.5),
         projmethod = "PCAproj")
dev.off()

pdf(file.path(out_dir, "32_hdr_bivariate_005.pdf"), width = 8, height = 8)
fboxplot(fds_obj, plot.type = "bivariate",
         type = "hdr", alpha = c(0.05, 0.5),
         projmethod = "PCAproj")
dev.off()

pdf(file.path(out_dir, "33_hdr_functional_007.pdf"), width = 10, height = 6)
fboxplot(fds_obj, plot.type = "functional",
         type = "hdr", alpha = c(0.07, 0.5),
         projmethod = "PCAproj")
dev.off()

pdf(file.path(out_dir, "34_hdr_functional_005.pdf"), width = 10, height = 6)
fboxplot(fds_obj, plot.type = "functional",
         type = "hdr", alpha = c(0.05, 0.5),
         projmethod = "PCAproj")
dev.off()


# ==============================================================================
# PART 12: FUNCTIONAL OUTLIER DETECTION
# (from 5_Rainbow2026.R lines 162-166)
# ==============================================================================

sink(file.path(out_dir, "35_foutliers_results.txt"))
cat("=== robMah ===\n")
print(foutliers(fds_obj, method = "robMah"))
cat("\n=== lrt ===\n")
print(foutliers(fds_obj, method = "lrt"))
cat("\n=== depth.trim ===\n")
print(foutliers(fds_obj, method = "depth.trim"))
cat("\n=== depth.pond ===\n")
print(foutliers(fds_obj, method = "depth.pond"))
cat("\n=== HUoutliers ===\n")
print(foutliers(fds_obj, method = "HUoutliers"))
sink()

cat("\n*** All EDA outputs saved (PDFs 04-34 + text file 35). Done! ***\n")






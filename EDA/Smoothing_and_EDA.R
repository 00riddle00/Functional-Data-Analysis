# vim: set ft=r tw=88 nu ai et ts=2 sw=2:
# ==============================================================================
#
#  Smoothing and Exploratory Data Analysis
#  Flanker Task EEG — Channel FC1 — Stimulus S2
#
#  Input:  Flanker_stimulus_FC1_channel.csv (from assemble_subject_data.R)
#  Output: PDF plots + fd_smooth.rds + text results in outputs/
#
# ==============================================================================

# Create output directory if it doesn't exist
out_dir <- file.path("./EDA/outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Clean previous outputs to avoid corrupt PDFs from incomplete overwrites
old_files <- list.files(out_dir, pattern = "\\.(pdf|rds|txt)$", full.names = TRUE)
if (length(old_files) > 0) {
  file.remove(old_files)
  cat("Cleaned", length(old_files), "old output files.\n")
}

# ==============================================================================
# Step 0: Plot raw curves before smoothing
# ==============================================================================

library(fda)
set.seed(42)

dat      <- read.csv("./EDA/Flanker_stimulus_FC1_channel.csv", check.names = FALSE)
time_vec <- dat$time
Y_mat    <- as.matrix(dat[, -1])
n_subj   <- ncol(Y_mat)

cairo_pdf(file.path(out_dir, "00_all_curves.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
matplot(time_vec, Y_mat, type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.4),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = paste0("Flanker Task — Channel FC1 — All ", n_subj, " Subjects"))
abline(h = 0, lty = 2, col = "grey30", lwd = 1.5)
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Individual subjects", "Stimulus onset"),
       col = c("steelblue", "red"),
       lty = c(1, 2), lwd = c(2, 2), cex = 0.8)
dev.off()
cat("Saved: 00_all_curves.pdf\n")

################################################################################
# Step 1: Smoothing
################################################################################

library(fda)

out_dir <- file.path("./EDA/outputs")

# ==============================================================================
# Load data
# ==============================================================================

dat      <- read.csv("./EDA/Flanker_stimulus_FC1_channel.csv", check.names = FALSE)
time_vec <- dat$time
Y_mat    <- as.matrix(dat[, -1])
n_subj   <- ncol(Y_mat)

cat("Time points:", length(time_vec), "\n")
cat("Time range: ", range(time_vec), "\n")
cat("Subjects:   ", n_subj, "\n")

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
cairo_pdf(file.path(out_dir, "01_gcv_lambda.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
plot(loglam, gcvsave, type = "b", lwd = 2,
     ylab = "GCV Criterion",
     xlab = expression(log[10](lambda)),
     main = "GCV for Smoothing Parameter Selection")
abline(v = loglam[best_idx],
       lty = 2, col = "forestgreen", lwd = 2.5)
y_pos <- min(gcvsave) + 0.6 * (max(gcvsave) - min(gcvsave))
text(loglam[best_idx], y_pos,
     labels = paste0("optimal: ", loglam[best_idx]),
     col = "forestgreen", pos = 4)
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

cairo_pdf(file.path(out_dir, "02_before_after_smoothing.pdf"), width = 12, height = 6,
          bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(1, 2))

matplot(time_vec, Y_mat, type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.4),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "Before Smoothing (raw)")
abline(h = 0, lty = 2, col = "grey30", lwd = 1.5)
abline(v = 0, lty = 2, col = "red")

plot(fd_smooth, col = adjustcolor("steelblue", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = bquote(bold("After Smoothing (" * lambda == 10^.(loglam[best_idx]) * ")")))
abline(v = 0, lty = 2, col = "red")

par(opar)
dev.off()
cat("Saved: 02_before_after_smoothing.pdf\n")

# ==============================================================================
# Plot: one subject close-up (raw vs smoothed)
# ==============================================================================

cairo_pdf(file.path(out_dir, "03_single_subject_fit.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
sub_idx <- 1
y_raw   <- Y_mat[, sub_idx]
y_fit   <- eval.fd(time_vec, fd_smooth[sub_idx])

plot(time_vec, y_raw, type = "l", col = "grey60",
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = paste0("Subject ", colnames(Y_mat)[sub_idx], ": Raw vs Smoothed"))
lines(time_vec, y_fit, col = "steelblue", lwd = 2)
abline(h = 0, lty = 3, col = "grey45", lwd = 1.5)
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Raw", "Smoothed"),
       col = c("grey60", "steelblue"),
       lty = 1, lwd = c(2, 3), cex = 0.8)
dev.off()
cat("Saved: 03_single_subject_fit.pdf\n")

# ==============================================================================
# Save the fd object for later use in EDA
# ==============================================================================

saveRDS(fd_smooth, file.path(out_dir, "fd_smooth.rds"))
cat("\nSmoothed fd object saved to:", file.path(out_dir, "fd_smooth.rds"), "\n")
cat("Use fd_smooth <- readRDS(file.path(out_dir, 'fd_smooth.rds')) in the EDA script.\n")

################################################################################
# Step 2: Functional EDA
################################################################################

library(fda)
library(fda.usc)
library(fdaoutlier)
library(rainbow)
library(fields)

# Load the smoothed fd object from Step 1
fd_smooth <- readRDS(file.path(out_dir, "fd_smooth.rds"))
cat("\nLoaded smoothed fd object from:", file.path(out_dir, "fd_smooth.rds"), "\n")

# Evaluate on a fine grid (for fdaoutlier and rainbow functions)
t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 501)
Y_eval <- eval.fd(t_fine, fd_smooth)
n_subj <- ncol(Y_eval)

# ==============================================================================
# PART 1: MEAN, SD, COVARIANCE
# ==============================================================================

meanfd <- mean.fd(fd_smooth)
sdfd   <- std.fd(fd_smooth)

cairo_pdf(file.path(out_dir, "04_mean_sd.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
plot(fd_smooth, col = "gray", lty = 1,
     xlab = "Time (s)",  ylab = "Amplitude (µV)",
     main = "Mean and Standard Deviation")
lines(meanfd, lwd = 4, lty = 2, col = 2)
lines(sdfd,   lwd = 4, lty = 2, col = 4)
lines(meanfd - sdfd, lwd = 4, lty = 2, col = 6)
lines(meanfd + sdfd, lwd = 4, lty = 2, col = 6)
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Curves", "Mean", "SD", "Mean ± SD"),
       col = c("gray", 2, 4, 6), lty = c(1, 2, 2, 2),
       lwd = c(2, 4, 4, 4), cex = 0.7)
dev.off()
cat("Saved: 04_mean_sd.pdf\n")

# --- Covariance ---
# (from 2_EDA_for_FDA2026.R lines 57-90)

covbifd <- var.fd(fd_smooth)
t_grid  <- seq(fd_smooth$basis$rangeval[1],
               fd_smooth$basis$rangeval[2], length.out = 101)
cov_mat <- eval.bifd(t_grid, t_grid, covbifd)

cairo_pdf(file.path(out_dir, "05_covariance_persp.pdf"), width = 10, height = 8,
          bg = "transparent")
par(bg = NA)
persp(t_grid, t_grid, cov_mat,
      theta = -45, phi = 25, r = 3, expand = 0.5,
      ticktype = "detailed",
      xlab = "Time s (s)",
      ylab = "Time t (s)",
      zlab = "Cov(s, t)  (μV²)",
      main = "Covariance Surface")
dev.off()
cat("Saved: 05_covariance_persp.pdf\n")

cairo_pdf(file.path(out_dir, "06_covariance_contour.pdf"), width = 10, height = 8,
          bg = "transparent")
par(bg = NA)
contour(t_grid, t_grid, cov_mat,
        xlab = "Time s (s)", ylab = "Time t (s)",
        main = "Covariance Contour", lwd = 2)
dev.off()
cat("Saved: 06_covariance_contour.pdf\n")

cairo_pdf(file.path(out_dir, "07_covariance_image.pdf"), width = 10, height = 8,
          bg = "transparent")
par(bg = NA)
image.plot(t_grid, t_grid, cov_mat,
           xlab = "Time s (s)", ylab = "Time t (s)",
           main = "Covariance Surface")
contour(t_grid, t_grid, cov_mat, col = "white", add = TRUE)
dev.off()
cat("Saved: 07_covariance_image.pdf\n")

# ==============================================================================
# PART 2: CENTRALITY AND DISPERSION MEASURES
# ==============================================================================

fdataobj <- fdata(t(Y_eval), argvals = t_fine)

cairo_pdf(file.path(out_dir, "08_centrality_dispersion.pdf"),
          width = 14, height = 7, bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(1, 2), bg = NA)

yrng <- quantile(Y_eval, probs = c(0.01, 0.99))

plot(func.mean(fdataobj), ylim = yrng,
     main = "Centrality Measures",
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     col = "black", lwd = 2)
legend("topleft", cex = 0.7,
       lty = 1:5,
       col = c("black", "#FF8C00", "forestgreen", "steelblue3", "cyan3"),
       lwd = c(2, 2, 2, 2, 2),
       legend = c("mean", "trim.mode", "trim.RP",
                  "median.mode", "median.RP"))
lines(func.trim.mode(fdataobj, trim = 0.15),
      col = "#FF8C00", lty = 2, lwd = 2)
lines(func.trim.RP(fdataobj, trim = 0.15),
      col = "forestgreen", lty = 3, lwd = 2)
lines(func.med.mode(fdataobj, trim = 0.15),
      col = "steelblue3", lty = 4, lwd = 2)
lines(func.med.RP(fdataobj, trim = 0.15),
      col = "cyan3", lty = 5, lwd = 1.2)
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")

plot(func.var(fdataobj),
     main = "Dispersion Measures",
     xlab = "Time (s)", ylab = "Variance (µV²)",
     col = "black", lwd = 2)
legend("topleft", cex = 0.7,
       lty = 1:3,
       col = c("black", "#FF8C00", "forestgreen"),
       lwd = c(2, 2, 2),
       legend = c("var", "trimvar.mode", "trimvar.RP"))
lines(func.trimvar.mode(fdataobj, trim = 0.15),
      col = "#FF8C00", lty = 2, lwd = 2)
lines(func.trimvar.RP(fdataobj, trim = 0.15),
      col = "forestgreen", lty = 3, lwd = 2)
abline(v = 0, lty = 2, col = "red")

par(opar)
dev.off()
cat("Saved: 08_centrality_dispersion.pdf\n")

# ==============================================================================
# PART 3: FUNCTIONAL DEPTH
# ==============================================================================

out.FM <- depth.FM(fdataobj, trim = 0.1, draw = FALSE)
cairo_pdf(file.path(out_dir, "09_depth_FM.pdf"),
          width = 10, height = 8, bg = "transparent")
par(bg = NA)
plot(fdataobj,
     xlab = "Time (s)",
     ylab = "Amplitude (µV)",
     main = "Functional Depth: Fraiman–Muniz (trim = 10%)",
     col = "grey60")
lines(fdataobj[out.FM$ltrim],
      col = "steelblue3", lwd = 1.3)
lines(out.FM$mtrim,
      col = "#FF8C00", lwd = 3, lty = 2)
lines(out.FM$median,
      col = "firebrick2", lwd = 2)
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("All curves", "Central curves",
                  "Trimmed mean (10%)", "Median"),
       col = c("grey60", "steelblue3",
               "orange", "firebrick2"),
       lty = c(1, 1, 2, 1),
       lwd = c(1, 1.3, 3, 2),
       cex = 0.8)
dev.off()
cat("Saved: 09_depth_FM.pdf\n")

out.mode <- depth.mode(fdataobj, trim = 0.1, draw = FALSE)
cairo_pdf(file.path(out_dir, "10_depth_mode.pdf"),
          width = 10, height = 8, bg = "transparent")
par(bg = NA)
plot(fdataobj,
     xlab = "Time (s)",
     ylab = "Amplitude (µV)",
     main = "Functional Depth: Mode (trim = 10%)",
     col = "grey60")
lines(fdataobj[out.mode$ltrim],
      col = "steelblue3", lwd = 1.3)
lines(out.mode$mtrim,
      col = "#FF8C00", lwd = 3, lty = 2)
lines(out.mode$median,
      col = "firebrick2", lwd = 2)
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("All curves", "Central curves",
                  "Trimmed mean (10%)", "Median"),
       col = c("grey60", "steelblue3",
               "#FF8C00", "firebrick2"),
       lty = c(1, 1, 2, 1),
       lwd = c(1, 1.3, 3, 2),
       cex = 0.8)
dev.off()
cat("Saved: 10_depth_mode.pdf\n")

out.RP <- depth.RP(fdataobj, trim = 0.1, draw = FALSE)
cairo_pdf(file.path(out_dir, "11_depth_RP.pdf"),
          width = 10, height = 8, bg = "transparent")
par(bg = NA)
plot(fdataobj,
     xlab = "Time (s)",
     ylab = "Amplitude (µV)",
     main = "Functional Depth: Random Projection (trim = 10%)",
     col = "grey60")
lines(fdataobj[out.RP$ltrim],
      col = "steelblue3", lwd = 1.3)
lines(out.RP$mtrim,
      col = "#FF8C00", lwd = 3, lty = 2)
lines(out.RP$median,
      col = "firebrick2", lwd = 2)
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("All curves", "Central curves",
                  "Trimmed mean (10%)", "Median"),
       col = c("grey60", "steelblue3",
               "#FF8C00", "firebrick2"),
       lty = c(1, 1, 2, 1),
       lwd = c(1, 1.3, 3, 2),
       cex = 0.8)
dev.off()
cat("Saved: 11_depth_RP.pdf\n")

# ==============================================================================
# PART 4: PCA
# ==============================================================================

nharm   <- 4
pcalist <- pca.fd(fd_smooth, nharm, centerfns = TRUE)

cat("\nVariance proportions:\n")
print(pcalist$varprop)
cat("Cumulative:", cumsum(pcalist$varprop), "\n")

harmonics_mat <- eval.fd(time_vec, pcalist$harmonics)
mean_vec      <- as.vector(eval.fd(time_vec, pcalist$meanfd))
cairo_pdf(file.path(out_dir, "12_pca_plot.pdf"),
          width = 10, height = 8, bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(2, 2), bg = NA)
for (i in 1:nharm) {
  fac    <- sqrt(pcalist$values[i])
  pc_plus  <- mean_vec + fac * harmonics_mat[, i]
  pc_minus <- mean_vec - fac * harmonics_mat[, i]
  ylim_i <- range(c(mean_vec, pc_plus, pc_minus))

  plot(time_vec, mean_vec, type = "l",
       col = "black", lwd = 2,
       xlab = "Time (s)",
       ylab = "Amplitude (µV)",
       ylim = ylim_i,
       main = sprintf("PCA fn. %d (Harmonic %d, Expl. variance: %.1f%%)",
                      i, i, 100 * pcalist$varprop[i]))
  lines(time_vec, pc_plus,  col = "steelblue3", lty = 2, lwd = 2)
  lines(time_vec, pc_minus, col = "#FF8C00",   lty = 3, lwd = 2)
  abline(h = 0, lty = 2, col = "grey30")
  abline(v = 0, col = "red", lty = 2)
  legend("topleft",
         legend = c("Mean", "+ variation", "− variation"),
         col = c("black", "steelblue3", "#FF8C00"),
         lty = c(1, 2, 3), lwd = c(1.5, 1.5, 1.5), cex = 0.8,
         bg = adjustcolor("white", alpha.f = 0.9))
}
par(opar)
dev.off()
cat("Saved: 12_pca_plot.pdf\n")

cairo_pdf(file.path(out_dir, "13_pca_harmonics.pdf"),
          width = 10, height = 6, bg = "transparent")
par(bg = NA)
matplot(time_vec, cbind(mean_vec, harmonics_mat),
        type = "l", lty = c(1, 2, 3, 4, 5),
        col = c("black", "steelblue3", "#FF8C00", "forestgreen", "purple3"),
        lwd = c(1.5, 2, 2, 2, 2),
        xlab = "Time (s)",
        ylab = "Eigenfunction value",
        main = "PCA Harmonics and Mean Function")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Mean", "PC1", "PC2", "PC3", "PC4"),
       col = c("black", "steelblue3", "#FF8C00", "forestgreen", "purple3"),
       lty = c(1, 2, 3, 4, 5),
       lwd = c(1.5, 1.5, 1.5, 1.5, 1.5),
       cex = 0.8)
dev.off()
cat("Saved: 13_pca_harmonics.pdf\n")

cairo_pdf(file.path(out_dir, "14_pca_scores.pdf"), width = 8, height = 6,
          bg = "transparent")
par(bg = NA)
plot(pcalist$scores[, 1], pcalist$scores[, 2],
     pch = 19, col = adjustcolor("steelblue", 0.6), cex = 1.2,
     xlab = paste0("PC1 (", round(pcalist$varprop[1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(pcalist$varprop[2] * 100, 1), "%)"),
     main = "fPCA Scores: PC1 vs PC2")
abline(h = 0, v = 0, lty = 3, col = "grey50", lwd = 2)
dev.off()
cat("Saved: 14_pca_scores.pdf\n")

# --- Perturbation plots ---
# (from 3_PCA_for_FDA2026.R lines 24-39)

c <- 2
mn <- pcalist$meanfd

cairo_pdf(file.path(out_dir, "15_pca_perturbation.pdf"), width = 12, height = 8,
          bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(2, 2))
for (k in 1:4) {
  phi    <- pcalist$harmonics[k]
  lambda <- pcalist$values[k]
  f1 <- mn - c * sqrt(lambda) * phi
  f2 <- mn + c * sqrt(lambda) * phi

  plot(mn, ylim = range(eval.fd(t_fine, f1), eval.fd(t_fine, f2)),
       lwd = 2, xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = paste0("PC", k, " Perturbation (Expl. variance: ",
                     round(pcalist$varprop[k] * 100, 1), "%)"))
  lines(f1, col = 2)
  lines(f2, col = 3)
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Mean", "+ variation", "− variation"),
         col = c("black", "green3", "red"),
         lty = c(1, 1, 1),
         lwd = c(1.5, 1.5, 1.5),
         cex = 0.8)
}
par(opar)
dev.off()
cat("Saved: 15_pca_perturbation.pdf\n")

# --- PCA reconstruction of first 5 subjects ---
# (from 3_PCA_for_FDA2026.R lines 42-80)

cairo_pdf(file.path(out_dir, "16_pca_reconstruction.pdf"), width = 12, height = 10,
          bg = "transparent")
par(bg = NA)
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

  orig  <- eval.fd(t_fine, fd_smooth[i])
  recon <- c(
    eval.fd(t_fine, fd.pca1),
    eval.fd(t_fine, fd.pca2),
    eval.fd(t_fine, fd.pca3),
    eval.fd(t_fine, fd.pca4)
  )

  yrng <- range(c(orig, recon))

  plot(fd.pca1, lwd = 1.5, ylim = yrng, xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = paste0("Subject ", colnames(fd_smooth$coefs)[i], " — 1 PC"))
  lines(fd_smooth[i], col = 2)
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Reconstruction", "Original"),
         col = c("black", "red"),
         lty = 1,
         lwd = c(2, 1.5),
         bg = adjustcolor("white", alpha.f = 0.9)
  )

  plot(fd.pca2, lwd = 1.5, ylim = yrng, xlab = "Time (s)", ylab = "Amplitude (µV)",
       # main = "2 PCs")
       main = paste0("Subject ", colnames(fd_smooth$coefs)[i], " — 2 PCs"))
  lines(fd_smooth[i], col = 2)
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Reconstruction", "Original"),
         col = c("black", "red"),
         lty = 1,
         lwd = c(2, 1.5),
         bg = adjustcolor("white", alpha.f = 0.9)
  )

  plot(fd.pca3, lwd = 1.5, ylim = yrng, xlab = "Time (s)", ylab = "Amplitude (µV)",
       # main = "3 PCs")
       main = paste0("Subject ", colnames(fd_smooth$coefs)[i], " — 3 PCs"))
  lines(fd_smooth[i], col = 2)
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Reconstruction", "Original"),
         col = c("black", "red"),
         lty = 1,
         lwd = c(2, 1.5),
         bg = adjustcolor("white", alpha.f = 0.9)
  )

  plot(fd.pca4, lwd = 1.5, ylim = yrng, xlab = "Time (s)", ylab = "Amplitude (µV)",
       # main = "4 PCs")
       main = paste0("Subject ", colnames(fd_smooth$coefs)[i], " — 4 PCs"))
  lines(fd_smooth[i], col = 2)
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Reconstruction", "Original"),
         col = c("black", "red"),
         lty = 1,
         lwd = c(2, 1.5),
         bg = adjustcolor("white", alpha.f = 0.9)
  )
}
par(opar)
dev.off()
cat("Saved: 16_pca_reconstruction.pdf\n")

# ==============================================================================
# PART 5: VARIMAX ROTATION
# ==============================================================================

varmx <- varmx.pca.fd(pcalist)

cat("\nVARIMAX variance proportions:\n")
print(varmx$varprop)

# Custom VARIMAX plot with correct "Rotated PC" labels
cairo_pdf(file.path(out_dir, "17_varimax_plot.pdf"), width = 12, height = 8,
          bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(2, 2))
for (k in 1:4) {
  phi    <- varmx$harmonics[k]
  lambda <- varmx$values[k]
  f1 <- mn - c * sqrt(lambda) * phi
  f2 <- mn + c * sqrt(lambda) * phi

  plot(mn, ylim = range(eval.fd(t_fine, f1), eval.fd(t_fine, f2)),
       lwd = 2, xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = paste0("Rotated PC", k, " (Expl. variance: ",
                      round(varmx$varprop[k] * 100, 1), "%)"))
  lines(f1, col = 2)
  lines(f2, col = 3)
  abline(h = 0, lty = 2, col = "grey30")
  abline(v = 0, lty = 2, col = "red")
  legend("topleft",
         legend = c("Mean", "+ variation", "− variation"),
         col = c("black", "green3", "red"),
         lty = c(1, 1, 1),
         lwd = c(1.5, 1.5, 1.5),
         cex = 0.8,
         bg = adjustcolor("white", alpha.f = 0.9))
}
par(opar)
dev.off()
cat("Saved: 17_varimax_plot.pdf\n")

cairo_pdf(file.path(out_dir, "18_varimax_harmonics.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
varmx_harm_mat <- eval.fd(time_vec, varmx$harmonics)
matplot(time_vec, varmx_harm_mat,
        type = "l",
        lty = c(1, 2, 3, 4),
        col = c("black", "red", "green3", "dodgerblue2"),
        lwd = 1.5,
        xlab = "Time (s)",
        ylab = "Eigenfunction value",
        main = "VARIMAX Rotated Harmonics")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Rotated PC1", "Rotated PC2", "Rotated PC3", "Rotated PC4"),
       col = c("black", "red", "green3", "dodgerblue2"),
       lty = c(1, 2, 3, 4),
       lwd = c(1.5, 1.5, 1.5, 1.5),
       cex = 0.8)
dev.off()
cat("Saved: 18_varimax_harmonics.pdf\n")

cairo_pdf(file.path(out_dir, "19_varimax_scores.pdf"), width = 8, height = 6,
          bg = "transparent")
par(bg = NA)
plot(varmx$scores[, 1], varmx$scores[, 2],
     pch = 19, col = adjustcolor("steelblue", 0.6), cex = 1.2,
     xlab = paste0("Rot. PC1 (", round(varmx$varprop[1] * 100, 1), "%)"),
     ylab = paste0("Rot. PC2 (", round(varmx$varprop[2] * 100, 1), "%)"),
     main = "VARIMAX Rotated Scores")
abline(h = 0, v = 0, lty = 3, col = "grey50", lwd = 2)
dev.off()
cat("Saved: 19_varimax_scores.pdf\n")

# ==============================================================================
# PART 6: FUNCTIONAL BOXPLOTS
# ==============================================================================

cairo_pdf(file.path(out_dir, "20_boxplot_fd.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
boxplot(fd_smooth,
        xlab = "Time (s)",
        ylab = "Amplitude (µV)",
        main = "Functional Boxplot")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
legend("topleft",
       legend = c("Median", "Central region", "Envelope"),
       col = c("black", "magenta", "dodgerblue"),
       lty = c(1, NA, 1),
       lwd = c(2.5, NA, 2.5),
       pch = c(NA, 15, NA),
       pt.cex = 2,
       cex = 0.8
)
dev.off()
cat("Saved: 20_boxplot_fd.pdf\n")

# fbplot: manually set x-axis to time values
# Note: fbplot(Y_eval) uses column indices for x by default.
# We suppress the default x-axis and draw our own with time labels.

n_t <- nrow(Y_eval)
tick_pos <- seq(1, n_t, length.out = 6)
tick_lab <- round(seq(t_fine[1], t_fine[n_t], length.out = 6), 2)
v0 <- which.min(abs(t_fine - 0))

cairo_pdf(file.path(out_dir, "21_fbplot_MBD.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
fbplot(Y_eval, method = "MBD", xaxt = "n",
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = "Functional Boxplot (MBD)")
axis(1, at = tick_pos, labels = tick_lab)
abline(h = 0, lty = 2, col = "grey30")
abline(v = v0, lty = 2, col = "grey30")
legend("topleft",
       legend = c("Median", "Central region", "Envelope", "Outlier cutoff"),
       col = c("black", "magenta", "dodgerblue", "red"),
       lty = c(1, NA, 1, 2),
       lwd = c(2.5, NA, 2.5, 2.5),
       pch = c(NA, 15, NA, NA),
       pt.cex = 2,
       cex = 0.8
)
dev.off()
cat("Saved: 21_fbplot_MBD.pdf\n")

cairo_pdf(file.path(out_dir, "22_fbplot_BD2.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
fbplot(Y_eval, method = "BD2", xaxt = "n",
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = "Functional Boxplot (BD2)")
axis(1, at = tick_pos, labels = tick_lab)
abline(h = 0, lty = 2, col = "grey30")
abline(v = v0, lty = 2, col = "red")
legend("topleft",
       legend = c("Median", "Central region", "Envelope"),
       col = c("black", "magenta", "dodgerblue"),
       lty = c(1, NA, 1),
       lwd = c(2.5, NA, 2.5),
       pch = c(NA, 15, NA),
       pt.cex = 2,
       cex = 0.8
)
dev.off()
cat("Saved: 22_fbplot_BD2.pdf\n")

# ==============================================================================
# PART 7: BAND DEPTH AND MODIFIED BAND DEPTH
# ==============================================================================

bd <- band_depth(dt = t(Y_eval))
names(bd) <- colnames(fd_smooth$coefs)
cat("\nBand Depth:\n")
print(bd)

mbd <- modified_band_depth(t(Y_eval))
names(mbd) <- colnames(fd_smooth$coefs)
cat("\nModified Band Depth:\n")
print(mbd)

cairo_pdf(file.path(out_dir, "23_band_depths.pdf"), width = 12, height = 6,
          bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(1, 2))
plot(bd,  type = "l", lwd = 1.5, main = "Band Depth (BD)",
     xlab = "Subject index", ylab = "BD")
plot(mbd, type = "l", lwd = 1.5, main = "Modified Band Depth (MBD)",
     xlab = "Subject index", ylab = "MBD")
par(opar)
dev.off()
cat("Saved: 23_band_depths.pdf\n")

# ==============================================================================
# PART 8: MUOD OUTLIER DETECTION
# ==============================================================================

m <- muod(t(Y_eval), cut_method = "boxplot")
cat("\nMUOD outliers:\n")
print(m$outliers)

cairo_pdf(file.path(out_dir, "24_muod.pdf"), width = 12, height = 4,
          bg = "transparent")
par(bg = NA)
opar <- par(mfrow = c(1, 3))
plot(m$indices$shape, type = "h", lwd = 1.5, main = "Shape Index (IS)",
     xlab = "Subject index", ylab = "IS")
plot(m$indices$magnitude, type = "h", lwd = 1.5, main = "Magnitude Index (IM)",
     xlab = "Subject index", ylab = "IM")
plot(m$indices$amplitude, type = "h", lwd = 1.5, main = "Amplitude Index (IA)",
     xlab = "Subject index", ylab = "IA")
par(opar)
dev.off()
cat("Saved: 24_muod.pdf\n")

# ==============================================================================
# PART 9: RAINBOW PLOTS
# ==============================================================================

fds_obj <- fds(x = t_fine, y = Y_eval,
               xname = "Time (s)", yname = "Amplitude (µV)")

cairo_pdf(file.path(out_dir, "25_rainbow_functions.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
plot(fds_obj, plot.type = "functions", plotlegend = TRUE,
     main = "Rainbow Plot — Functions")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
dev.off()
cat("Saved: 25_rainbow_functions.pdf\n")

cairo_pdf(file.path(out_dir, "26_rainbow_depth.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
plot(fds_obj, plot.type = "depth", plotlegend = TRUE,
     main = "Rainbow Plot — Depth Ordering")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
dev.off()
cat("Saved: 26_rainbow_depth.pdf\n")

# ==============================================================================
# PART 10: FUNCTIONAL BAGPLOT
# ==============================================================================

cairo_pdf(file.path(out_dir, "27_bagplot_bivariate.pdf"), width = 8, height = 8,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "bivariate", type = "bag",
         projmethod = "PCAproj", main = "Bagplot — PCA Scores")
dev.off()
cat("Saved: 27_bagplot_bivariate.pdf\n")

cairo_pdf(file.path(out_dir, "28_bagplot_functional.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "functional",
         type = "bag", projmethod = "PCAproj")
title(main = "Functional Bagplot — PCA Projection")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
dev.off()
cat("Saved: 28_bagplot_functional.pdf\n")

# ==============================================================================
# PART 11: HDR BOXPLOT
# ==============================================================================

cairo_pdf(file.path(out_dir, "29_hdr_bivariate_007.pdf"), width = 8, height = 8,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "bivariate", type = "hdr",
         alpha = c(0.07, 0.5), projmethod = "PCAproj",
         main = "HDR Boxplot — Bivariate PCA Scores, α = (0.07, 0.50)")
dev.off()
cat("Saved: 29_hdr_bivariate_007.pdf\n")

cairo_pdf(file.path(out_dir, "30_hdr_bivariate_005.pdf"), width = 8, height = 8,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "bivariate", type = "hdr",
         alpha = c(0.05, 0.5), projmethod = "PCAproj",
         main = "HDR Boxplot — Bivariate PCA Scores, α = (0.05, 0.50)")
dev.off()
cat("Saved: 30_hdr_bivariate_005.pdf\n")

cairo_pdf(file.path(out_dir, "31_hdr_functional_007.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "functional",
         type = "hdr", alpha = c(0.07, 0.5),
         projmethod = "PCAproj")
title(main = "HDR Functional Boxplot, α = (0.07, 0.50)")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
dev.off()
cat("Saved: 31_hdr_functional_007.pdf\n")

cairo_pdf(file.path(out_dir, "32_hdr_functional_005.pdf"), width = 10, height = 6,
          bg = "transparent")
par(bg = NA)
fboxplot(fds_obj, plot.type = "functional",
         type = "hdr", alpha = c(0.05, 0.5),
         projmethod = "PCAproj")
title(main = "HDR Functional Boxplot, α = (0.05, 0.50)")
abline(h = 0, lty = 2, col = "grey30")
abline(v = 0, lty = 2, col = "red")
dev.off()
cat("Saved: 32_hdr_functional_005.pdf\n")

# ==============================================================================
# PART 12: FUNCTIONAL OUTLIER DETECTION
# ==============================================================================

sink(file.path(out_dir, "33_foutliers_results.txt"), split = TRUE)
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
cat("Saved: 33_foutliers_results.txt\n")

cat("\n*** All EDA outputs saved. Done! ***\n")

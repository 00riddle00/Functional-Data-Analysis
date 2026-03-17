# ==============================================================================
# 0. Setup & Load Data
# ==============================================================================


install.packages(c("fda", "ggplot2", "reshape2", "viridis", "gridExtra"))

library(fda)
library(ggplot2)
library(reshape2)
library(viridis)
library(gridExtra)

# Load the functional data object
fd_obj <- readRDS("C:/Users/berzi/Downloads/task-flanker_stimulus_s11.rds")


# Just look at it
str(fd_obj)


# ==============================================================================
# 1. Inspect the Data Structure
# ==============================================================================

cat("========== DATA STRUCTURE ==========\n")
cat("Class:", class(fd_obj), "\n")
cat("Basis type:", fd_obj$basis$type, "\n")
cat("Number of basis functions:", fd_obj$basis$nbasis, "\n")
cat("Basis order (degree+1):", 4, "(cubic B-spline)\n")
cat("Time range:", fd_obj$basis$rangeval, "seconds\n")
cat("Coefficient matrix dimensions:", dim(fd_obj$coefs), "\n")
cat("  -> Rows = basis functions:", nrow(fd_obj$coefs), "\n")
cat("  -> Columns = trials/curves:", ncol(fd_obj$coefs), "\n")
cat("fdnames:", fd_obj$fdnames[[1]], ",", fd_obj$fdnames[[2]], ",",
    fd_obj$fdnames[[3]], "\n")
cat("Interior knots:", length(fd_obj$basis$params), "\n")
cat("Knot positions (first 10):", head(fd_obj$basis$params, 10), "\n")
cat("====================================\n\n")

# Define time grid for evaluation
time_range <- fd_obj$basis$rangeval
t_fine     <- seq(time_range[1], time_range[2], length.out = 501)
n_curves   <- ncol(fd_obj$coefs)

# Evaluate the fd object on the fine grid
Y_mat <- eval.fd(t_fine, fd_obj)   # 501 x n_curves matrix

cat("Evaluated data matrix dimensions:", dim(Y_mat), "\n")
cat("  -> 501 time points x", n_curves, "trials\n\n")

# ==============================================================================
# 2. EXPLORATORY DATA ANALYSIS (EDA) - Raw Functional Curves
# ==============================================================================

# --- 2a. Plot ALL individual curves (spaghetti plot) ---

pdf("01_all_curves_spaghetti.pdf", width = 10, height = 6)
matplot(t_fine, Y_mat, type = "l", lty = 1, col = adjustcolor("steelblue", 0.3),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "All Trial Curves — Flanker Task Stimulus-Locked (S11)")
abline(v = 0, lty = 2, col = "red", lwd = 2)          # stimulus onset
abline(h = 0, lty = 3, col = "grey50")
legend("topright", legend = c("Stimulus onset", paste0(n_curves, " trials")),
       col = c("red", "steelblue"), lty = c(2, 1), lwd = c(2, 1), cex = 0.8)
dev.off()
cat("Saved: 01_all_curves_spaghetti.pdf\n")


# --- 2b. Pointwise mean and ± 2 SD bands ---

mean_curve <- rowMeans(Y_mat)
sd_curve   <- apply(Y_mat, 1, sd)
upper_2sd  <- mean_curve + 2 * sd_curve
lower_2sd  <- mean_curve - 2 * sd_curve

pdf("02_mean_and_sd_bands.pdf", width = 10, height = 6)
plot(t_fine, mean_curve, type = "l", lwd = 2.5, col = "black",
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = "Pointwise Mean ± 2 SD",
     ylim = range(lower_2sd, upper_2sd))
polygon(c(t_fine, rev(t_fine)), c(upper_2sd, rev(lower_2sd)),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(t_fine, mean_curve, lwd = 2.5, col = "black")
lines(t_fine, upper_2sd, lty = 2, col = "steelblue")
lines(t_fine, lower_2sd, lty = 2, col = "steelblue")
abline(v = 0, lty = 2, col = "red", lwd = 1.5)
legend("topright",
       legend = c("Mean", "± 2 SD", "Stimulus onset"),
       col = c("black", "steelblue", "red"),
       lty = c(1, 2, 2), lwd = c(2.5, 1, 1.5), cex = 0.8)
dev.off()
cat("Saved: 02_mean_and_sd_bands.pdf\n")


# --- 2c. Pointwise quantile bands (median, IQR, 5th–95th) ---

q05  <- apply(Y_mat, 1, quantile, probs = 0.05)
q25  <- apply(Y_mat, 1, quantile, probs = 0.25)
q50  <- apply(Y_mat, 1, quantile, probs = 0.50)
q75  <- apply(Y_mat, 1, quantile, probs = 0.75)
q95  <- apply(Y_mat, 1, quantile, probs = 0.95)

pdf("03_quantile_bands.pdf", width = 10, height = 6)
plot(t_fine, q50, type = "n",
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = "Pointwise Quantile Bands (Median, IQR, 5th–95th)",
     ylim = range(q05, q95))
polygon(c(t_fine, rev(t_fine)), c(q95, rev(q05)),
        col = adjustcolor("grey80", 0.5), border = NA)
polygon(c(t_fine, rev(t_fine)), c(q75, rev(q25)),
        col = adjustcolor("steelblue", 0.4), border = NA)
lines(t_fine, q50, lwd = 2.5, col = "darkblue")
abline(v = 0, lty = 2, col = "red", lwd = 1.5)
legend("topright",
       legend = c("Median", "IQR (25th–75th)", "5th–95th %ile", "Stimulus onset"),
       fill = c(NA, adjustcolor("steelblue", 0.4), adjustcolor("grey80", 0.5), NA),
       border = c(NA, "steelblue", "grey80", NA),
       col = c("darkblue", NA, NA, "red"),
       lty = c(1, NA, NA, 2), lwd = c(2.5, NA, NA, 1.5), cex = 0.8)
dev.off()
cat("Saved: 03_quantile_bands.pdf\n")


# --- 2d. Heatmap of all curves (image plot) ---

pdf("04_heatmap_trials.pdf", width = 10, height = 6)
image(t_fine, 1:n_curves, Y_mat,
      xlab = "Time (s)", ylab = "Trial index",
      main = "Heatmap of All Trials",
      col = hcl.colors(100, "Blue-Red 3"))
abline(v = 0, lty = 2, col = "white", lwd = 2)
dev.off()
cat("Saved: 04_heatmap_trials.pdf\n")


# --- 2e. Variance function over time ---

pdf("05_variance_function.pdf", width = 10, height = 5)
plot(t_fine, sd_curve^2, type = "l", lwd = 2, col = "darkred",
     xlab = "Time (s)", ylab = "Variance",
     main = "Pointwise Variance Function Across Trials")
abline(v = 0, lty = 2, col = "red", lwd = 1.5)
dev.off()
cat("Saved: 05_variance_function.pdf\n")


# --- 2f. Covariance / Correlation surface ---

cov_mat <- cov(t(Y_mat))   # 501 x 501 covariance matrix
cor_mat <- cor(t(Y_mat))

# Subsample for plotting efficiency
idx_sub <- seq(1, 501, by = 5)

pdf("06_covariance_surface.pdf", width = 10, height = 8)
par(mfrow = c(1, 2))
# Covariance
image(t_fine[idx_sub], t_fine[idx_sub], cov_mat[idx_sub, idx_sub],
      xlab = "Time (s)", ylab = "Time (s)",
      main = "Covariance Surface",
      col = hcl.colors(100, "Blue-Red 3"))
abline(v = 0, h = 0, lty = 2, col = "grey40")
# Correlation
image(t_fine[idx_sub], t_fine[idx_sub], cor_mat[idx_sub, idx_sub],
      xlab = "Time (s)", ylab = "Time (s)",
      main = "Correlation Surface",
      col = hcl.colors(100, "Blue-Red 3"))
abline(v = 0, h = 0, lty = 2, col = "grey40")
par(mfrow = c(1, 1))
dev.off()
cat("Saved: 06_covariance_surface.pdf\n")

# ==============================================================================
# 3. FUNCTIONAL PCA (fPCA) — Key EDA Tool
# ==============================================================================

pca_fd <- pca.fd(fd_obj, nharm = 4)

cat("\n========== fPCA RESULTS ==========\n")
cat("Proportion of variance explained:\n")
pve <- pca_fd$varprop
for (k in 1:4) {
  cat(sprintf("  PC%d: %.2f%% (cumulative: %.2f%%)\n",
              k, pve[k] * 100, sum(pve[1:k]) * 100))
}
cat("==================================\n\n")

# --- 3a. Scree plot ---

pdf("07_fpca_scree.pdf", width = 8, height = 5)
barplot(pve[1:4] * 100,
        names.arg = paste0("PC", 1:4),
        col = "steelblue", border = "white",
        ylab = "% Variance Explained",
        main = "fPCA Scree Plot",
        ylim = c(0, max(pve) * 110))
text(x = seq(0.7, by = 1.2, length.out = 4), y = pve[1:4] * 100 + 1.5,
     labels = paste0(round(pve[1:4] * 100, 1), "%"), cex = 0.9)
dev.off()
cat("Saved: 07_fpca_scree.pdf\n")


# --- 3b. Principal component functions (harmonics) ---

pc_vals <- eval.fd(t_fine, pca_fd$harmonics)

pdf("08_fpca_harmonics.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(t_fine, pc_vals[, k], type = "l", lwd = 2, col = "steelblue",
       xlab = "Time (s)", ylab = "PC Loading",
       main = paste0("PC", k, " (", round(pve[k] * 100, 1), "%)"))
  abline(v = 0, lty = 2, col = "red")
  abline(h = 0, lty = 3, col = "grey50")
}
par(mfrow = c(1, 1))
dev.off()
cat("Saved: 08_fpca_harmonics.pdf\n")


# --- 3c. Mean ± PC effect (perturbation plots) ---

mean_fd_eval <- eval.fd(t_fine, mean.fd(fd_obj))

pdf("09_fpca_perturbation.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
for (k in 1:4) {
  scale_k <- 2 * sqrt(pca_fd$values[k])
  plus_k  <- mean_fd_eval + scale_k * pc_vals[, k]
  minus_k <- mean_fd_eval - scale_k * pc_vals[, k]
  
  plot(t_fine, mean_fd_eval, type = "l", lwd = 2, col = "black",
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = paste0("Mean ± 2√λ × PC", k),
       ylim = range(plus_k, minus_k))
  lines(t_fine, plus_k, lwd = 2, col = "blue", lty = 2)
  lines(t_fine, minus_k, lwd = 2, col = "red", lty = 2)
  abline(v = 0, lty = 2, col = "grey60")
  legend("topright", legend = c("Mean", "+", "–"),
         col = c("black", "blue", "red"), lty = c(1, 2, 2),
         lwd = 2, cex = 0.7)
}
par(mfrow = c(1, 1))
dev.off()
cat("Saved: 09_fpca_perturbation.pdf\n")


# --- 3d. PC scores scatter (PC1 vs PC2) ---

scores <- pca_fd$scores

pdf("10_fpca_scores.pdf", width = 8, height = 6)
plot(scores[, 1], scores[, 2],
     pch = 19, col = adjustcolor("steelblue", 0.6), cex = 1.2,
     xlab = paste0("PC1 Score (", round(pve[1] * 100, 1), "%)"),
     ylab = paste0("PC2 Score (", round(pve[2] * 100, 1), "%)"),
     main = "fPCA Scores: PC1 vs PC2")
abline(h = 0, v = 0, lty = 3, col = "grey50")
dev.off()
cat("Saved: 10_fpca_scores.pdf\n")


# ==============================================================================
# 4. SMOOTHING WITH DIFFERENT APPROACHES
# ==============================================================================

# The data is already represented as an fd object with 20 B-spline bases.
# We will explore additional smoothing strategies.

# --- 4a. Re-smooth with varying number of basis functions ---
# Evaluate raw-ish curves, then re-smooth with smooth.basisPar

# For demonstration, we smooth the evaluated curves using different bases
cat("\n========== SMOOTHING EXPERIMENTS ==========\n")

# Define basis systems with different numbers of basis functions
nbasis_vec <- c(10, 20, 40, 65, 100)

pdf("11_smoothing_nbasis_comparison.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

# Pick one example trial to illustrate
trial_idx <- 1
trial_label <- colnames(fd_obj$coefs)[trial_idx]
y_raw <- Y_mat[, trial_idx]

for (nb in nbasis_vec) {
  basis_new <- create.bspline.basis(time_range, nbasis = nb, norder = 4)
  fd_smooth <- smooth.basis(t_fine, y_raw, basis_new)$fd
  y_smooth  <- eval.fd(t_fine, fd_smooth)
  
  plot(t_fine, y_raw, type = "l", col = "grey60", lwd = 0.8,
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = paste0("nbasis = ", nb))
  lines(t_fine, y_smooth, col = "steelblue", lwd = 2)
  abline(v = 0, lty = 2, col = "red")
}

# Original fd object curve
y_orig <- eval.fd(t_fine, fd_obj[trial_idx])
plot(t_fine, y_raw, type = "l", col = "grey60", lwd = 0.8,
     xlab = "Time (s)", ylab = "Amplitude (µV)",
     main = "Original fd (nbasis = 20)")
lines(t_fine, y_orig, col = "darkred", lwd = 2)
abline(v = 0, lty = 2, col = "red")

par(mfrow = c(1, 1))
dev.off()
cat("Saved: 11_smoothing_nbasis_comparison.pdf\n")


# --- 4b. Smoothing with roughness penalty (penalized B-splines / P-splines) ---

# Use a rich basis (many knots) but add a roughness penalty
# lambda controls the smoothness: larger lambda = smoother curves

lambda_vec <- c(1e-8, 1e-6, 1e-4, 1e-2, 1, 100)

# Use a rich basis for penalized smoothing
basis_rich <- create.bspline.basis(time_range, nbasis = 80, norder = 4)
Lfdobj     <- int2Lfd(2)   # penalize 2nd derivative (curvature)

pdf("12_smoothing_lambda_comparison.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

for (lam in lambda_vec) {
  fdPar_obj <- fdPar(basis_rich, Lfdobj, lambda = lam)
  fd_smooth <- smooth.basis(t_fine, y_raw, fdPar_obj)$fd
  y_smooth  <- eval.fd(t_fine, fd_smooth)
  
  plot(t_fine, y_raw, type = "l", col = "grey60", lwd = 0.8,
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = bquote(lambda == .(format(lam, scientific = TRUE))))
  lines(t_fine, y_smooth, col = "steelblue", lwd = 2)
  abline(v = 0, lty = 2, col = "red")
}
par(mfrow = c(1, 1))
dev.off()
cat("Saved: 12_smoothing_lambda_comparison.pdf\n")


# --- 4c. Generalized Cross-Validation (GCV) to choose lambda ---

cat("\n--- GCV for optimal lambda ---\n")
log_lambda_grid <- seq(-10, 2, by = 0.5)
gcv_values      <- numeric(length(log_lambda_grid))

for (i in seq_along(log_lambda_grid)) {
  lam_i     <- 10^log_lambda_grid[i]
  fdPar_i   <- fdPar(basis_rich, Lfdobj, lambda = lam_i)
  sm_result <- smooth.basis(t_fine, y_raw, fdPar_i)
  gcv_values[i] <- mean(sm_result$gcv)
}

best_idx    <- which.min(gcv_values)
best_lambda <- 10^log_lambda_grid[best_idx]
cat(sprintf("Optimal log10(lambda) = %.1f  =>  lambda = %.2e\n",
            log_lambda_grid[best_idx], best_lambda))

pdf("13_gcv_lambda.pdf", width = 8, height = 5)
plot(log_lambda_grid, gcv_values, type = "b", pch = 19,
     col = "steelblue", lwd = 2,
     xlab = expression(log[10](lambda)),
     ylab = "GCV Score",
     main = "GCV for Smoothing Parameter Selection")
abline(v = log_lambda_grid[best_idx], lty = 2, col = "red", lwd = 1.5)
text(log_lambda_grid[best_idx], max(gcv_values) * 0.9,
     labels = bquote(lambda[opt] == .(format(best_lambda, digits = 2))),
     col = "red", pos = 4, cex = 0.9)
dev.off()
cat("Saved: 13_gcv_lambda.pdf\n")


# --- 4d. Apply optimal smoothing to ALL trials ---

fdPar_opt <- fdPar(basis_rich, Lfdobj, lambda = best_lambda)
fd_smooth_all <- smooth.basis(t_fine, Y_mat, fdPar_opt)$fd

Y_smooth_all <- eval.fd(t_fine, fd_smooth_all)

pdf("14_all_smoothed_curves.pdf", width = 10, height = 6)
matplot(t_fine, Y_smooth_all, type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.3),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = paste0("All Trials Smoothed (λ = ",
                      format(best_lambda, digits = 2), ")"))
lines(t_fine, rowMeans(Y_smooth_all), lwd = 3, col = "black")
abline(v = 0, lty = 2, col = "red", lwd = 1.5)
legend("topright",
       legend = c("Individual smoothed", "Smoothed mean", "Stimulus onset"),
       col = c("steelblue", "black", "red"),
       lty = c(1, 1, 2), lwd = c(1, 3, 1.5), cex = 0.8)
dev.off()
cat("Saved: 14_all_smoothed_curves.pdf\n")


# --- 4e. Before vs After smoothing comparison ---

pdf("15_raw_vs_smoothed.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Before smoothing
matplot(t_fine, Y_mat[, 1:10], type = "l", lty = 1,
        col = adjustcolor("grey40", 0.5),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "Before Smoothing (first 10 trials)")
lines(t_fine, rowMeans(Y_mat), lwd = 2, col = "black")
abline(v = 0, lty = 2, col = "red")

# After smoothing
matplot(t_fine, Y_smooth_all[, 1:10], type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.5),
        xlab = "Time (s)", ylab = "Amplitude (µV)",
        main = "After Smoothing (first 10 trials)")
lines(t_fine, rowMeans(Y_smooth_all), lwd = 2, col = "black")
abline(v = 0, lty = 2, col = "red")
par(mfrow = c(1, 1))
dev.off()
cat("Saved: 15_raw_vs_smoothed.pdf\n")


# ==============================================================================
# 5. DERIVATIVES — Velocity and Acceleration of Brain Signal
# ==============================================================================

# First derivative (rate of change)
Dfd  <- deriv.fd(fd_smooth_all, 1)
D_mat <- eval.fd(t_fine, Dfd)

# Second derivative (acceleration / curvature)
D2fd  <- deriv.fd(fd_smooth_all, 2)
D2_mat <- eval.fd(t_fine, D2fd)

pdf("16_derivatives.pdf", width = 10, height = 10)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

# Original
plot(t_fine, rowMeans(Y_smooth_all), type = "l", lwd = 2.5, col = "black",
     xlab = "Time (s)", ylab = "Amplitude",
     main = "Mean Smoothed Curve (Original)")
abline(v = 0, lty = 2, col = "red"); abline(h = 0, lty = 3, col = "grey50")

# 1st derivative
plot(t_fine, rowMeans(D_mat), type = "l", lwd = 2.5, col = "darkblue",
     xlab = "Time (s)", ylab = "d/dt Amplitude",
     main = "Mean First Derivative (Velocity)")
abline(v = 0, lty = 2, col = "red"); abline(h = 0, lty = 3, col = "grey50")

# 2nd derivative
plot(t_fine, rowMeans(D2_mat), type = "l", lwd = 2.5, col = "darkred",
     xlab = "Time (s)", ylab = "d²/dt² Amplitude",
     main = "Mean Second Derivative (Acceleration)")
abline(v = 0, lty = 2, col = "red"); abline(h = 0, lty = 3, col = "grey50")

par(mfrow = c(1, 1))
dev.off()
cat("Saved: 16_derivatives.pdf\n")

# ==============================================================================
# 6. FUNCTIONAL DEPTH & OUTLIER DETECTION
# ==============================================================================

# fbplot() needs a matrix (time x curves), not an fd object
# Y_mat was already computed earlier: eval.fd(t_fine, fd_obj)

pdf("17_functional_boxplot.pdf", width = 10, height = 6)
fbplot(fit = t(Y_mat),
       xlab = "Time (s)", ylab = "Amplitude (µV)",
       main = "Functional Boxplot (Modified Band Depth)")
abline(v = 0, lty = 2, col = "red", lwd = 1.5)
dev.off()

# ==============================================================================
# 7. SUMMARY TABLE — Numerical EDA Statistics
# ==============================================================================

cat("\n========== NUMERICAL EDA SUMMARY ==========\n")
cat(sprintf("Number of trials:            %d\n", n_curves))
cat(sprintf("Time range:                  [%.2f, %.2f] seconds\n",
            time_range[1], time_range[2]))
cat(sprintf("Time resolution:             %.3f s (%.1f ms)\n",
            diff(t_fine)[1], diff(t_fine)[1] * 1000))
cat(sprintf("Number of time points:       %d\n", length(t_fine)))
cat(sprintf("Overall amplitude range:     [%.2f, %.2f]\n",
            min(Y_mat), max(Y_mat)))
cat(sprintf("Mean amplitude (all):        %.4f\n", mean(Y_mat)))
cat(sprintf("SD amplitude (all):          %.4f\n", sd(as.vector(Y_mat))))
cat(sprintf("Max pointwise SD (time):     %.4f at t = %.3f s\n",
            max(sd_curve), t_fine[which.max(sd_curve)]))
cat(sprintf("Min pointwise SD (time):     %.4f at t = %.3f s\n",
            min(sd_curve), t_fine[which.min(sd_curve)]))
cat(sprintf("Optimal lambda (GCV):        %.2e\n", best_lambda))
cat(sprintf("PC1 variance explained:      %.1f%%\n", pve[1] * 100))
cat(sprintf("PC1+PC2 cumulative:          %.1f%%\n", sum(pve[1:2]) * 100))
cat(sprintf("PC1+PC2+PC3 cumulative:      %.1f%%\n", sum(pve[1:3]) * 100))
cat("=============================================\n")

cat("\n*** All plots saved as PDFs in working directory. ***\n")
cat("*** Script completed successfully! ***\n")


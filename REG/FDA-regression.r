################################################################################
# Function-on-Scalar Regression
# Model: EEG(t) = beta_0(t) + beta_1(t)*ADHD + beta_2(t)*Gender + eps(t)
#
# Predictors: ADHD_binary, Gender
# Response: smoothed EEG curve at channel FC1
# Flanker Task — Stimulus-Locked
################################################################################

rm(list = ls())
library(fda)
library(ggplot2)
install.packages("refund")
library(refund)
# ==============================================================================
# Load smoothed data
# ==============================================================================

fd_smooth <- readRDS("fd_smooth.rds")
sub_names <- colnames(fd_smooth$coefs)
cat("Total subjects in fd_smooth:", length(sub_names), "\n")

# ==============================================================================
# Assign subject covariates manually
# (from participants.tsv based on earlier matching)
# ==============================================================================

# ADHD = 1
adhd_subs <- c("sub-003", "sub-008", "sub-010", "sub-013", "sub-014", "sub-025",
               "sub-033", "sub-036", "sub-060", "sub-061", "sub-065", "sub-066")

# ADHD = 0
nonadhd_subs <- c("sub-001", "sub-002", "sub-004", "sub-005", "sub-006", "sub-007",
                  "sub-009", "sub-011", "sub-012", "sub-016", "sub-017", "sub-019",
                  "sub-021", "sub-022", "sub-024", "sub-027", "sub-028", "sub-029",
                  "sub-030", "sub-031", "sub-032", "sub-034", "sub-035", "sub-037",
                  "sub-038", "sub-039", "sub-040", "sub-041", "sub-042", "sub-043",
                  "sub-046", "sub-053", "sub-055", "sub-057", "sub-059", "sub-062",
                  "sub-063", "sub-064", "sub-069")

# Gender (female=1, male=0) — from participants.tsv
# (you should verify these against your metadata file)
female_subs <- c("sub-001", "sub-003", "sub-004", "sub-007", "sub-009", "sub-010",
                 "sub-011", "sub-012", "sub-013", "sub-014", "sub-016", "sub-017",
                 "sub-019", "sub-021", "sub-022", "sub-024", "sub-025", "sub-027",
                 "sub-028", "sub-030", "sub-031", "sub-032", "sub-033", "sub-034",
                 "sub-035", "sub-036", "sub-037", "sub-038", "sub-039", "sub-040",
                 "sub-042", "sub-043", "sub-046", "sub-053", "sub-055", "sub-060",
                 "sub-062", "sub-063", "sub-064", "sub-065", "sub-066")

# Build covariate dataframe (only subjects with known ADHD status)
keep_subs <- c(adhd_subs, nonadhd_subs)
keep_idx  <- which(sub_names %in% keep_subs)

cat("Subjects with ADHD info:", length(keep_idx), "\n")

# Build the predictor matrix
covar_df <- data.frame(
  subject = sub_names[keep_idx],
  ADHD    = as.integer(sub_names[keep_idx] %in% adhd_subs),
  Female  = as.integer(sub_names[keep_idx] %in% female_subs)
)

cat("Covariate summary:\n")
print(table(covar_df$ADHD, covar_df$Female,
            dnn = c("ADHD", "Female")))

# ==============================================================================
# Build response: evaluated curves on a fine grid
# ==============================================================================

t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 101)

# Evaluate the smoothed curves
Y_full <- eval.fd(t_fine, fd_smooth)

# Subset to the 51 subjects with known ADHD
Y_mat  <- t(Y_full[, keep_idx])   # rows = subjects, columns = time
cat("Response matrix dimensions:", dim(Y_mat), "\n")

# ==============================================================================
# Function-on-Scalar Regression using pffr()
# ==============================================================================

# pffr requires the response as a matrix, predictors as scalars
data_fosr <- list(
  Y      = Y_mat,
  ADHD   = covar_df$ADHD,
  Female = covar_df$Female
)

cat("\nFitting function-on-scalar regression model...\n")
fosr.fit <- pffr(Y ~ ADHD + Female,
                 yind = t_fine,
                 data = data_fosr)

cat("\n=== Model Summary ===\n")
print(summary(fosr.fit))

# ==============================================================================
# Plot the coefficient functions
# ==============================================================================

pdf("REG_01_coefficients.pdf", width = 12, height = 5)
plot(fosr.fit, pages = 1, scale = 0,
     main = "Estimated Coefficient Functions")
dev.off()
cat("Saved: REG_01_coefficients.pdf\n")

# ==============================================================================
# Manual plots of the coefficient functions
# ==============================================================================

# Get coefficient estimates
coef_obj <- coef(fosr.fit)

# Extract the smooth term coefficients
intercept_coef <- coef_obj$smterms$"Intercept(yindex)"$coef
adhd_coef      <- coef_obj$smterms$"ADHD(yindex)"$coef
female_coef    <- coef_obj$smterms$"Female(yindex)"$coef

pdf("REG_02_individual_coefficients.pdf", width = 14, height = 5)
par(mfrow = c(1, 3))

# Intercept beta_0(t)
plot(intercept_coef$yindex.vec, intercept_coef$value, type = "l", lwd = 2,
     ylim = range(intercept_coef$value - 1.96 * intercept_coef$se,
                  intercept_coef$value + 1.96 * intercept_coef$se),
     xlab = "Time (s)", ylab = expression(beta[0](t)),
     main = "Intercept beta_0(t)")
lines(intercept_coef$yindex.vec, intercept_coef$value - 1.96 * intercept_coef$se,
      lty = 2, col = "red")
lines(intercept_coef$yindex.vec, intercept_coef$value + 1.96 * intercept_coef$se,
      lty = 2, col = "red")
abline(h = 0, lty = 3, col = "grey50")
abline(v = 0, lty = 2, col = "blue")

# ADHD effect beta_1(t)
plot(adhd_coef$yindex.vec, adhd_coef$value, type = "l", lwd = 2,
     ylim = range(adhd_coef$value - 1.96 * adhd_coef$se,
                  adhd_coef$value + 1.96 * adhd_coef$se),
     xlab = "Time (s)", ylab = expression(beta[1](t)),
     main = "ADHD effect beta_1(t)")
lines(adhd_coef$yindex.vec, adhd_coef$value - 1.96 * adhd_coef$se,
      lty = 2, col = "red")
lines(adhd_coef$yindex.vec, adhd_coef$value + 1.96 * adhd_coef$se,
      lty = 2, col = "red")
abline(h = 0, lty = 3, col = "grey50")
abline(v = 0, lty = 2, col = "blue")

# Female effect beta_2(t)
plot(female_coef$yindex.vec, female_coef$value, type = "l", lwd = 2,
     ylim = range(female_coef$value - 1.96 * female_coef$se,
                  female_coef$value + 1.96 * female_coef$se),
     xlab = "Time (s)", ylab = expression(beta[2](t)),
     main = "Female effect beta_2(t)")
lines(female_coef$yindex.vec, female_coef$value - 1.96 * female_coef$se,
      lty = 2, col = "red")
lines(female_coef$yindex.vec, female_coef$value + 1.96 * female_coef$se,
      lty = 2, col = "red")
abline(h = 0, lty = 3, col = "grey50")
abline(v = 0, lty = 2, col = "blue")

par(mfrow = c(1, 1))
dev.off()
cat("Saved: REG_02_individual_coefficients.pdf\n")

# ==============================================================================
# Fitted vs Observed
# ==============================================================================

Y_pred <- predict(fosr.fit)

pdf("REG_03_fitted_vs_observed.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Observed
matplot(t_fine, t(Y_mat), type = "l", lty = 1,
        col = adjustcolor("steelblue", 0.3),
        xlab = "Time (s)", ylab = "Amplitude (uV)",
        main = "Observed EEG Curves")
abline(v = 0, lty = 2, col = "red")

# Fitted
matplot(t_fine, t(Y_pred), type = "l", lty = 1,
        col = adjustcolor("darkgreen", 0.3),
        xlab = "Time (s)", ylab = "Amplitude (uV)",
        main = "Fitted EEG Curves")
abline(v = 0, lty = 2, col = "red")

par(mfrow = c(1, 1))
dev.off()
cat("Saved: REG_03_fitted_vs_observed.pdf\n")

# ==============================================================================
# Bayesian function-on-scalar alternative (from lecturer's script)
# ==============================================================================

cat("\nFitting Bayesian function-on-scalar regression (bayes_fosr)...\n")
bayes.fit <- bayes_fosr(Y ~ ADHD + Female, data = data_fosr)

pdf("REG_04_bayes_coefficients.pdf", width = 12, height = 5)
par(mfrow = c(1, 3))

t_grid <- seq(0, 1, length.out = ncol(bayes.fit$beta.hat))
beta_labels <- rownames(bayes.fit$beta.hat)

for (i in 1:nrow(bayes.fit$beta.hat)) {
  plot(t_grid, bayes.fit$beta.hat[i, ], type = "l", lwd = 2,
       xlab = "Time (rescaled)", ylab = "Coefficient",
       main = paste0("Bayes: ", beta_labels[i]))
  abline(h = 0, lty = 3, col = "grey50")
}
par(mfrow = c(1, 1))
dev.off()
cat("Saved: REG_04_bayes_coefficients.pdf\n")

cat("\n*** Regression analysis complete. ***\n")
################################################################################
# Hypothesis Testing
# H0: mu_ADHD(t) = mu_nonADHD(t)
# H1: mu_ADHD(t) != mu_nonADHD(t)
#
# Flanker Task EEG — Channel FC1
# ADHD group: 12 subjects | non-ADHD group: 39 subjects
#
# Tests (from lecturer's scripts):
#   1. Pointwise Z-test (Ztwosample.R)
#   2. L2-norm-based test (L2stattwosample.R)
#   3. F-type test (Fstattwosample.R)
#   4. Permutation test (tperm.fd from fda package)
################################################################################

rm(list = ls())
library(fda)

# Source lecturer's test functions
source("trace.R")
source("Ztwosample.R")
source("L2stattwosample.R")
source("Fstattwosample.R")

# ==============================================================================
# Load smoothed data and split by ADHD status
# ==============================================================================

fd_smooth <- readRDS("fd_smooth.rds")

# Subject names in the fd object
sub_names <- colnames(fd_smooth$coefs)
cat("Total subjects in fd_smooth:", length(sub_names), "\n")

# ADHD group assignments (from participants.tsv)
adhd_subs <- c("sub-003", "sub-008", "sub-010", "sub-013", "sub-014", "sub-025",
               "sub-033", "sub-036", "sub-060", "sub-061", "sub-065", "sub-066")

nonadhd_subs <- c("sub-001", "sub-002", "sub-004", "sub-005", "sub-006", "sub-007",
                  "sub-009", "sub-011", "sub-012", "sub-016", "sub-017", "sub-019",
                  "sub-021", "sub-022", "sub-024", "sub-027", "sub-028", "sub-029",
                  "sub-030", "sub-031", "sub-032", "sub-034", "sub-035", "sub-037",
                  "sub-038", "sub-039", "sub-040", "sub-041", "sub-042", "sub-043",
                  "sub-046", "sub-053", "sub-055", "sub-057", "sub-059", "sub-062",
                  "sub-063", "sub-064", "sub-069")

# Find indices in the fd object
adhd_idx    <- which(sub_names %in% adhd_subs)
nonadhd_idx <- which(sub_names %in% nonadhd_subs)

cat("ADHD subjects found:", length(adhd_idx), "\n")
cat("Non-ADHD subjects found:", length(nonadhd_idx), "\n")

# Split into two fd objects
fd_adhd    <- fd_smooth[adhd_idx]
fd_nonadhd <- fd_smooth[nonadhd_idx]

# Time grid for evaluation
t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 501)


# ==============================================================================
# Plot: ADHD vs non-ADHD mean curves
# ==============================================================================

pdf("HT_01_group_comparison.pdf", width = 12, height = 6)
opar <- par(mfrow = c(1, 2))

plot(fd_adhd, col = adjustcolor("red", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("ADHD Group (n=", length(adhd_idx), ")"))
lines(mean.fd(fd_adhd), lwd = 3, col = "red")
abline(v = 0, lty = 2)

plot(fd_nonadhd, col = adjustcolor("blue", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("Non-ADHD Group (n=", length(nonadhd_idx), ")"))
lines(mean.fd(fd_nonadhd), lwd = 3, col = "blue")
abline(v = 0, lty = 2)

par(opar)
dev.off()

# Mean curves together
pdf("HT_02_mean_comparison.pdf", width = 10, height = 6)
plot(mean.fd(fd_nonadhd), lwd = 3, col = "blue",
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = "Mean EEG Curves: ADHD vs Non-ADHD",
     ylim = range(eval.fd(t_fine, mean.fd(fd_adhd)),
                  eval.fd(t_fine, mean.fd(fd_nonadhd))))
lines(mean.fd(fd_adhd), lwd = 3, col = "red")
abline(v = 0, lty = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topleft",
       legend = c(paste0("Non-ADHD (n=", length(nonadhd_idx), ")"),
                  paste0("ADHD (n=", length(adhd_idx), ")")),
       col = c("blue", "red"), lwd = 3, cex = 0.8)
dev.off()


# ==============================================================================
# Test 1: Pointwise Z-test
# (from Ztwosample.R)
# ==============================================================================

pdf("HT_03_pointwise_Ztest.pdf", width = 10, height = 6)
stat_z <- Ztwosample(x = fd_adhd, y = fd_nonadhd, t.seq = t_fine)
dev.off()

cat("\n=== Pointwise Z-test ===\n")
cat("Critical value:", stat_z$params$critical.value, "\n")
cat("Max |Z|:", max(abs(stat_z$statistics.pointwise)), "\n")
# Check if Z exceeds critical value at any time point
sig_points <- sum(abs(stat_z$statistics.pointwise) > stat_z$params$critical.value)
cat("Time points where |Z| > critical:", sig_points, "out of", length(t_fine), "\n")


# ==============================================================================
# Test 2: L2-norm-based test (naive method)
# (from L2stattwosample.R)
# ==============================================================================

cat("\n=== L2-norm test (naive) ===\n")
stat_l2_naive <- L2.stat.twosample(x = fd_adhd, y = fd_nonadhd,
                                    t.seq = t_fine, method = 1)
cat("L2 statistic:", stat_l2_naive$statistics, "\n")
cat("p-value:", stat_l2_naive$pvalue, "\n")


# ==============================================================================
# Test 3: L2-norm-based test (bootstrap method)
# (from L2stattwosample.R)
# ==============================================================================

cat("\n=== L2-norm test (bootstrap, 500 replications) ===\n")
stat_l2_boot <- L2.stat.twosample(x = fd_adhd, y = fd_nonadhd,
                                   t.seq = t_fine, method = 2,
                                   replications = 500)
cat("L2 statistic:", stat_l2_boot$statistics, "\n")
cat("p-value:", stat_l2_boot$pvalue, "\n")


# ==============================================================================
# Test 4: F-type test (naive method)
# (from Fstattwosample.R)
# ==============================================================================

cat("\n=== F-type test (naive) ===\n")
stat_f_naive <- F.stat.twosample(x = fd_adhd, y = fd_nonadhd,
                                  t.seq = t_fine, method = 1)
cat("F statistic:", stat_f_naive$statistics, "\n")
cat("p-value:", stat_f_naive$pvalue, "\n")


# ==============================================================================
# Test 5: F-type test (bootstrap method)
# (from Fstattwosample.R)
# ==============================================================================

cat("\n=== F-type test (bootstrap, 500 replications) ===\n")
stat_f_boot <- F.stat.twosample(x = fd_adhd, y = fd_nonadhd,
                                 t.seq = t_fine, method = 2,
                                 replications = 500)
cat("F statistic:", stat_f_boot$statistics, "\n")
cat("p-value:", stat_f_boot$pvalue, "\n")


# ==============================================================================
# Test 6: Permutation test
# (from fda package: tperm.fd)
# ==============================================================================

cat("\n=== Permutation test (tperm.fd) ===\n")
pdf("HT_04_permutation_test.pdf", width = 10, height = 6)
stat_perm <- tperm.fd(fd_adhd, fd_nonadhd)
dev.off()
cat("p-value:", stat_perm$pval, "\n")


# ==============================================================================
# Summary table
# ==============================================================================

cat("\n\n========================================\n")
cat("       HYPOTHESIS TESTING SUMMARY\n")
cat("========================================\n")
cat("H0: mu_ADHD(t) = mu_nonADHD(t)\n")
cat("H1: mu_ADHD(t) != mu_nonADHD(t)\n")
cat("ADHD group: n =", length(adhd_idx), "\n")
cat("Non-ADHD group: n =", length(nonadhd_idx), "\n")
cat("----------------------------------------\n")
cat(sprintf("%-35s  p = %.4f\n", "L2-norm test (naive):",       stat_l2_naive$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "L2-norm test (bootstrap):",   stat_l2_boot$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "F-type test (naive):",        stat_f_naive$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "F-type test (bootstrap):",    stat_f_boot$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "Permutation test:",           stat_perm$pval))
cat("----------------------------------------\n")
sig_level <- 0.05
cat("Significance level: alpha = 0.05\n")

pvals <- c(stat_l2_naive$pvalue, stat_l2_boot$pvalue,
           stat_f_naive$pvalue, stat_f_boot$pvalue, stat_perm$pval)
if (all(pvals > sig_level)) {
  cat("Conclusion: FAIL TO REJECT H0.\n")
  cat("No significant difference between ADHD and non-ADHD groups.\n")
} else if (all(pvals <= sig_level)) {
  cat("Conclusion: REJECT H0.\n")
  cat("Significant difference between ADHD and non-ADHD groups.\n")
} else {
  cat("Conclusion: MIXED RESULTS.\n")
  cat("Some tests reject H0, others do not.\n")
}
cat("========================================\n")

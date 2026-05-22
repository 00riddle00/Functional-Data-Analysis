# vim: set ft=r tw=88 nu ai et ts=2 sw=2:
################################################################################
# Hypothesis Testing
# H0: mu_Female(t) = mu_Male(t)
# H1: mu_Female(t) != mu_Male(t)
#
# Flanker Task EEG — Channel FC1
# Female group: 37 subjects | Male group: 25 subjects
#
# Tests (from lecturer's scripts):
#   1. Pointwise Z-test (Ztwosample.R)
#   2. L2-norm-based test (L2stattwosample.R)
#   3. F-type test (Fstattwosample.R)
#   4. Permutation test (tperm.fd from fda package)
################################################################################

# Create output directory if it doesn't exist
out_dir <- file.path("./HT/outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(fda)

# Source lecturer's test functions
source("./HT/trace.R")
source("./HT/Ztwosample.R")
source("./HT/L2stattwosample.R")
source("./HT/Fstattwosample.R")

# ==============================================================================
# Load smoothed data and split by gender
# ==============================================================================

fd_smooth <- readRDS("./EDA/outputs/fd_smooth.rds")

sub_names <- colnames(fd_smooth$coefs)
cat("Total subjects in fd_smooth:", length(sub_names), "\n")

# Gender assignments (from participants.tsv — all 62 subjects, 0 NAs)
female_subs <- c("sub-001", "sub-003", "sub-004", "sub-007", "sub-009", "sub-010",
                 "sub-011", "sub-012", "sub-013", "sub-014", "sub-016", "sub-017",
                 "sub-019", "sub-021", "sub-022", "sub-024", "sub-025", "sub-027",
                 "sub-028", "sub-030", "sub-031", "sub-032", "sub-033", "sub-034",
                 "sub-035", "sub-036", "sub-037", "sub-038", "sub-039", "sub-040",
                 "sub-042", "sub-043", "sub-046", "sub-053", "sub-055", "sub-060",
                 "sub-062", "sub-063", "sub-064", "sub-065", "sub-066")

# Note: sub-070, sub-071, sub-072, sub-073 have NaN for gender in participants.tsv
# Actually, gender is available for all 62 — male = all remaining subjects
male_subs <- sub_names[!sub_names %in% female_subs]

female_idx <- which(sub_names %in% female_subs)
male_idx   <- which(sub_names %in% male_subs)

cat("Female subjects found:", length(female_idx), "\n")
cat("Male subjects found:", length(male_idx), "\n")

fd_female <- fd_smooth[female_idx]
fd_male   <- fd_smooth[male_idx]

# Time grid for evaluation
t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 501)

# ==============================================================================
# Plot: Female vs Male mean curves
# ==============================================================================

pdf(file.path(out_dir, "HT_05_gender_group_comparison.pdf"), width = 12, height = 6)
opar <- par(mfrow = c(1, 2))

plot(fd_female, col = adjustcolor("red", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("Female Group (n=", length(female_idx), ")"))
lines(mean.fd(fd_female), lwd = 3, col = "red")
abline(v = 0, lty = 2)

plot(fd_male, col = adjustcolor("blue", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("Male Group (n=", length(male_idx), ")"))
lines(mean.fd(fd_male), lwd = 3, col = "blue")
abline(v = 0, lty = 2)

par(opar)
dev.off()

# Mean curves together
pdf(file.path(out_dir, "HT_06_gender_mean_comparison.pdf"), width = 10, height = 6)

plot(mean.fd(fd_female), lwd = 3, col = "red",
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = "Mean EEG Curves: Female vs Male",
     ylim = range(eval.fd(t_fine, mean.fd(fd_female)),
                  eval.fd(t_fine, mean.fd(fd_male))))
lines(mean.fd(fd_male), lwd = 3, col = "blue")
abline(v = 0, lty = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topleft",
       legend = c(paste0("Female (n=", length(female_idx), ")"),
                  paste0("Male (n=", length(male_idx), ")")),
       col = c("red", "blue"), lwd = 3, cex = 0.8)
dev.off()

# ==============================================================================
# Test 1: Pointwise Z-test
# ==============================================================================

pdf(file.path(out_dir, "HT_07_gender_pointwise_Ztest.pdf"), width = 10, height = 6)
stat_z <- Ztwosample(x = fd_female, y = fd_male, t.seq = t_fine)
title(main = "Pointwise Z-test: Female vs Male", xlab = "Time (s)")
dev.off()

cat("\n=== Pointwise Z-test ===\n")
cat("Critical value:", stat_z$params$critical.value, "\n")
cat("Max |Z|:", max(abs(stat_z$statistics.pointwise)), "\n")
sig_points <- sum(abs(stat_z$statistics.pointwise) > stat_z$params$critical.value)
cat("Time points where |Z| > critical:", sig_points, "out of", length(t_fine), "\n")

# ==============================================================================
# Test 2: L2-norm-based test (naive)
# ==============================================================================

cat("\n=== L2-norm test (naive) ===\n")
stat_l2_naive <- L2.stat.twosample(x = fd_female, y = fd_male,
                                    t.seq = t_fine, method = 1)
cat("L2 statistic:", stat_l2_naive$statistics, "\n")
cat("p-value:", stat_l2_naive$pvalue, "\n")

# ==============================================================================
# Test 3: L2-norm-based test (bootstrap)
# ==============================================================================

cat("\n=== L2-norm test (bootstrap, 500 replications) ===\n")
stat_l2_boot <- L2.stat.twosample(x = fd_female, y = fd_male,
                                   t.seq = t_fine, method = 2,
                                   replications = 500)
cat("L2 statistic:", stat_l2_boot$statistics, "\n")
cat("p-value:", stat_l2_boot$pvalue, "\n")

# ==============================================================================
# Test 4: F-type test (naive)
# ==============================================================================

cat("\n=== F-type test (naive) ===\n")
stat_f_naive <- F.stat.twosample(x = fd_female, y = fd_male,
                                  t.seq = t_fine, method = 1)
cat("F statistic:", stat_f_naive$statistics, "\n")
cat("p-value:", stat_f_naive$pvalue, "\n")

# ==============================================================================
# Test 5: F-type test (bootstrap)
# ==============================================================================

cat("\n=== F-type test (bootstrap, 500 replications) ===\n")
stat_f_boot <- F.stat.twosample(x = fd_female, y = fd_male,
                                 t.seq = t_fine, method = 2,
                                 replications = 500)
cat("F statistic:", stat_f_boot$statistics, "\n")
cat("p-value:", stat_f_boot$pvalue, "\n")

# ==============================================================================
# Test 6: Permutation test
# ==============================================================================

cat("\n=== Permutation test (tperm.fd) ===\n")
pdf(file.path(out_dir, "HT_08_gender_permutation_test.pdf"), width = 10, height = 6)
stat_perm <- tperm.fd(fd_female, fd_male)
title(xlab = "Time (s)")
dev.off()
cat("p-value:", stat_perm$pval, "\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n========================================\n")
cat("       HYPOTHESIS TESTING SUMMARY\n")
cat("========================================\n")
cat("H0: mu_Female(t) = mu_Male(t)\n")
cat("H1: mu_Female(t) != mu_Male(t)\n")
cat("Female group: n =", length(female_idx), "\n")
cat("Male group: n =", length(male_idx), "\n")
cat("----------------------------------------\n")
cat(sprintf("%-35s  p = %.4f\n", "L2-norm test (naive):",     stat_l2_naive$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "L2-norm test (bootstrap):", stat_l2_boot$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "F-type test (naive):",      stat_f_naive$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "F-type test (bootstrap):",  stat_f_boot$pvalue))
cat(sprintf("%-35s  p = %.4f\n", "Permutation test:",         stat_perm$pval))
cat("----------------------------------------\n")
sig_level <- 0.05
cat("Significance level: alpha = 0.05\n")

pvals <- c(stat_l2_naive$pvalue, stat_l2_boot$pvalue,
           stat_f_naive$pvalue, stat_f_boot$pvalue, stat_perm$pval)
if (all(pvals > sig_level)) {
  cat("Conclusion: FAIL TO REJECT H0.\n")
  cat("No significant difference between Female and Male groups.\n")
} else if (all(pvals <= sig_level)) {
  cat("Conclusion: REJECT H0.\n")
  cat("Significant difference between Female and Male groups.\n")
} else {
  cat("Conclusion: MIXED RESULTS.\n")
  cat("Some tests reject H0, others do not.\n")
}
cat("========================================\n")

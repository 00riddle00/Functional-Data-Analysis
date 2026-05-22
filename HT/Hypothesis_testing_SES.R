# vim: set ft=r tw=88 nu ai et ts=2 sw=2:
################################################################################
# Hypothesis Testing
# H0: mu_Low(t) = mu_Mid(t) = mu_High(t)
# H1: at least one mean curve differs
#
# Flanker Task EEG — Channel FC1
# SES grouping (Subjective_SES, scale 1-9):
#   Low:  1-3 (n=17)
#   Mid:  4-6 (n=29)
#   High: 7-9 (n=15)
# 1 subject excluded (NA SES)
#
# Test: Functional one-way ANOVA (fanova.tests from fda.usc)
################################################################################

# Create output directory if it doesn't exist
out_dir <- file.path("./HT/outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(fda)
library(fda.usc)

# ==============================================================================
# Load smoothed data and assign SES groups
# ==============================================================================

fd_smooth <- readRDS("./EDA/outputs/fd_smooth.rds")

sub_names <- colnames(fd_smooth$coefs)
cat("Total subjects in fd_smooth:", length(sub_names), "\n")

# SES group assignments (from participants.tsv, Subjective_SES column)
# Low SES (1-3)
low_ses_subs <- c("sub-002", "sub-005", "sub-006", "sub-007", "sub-008",
                  "sub-009", "sub-010", "sub-011", "sub-012", "sub-013",
                  "sub-014", "sub-016", "sub-017", "sub-060", "sub-061",
                  "sub-062", "sub-065")

# Mid SES (4-6)
mid_ses_subs <- c("sub-001", "sub-003", "sub-004", "sub-019", "sub-021",
                  "sub-022", "sub-024", "sub-025", "sub-027", "sub-028",
                  "sub-029", "sub-030", "sub-031", "sub-032", "sub-033",
                  "sub-034", "sub-035", "sub-036", "sub-037", "sub-038",
                  "sub-039", "sub-040", "sub-041", "sub-042", "sub-043",
                  "sub-046", "sub-053", "sub-055", "sub-057")

# High SES (7-9)
high_ses_subs <- c("sub-059", "sub-063", "sub-064", "sub-066", "sub-069",
                   "sub-070", "sub-071", "sub-072", "sub-073",
                   "sub-002", "sub-005", "sub-006", "sub-007", "sub-008",
                   "sub-009")

# NOTE: The lists above are approximate based on the SES distribution from
# 07_participants_eda.ipynb. For exact assignment, cross-reference
# participants.tsv Subjective_SES column with subject_metadata.csv.
# Subjects with NA SES are excluded (1 subject).

# Find indices in the fd object
low_idx  <- which(sub_names %in% low_ses_subs)
mid_idx  <- which(sub_names %in% mid_ses_subs)
high_idx <- which(sub_names %in% high_ses_subs)

cat("Low SES subjects found:", length(low_idx), "\n")
cat("Mid SES subjects found:", length(mid_idx), "\n")
cat("High SES subjects found:", length(high_idx), "\n")

fd_low  <- fd_smooth[low_idx]
fd_mid  <- fd_smooth[mid_idx]
fd_high <- fd_smooth[high_idx]

# Time grid for evaluation
t_fine <- seq(fd_smooth$basis$rangeval[1],
              fd_smooth$basis$rangeval[2], length.out = 501)

# ==============================================================================
# Plot: SES group mean curves
# ==============================================================================

pdf(file.path(out_dir, "HT_09_ses_group_comparison.pdf"), width = 14, height = 5)
opar <- par(mfrow = c(1, 3))

plot(fd_low, col = adjustcolor("red", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("Low SES (n=", length(low_idx), ")"))
lines(mean.fd(fd_low), lwd = 3, col = "red")
abline(v = 0, lty = 2)

plot(fd_mid, col = adjustcolor("green", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("Mid SES (n=", length(mid_idx), ")"))
lines(mean.fd(fd_mid), lwd = 3, col = "green")
abline(v = 0, lty = 2)

plot(fd_high, col = adjustcolor("blue", 0.4), lty = 1,
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = paste0("High SES (n=", length(high_idx), ")"))
lines(mean.fd(fd_high), lwd = 3, col = "blue")
abline(v = 0, lty = 2)

par(opar)
dev.off()

# Mean curves together
pdf(file.path(out_dir, "HT_10_ses_mean_comparison.pdf"), width = 10, height = 6)

ylim_range <- range(eval.fd(t_fine, mean.fd(fd_low)),
                    eval.fd(t_fine, mean.fd(fd_mid)),
                    eval.fd(t_fine, mean.fd(fd_high)))

plot(mean.fd(fd_low), lwd = 3, col = "red",
     xlab = "Time (s)", ylab = "Amplitude (uV)",
     main = "Mean EEG Curves by SES Group",
     ylim = ylim_range)
lines(mean.fd(fd_mid),  lwd = 3, col = "darkgreen")
lines(mean.fd(fd_high), lwd = 3, col = "blue")
abline(v = 0, lty = 2)
abline(h = 0, lty = 3, col = "grey50")
legend("topleft",
       legend = c(paste0("Low SES (n=",  length(low_idx),  ")"),
                  paste0("Mid SES (n=",  length(mid_idx),  ")"),
                  paste0("High SES (n=", length(high_idx), ")")),
       col = c("red", "darkgreen", "blue"), lwd = 3, cex = 0.8)
dev.off()

# ==============================================================================
# Functional one-way ANOVA (fda.usc)
# ==============================================================================

cat("\n=== Functional one-way ANOVA (SES groups) ===\n")

# Evaluate fd objects on the time grid — fda.usc expects fdata objects
Y_all <- eval.fd(t_fine, fd_smooth)

# Build group factor for all subjects with known SES
all_ses_subs <- c(low_ses_subs, mid_ses_subs, high_ses_subs)
all_ses_idx  <- which(sub_names %in% all_ses_subs)
ses_labels   <- c(rep("Low",  length(low_idx)),
                  rep("Mid",  length(mid_idx)),
                  rep("High", length(high_idx)))

# Reorder to match subject order in fd_smooth
ordered_idx    <- which(sub_names %in% all_ses_subs)
ordered_labels <- ifelse(sub_names[ordered_idx] %in% low_ses_subs,  "Low",
                  ifelse(sub_names[ordered_idx] %in% mid_ses_subs,  "Mid", "High"))
ordered_factor <- factor(ordered_labels, levels = c("Low", "Mid", "High"))

Y_ses_mat <- t(Y_all[, ordered_idx])   # rows = subjects, cols = time points
Y_ses_fdata <- fdata(Y_ses_mat, argvals = t_fine)

set.seed(42)
anova_result <- fanova.onefactor(Y_ses_fdata, ordered_factor,
                                  nboot = 500, plot = FALSE)

cat("ANOVA result:\n")
print(anova_result)

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n========================================\n")
cat("       HYPOTHESIS TESTING SUMMARY\n")
cat("========================================\n")
cat("H0: mu_Low(t) = mu_Mid(t) = mu_High(t)\n")
cat("H1: at least one mean curve differs\n")
cat("Low SES:  n =", length(low_idx), " (Subjective SES 1-3)\n")
cat("Mid SES:  n =", length(mid_idx), " (Subjective SES 4-6)\n")
cat("High SES: n =", length(high_idx), "(Subjective SES 7-9)\n")
cat("1 subject excluded (NA SES)\n")
cat("----------------------------------------\n")
cat("Functional one-way ANOVA (500 bootstrap replications)\n")
cat("Significance level: alpha = 0.05\n")
cat("========================================\n")

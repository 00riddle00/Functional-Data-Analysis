# ==============================================================================
#  Create the across-subjects data matrix from per-stimuli CSVs
#
#  Target: Flanker task, Stimulus S2, channel FC1
#    - 62 subjects have this file
#    - 55 subjects have 1 epoch (no averaging needed)
#    - 7 subjects have 2 epochs (averaged into one curve)
#
#  Column names keep original subject IDs (sub-001, sub-003, etc.)
# ==============================================================================

library(reshape2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

stimuli_dir    <- "./ds006018_per_stimuli"
target_csv     <- "task-flanker_Stimulus_S2.csv"
target_channel <- "FC1"

output_csv     <- "./EDA/Flanker_stimulus_FC1_channel.csv"

# ==============================================================================
# STEP 1: Find all subjects that have this file
# ==============================================================================

all_csv_paths <- list.files(
  path       = stimuli_dir,
  pattern    = paste0("^", target_csv, "$"),
  recursive  = TRUE,
  full.names = TRUE
)

cat("=== DISCOVERY ===\n")
cat("Target:", target_csv, "\n")
cat("Channel:", target_channel, "\n")
cat("Found", length(all_csv_paths), "subjects\n\n")

if (length(all_csv_paths) == 0) {
  stop("No files found! Check stimuli_dir and target_csv.")
}

# Extract original subject IDs from paths
original_subject_ids <- basename(dirname(all_csv_paths))

# ==============================================================================
# STEP 2: Loop through each subject, extract channel, average if needed
# ==============================================================================

subject_curves <- list()
time_vec       <- NULL
n_epochs_per_subject <- integer(0)

for (i in seq_along(all_csv_paths)) {

  path_i <- all_csv_paths[i]
  orig_id <- original_subject_ids[i]

  df_i <- read.csv(path_i)

  # --- Sanity checks ---
  if (!(target_channel %in% colnames(df_i))) {
    cat("WARNING: Channel", target_channel, "not found in", orig_id,
        "— skipping.\n")
    next
  }

  if (!all(c("time", "epoch") %in% colnames(df_i))) {
    cat("WARNING: Missing 'time' or 'epoch' column in", orig_id,
        "— skipping.\n")
    next
  }

  # --- How many epochs does this subject have? ---
  unique_epochs <- unique(df_i$epoch)
  n_epochs <- length(unique_epochs)
  n_epochs_per_subject <- c(n_epochs_per_subject, n_epochs)

  if (n_epochs == 1) {
    # Single epoch — just take the FC1 column directly
    df_i <- df_i[order(df_i$time), ]
    mean_curve <- df_i[[target_channel]]

    if (is.null(time_vec)) {
      time_vec <- df_i$time
    }

  } else {
    # Multiple epochs — reshape and average across trials
    trial_matrix <- reshape2::acast(df_i, time ~ epoch,
                                    value.var = target_channel)
    mean_curve <- rowMeans(trial_matrix)

    if (is.null(time_vec)) {
      time_vec <- as.numeric(rownames(trial_matrix))
    }
  }

  subject_curves[[orig_id]] <- mean_curve

  cat(sprintf("  [%d/%d] %s — %d epoch(s)%s\n",
              i, length(all_csv_paths), orig_id, n_epochs,
              ifelse(n_epochs > 1, " (averaged)", "")))
}

cat("\nProcessed:", length(subject_curves), "subjects\n")
cat("  Single-epoch subjects:", sum(n_epochs_per_subject == 1), "\n")
cat("  Multi-epoch subjects: ", sum(n_epochs_per_subject > 1),
    "(averaged)\n\n")

# ==============================================================================
# STEP 3: Combine subjects (keep original IDs)
# ==============================================================================

Y_mat <- do.call(cbind, subject_curves)

# Column names are already the original sub-XXX IDs
# from names(subject_curves), no renaming needed
n_subjects <- ncol(Y_mat)

cat("=== SUBJECTS INCLUDED ===\n")
cat(paste(colnames(Y_mat), collapse = ", "), "\n")
cat("Total:", n_subjects, "\n")
cat("=========================\n\n")

# Subject metadata — records how many epochs each had
metadata_df <- data.frame(
  subject_id = names(subject_curves),
  n_epochs   = n_epochs_per_subject[seq_len(n_subjects)],
  stringsAsFactors = FALSE
)

cat("=== EPOCH COUNTS ===\n")
print(metadata_df)
cat("====================\n\n")

# ==============================================================================
# STEP 4: Export
# ==============================================================================

final_df <- data.frame(time = time_vec, Y_mat, check.names = FALSE)

cat("Final matrix:", nrow(final_df), "time points ×",
    n_subjects, "subjects\n\n")

write.csv(final_df, output_csv, row.names = FALSE)
cat("Saved:", output_csv, "\n")

write.csv(metadata_df, "./EDA/subject_metadata.csv", row.names = FALSE)
cat("Saved: ./EDA/subject_metadata.csv\n")

cat("\n*** Done! All", n_subjects, "subjects with S2 flanker data. ***\n")

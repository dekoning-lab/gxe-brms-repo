#!/usr/bin/env Rscript
############################################################
# 45_augment_imputed_data.R
# Add 'pass' and 'diff' parenting variables to existing
# imputed datasets, preserving all original imputations.
#
# pass and diff have only 1 missing value (row 156, same row
# missing all parenting variables). For this single row, we
# impute pass/diff from the conditional distribution given
# the other (already imputed) parenting variables.
#
# Input:  data/imputed_datasets_for_brms_m100.rds     (original, PRESERVED)
#         data/Fulldata_with_PCs_and_maternal_PCs.txt  (raw data)
# Output: data/imputed_datasets_for_brms_m100_v2.rds  (augmented)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

cat("=================================================================\n")
cat("AUGMENT IMPUTED DATA: add pass and diff\n")
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

# ── Load original imputed datasets ────────────────────────
orig_file <- "data/imputed_datasets_for_brms_m100.rds"
imputed_list <- readRDS(orig_file)
cat(sprintf("Loaded %d imputed datasets from %s\n", length(imputed_list), orig_file))
cat(sprintf("  n = %d, p = %d\n", nrow(imputed_list[[1]]), ncol(imputed_list[[1]])))
cat(sprintf("  Variables: %s\n\n", paste(names(imputed_list[[1]]), collapse = ", ")))

# ── Load raw data ─────────────────────────────────────────
raw_file <- "data/Fulldata_with_PCs_and_maternal_PCs.txt"
raw <- read.table(raw_file, header = TRUE, sep = "\t",
                  check.names = FALSE, encoding = "latin-1")
names(raw) <- gsub("-", ".", names(raw))
names(raw) <- make.names(names(raw), unique = TRUE)

# Filter to rows with outcome (matching imputed data)
raw <- raw[!is.na(raw$Final7PointLikertScale), ]
cat(sprintf("Raw data filtered to n = %d (with outcome)\n", nrow(raw)))

# Verify row alignment: check that the raw outcome matches imputed
stopifnot(nrow(raw) == nrow(imputed_list[[1]]))
stopifnot(all(raw$Final7PointLikertScale == imputed_list[[1]]$Final7PointLikertScale))
cat("  Row alignment verified (outcome values match)\n\n")

# ── Check missingness ────────────────────────────────────
miss_pass <- which(is.na(raw$pass))
miss_diff <- which(is.na(raw$diff))
cat(sprintf("Missing pass: %d rows (indices: %s)\n", length(miss_pass),
            paste(miss_pass, collapse = ", ")))
cat(sprintf("Missing diff: %d rows (indices: %s)\n", length(miss_diff),
            paste(miss_diff, collapse = ", ")))

# ── Scale pass and diff using raw observed values ────────
# Compute mean and sd from observed values only
pass_obs <- raw$pass[!is.na(raw$pass)]
diff_obs <- raw$diff[!is.na(raw$diff)]

pass_mean <- mean(pass_obs)
pass_sd   <- sd(pass_obs)
diff_mean <- mean(diff_obs)
diff_sd   <- sd(diff_obs)

cat(sprintf("\nScaling parameters (from observed values):\n"))
cat(sprintf("  pass: mean = %.4f, sd = %.4f\n", pass_mean, pass_sd))
cat(sprintf("  diff: mean = %.4f, sd = %.4f\n", diff_mean, diff_sd))

# ── Augment each imputed dataset ─────────────────────────
cat(sprintf("\nAugmenting %d imputed datasets...\n", length(imputed_list)))

augmented_list <- map(seq_along(imputed_list), function(m) {
  df <- imputed_list[[m]]

  # Start with observed raw values (scaled)
  df$pass <- (raw$pass - pass_mean) / pass_sd
  df$diff <- (raw$diff - diff_mean) / diff_sd

  # For the missing row(s), impute from conditional on other parenting vars
  if (length(miss_pass) > 0 || length(miss_diff) > 0) {
    # Use the already-imputed parenting variables to predict pass/diff
    # via linear regression on the observed data (donor approach)
    obs_idx <- setdiff(1:nrow(df), union(miss_pass, miss_diff))

    # Fit pass ~ sens + cont + unre on observed
    lm_pass <- lm(pass ~ sens + cont + unre, data = df[obs_idx, ])
    lm_diff <- lm(diff ~ sens + cont + unre, data = df[obs_idx, ])

    # Predict + add noise (proper imputation, not just point prediction)
    for (idx in miss_pass) {
      pred_pass <- predict(lm_pass, newdata = df[idx, , drop = FALSE])
      resid_sd <- sigma(lm_pass)
      df$pass[idx] <- pred_pass + rnorm(1, 0, resid_sd)
    }
    for (idx in miss_diff) {
      pred_diff <- predict(lm_diff, newdata = df[idx, , drop = FALSE])
      resid_sd <- sigma(lm_diff)
      df$diff[idx] <- pred_diff + rnorm(1, 0, resid_sd)
    }
  }

  return(df)
})

# Set different random seed per imputation for reproducibility
set.seed(20260308)
augmented_list <- map(seq_along(imputed_list), function(m) {
  set.seed(20260308 + m)
  df <- imputed_list[[m]]

  df$pass <- (raw$pass - pass_mean) / pass_sd
  df$diff <- (raw$diff - diff_mean) / diff_sd

  if (length(miss_pass) > 0 || length(miss_diff) > 0) {
    obs_idx <- setdiff(1:nrow(df), union(miss_pass, miss_diff))
    lm_pass <- lm(pass ~ sens + cont + unre, data = df[obs_idx, ])
    lm_diff <- lm(diff ~ sens + cont + unre, data = df[obs_idx, ])

    for (idx in miss_pass) {
      pred_pass <- predict(lm_pass, newdata = df[idx, , drop = FALSE])
      df$pass[idx] <- pred_pass + rnorm(1, 0, sigma(lm_pass))
    }
    for (idx in miss_diff) {
      pred_diff <- predict(lm_diff, newdata = df[idx, , drop = FALSE])
      df$diff[idx] <- pred_diff + rnorm(1, 0, sigma(lm_diff))
    }
  }

  return(df)
})

# ── Verify ───────────────────────────────────────────────
cat(sprintf("\nVerification:\n"))
cat(sprintf("  Original variables preserved: %s\n",
            all(sapply(seq_along(augmented_list), function(m) {
              identical(augmented_list[[m]][, names(imputed_list[[m]])],
                        imputed_list[[m]])
            }))))
cat(sprintf("  New variables added: pass, diff\n"))
cat(sprintf("  Any NA in pass: %s\n",
            any(sapply(augmented_list, function(df) any(is.na(df$pass))))))
cat(sprintf("  Any NA in diff: %s\n",
            any(sapply(augmented_list, function(df) any(is.na(df$diff))))))

# Check distributions
pass_vals <- sapply(augmented_list, function(df) mean(df$pass))
diff_vals <- sapply(augmented_list, function(df) mean(df$diff))
cat(sprintf("  Mean of pass means: %.4f (should be ~0)\n", mean(pass_vals)))
cat(sprintf("  Mean of diff means: %.4f (should be ~0)\n", mean(diff_vals)))

# ── Save augmented data ──────────────────────────────────
out_file <- "data/imputed_datasets_for_brms_m100_v2.rds"
saveRDS(augmented_list, out_file)
cat(sprintf("\nSaved augmented data: %s\n", out_file))
cat(sprintf("  Original preserved at: %s\n", orig_file))

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=== DONE ===\n")

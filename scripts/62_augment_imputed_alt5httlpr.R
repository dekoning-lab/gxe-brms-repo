#!/usr/bin/env Rscript
############################################################
# 62_augment_imputed_alt5httlpr.R
#
# ALTERNATE-CODING RUN of the augmentation step (adds 'pass' and
# 'diff' to each imputed dataset). Near-copy of
# 45_augment_imputed_data.R with paths swapped to the alt dataset.
#
# Same augmentation seed (20260308 + m) so the pass/diff sampling is
# byte-identical between original and alternate runs.
#
# Input:  data_alt_5httlpr/imputed_datasets_for_brms_m100.rds
#         data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt
# Output: data_alt_5httlpr/imputed_datasets_for_brms_m100_v2.rds
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(purrr)
})

cat("=================================================================\n")
cat("AUGMENT alt-coded IMPUTED DATA: add pass and diff\n")
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

orig_file <- "data_alt_5httlpr/imputed_datasets_for_brms_m100.rds"
imputed_list <- readRDS(orig_file)
cat(sprintf("Loaded %d imputed datasets from %s\n",
            length(imputed_list), orig_file))
cat(sprintf("  n = %d, p = %d\n", nrow(imputed_list[[1]]),
            ncol(imputed_list[[1]])))

raw_file <- "data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt"
raw <- read.table(raw_file, header = TRUE, sep = "\t",
                  check.names = FALSE, encoding = "latin-1")
names(raw) <- gsub("-", ".", names(raw))
names(raw) <- make.names(names(raw), unique = TRUE)
raw <- raw[!is.na(raw$Final7PointLikertScale), ]
cat(sprintf("Raw data filtered to n = %d (with outcome)\n", nrow(raw)))

stopifnot(nrow(raw) == nrow(imputed_list[[1]]))
stopifnot(all(raw$Final7PointLikertScale == imputed_list[[1]]$Final7PointLikertScale))
cat("  Row alignment verified.\n\n")

miss_pass <- which(is.na(raw$pass))
miss_diff <- which(is.na(raw$diff))
cat(sprintf("Missing pass: %d rows\nMissing diff: %d rows\n",
            length(miss_pass), length(miss_diff)))

pass_obs <- raw$pass[!is.na(raw$pass)]
diff_obs <- raw$diff[!is.na(raw$diff)]
pass_mean <- mean(pass_obs); pass_sd <- sd(pass_obs)
diff_mean <- mean(diff_obs); diff_sd <- sd(diff_obs)
cat(sprintf("pass: mean = %.4f sd = %.4f\n", pass_mean, pass_sd))
cat(sprintf("diff: mean = %.4f sd = %.4f\n\n", diff_mean, diff_sd))

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
      df$pass[idx] <- predict(lm_pass, newdata = df[idx, , drop = FALSE]) +
                      rnorm(1, 0, sigma(lm_pass))
    }
    for (idx in miss_diff) {
      df$diff[idx] <- predict(lm_diff, newdata = df[idx, , drop = FALSE]) +
                      rnorm(1, 0, sigma(lm_diff))
    }
  }
  df
})

cat("Verification:\n")
cat(sprintf("  Any NA in pass: %s\n",
            any(sapply(augmented_list, function(df) any(is.na(df$pass))))))
cat(sprintf("  Any NA in diff: %s\n",
            any(sapply(augmented_list, function(df) any(is.na(df$diff))))))

out_file <- "data_alt_5httlpr/imputed_datasets_for_brms_m100_v2.rds"
saveRDS(augmented_list, out_file)
cat(sprintf("\nSaved augmented data: %s\n", out_file))
cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

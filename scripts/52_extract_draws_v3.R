#!/usr/bin/env Rscript
############################################################
# 52_extract_draws_v3.R
# Extract posterior draws for key parameters from M4 sens fits
# to create actual histograms (not normal approximations).
#
# Loads all 100 imputation fits, extracts draws for:
#   - b_gene_SLC6A3.10R  (gene main effect)
#   - b_int_sens:CNR1.77  (interaction)
#   - b_int_sens:X5HTTLPR (interaction)
# Pools draws across imputations (draw-level pooling).
#
# Run from: the repository root.
# Usage: Rscript scripts/52_extract_draws_v3.R
############################################################

suppressPackageStartupMessages(library(brms))

cat("=================================================================\n")
cat("EXTRACTING POSTERIOR DRAWS FOR HISTOGRAM\n")
cat("=================================================================\n\n")

fit_dir <- "results/v3_m4/sens"
out_file <- "results/v3_figures/posterior_draws_sens.csv"

fit_files <- sort(list.files(fit_dir, pattern = "^imp_\\d+\\.rds$", full.names = TRUE))
cat(sprintf("Found %d fit files in %s\n\n", length(fit_files), fit_dir))

# Target columns (brms internal naming)
target_cols <- c("b_gene_SLC6A3.10R", "b_int_sens:CNR1.77", "b_int_sens:X5HTTLPR")
# Display names
display_names <- c("SLC6A3.10R (gene main)", "Sens x CNR1.77", "Sens x 5-HTTLPR")

all_draws <- list()

for (i in seq_along(fit_files)) {
  f <- fit_files[i]
  imp_num <- as.integer(sub(".*imp_(\\d+)\\.rds$", "\\1", basename(f)))
  cat(sprintf("  Loading imp %03d (%d/%d)... ", imp_num, i, length(fit_files)))

  fit <- readRDS(f)
  draws <- as.data.frame(fit)

  # Check column names exist
  found <- target_cols %in% names(draws)
  if (!all(found)) {
    # Try alternate naming (brms sometimes uses . instead of :)
    alt_cols <- gsub(":", ".", target_cols)
    found2 <- alt_cols %in% names(draws)
    if (all(found2)) {
      target_cols_use <- alt_cols
    } else {
      cat(sprintf("WARNING: missing columns: %s\n",
                  paste(target_cols[!found], collapse = ", ")))
      # Print available columns for debugging
      if (i == 1) {
        int_cols <- grep("b_int", names(draws), value = TRUE)
        gene_cols <- grep("b_gene", names(draws), value = TRUE)
        cat("  Available gene cols:", paste(gene_cols, collapse = ", "), "\n")
        cat("  Available int cols:", paste(int_cols, collapse = ", "), "\n")
      }
      next
    }
  } else {
    target_cols_use <- target_cols
  }

  n_draws <- nrow(draws)
  for (j in seq_along(target_cols_use)) {
    all_draws[[length(all_draws) + 1]] <- data.frame(
      parameter = display_names[j],
      value = draws[[target_cols_use[j]]],
      imputation = imp_num,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("%d draws extracted\n", n_draws))
  rm(fit, draws); gc(verbose = FALSE)
}

cat("\nCombining draws...\n")
draws_df <- do.call(rbind, all_draws)
cat(sprintf("Total: %d draws across %d parameters\n", nrow(draws_df), length(display_names)))

# Summary
for (p in display_names) {
  sub <- draws_df[draws_df$parameter == p, ]
  cat(sprintf("  %s: mean=%.3f, SD=%.3f, n=%d draws\n",
              p, mean(sub$value), sd(sub$value), nrow(sub)))
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(draws_df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s (%.1f MB)\n", out_file, file.size(out_file) / 1e6))
cat("=== DONE ===\n")

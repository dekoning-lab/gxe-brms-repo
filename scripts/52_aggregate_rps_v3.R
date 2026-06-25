#!/usr/bin/env Rscript
############################################################
# 52_aggregate_rps_v3.R
# Aggregate three-scoring-rule results for v3 models
#
# Run from: the repository root.
# Usage: Rscript scripts/52_aggregate_rps_v3.R
############################################################

cat("=================================================================\n")
cat("AGGREGATING V3 THREE-SCORING-RULE RESULTS\n")
cat("(pass + diff as child behavior covariates in all models)\n")
cat("=================================================================\n\n")

pv_list <- c("sens", "cont", "unre")
pv_labels <- c(sens = "Sensitivity", cont = "Controlling", unre = "Unresponsiveness")

# ── Load all RPS results ──────────────────────────────────
all_rps <- list()

for (pv in pv_list) {
  rps_dir <- file.path("results/v3_loo_rps", pv)
  rps_files <- list.files(rps_dir, pattern = "^rps_\\d+\\.csv$", full.names = TRUE)
  if (length(rps_files) > 0) {
    pv_data <- do.call(rbind, lapply(rps_files, read.csv, stringsAsFactors = FALSE))
    pv_data$parenting_label <- pv_labels[pv]
    all_rps[[pv]] <- pv_data
    cat(sprintf("  %s: loaded %d files (%d rows)\n", pv, length(rps_files), nrow(pv_data)))
  }
}

if (length(all_rps) == 0) {
  cat("No RPS results found.\n")
  quit(save = "no")
}

rps_df <- do.call(rbind, all_rps)
rownames(rps_df) <- NULL

write.csv(rps_df, "results/v3_loo_rps/all_rps_combined.csv", row.names = FALSE)
cat(sprintf("\nTotal: %d rows saved to results/v3_loo_rps/all_rps_combined.csv\n\n",
            nrow(rps_df)))

# ── Summary by model ─────────────────────────────────────
cat("=== MODEL SUMMARY (AVERAGED ACROSS IMPUTATIONS) ===\n\n")

cat(sprintf("%-8s %-16s %8s %8s %8s %8s %6s\n",
            "Parent.", "Model", "ELPD", "RPS", "Binary", "p_loo", "bad_k"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (pv in pv_list) {
  pv_data <- rps_df[rps_df$parenting == pv, ]
  for (m in c("M1", "M2", "M3", "M4")) {
    sub <- pv_data[pv_data$model == m, ]
    if (nrow(sub) > 0) {
      cat(sprintf("%-8s %-16s %8.1f %8.4f %8.1f %8.1f %6.1f\n",
                  pv, m, mean(sub$loo_elpd), mean(sub$mean_rps),
                  mean(sub$binary_elpd), mean(sub$p_loo), mean(sub$n_high_k)))
    }
  }
  cat("\n")
}

# ── Pairwise comparisons ─────────────────────────────────
cat("=== PAIRWISE COMPARISONS ===\n\n")

pairwise_rows <- list()

comparisons <- list(
  c("M2", "M1", "Does parenting improve on demographics?"),
  c("M3", "M1", "Do genes improve on demographics?"),
  c("M4", "M2", "Do genes + GxE improve on parenting?"),
  c("M4", "M1", "Does full model improve on demographics?")
)

for (pv in pv_list) {
  pv_data <- rps_df[rps_df$parenting == pv, ]
  cat(sprintf("--- %s ---\n", pv_labels[pv]))

  for (comp in comparisons) {
    complex <- comp[1]; simple <- comp[2]; question <- comp[3]

    base <- pv_data[pv_data$model == simple, ]
    test <- pv_data[pv_data$model == complex, ]

    shared <- intersect(base$imputation, test$imputation)
    if (length(shared) == 0) next

    base_m <- base[match(shared, base$imputation), ]
    test_m <- test[match(shared, test$imputation), ]

    d_elpd <- test_m$loo_elpd - base_m$loo_elpd
    # For RPS: lower is better, so d_rps = base - test (positive = complex better)
    d_rps  <- base_m$mean_rps - test_m$mean_rps
    d_bin  <- test_m$binary_elpd - base_m$binary_elpd

    cat(sprintf("  %s vs %s: %s\n", complex, simple, question))
    cat(sprintf("    dELPD:   %+.2f (SD %.2f) [%0.f%% positive]\n",
                mean(d_elpd), sd(d_elpd), mean(d_elpd > 0) * 100))
    cat(sprintf("    dRPS:    %+.4f (SD %.4f) [%0.f%% positive]\n",
                mean(d_rps), sd(d_rps), mean(d_rps > 0) * 100))
    cat(sprintf("    dBinary: %+.2f (SD %.2f) [%0.f%% positive]\n",
                mean(d_bin), sd(d_bin), mean(d_bin > 0) * 100))

    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      parenting = pv, comparison = paste(complex, "vs", simple),
      mean_d_elpd = mean(d_elpd), sd_d_elpd = sd(d_elpd),
      mean_d_rps = mean(d_rps), sd_d_rps = sd(d_rps),
      mean_d_binary = mean(d_bin), sd_d_binary = sd(d_bin),
      pct_complex_better_elpd = mean(d_elpd > 0) * 100,
      pct_complex_better_rps = mean(d_rps > 0) * 100,
      pct_complex_better_binary = mean(d_bin > 0) * 100,
      n_imps = length(shared),
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

if (length(pairwise_rows) > 0) {
  pw_df <- do.call(rbind, pairwise_rows)
  write.csv(pw_df, "results/v3_loo_rps/pairwise_summary.csv", row.names = FALSE)
  cat("Pairwise results saved to results/v3_loo_rps/pairwise_summary.csv\n")
}

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

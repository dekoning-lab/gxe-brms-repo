#!/usr/bin/env Rscript
############################################################
# 52_aggregate_r2_v3.R
# Aggregate Bayesian R² and McKelvey-Zavoina R² results
#
# Run from: the repository root.
# Usage: Rscript scripts/52_aggregate_r2_v3.R
############################################################

cat("=================================================================\n")
cat("V3 R² AGGREGATION\n")
cat("=================================================================\n\n")

parenting_vars  <- c("sens", "cont", "unre")
parenting_label <- c(sens = "Sensitivity", cont = "Controlling",
                     unre = "Unresponsiveness")

# ── Read all per-imputation R² files ────────────────────
all_r2 <- list()
for (pv in parenting_vars) {
  dir_path <- file.path("results", "v3_r2", pv)
  files <- list.files(dir_path, pattern = "^r2_\\d+\\.csv$", full.names = TRUE)
  cat(sprintf("  %s: %d files\n", pv, length(files)))

  for (f in files) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    all_r2[[length(all_r2) + 1]] <- d
  }
}

r2_all <- do.call(rbind, all_r2)
cat(sprintf("\nTotal: %d rows\n", nrow(r2_all)))

# ── Summary: mean across imputations ───────────────────
cat("\n=== R² SUMMARY ===\n\n")

r2_summary <- aggregate(
  cbind(bayes_r2_mean, bayes_r2_sd, bayes_r2_q025, bayes_r2_q975, bayes_r2_median,
        mz_r2_mean, mz_r2_sd, mz_r2_q025, mz_r2_q975, mz_r2_median) ~
    model + parenting,
  data = r2_all,
  FUN = mean
)

# Also compute SD across imputations for the means
r2_imp_sd <- aggregate(
  cbind(bayes_r2_mean, mz_r2_mean) ~ model + parenting,
  data = r2_all,
  FUN = sd
)
names(r2_imp_sd)[3:4] <- c("bayes_r2_imp_sd", "mz_r2_imp_sd")

r2_summary <- merge(r2_summary, r2_imp_sd, by = c("model", "parenting"))
r2_summary <- r2_summary[order(r2_summary$parenting, r2_summary$model), ]

cat(sprintf("%-8s %-5s %8s %8s %8s %8s %8s %8s\n",
            "Parent.", "Model",
            "BayR2", "BR2_SD", "BR2_025", "BR2_975",
            "MZ_R2", "MZ_025"))
cat(paste(rep("-", 72), collapse = ""), "\n")

for (i in seq_len(nrow(r2_summary))) {
  r <- r2_summary[i, ]
  cat(sprintf("%-8s %-5s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
              r$parenting, r$model,
              r$bayes_r2_mean, r$bayes_r2_sd,
              r$bayes_r2_q025, r$bayes_r2_q975,
              r$mz_r2_mean, r$mz_r2_q025))
}

# ── Save ──────────────────────────────────────────────
dir.create("results/v3_r2", recursive = TRUE, showWarnings = FALSE)
write.csv(r2_summary, "results/v3_r2/r2_summary.csv", row.names = FALSE)
write.csv(r2_all, "results/v3_r2/r2_all_imputations.csv", row.names = FALSE)
cat("\nSaved: results/v3_r2/r2_summary.csv\n")
cat("Saved: results/v3_r2/r2_all_imputations.csv\n")
cat("\n=== DONE ===\n")

#!/usr/bin/env Rscript
############################################################
# Step 54: Outlier-exclusion sensitivity analysis
#
# Re-runs the draw-level pooling used to produce
# results/v3_m4/<parenting>/summary.csv (see 52_fit_m4_v3.R
# lines 234–298), but EXCLUDING the catastrophic-pathology
# imputation for the requested parenting model:
#
#   sens: exclude imp 3
#   cont: exclude imp 1
#   unre: exclude imp 38
#
# These are the three imputed datasets identified in
# scripts/53_diagnostics_v3_m4.R as having severe sampling
# pathology (Rhat ~ 1.55–1.59, ~3000–4000 divergences,
# min ESS = 7, BFMI < 0.27).
#
# The output is a parallel summary.csv computed from the
# remaining 99 fits, with identical structure. A separate
# join script (or the manuscript update) produces the
# original-vs-excluded comparison.
#
# Run from: the repository root (locally, or on a cluster after sourcing config.sh; see README).
# Usage:    Rscript scripts/54_sensitivity_excl_outlier.R <parenting>
#
# Input:  results/v3_m4/<parenting>/imp_*.rds  (skips one)
# Output: results/v3_m4_excl_outlier/<parenting>/summary.csv
#         results/v3_m4_excl_outlier/<parenting>/excluded_imp.txt
############################################################

suppressPackageStartupMessages({
  library(brms)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/54_sensitivity_excl_outlier.R <parenting>")
}
parenting_var <- args[1]
stopifnot(parenting_var %in% c("sens", "cont", "unre"))

# Catastrophic-pathology imputation per parenting model
# (sourced from 53_diagnostics_v3_m4.R / Table S6 of manuscript)
excluded_imp <- switch(parenting_var,
  sens = 3L,
  cont = 1L,
  unre = 38L)

cat("=================================================================\n")
cat(sprintf("OUTLIER-EXCLUSION SENSITIVITY: %s (excluding imp %d)\n",
            parenting_var, excluded_imp))
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

fit_dir <- file.path("results", "v3_m4", parenting_var)
out_dir <- file.path("results", "v3_m4_excl_outlier", parenting_var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_files <- sort(list.files(fit_dir, pattern = "^imp_\\d+\\.rds$",
                             full.names = TRUE))
stopifnot(length(all_files) >= 99)

# Identify and remove the excluded imputation file
excl_basename <- sprintf("imp_%03d.rds", excluded_imp)
excluded_path <- file.path(fit_dir, excl_basename)
stopifnot(excluded_path %in% all_files)
kept_files <- all_files[all_files != excluded_path]
cat(sprintf("Total fits: %d; excluded: %s; kept: %d\n",
            length(all_files), excl_basename, length(kept_files)))

# Write a small record of what we did
writeLines(c(
  sprintf("parenting_var = %s", parenting_var),
  sprintf("excluded_imp = %d", excluded_imp),
  sprintf("excluded_file = %s", excluded_path),
  sprintf("kept_n = %d", length(kept_files)),
  sprintf("rationale = catastrophic sampling pathology (see Table S6 of manuscript)"),
  sprintf("script = scripts/54_sensitivity_excl_outlier.R"),
  sprintf("started = %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
), con = file.path(out_dir, "excluded_imp.txt"))

# ── Pool draws across the kept imputations ─────────────────
# Use the same procedure as the production aggregator
# (52_fit_m4_v3.R lines 240–245).
all_draws <- data.frame()
t0 <- Sys.time()
for (i in seq_along(kept_files)) {
  f <- kept_files[i]
  imp <- as.integer(sub(".*imp_(\\d+)\\.rds", "\\1", basename(f)))
  fit <- readRDS(f)
  draws <- as.data.frame(fit)
  all_draws <- rbind(all_draws, draws)
  rm(fit); gc(verbose = FALSE)
  if (i %% 10 == 0 || i == length(kept_files)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  loaded %d/%d fits (imp %d); pooled draws so far: %d  [%.0fs]\n",
                i, length(kept_files), imp, nrow(all_draws), el))
  }
}

cat(sprintf("\nPooled draws total: %d rows × %d columns\n",
            nrow(all_draws), ncol(all_draws)))

# ── Compute pooled summary (identical structure to summary.csv) ─
cov_cols  <- grep("^b_cov_", names(all_draws), value = TRUE)
gene_cols <- grep("^b_gene_", names(all_draws), value = TRUE)
int_cols  <- grep("^b_int_", names(all_draws), value = TRUE)
fixef_cols <- c(cov_cols, gene_cols, int_cols)

term_names <- gsub("^b_(cov|gene|int)_", "", fixef_cols)
term_names <- sapply(term_names, function(t) {
  if (grepl(":", t)) {
    parts <- strsplit(t, ":")[[1]]
    if (parts[2] == parenting_var) paste(parts[2], parts[1], sep = ":") else t
  } else t
}, USE.NAMES = FALSE)

effect_types <- ifelse(fixef_cols %in% cov_cols, "covariate",
                ifelse(fixef_cols %in% gene_cols, "gene_main",
                ifelse(fixef_cols %in% int_cols, "interaction", "unknown")))

results <- data.frame(
  parenting = parenting_var,
  term = term_names,
  estimate = sapply(fixef_cols, function(col) mean(all_draws[[col]])),
  sd = sapply(fixef_cols, function(col) sd(all_draws[[col]])),
  q025 = sapply(fixef_cols, function(col) quantile(all_draws[[col]], 0.025)),
  q975 = sapply(fixef_cols, function(col) quantile(all_draws[[col]], 0.975)),
  post_prob_pos = sapply(fixef_cols, function(col) mean(all_draws[[col]] > 0)),
  post_prob_neg = sapply(fixef_cols, function(col) mean(all_draws[[col]] < 0)),
  stringsAsFactors = FALSE, row.names = NULL
)
results$post_prob <- pmax(results$post_prob_pos, results$post_prob_neg)
results$direction <- ifelse(results$post_prob_pos > results$post_prob_neg,
                            "positive", "negative")
results$effect_type <- effect_types

# Hyperparameters
for (hp in c("sigma_gene", "sigma_int")) {
  if (hp %in% names(all_draws)) {
    results <- rbind(results, data.frame(
      parenting = parenting_var, term = hp,
      estimate = mean(all_draws[[hp]]), sd = sd(all_draws[[hp]]),
      q025 = quantile(all_draws[[hp]], 0.025),
      q975 = quantile(all_draws[[hp]], 0.975),
      post_prob_pos = NA, post_prob_neg = NA,
      post_prob = NA, direction = NA,
      effect_type = "hyperparameter",
      stringsAsFactors = FALSE, row.names = NULL
    ))
  }
}

out_csv <- file.path(out_dir, "summary.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)\n", out_csv, nrow(results)))
cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

#!/usr/bin/env Rscript
############################################################
# 52_compute_r2_v3.R
# Compute Bayesian R² and McKelvey-Zavoina R² for v3 models
#
# For each imputation, loads M1, M2, M3, M4 and computes:
#   1. Bayesian R² (Gelman et al., 2019) via bayes_R2()
#   2. McKelvey-Zavoina R² = Var(Xβ) / [Var(Xβ) + π²/3]
#
# Args: <parenting> <start_imp> <end_imp>
# Run from: the repository root (locally or on a cluster).
############################################################

suppressPackageStartupMessages({
  library(brms)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 52_compute_r2_v3.R <parenting> [start] [end]")

parenting_var <- args[1]
imp_start <- if (length(args) >= 2) as.integer(args[2]) else 1L
imp_end   <- if (length(args) >= 3) as.integer(args[3]) else 100L

stopifnot(parenting_var %in% c("sens", "cont", "unre"))

cat(sprintf("=== V3 R² COMPUTATION: %s, imputations %d-%d ===\n\n",
            parenting_var, imp_start, imp_end))

out_dir <- sprintf("results/v3_r2/%s", parenting_var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Logistic residual variance (for McKelvey-Zavoina)
logistic_var <- pi^2 / 3

# ── Function to compute R² measures from a brms fit ──────
compute_r2 <- function(fit, model_name) {
  # 1. Bayesian R² (posterior distribution)
  br2 <- tryCatch({
    r2_draws <- bayes_R2(fit, summary = FALSE)  # [S draws]
    data.frame(
      mean  = mean(r2_draws),
      sd    = sd(r2_draws),
      q025  = quantile(r2_draws, 0.025),
      q975  = quantile(r2_draws, 0.975),
      median = median(r2_draws)
    )
  }, error = function(e) {
    cat(sprintf("    bayes_R2 error for %s: %s\n", model_name, e$message))
    data.frame(mean = NA, sd = NA, q025 = NA, q975 = NA, median = NA)
  })

  # 2. McKelvey-Zavoina R²
  # Compute Var(linear predictor) across observations for each draw
  mz_r2 <- tryCatch({
    # posterior_linpred gives [S draws × N observations]
    linpred <- posterior_linpred(fit)  # matrix: S × N
    # For each draw, compute variance across observations
    var_xb_per_draw <- apply(linpred, 1, var)
    # R²_MZ = Var(Xβ) / [Var(Xβ) + π²/3]
    mz_draws <- var_xb_per_draw / (var_xb_per_draw + logistic_var)
    data.frame(
      mean  = mean(mz_draws),
      sd    = sd(mz_draws),
      q025  = quantile(mz_draws, 0.025),
      q975  = quantile(mz_draws, 0.975),
      median = median(mz_draws)
    )
  }, error = function(e) {
    cat(sprintf("    MZ R² error for %s: %s\n", model_name, e$message))
    data.frame(mean = NA, sd = NA, q025 = NA, q975 = NA, median = NA)
  })

  data.frame(
    model = model_name,
    bayes_r2_mean = br2$mean, bayes_r2_sd = br2$sd,
    bayes_r2_q025 = br2$q025, bayes_r2_q975 = br2$q975,
    bayes_r2_median = br2$median,
    mz_r2_mean = mz_r2$mean, mz_r2_sd = mz_r2$sd,
    mz_r2_q025 = mz_r2$q025, mz_r2_q975 = mz_r2$q975,
    mz_r2_median = mz_r2$median,
    stringsAsFactors = FALSE
  )
}

# ── Loop over imputations ──────────────────────────────────
for (imp in imp_start:imp_end) {
  out_file <- file.path(out_dir, sprintf("r2_%03d.csv", imp))
  if (file.exists(out_file)) {
    cat(sprintf("  Imputation %d: already done, skipping.\n", imp))
    next
  }

  cat(sprintf("\n--- Imputation %d ---\n", imp))
  imp_str <- sprintf("%03d", imp)

  m1_file <- sprintf("results/v3_m1/imp_%s.rds", imp_str)
  m2_file <- sprintf("results/v3_m2/%s/imp_%s.rds", parenting_var, imp_str)
  m3_file <- sprintf("results/v3_m3/imp_%s.rds", imp_str)
  m4_file <- sprintf("results/v3_m4/%s/imp_%s.rds", parenting_var, imp_str)

  # Check all files exist
  missing <- character()
  for (f in c(m1_file, m2_file, m3_file, m4_file)) {
    if (!file.exists(f)) missing <- c(missing, f)
  }
  if (length(missing) > 0) {
    cat(sprintf("  MISSING: %s, skipping.\n", paste(missing, collapse = ", ")))
    next
  }

  results <- list()

  for (spec in list(
    list(file = m1_file, name = "M1"),
    list(file = m2_file, name = "M2"),
    list(file = m3_file, name = "M3"),
    list(file = m4_file, name = "M4")
  )) {
    cat(sprintf("  Loading & computing R² for %s...\n", spec$name))
    t0 <- Sys.time()
    fit <- readRDS(spec$file)
    r2 <- compute_r2(fit, spec$name)
    rm(fit); gc(verbose = FALSE)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("    Bayes R² = %.3f [%.3f, %.3f], MZ R² = %.3f [%.3f, %.3f] (%.0fs)\n",
                r2$bayes_r2_mean, r2$bayes_r2_q025, r2$bayes_r2_q975,
                r2$mz_r2_mean, r2$mz_r2_q025, r2$mz_r2_q975, elapsed))
    results[[length(results) + 1]] <- r2
  }

  result_df <- do.call(rbind, results)
  result_df$parenting <- parenting_var
  result_df$imputation <- imp
  write.csv(result_df, out_file, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\n=== DONE ===\n")

#!/usr/bin/env Rscript
############################################################
# 52_loo_rps_v3.R
# Three-scoring-rule LOO comparison for v3 models
# (pass + diff as child behavior covariates in all models)
#
# Computes ELPD, RPS, and binary ELPD for M1, M2, M3, M4
# by loading pre-fitted model objects.
#
# Args: <parenting> <start_imp> <end_imp>
# Run from: the repository root (locally or on a cluster).
# Usage:
#   Rscript scripts/52_loo_rps_v3.R sens 1 1
#   Rscript scripts/52_loo_rps_v3.R cont 50 50
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 52_loo_rps_v3.R <parenting> [start] [end]")

parenting_var <- args[1]
stopifnot(parenting_var %in% c("sens", "cont", "unre"))

imp_start <- if (length(args) >= 2) as.integer(args[2]) else 1L
imp_end   <- if (length(args) >= 3) as.integer(args[3]) else 100L

cat(sprintf("=== V3 LOO-RPS: %s, imputations %d-%d ===\n\n",
            parenting_var, imp_start, imp_end))

# Output directory
out_dir <- sprintf("results/v3_loo_rps/%s", parenting_var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load imputed datasets (v2 with pass+diff)
cat("Loading imputed data...\n")
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")

# ── Scoring function ──────────────────────────────────────
compute_scores <- function(fit, y_num, K, model_name) {
  n <- length(y_num)

  # Posterior predicted category probabilities [draws, obs, K]
  pp <- posterior_epred(fit)

  # LOO with saved PSIS object
  loo_obj <- loo(fit, save_psis = TRUE, cores = 1)

  elpd     <- loo_obj$estimates["elpd_loo", "Estimate"]
  elpd_se  <- loo_obj$estimates["elpd_loo", "SE"]
  p_loo_val <- loo_obj$estimates["p_loo", "Estimate"]
  n_high_k <- sum(loo_obj$diagnostics$pareto_k > 0.7)
  max_k    <- max(loo_obj$diagnostics$pareto_k)

  # Normalized PSIS weights [S, n]
  lw <- weights(loo_obj$psis_object, normalize = TRUE, log = FALSE)

  # Compute RPS and binary LOO per observation
  rps <- numeric(n)
  binary_ll <- numeric(n)
  threshold <- 4L  # Category 4+ = "secure"

  for (i in 1:n) {
    p_loo <- as.numeric(crossprod(lw[, i], pp[, i, ]))
    # RPS: squared difference of cumulative distributions
    F_loo <- cumsum(p_loo)
    obs_cum <- as.numeric(y_num[i] <= seq_len(K))
    rps[i] <- sum((F_loo - obs_cum)^2) / (K - 1)
    # Binary threshold: P(Y >= threshold)
    p_ge <- sum(p_loo[threshold:K])
    p_ge <- max(min(p_ge, 1 - 1e-10), 1e-10)
    binary_ll[i] <- if (y_num[i] >= threshold) log(p_ge) else log(1 - p_ge)
  }

  data.frame(
    model       = model_name,
    loo_elpd    = elpd,
    loo_elpd_se = elpd_se,
    p_loo       = p_loo_val,
    mean_rps    = mean(rps),
    total_rps   = sum(rps),
    binary_elpd = sum(binary_ll),
    n_high_k    = n_high_k,
    max_k       = max_k,
    stringsAsFactors = FALSE
  )
}

# ── Loop over imputations ──────────────────────────────────
for (imp in imp_start:imp_end) {
  out_file <- file.path(out_dir, sprintf("rps_%03d.csv", imp))
  if (file.exists(out_file)) {
    cat(sprintf("  Imputation %d: already done, skipping.\n", imp))
    next
  }

  cat(sprintf("\n--- Imputation %d ---\n", imp))

  # Get outcome vector
  dat_raw <- imputed_data[[imp]]
  dat <- dat_raw
  if (!"Final7PointLikertScaleFactor" %in% names(dat)) {
    dat$Final7PointLikertScaleFactor <- factor(dat$Final7PointLikertScale)
  }
  dat <- dat[!is.na(dat$Final7PointLikertScaleFactor), ]
  y_num <- as.numeric(dat$Final7PointLikertScaleFactor)
  K <- max(y_num)

  # ── Load pre-fitted models ──
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

  # Load and compute scores for each model
  cat("  Loading & scoring M1...\n")
  t0 <- Sys.time()
  fit_m1 <- readRDS(m1_file)
  sc_m1 <- compute_scores(fit_m1, y_num, K, "M1")
  rm(fit_m1); gc(verbose = FALSE)
  cat(sprintf("    done (%.1f min)\n", difftime(Sys.time(), t0, units = "mins")))

  cat("  Loading & scoring M2...\n")
  t0 <- Sys.time()
  fit_m2 <- readRDS(m2_file)
  sc_m2 <- compute_scores(fit_m2, y_num, K, "M2")
  rm(fit_m2); gc(verbose = FALSE)
  cat(sprintf("    done (%.1f min)\n", difftime(Sys.time(), t0, units = "mins")))

  cat("  Loading & scoring M3...\n")
  t0 <- Sys.time()
  fit_m3 <- readRDS(m3_file)
  sc_m3 <- compute_scores(fit_m3, y_num, K, "M3")
  rm(fit_m3); gc(verbose = FALSE)
  cat(sprintf("    done (%.1f min)\n", difftime(Sys.time(), t0, units = "mins")))

  cat("  Loading & scoring M4...\n")
  t0 <- Sys.time()
  fit_m4 <- readRDS(m4_file)
  sc_m4 <- compute_scores(fit_m4, y_num, K, "M4")
  rm(fit_m4); gc(verbose = FALSE)
  cat(sprintf("    done (%.1f min)\n", difftime(Sys.time(), t0, units = "mins")))

  # Combine and save
  result <- rbind(sc_m1, sc_m2, sc_m3, sc_m4)
  result$parenting <- parenting_var
  result$imputation <- imp
  write.csv(result, out_file, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", out_file))

  # Print summary
  cat(sprintf("  ELPD:   M1=%.1f  M2=%.1f  M3=%.1f  M4=%.1f\n",
              sc_m1$loo_elpd, sc_m2$loo_elpd, sc_m3$loo_elpd, sc_m4$loo_elpd))
  cat(sprintf("  RPS:    M1=%.4f  M2=%.4f  M3=%.4f  M4=%.4f\n",
              sc_m1$mean_rps, sc_m2$mean_rps, sc_m3$mean_rps, sc_m4$mean_rps))
  cat(sprintf("  Binary: M1=%.1f  M2=%.1f  M3=%.1f  M4=%.1f\n",
              sc_m1$binary_elpd, sc_m2$binary_elpd, sc_m3$binary_elpd, sc_m4$binary_elpd))
}

cat("\n=== DONE ===\n")

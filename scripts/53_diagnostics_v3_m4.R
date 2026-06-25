#!/usr/bin/env Rscript
############################################################
# Step 53: MCMC Diagnostics for v3 M4 (Production) Fits
#
# Extracts comprehensive MCMC diagnostics from each fitted
# production v3 M4 brms model: ESS (bulk & tail), Rhat,
# divergent transitions, max treedepth hits, and BFMI.
#
# Adapted from 30_diagnostics_hierarchical.R to read the
# production v3 fit objects (with pass + diff covariates)
# saved at results/v3_m4/<parenting>/imp_NNN.rds.
#
# Run from: the repository root.
# Usage:    Rscript scripts/53_diagnostics_v3_m4.R <parenting>
#
# Input:  results/v3_m4/<parenting>/imp_*.rds
# Output: results/diagnostics_v3/<parenting>/diagnostics_per_fit.csv
#         results/diagnostics_v3/<parenting>/diagnostics_per_param.csv
#         results/diagnostics_v3/<parenting>/diagnostics_param_summary.csv
#         results/diagnostics_v3/<parenting>/diagnostics_summary.txt
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/53_diagnostics_v3_m4.R <parenting>")
}

parenting_var <- args[1]
stopifnot(parenting_var %in% c("sens", "cont", "unre"))

cat(sprintf("=================================================================\n"))
cat(sprintf("MCMC DIAGNOSTICS (v3 M4 PRODUCTION): %s\n", parenting_var))
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("=================================================================\n\n"))

fit_dir <- file.path("results", "v3_m4", parenting_var)
out_dir <- file.path("results", "diagnostics_v3", parenting_var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Find fit files
fit_files <- sort(list.files(fit_dir, pattern = "^imp_\\d+\\.rds$",
                             full.names = TRUE))
fit_files <- fit_files[!grepl("test", fit_files)]

if (length(fit_files) == 0) {
  stop(sprintf("No fit files found in %s", fit_dir))
}

cat(sprintf("Found %d fit files in %s\n\n", length(fit_files), fit_dir))

# Open summary text file
summary_file <- file.path(out_dir, "diagnostics_summary.txt")
sink(summary_file, split = TRUE)

# Storage
fit_rows <- list()
param_rows <- list()

for (i in seq_along(fit_files)) {
  f <- fit_files[i]
  imp_num <- as.integer(sub(".*imp_(\\d+)\\.rds$", "\\1", basename(f)))

  cat(sprintf("[%d/%d] imp_%03d: ", i, length(fit_files), imp_num))
  t0 <- Sys.time()

  fit <- readRDS(f)

  # --- NUTS diagnostics ---
  np <- nuts_params(fit)

  # Divergences
  n_divergent <- sum(np$Value[np$Parameter == "divergent__"])

  # Max treedepth hits (threshold = 15 as set in 52_fit_m4_v3.R)
  n_max_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= 15)

  # BFMI per chain
  energy_df <- np[np$Parameter == "energy__", ]
  bfmi_vals <- tapply(energy_df$Value, energy_df$Chain, function(e) {
    if (length(e) < 2) return(NA_real_)
    var(diff(e)) / var(e)
  })
  min_bfmi <- min(bfmi_vals, na.rm = TRUE)
  mean_bfmi <- mean(bfmi_vals, na.rm = TRUE)
  n_low_bfmi <- sum(bfmi_vals < 0.3, na.rm = TRUE)

  # --- Per-parameter ESS and Rhat ---
  draws <- as_draws_df(fit)

  meta_cols <- c(".chain", ".iteration", ".draw", "lp__", "lprior")
  param_cols <- setdiff(names(draws), meta_cols)

  param_summary <- summarise_draws(
    draws[, c(param_cols, ".chain", ".iteration", ".draw")],
    ess_bulk, ess_tail, rhat
  )

  for (j in seq_len(nrow(param_summary))) {
    param_rows[[length(param_rows) + 1]] <- data.frame(
      parenting = parenting_var,
      imputation = imp_num,
      parameter = param_summary$variable[j],
      ess_bulk = param_summary$ess_bulk[j],
      ess_tail = param_summary$ess_tail[j],
      rhat = param_summary$rhat[j],
      stringsAsFactors = FALSE
    )
  }

  # Per-fit summary stats
  min_ess_bulk <- min(param_summary$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(param_summary$ess_tail, na.rm = TRUE)
  median_ess_bulk <- median(param_summary$ess_bulk, na.rm = TRUE)
  max_rhat <- max(param_summary$rhat, na.rm = TRUE)
  worst_rhat_param <- param_summary$variable[which.max(param_summary$rhat)]
  n_rhat_above <- sum(param_summary$rhat > 1.01, na.rm = TRUE)

  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    parenting = parenting_var,
    imputation = imp_num,
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth,
    min_bfmi = round(min_bfmi, 4),
    mean_bfmi = round(mean_bfmi, 4),
    n_low_bfmi = n_low_bfmi,
    min_ess_bulk = round(min_ess_bulk, 0),
    min_ess_tail = round(min_ess_tail, 0),
    median_ess_bulk = round(median_ess_bulk, 0),
    max_rhat = round(max_rhat, 4),
    worst_rhat_param = worst_rhat_param,
    n_rhat_above_1_01 = n_rhat_above,
    stringsAsFactors = FALSE
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("div=%d td=%d bfmi=%.3f minESS=%.0f maxRhat=%.4f [%.0fs]\n",
              n_divergent, n_max_treedepth, min_bfmi,
              min_ess_bulk, max_rhat, elapsed))
}

# ── Save per-fit diagnostics ──────────────────────────────────────────
fit_df <- do.call(rbind, fit_rows)
fit_out <- file.path(out_dir, "diagnostics_per_fit.csv")
write.csv(fit_df, fit_out, row.names = FALSE)
cat(sprintf("\nSaved per-fit diagnostics: %s (%d fits)\n", fit_out, nrow(fit_df)))

# ── Save per-parameter diagnostics ────────────────────────────────────
param_df <- do.call(rbind, param_rows)
param_out <- file.path(out_dir, "diagnostics_per_param.csv")
write.csv(param_df, param_out, row.names = FALSE)
cat(sprintf("Saved per-parameter diagnostics: %s (%d rows)\n\n",
            param_out, nrow(param_df)))

# ── Overall summary ──────────────────────────────────────────────────
cat("=================================================================\n")
cat(sprintf("DIAGNOSTICS SUMMARY: %s (%d fits)\n", parenting_var, nrow(fit_df)))
cat("=================================================================\n\n")

cat("DIVERGENT TRANSITIONS:\n")
cat(sprintf("  Total: %d across %d fits\n", sum(fit_df$n_divergent), nrow(fit_df)))
cat(sprintf("  Mean per fit: %.1f\n", mean(fit_df$n_divergent)))
cat(sprintf("  Median per fit: %.0f\n", median(fit_df$n_divergent)))
cat(sprintf("  Max in single fit: %d (imp %d)\n",
            max(fit_df$n_divergent),
            fit_df$imputation[which.max(fit_df$n_divergent)]))
cat(sprintf("  Fits with 0 divergences: %d/%d (%.0f%%)\n",
            sum(fit_df$n_divergent == 0), nrow(fit_df),
            100 * mean(fit_df$n_divergent == 0)))
cat(sprintf("  Fits with >50 divergences: %d/%d\n\n",
            sum(fit_df$n_divergent > 50), nrow(fit_df)))

cat("MAX TREEDEPTH:\n")
cat(sprintf("  Total hits: %d across %d fits\n",
            sum(fit_df$n_max_treedepth), nrow(fit_df)))
cat(sprintf("  Fits with treedepth hits: %d/%d\n\n",
            sum(fit_df$n_max_treedepth > 0), nrow(fit_df)))

cat("BFMI (Bayesian Fraction of Missing Information):\n")
cat(sprintf("  Mean across fits: %.3f\n", mean(fit_df$mean_bfmi)))
cat(sprintf("  Min across fits: %.3f (imp %d)\n",
            min(fit_df$min_bfmi),
            fit_df$imputation[which.min(fit_df$min_bfmi)]))
cat(sprintf("  Fits with any low-BFMI chain (<0.3): %d/%d\n\n",
            sum(fit_df$n_low_bfmi > 0), nrow(fit_df)))

cat("EFFECTIVE SAMPLE SIZE (bulk):\n")
cat(sprintf("  Median of min-ESS across fits: %.0f\n", median(fit_df$min_ess_bulk)))
cat(sprintf("  Worst min-ESS: %.0f (imp %d)\n",
            min(fit_df$min_ess_bulk),
            fit_df$imputation[which.min(fit_df$min_ess_bulk)]))
cat(sprintf("  Median of median-ESS: %.0f\n", median(fit_df$median_ess_bulk)))
cat(sprintf("  Fits with min-ESS < 400: %d/%d\n\n",
            sum(fit_df$min_ess_bulk < 400), nrow(fit_df)))

cat("EFFECTIVE SAMPLE SIZE (tail):\n")
cat(sprintf("  Median of min-tail-ESS: %.0f\n", median(fit_df$min_ess_tail)))
cat(sprintf("  Worst min-tail-ESS: %.0f (imp %d)\n",
            min(fit_df$min_ess_tail),
            fit_df$imputation[which.min(fit_df$min_ess_tail)]))
cat(sprintf("  Fits with min-tail-ESS < 400: %d/%d\n\n",
            sum(fit_df$min_ess_tail < 400), nrow(fit_df)))

cat("RHAT CONVERGENCE:\n")
cat(sprintf("  Mean of max-Rhat across fits: %.4f\n", mean(fit_df$max_rhat)))
cat(sprintf("  Worst max-Rhat: %.4f (imp %d)\n",
            max(fit_df$max_rhat),
            fit_df$imputation[which.max(fit_df$max_rhat)]))
cat(sprintf("  Fits with any Rhat > 1.01: %d/%d\n",
            sum(fit_df$n_rhat_above_1_01 > 0), nrow(fit_df)))
cat(sprintf("  Fits with any Rhat > 1.05: %d/%d\n\n",
            sum(fit_df$max_rhat > 1.05), nrow(fit_df)))

# ── Per-parameter summary across imputations ─────────────────────────
cat("=================================================================\n")
cat(sprintf("PER-PARAMETER SUMMARY (across %d imputations)\n", nrow(fit_df)))
cat("=================================================================\n\n")

# Classify parameters
param_df$param_type <- ifelse(
  param_df$parameter %in% c("sigma_gene", "sigma_int"), "hyperparameter",
  ifelse(grepl("^b_cov_", param_df$parameter), "covariate",
  ifelse(grepl("^b_gene_", param_df$parameter), "gene_main",
  ifelse(grepl("^b_int_", param_df$parameter), "interaction",
  ifelse(grepl("^b_Intercept", param_df$parameter), "threshold", "other")))))

param_agg <- do.call(rbind, lapply(split(param_df, param_df$parameter), function(d) {
  data.frame(
    parameter = d$parameter[1],
    param_type = d$param_type[1],
    n_imps = nrow(d),
    median_ess_bulk = median(d$ess_bulk, na.rm = TRUE),
    min_ess_bulk = min(d$ess_bulk, na.rm = TRUE),
    q10_ess_bulk = quantile(d$ess_bulk, 0.10, na.rm = TRUE),
    median_ess_tail = median(d$ess_tail, na.rm = TRUE),
    min_ess_tail = min(d$ess_tail, na.rm = TRUE),
    median_rhat = median(d$rhat, na.rm = TRUE),
    max_rhat = max(d$rhat, na.rm = TRUE),
    n_rhat_above_1_01 = sum(d$rhat > 1.01, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
row.names(param_agg) <- NULL

for (ptype in c("hyperparameter", "covariate", "gene_main",
                "interaction", "threshold")) {
  pdata <- param_agg[param_agg$param_type == ptype, ]
  if (nrow(pdata) == 0) next

  ptype_label <- switch(ptype,
    hyperparameter = "HYPERPARAMETERS",
    covariate = "COVARIATES",
    gene_main = "GENE MAIN EFFECTS",
    interaction = "INTERACTIONS",
    threshold = "THRESHOLDS"
  )
  cat(sprintf("\n%s:\n", ptype_label))
  cat(sprintf("  %-25s %8s %8s %8s %8s %8s %8s\n",
              "Parameter", "MedESS", "MinESS", "MedTail", "MinTail",
              "MedRhat", "MaxRhat"))
  cat(sprintf("  %s\n", paste(rep("-", 83), collapse = "")))

  pdata <- pdata[order(pdata$min_ess_bulk), ]

  for (j in seq_len(nrow(pdata))) {
    r <- pdata[j, ]
    flag <- ""
    if (r$min_ess_bulk < 400) flag <- paste0(flag, " *LOW_ESS*")
    if (r$max_rhat > 1.01) flag <- paste0(flag, " *HIGH_RHAT*")
    cat(sprintf("  %-25s %8.0f %8.0f %8.0f %8.0f %8.4f %8.4f%s\n",
                r$parameter, r$median_ess_bulk, r$min_ess_bulk,
                r$median_ess_tail, r$min_ess_tail,
                r$median_rhat, r$max_rhat, flag))
  }
}

param_agg_out <- file.path(out_dir, "diagnostics_param_summary.csv")
write.csv(param_agg, param_agg_out, row.names = FALSE)
cat(sprintf("\n\nSaved parameter summary: %s\n", param_agg_out))

# ── Publication-ready summary ─────────────────────────────────────────
cat("\n=================================================================\n")
cat("PUBLICATION SUMMARY\n")
cat("=================================================================\n\n")

# v3 M4 fits use 8000 iter, 4000 warmup => 4000 post-warmup per chain
post_warmup_per_chain <- 4000L
chains <- 4L
draws_per_fit <- post_warmup_per_chain * chains

cat(sprintf("Model: %s (v3 M4 production, hierarchical partial pooling)\n", parenting_var))
cat(sprintf("Imputations analyzed: %d\n", nrow(fit_df)))
cat(sprintf("Chains per fit: %d, Post-warmup draws per chain: %d\n",
            chains, post_warmup_per_chain))
cat(sprintf("Total posterior draws per imputation: %d\n", draws_per_fit))
cat(sprintf("adapt_delta: 0.99, max_treedepth: 15, init: 0.5, QR decomposition: enabled\n\n"))

total_draws <- nrow(fit_df) * draws_per_fit
cat(sprintf("Across %d fits (%d total posterior draws):\n", nrow(fit_df), total_draws))
cat(sprintf("  Divergent transitions: %d/%d (%.3f%%)\n",
            sum(fit_df$n_divergent),
            total_draws,
            100 * sum(fit_df$n_divergent) / total_draws))
cat(sprintf("  Max treedepth hits: %d\n", sum(fit_df$n_max_treedepth)))
cat(sprintf("  Low BFMI chains: %d (across %d fits)\n",
            sum(fit_df$n_low_bfmi), sum(fit_df$n_low_bfmi > 0)))
cat(sprintf("  Min bulk ESS (any param, any imp): %.0f\n", min(fit_df$min_ess_bulk)))
cat(sprintf("  Min tail ESS (any param, any imp): %.0f\n", min(fit_df$min_ess_tail)))
cat(sprintf("  Max Rhat (any param, any imp): %.4f\n", max(fit_df$max_rhat)))

problem_fits <- fit_df[fit_df$n_divergent > 100 |
                       fit_df$min_ess_bulk < 200 |
                       fit_df$max_rhat > 1.05, ]
if (nrow(problem_fits) > 0) {
  cat(sprintf("\nPROBLEMATIC FITS (%d):\n", nrow(problem_fits)))
  for (j in seq_len(nrow(problem_fits))) {
    r <- problem_fits[j, ]
    cat(sprintf("  imp_%03d: div=%d minESS=%.0f maxRhat=%.4f bfmi=%.3f\n",
                r$imputation, r$n_divergent, r$min_ess_bulk,
                r$max_rhat, r$min_bfmi))
  }
} else {
  cat("\nNo severely problematic fits detected.\n")
}

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

sink()
cat(sprintf("Diagnostics summary saved to: %s\n", summary_file))

#!/usr/bin/env Rscript
############################################################
# 52_fit_m1_v3.R
# M1 (v3): Demographics + child behavior indices
#
# Model:
#   Attachment ~ PC1 + PC2 + PC3 + Sex + infant_age + pass + diff
#
# Key change from 41_fit_demographics_only.R:
#   - Adds passivity and difficultness as child behavior
#     covariates (they are child indices from CARE, not
#     parenting variables)
#   - Uses augmented v2 imputed data (has pass/diff)
#   - Parenting-independent: 100 fits serve all comparisons
#
# Run from: the repository root (locally, or on a cluster after sourcing config.sh; see README).
# Usage:    Rscript scripts/52_fit_m1_v3.R [start] [end]
#
# Output:   results/v3_m1/imp_NNN.rds
#           results/v3_m1/loo_NNN.rds
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
imp_start <- if (length(args) >= 1) as.integer(args[1]) else 1L
imp_end   <- if (length(args) >= 2) as.integer(args[2]) else 100L

cat("=================================================================\n")
cat(sprintf("M1 (v3): DEMOGRAPHICS + CHILD BEHAVIOR: imp %d-%d\n", imp_start, imp_end))
cat("  Formula: Attachment ~ PC1 + PC2 + PC3 + Sex + infant_age + pass + diff\n")
cat("  Priors:  Normal(0, 2) on all coefficients\n")
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

# ── Load augmented imputed data (v2, has pass and diff) ─────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets (n=%d each)\n",
            length(imputed_data), nrow(imputed_data[[1]])))

# Verify pass and diff are present
stopifnot("pass" %in% names(imputed_data[[1]]))
stopifnot("diff" %in% names(imputed_data[[1]]))
cat("  Confirmed: pass and diff present in data\n\n")

# ── Output directory ─────────────────────────────────────
out_dir <- "results/v3_m1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Model formula ────────────────────────────────────────
bf_formula <- bf(Final7PointLikertScaleFactor ~ PC1 + PC2 + PC3 + Sex + infant_age + pass + diff)

prior_spec <- c(
  prior(normal(0, 2), class = "b")
)

# ── Fit + LOO loop ───────────────────────────────────────
loo_rows <- list()

for (imp in imp_start:imp_end) {
  out_file_rds <- file.path(out_dir, sprintf("imp_%03d.rds", imp))
  out_file_loo <- file.path(out_dir, sprintf("loo_%03d.rds", imp))

  # Skip if LOO already computed
  if (file.exists(out_file_loo)) {
    cat(sprintf("[imp %03d] LOO already exists, loading...\n", imp))
    loo_result <- readRDS(out_file_loo)
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      model = "m1_v3",
      imputation = imp,
      elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
      se_elpd = loo_result$estimates["elpd_loo", "SE"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      n_high_pareto_k = sum(loo_result$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
    next
  }

  cat(sprintf("[imp %03d] Fitting M1 (v3)... ", imp))
  t0 <- Sys.time()

  dat <- imputed_data[[imp]]

  fit <- tryCatch(
    brm(bf_formula,
        data = dat,
        family = cumulative("logit"),
        prior = prior_spec,
        chains = 4, iter = 4000, warmup = 2000,
        cores = 2,
        init = 0.5,
        control = list(adapt_delta = 0.95),
        seed = 20250310 + imp,
        silent = 2, refresh = 0),
    error = function(e) {
      cat(sprintf("FIT ERROR: %s\n", e$message))
      return(NULL)
    }
  )

  if (is.null(fit)) next

  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("fit %.0fs, ", t_fit))

  saveRDS(fit, out_file_rds)

  # Compute LOO
  loo_result <- tryCatch(
    loo(fit, cores = 1),
    error = function(e) {
      cat(sprintf("LOO ERROR: %s\n", e$message))
      return(NULL)
    }
  )

  if (!is.null(loo_result)) {
    saveRDS(loo_result, out_file_loo)
    elpd <- loo_result$estimates["elpd_loo", "Estimate"]
    p_loo <- loo_result$estimates["p_loo", "Estimate"]
    n_bad <- sum(loo_result$diagnostics$pareto_k > 0.7)

    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      model = "m1_v3",
      imputation = imp,
      elpd_loo = elpd,
      se_elpd = loo_result$estimates["elpd_loo", "SE"],
      p_loo = p_loo,
      n_high_pareto_k = n_bad,
      stringsAsFactors = FALSE
    )
    cat(sprintf("LOO ELPD=%.1f (p_loo=%.1f, bad_k=%d)\n", elpd, p_loo, n_bad))
  }

  rm(fit); gc(verbose = FALSE)
}

# ── Save combined LOO results ────────────────────────────
if (length(loo_rows) > 0) {
  loo_df <- do.call(rbind, loo_rows)
  out_csv <- file.path(out_dir, "loo_results.csv")
  write.csv(loo_df, out_csv, row.names = FALSE)
  cat(sprintf("\n=== SUMMARY ===\n"))
  cat(sprintf("Saved %d LOO results to %s\n", nrow(loo_df), out_csv))
  cat(sprintf("Mean ELPD: %.1f (sd=%.1f)\n", mean(loo_df$elpd_loo), sd(loo_df$elpd_loo)))
}

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

#!/usr/bin/env Rscript
############################################################
# 99_test_proportional_odds.R
# Quick frequentist test of the proportional odds assumption
# using ordinal::clm() nominal_test() and scale_test()
#
# Tests on first 5 imputed datasets to check stability.
# Uses the sensitivity model (M4) since that's where the
# strongest effects were found.
#
# Run from: the repository root.
# Usage: Rscript scripts/99_test_proportional_odds.R
############################################################

suppressPackageStartupMessages({
  library(ordinal)
})

cat("=================================================================\n")
cat("PROPORTIONAL ODDS ASSUMPTION TEST\n")
cat("Model: M4 (sensitivity), ordinal::clm()\n")
cat("=================================================================\n\n")

# ── Load imputed data ────────────────────────────────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets (n=%d each)\n\n",
            length(imputed_data), nrow(imputed_data[[1]])))

# ── Gene list ────────────────────────────────────────────
genes <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
           "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

# ── Formula (M4 with sensitivity) ────────────────────────
# Matches the brms v3 specification: cov + gene + int
f <- as.formula(paste0(
  "Final7PointLikertScaleFactor ~ PC1 + PC2 + PC3 + Sex + infant_age + ",
  "pass + diff + sens + ",
  paste(genes, collapse = " + "), " + ",
  paste(paste0("sens:", genes), collapse = " + ")
))
cat("Formula:\n")
print(f)
cat("\n")

# ── Test on multiple imputations ─────────────────────────
n_imp <- 10  # test on first 10 imputations for stability

nominal_results <- list()
scale_results <- list()

for (i in 1:n_imp) {
  cat(sprintf("── Imputation %d ──\n", i))

  df <- imputed_data[[i]]

  # Fit proportional odds model
  fit <- tryCatch(
    clm(f, data = df, link = "logit"),
    error = function(e) {
      cat(sprintf("  clm() failed: %s\n", e$message))
      NULL
    }
  )

  if (is.null(fit)) next

  cat(sprintf("  Fitted: %d obs, %d coefficients\n",
              fit$nobs, length(coef(fit))))

  # Nominal test: tests whether coefficients vary across thresholds
  # (i.e., whether each predictor has a different effect at each cutpoint)
  nom <- tryCatch(
    nominal_test(fit),
    error = function(e) {
      cat(sprintf("  nominal_test() failed: %s\n", e$message))
      NULL
    }
  )

  if (!is.null(nom)) {
    nominal_results[[i]] <- nom
    cat("  nominal_test() results:\n")
    # Print the test table
    print(nom)
    cat("\n")
  }

  # Scale test: tests whether residual variance changes across thresholds
  sc <- tryCatch(
    scale_test(fit),
    error = function(e) {
      cat(sprintf("  scale_test() failed: %s\n", e$message))
      NULL
    }
  )

  if (!is.null(sc)) {
    scale_results[[i]] <- sc
    cat("  scale_test() results:\n")
    print(sc)
    cat("\n")
  }
}

# ── Summarize across imputations ─────────────────────────
cat("\n=================================================================\n")
cat("SUMMARY ACROSS IMPUTATIONS\n")
cat("=================================================================\n\n")

# Collect nominal test p-values across imputations
if (length(nominal_results) > 0) {
  cat("NOMINAL TEST (tests whether slopes vary across cutpoints)\n")
  cat("H0: proportional odds holds for each predictor\n\n")

  # Get predictor names from first result
  first_nom <- nominal_results[[which(!sapply(nominal_results, is.null))[1]]]
  pred_names <- rownames(first_nom)

  # Collect p-values
  pval_matrix <- matrix(NA, nrow = length(pred_names), ncol = length(nominal_results))
  rownames(pval_matrix) <- pred_names

  for (i in seq_along(nominal_results)) {
    if (!is.null(nominal_results[[i]])) {
      # Extract p-values (last column)
      pvals <- nominal_results[[i]][, ncol(nominal_results[[i]])]
      pval_matrix[, i] <- pvals
    }
  }

  # Report median p-value and proportion significant across imputations
  cat(sprintf("%-25s  Median p   Min p    Max p    %% sig (p<0.05)\n",
              "Predictor"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (j in seq_len(nrow(pval_matrix))) {
    pvals <- pval_matrix[j, ]
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) > 0) {
      cat(sprintf("%-25s  %.4f    %.4f   %.4f   %.0f%%\n",
                  pred_names[j],
                  median(pvals), min(pvals), max(pvals),
                  100 * mean(pvals < 0.05)))
    }
  }
  cat("\n")
}

if (length(scale_results) > 0) {
  cat("SCALE TEST (tests whether residual variance depends on predictors)\n")
  cat("H0: constant residual scale across predictor levels\n\n")

  first_sc <- scale_results[[which(!sapply(scale_results, is.null))[1]]]
  pred_names_sc <- rownames(first_sc)

  pval_matrix_sc <- matrix(NA, nrow = length(pred_names_sc), ncol = length(scale_results))
  rownames(pval_matrix_sc) <- pred_names_sc

  for (i in seq_along(scale_results)) {
    if (!is.null(scale_results[[i]])) {
      pvals <- scale_results[[i]][, ncol(scale_results[[i]])]
      pval_matrix_sc[, i] <- pvals
    }
  }

  cat(sprintf("%-25s  Median p   Min p    Max p    %% sig (p<0.05)\n",
              "Predictor"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (j in seq_len(nrow(pval_matrix_sc))) {
    pvals <- pval_matrix_sc[j, ]
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) > 0) {
      cat(sprintf("%-25s  %.4f    %.4f   %.4f   %.0f%%\n",
                  pred_names_sc[j],
                  median(pvals), min(pvals), max(pvals),
                  100 * mean(pvals < 0.05)))
    }
  }
  cat("\n")
}

cat("=== DONE ===\n")

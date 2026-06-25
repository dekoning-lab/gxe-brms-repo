#!/usr/bin/env Rscript
############################################################
# 52_aggregate_sensitivity_v3.R
# Aggregate prior sensitivity analysis results (s=0.25, 0.50, 1.00)
# for the M4 (full) model across 3 parenting variables.
#
# Run from: the repository root.
# Usage: Rscript scripts/52_aggregate_sensitivity_v3.R
############################################################

cat("=================================================================\n")
cat("V3 PRIOR SENSITIVITY AGGREGATION\n")
cat("=================================================================\n\n")

parenting_vars  <- c("sens", "cont", "unre")
parenting_label <- c(sens = "Sensitivity", cont = "Controlling",
                     unre = "Unresponsiveness")

scales <- c("s025" = 0.25, "s050" = 0.50, "s100" = 1.00)

# ── Load all summaries ──────────────────────────────────
all_sens <- list()

for (scale_label in names(scales)) {
  sigma_val <- scales[scale_label]

  # s050 = primary results in v3_m4/
  if (scale_label == "s050") {
    base_dir <- "results/v3_m4"
  } else {
    base_dir <- file.path("results", paste0("v3_m4_", scale_label))
  }

  for (pv in parenting_vars) {
    f <- file.path(base_dir, pv, "summary.csv")
    if (file.exists(f)) {
      d <- read.csv(f, stringsAsFactors = FALSE)
      d$sigma_scale <- sigma_val
      d$scale_label <- scale_label
      all_sens[[length(all_sens) + 1]] <- d
      cat(sprintf("  Loaded %s/%s: %d rows\n", scale_label, pv, nrow(d)))
    } else {
      cat(sprintf("  WARNING: %s not found\n", f))
    }
  }
}

sens_df <- do.call(rbind, all_sens)
rownames(sens_df) <- NULL
cat(sprintf("\nTotal: %d rows across %d scales\n\n", nrow(sens_df), length(scales)))

# ── Hyperparameter comparison ───────────────────────────
cat("=== HYPERPARAMETER POSTERIORS ACROSS SCALES ===\n\n")

hyper <- sens_df[sens_df$effect_type == "hyperparameter", ]

cat(sprintf("%-8s %-6s %-12s %8s %8s %8s %8s\n",
            "Parent.", "Scale", "Parameter", "Mean", "SD", "Q2.5", "Q97.5"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (pv in parenting_vars) {
  for (param in c("sigma_gene", "sigma_int")) {
    for (sl in names(scales)) {
      r <- hyper[hyper$parenting == pv & hyper$term == param & hyper$scale_label == sl, ]
      if (nrow(r) > 0) {
        cat(sprintf("%-8s %-6s %-12s %8.4f %8.4f %8.4f %8.4f\n",
                    pv, sl, param, r$estimate[1], r$sd[1], r$q025[1], r$q975[1]))
      }
    }
  }
  cat("\n")
}

# ── PP stability across scales ──────────────────────────
cat("=== POSTERIOR PROBABILITY STABILITY ACROSS SCALES ===\n\n")

effects <- sens_df[sens_df$effect_type %in% c("gene_main", "interaction"), ]

# Wide format: PP by scale
pp_wide <- reshape(
  effects[, c("parenting", "term", "effect_type", "post_prob", "estimate", "scale_label")],
  timevar = "scale_label",
  idvar = c("parenting", "term", "effect_type"),
  direction = "wide"
)

# Compute PP range across scales
pp_cols <- grep("^post_prob\\.", names(pp_wide), value = TRUE)
est_cols <- grep("^estimate\\.", names(pp_wide), value = TRUE)

pp_wide$pp_range <- apply(pp_wide[, pp_cols], 1, function(x) diff(range(x, na.rm = TRUE)))
pp_wide$max_pp <- apply(pp_wide[, pp_cols], 1, function(x) max(x, na.rm = TRUE))

# Sort by max PP range (most sensitive first)
pp_wide <- pp_wide[order(-pp_wide$pp_range), ]

cat(sprintf("%-8s %-20s %8s %8s %8s %8s\n",
            "Parent.", "Term", "PP s025", "PP s050", "PP s100", "PP range"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (i in seq_len(min(30, nrow(pp_wide)))) {
  r <- pp_wide[i, ]
  cat(sprintf("%-8s %-20s %7.1f%% %7.1f%% %7.1f%% %7.1f pp\n",
              r$parenting, r$term,
              r$post_prob.s025 * 100, r$post_prob.s050 * 100, r$post_prob.s100 * 100,
              r$pp_range * 100))
}

# ── LOO comparison across scales ──────────────────────────
cat("\n\n=== LOO-ELPD ACROSS SCALES ===\n\n")

for (sl in names(scales)) {
  if (sl == "s050") {
    base_dir <- "results/v3_m4"
  } else {
    base_dir <- file.path("results", paste0("v3_m4_", sl))
  }

  for (pv in parenting_vars) {
    f <- file.path(base_dir, pv, "loo_results.csv")
    if (file.exists(f)) {
      d <- read.csv(f, stringsAsFactors = FALSE)
      cat(sprintf("  %s/%s: mean ELPD = %.2f (n=%d imps)\n",
                  sl, pv, mean(d$elpd_loo), nrow(d)))
    }
  }
}

# ── Save comprehensive sensitivity results ────────────────
dir.create("results/v3_sensitivity", recursive = TRUE, showWarnings = FALSE)
write.csv(sens_df, "results/v3_sensitivity/all_sensitivity.csv", row.names = FALSE)
write.csv(pp_wide, "results/v3_sensitivity/pp_stability.csv", row.names = FALSE)

# Save hyperparameter summary in nice format
hyper_wide <- reshape(
  hyper[, c("parenting", "term", "estimate", "sd", "q025", "q975", "scale_label")],
  timevar = "scale_label",
  idvar = c("parenting", "term"),
  direction = "wide"
)
write.csv(hyper_wide, "results/v3_sensitivity/hyperparameter_comparison.csv", row.names = FALSE)

cat("\nSaved: results/v3_sensitivity/all_sensitivity.csv\n")
cat("Saved: results/v3_sensitivity/pp_stability.csv\n")
cat("Saved: results/v3_sensitivity/hyperparameter_comparison.csv\n")
cat("\n=== DONE ===\n")

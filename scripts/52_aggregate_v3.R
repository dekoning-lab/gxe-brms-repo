#!/usr/bin/env Rscript
############################################################
# 52_aggregate_v3.R
# Aggregate v3 model comparison results (with pass+diff
# as child behavior covariates in all models).
#
# Reads LOO results from M1-M4 and produces comparison tables.
#
# Run from: the repository root.
# Usage: Rscript scripts/52_aggregate_v3.R
############################################################

cat("=================================================================\n")
cat("AGGREGATING V3 MODEL COMPARISON RESULTS\n")
cat("(pass + diff as child behavior covariates in all models)\n")
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

pv_list <- c("sens", "cont", "unre")
pv_labels <- c(sens = "Sensitivity", cont = "Controlling", unre = "Unresponsiveness")

# ── Load M1 LOO ─────────────────────────────────────────
m1_dir <- "results/v3_m1"
m1_files <- list.files(m1_dir, pattern = "^loo_\\d+\\.rds$", full.names = TRUE)

if (length(m1_files) > 0) {
  m1_loo <- do.call(rbind, lapply(m1_files, function(f) {
    imp <- as.integer(sub(".*loo_(\\d+)\\.rds", "\\1", basename(f)))
    loo_obj <- readRDS(f)
    data.frame(
      model = "M1", parenting = "none", imputation = imp,
      elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
      se_elpd = loo_obj$estimates["elpd_loo", "SE"],
      p_loo = loo_obj$estimates["p_loo", "Estimate"],
      n_high_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
  }))
  cat(sprintf("M1: loaded %d LOO results, mean ELPD = %.1f\n",
              nrow(m1_loo), mean(m1_loo$elpd_loo)))
} else {
  cat("M1: no LOO results found\n")
  m1_loo <- data.frame()
}

# ── Load M3 LOO ─────────────────────────────────────────
m3_dir <- "results/v3_m3"
m3_files <- list.files(m3_dir, pattern = "^loo_\\d+\\.rds$", full.names = TRUE)

if (length(m3_files) > 0) {
  m3_loo <- do.call(rbind, lapply(m3_files, function(f) {
    imp <- as.integer(sub(".*loo_(\\d+)\\.rds", "\\1", basename(f)))
    loo_obj <- readRDS(f)
    data.frame(
      model = "M3", parenting = "none", imputation = imp,
      elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
      se_elpd = loo_obj$estimates["elpd_loo", "SE"],
      p_loo = loo_obj$estimates["p_loo", "Estimate"],
      n_high_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
  }))
  cat(sprintf("M3: loaded %d LOO results, mean ELPD = %.1f\n",
              nrow(m3_loo), mean(m3_loo$elpd_loo)))
} else {
  cat("M3: no LOO results found\n")
  m3_loo <- data.frame()
}

# ── Load M2 and M4 LOO (per parenting) ──────────────────
m2_loo <- data.frame()
m4_loo <- data.frame()

for (pv in pv_list) {
  # M2
  m2_files <- list.files(file.path("results/v3_m2", pv),
                          pattern = "^loo_\\d+\\.rds$", full.names = TRUE)
  if (length(m2_files) > 0) {
    m2_pv <- do.call(rbind, lapply(m2_files, function(f) {
      imp <- as.integer(sub(".*loo_(\\d+)\\.rds", "\\1", basename(f)))
      loo_obj <- readRDS(f)
      data.frame(
        model = "M2", parenting = pv, imputation = imp,
        elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
        se_elpd = loo_obj$estimates["elpd_loo", "SE"],
        p_loo = loo_obj$estimates["p_loo", "Estimate"],
        n_high_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
        stringsAsFactors = FALSE
      )
    }))
    m2_loo <- rbind(m2_loo, m2_pv)
    cat(sprintf("M2 %s: %d LOO results, mean ELPD = %.1f\n",
                pv, nrow(m2_pv), mean(m2_pv$elpd_loo)))
  }

  # M4
  m4_files <- list.files(file.path("results/v3_m4", pv),
                          pattern = "^loo_\\d+\\.rds$", full.names = TRUE)
  if (length(m4_files) > 0) {
    m4_pv <- do.call(rbind, lapply(m4_files, function(f) {
      imp <- as.integer(sub(".*loo_(\\d+)\\.rds", "\\1", basename(f)))
      loo_obj <- readRDS(f)
      data.frame(
        model = "M4", parenting = pv, imputation = imp,
        elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
        se_elpd = loo_obj$estimates["elpd_loo", "SE"],
        p_loo = loo_obj$estimates["p_loo", "Estimate"],
        n_high_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
        stringsAsFactors = FALSE
      )
    }))
    m4_loo <- rbind(m4_loo, m4_pv)
    cat(sprintf("M4 %s: %d LOO results, mean ELPD = %.1f\n",
                pv, nrow(m4_pv), mean(m4_pv$elpd_loo)))
  }
}

# ── Combine all ──────────────────────────────────────────
all_loo <- rbind(m1_loo, m2_loo, m3_loo, m4_loo)
write.csv(all_loo, "results/v3_all_loo.csv", row.names = FALSE)
cat(sprintf("\nTotal: %d LOO results saved to results/v3_all_loo.csv\n", nrow(all_loo)))

# ── Pairwise comparisons ────────────────────────────────
cat("\n=== PAIRWISE COMPARISONS ===\n\n")

pairwise_rows <- list()

for (pv in pv_list) {
  # Get matched imputations
  m1_sub <- m1_loo
  m2_sub <- m2_loo[m2_loo$parenting == pv, ]
  m3_sub <- m3_loo
  m4_sub <- m4_loo[m4_loo$parenting == pv, ]

  # M2 vs M1: does parenting improve on demographics + child behavior?
  if (nrow(m1_sub) > 0 && nrow(m2_sub) > 0) {
    shared <- intersect(m1_sub$imputation, m2_sub$imputation)
    d <- m2_sub$elpd_loo[match(shared, m2_sub$imputation)] -
         m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      parenting = pv, comparison = "M2 vs M1",
      mean_delta = mean(d), sd_delta = sd(d),
      pct_positive = mean(d > 0) * 100,
      n_imputations = length(shared),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s M2 vs M1: dELPD = %+.2f (%.0f%% positive, n=%d)\n",
                pv_labels[pv], mean(d), mean(d > 0) * 100, length(shared)))
  }

  # M3 vs M1: do genes improve on demographics + child behavior?
  if (nrow(m1_sub) > 0 && nrow(m3_sub) > 0) {
    shared <- intersect(m1_sub$imputation, m3_sub$imputation)
    d <- m3_sub$elpd_loo[match(shared, m3_sub$imputation)] -
         m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      parenting = pv, comparison = "M3 vs M1",
      mean_delta = mean(d), sd_delta = sd(d),
      pct_positive = mean(d > 0) * 100,
      n_imputations = length(shared),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s M3 vs M1: dELPD = %+.2f (%.0f%% positive, n=%d)\n",
                pv_labels[pv], mean(d), mean(d > 0) * 100, length(shared)))
  }

  # M4 vs M2: do genes + GxE improve on demographics + parenting + child behavior?
  if (nrow(m2_sub) > 0 && nrow(m4_sub) > 0) {
    shared <- intersect(m2_sub$imputation, m4_sub$imputation)
    d <- m4_sub$elpd_loo[match(shared, m4_sub$imputation)] -
         m2_sub$elpd_loo[match(shared, m2_sub$imputation)]
    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      parenting = pv, comparison = "M4 vs M2",
      mean_delta = mean(d), sd_delta = sd(d),
      pct_positive = mean(d > 0) * 100,
      n_imputations = length(shared),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s M4 vs M2: dELPD = %+.2f (%.0f%% positive, n=%d)\n",
                pv_labels[pv], mean(d), mean(d > 0) * 100, length(shared)))
  }

  # M4 vs M1: does the full model improve on demographics + child behavior?
  if (nrow(m1_sub) > 0 && nrow(m4_sub) > 0) {
    shared <- intersect(m1_sub$imputation, m4_sub$imputation)
    d <- m4_sub$elpd_loo[match(shared, m4_sub$imputation)] -
         m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      parenting = pv, comparison = "M4 vs M1",
      mean_delta = mean(d), sd_delta = sd(d),
      pct_positive = mean(d > 0) * 100,
      n_imputations = length(shared),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s M4 vs M1: dELPD = %+.2f (%.0f%% positive, n=%d)\n",
                pv_labels[pv], mean(d), mean(d > 0) * 100, length(shared)))
  }

  cat("\n")
}

if (length(pairwise_rows) > 0) {
  pairwise_df <- do.call(rbind, pairwise_rows)
  write.csv(pairwise_df, "results/v3_pairwise.csv", row.names = FALSE)
  cat("Pairwise results saved to results/v3_pairwise.csv\n")
}

# ── Summary table ────────────────────────────────────────
cat("\n=== LOO-ELPD SUMMARY TABLE ===\n\n")
cat(sprintf("%-14s  %-8s  %6s  %6s  %6s  %5s\n",
            "Model", "Parent.", "ELPD", "SD", "p_loo", "bad_k"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (m in c("M1", "M3")) {
  sub <- all_loo[all_loo$model == m, ]
  if (nrow(sub) > 0) {
    cat(sprintf("%-14s  %-8s  %6.1f  %6.2f  %6.1f  %5.1f\n",
                m, "—", mean(sub$elpd_loo), sd(sub$elpd_loo),
                mean(sub$p_loo), mean(sub$n_high_k)))
  }
}

for (pv in pv_list) {
  for (m in c("M2", "M4")) {
    sub <- all_loo[all_loo$model == m & all_loo$parenting == pv, ]
    if (nrow(sub) > 0) {
      cat(sprintf("%-14s  %-8s  %6.1f  %6.2f  %6.1f  %5.1f\n",
                  m, pv, mean(sub$elpd_loo), sd(sub$elpd_loo),
                  mean(sub$p_loo), mean(sub$n_high_k)))
    }
  }
}

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

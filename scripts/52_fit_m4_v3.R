#!/usr/bin/env Rscript
############################################################
# 52_fit_m4_v3.R
# M4 (v3): Full hierarchical model with child behavior covariates
#
# Model (for parenting variable P):
#   Final7PointLikertScaleFactor ~ cov + gene + int
#   cov  ~ 0 + PC1 + PC2 + PC3 + Sex + infant_age + pass + diff + P  (8 terms)
#   gene ~ 0 + BDNF + ... + X5HTTLPR                                  (9 terms)
#   int  ~ 0 + P:BDNF + ... + P:X5HTTLPR                              (9 terms)
#
# Priors:
#   cov:  Normal(0, 2)
#   gene: Normal(0, sigma_gene), sigma_gene ~ half-Normal(0, 0.5)
#   int:  Normal(0, sigma_int),  sigma_int  ~ half-Normal(0, 0.5)
#
# Key change from 23_fit_hierarchical.R:
#   - Adds pass + diff as child behavior covariates in cov block
#   - Uses augmented v2 imputed data
#   - Only sens/cont/unre as parenting variable
#   - Computes LOO inline
#
# Run from: the repository root (locally, or on a cluster after sourcing config.sh; see README).
# Usage:    Rscript scripts/52_fit_m4_v3.R <parenting> [start] [end]
#
# Output:   results/v3_m4/<parenting>/imp_NNN.rds
#           results/v3_m4/<parenting>/loo_NNN.rds
#           results/v3_m4/<parenting>/summary.csv
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

# ── Parse arguments ──────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/52_fit_m4_v3.R <parenting> [start] [end]")
}

parenting_var <- args[1]
start_idx <- if (length(args) >= 2) as.integer(args[2]) else 1L
end_idx   <- if (length(args) >= 3) as.integer(args[3]) else 100L
sigma_scale <- 0.5

stopifnot(parenting_var %in% c("sens", "cont", "unre"))

cat("=================================================================\n")
cat(sprintf("M4 (v3): FULL HIERARCHICAL: %s x genes\n", parenting_var))
cat("  cov  ~ PC1-3 + Sex + infant_age + pass + diff + P (8 terms)\n")
cat("  gene ~ 9 candidate genes\n")
cat("  int  ~ P × 9 genes\n")
cat(sprintf("  sigma_gene, sigma_int ~ half-Normal(0, %.2f)\n", sigma_scale))
cat(sprintf("Imputations: %d to %d\n", start_idx, end_idx))
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

# ── Output directory ─────────────────────────────────────
out_dir <- file.path("results", "v3_m4", parenting_var)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load augmented imputed data ─────────────────────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets (n=%d each)\n",
            length(imputed_data), nrow(imputed_data[[1]])))
stopifnot("pass" %in% names(imputed_data[[1]]))
stopifnot("diff" %in% names(imputed_data[[1]]))
cat("  Confirmed: pass and diff present in data\n\n")

# ── Build non-linear model formula (3 nlpar blocks) ─────
genes <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
           "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

cov_rhs  <- paste(c("0", "PC1", "PC2", "PC3", "Sex", "infant_age",
                     "pass", "diff", parenting_var), collapse = " + ")
gene_rhs <- paste(c("0", genes), collapse = " + ")
int_rhs  <- paste(c("0", paste0(parenting_var, ":", genes)), collapse = " + ")

model_formula <- bf(
  as.formula("Final7PointLikertScaleFactor ~ cov + gene + int"),
  as.formula(paste("cov ~", cov_rhs)),
  as.formula(paste("gene ~", gene_rhs)),
  as.formula(paste("int ~", int_rhs)),
  nl = TRUE,
  decomp = "QR"
)

cat("Non-linear formula:\n")
cat(sprintf("  cov  ~ %s\n", cov_rhs))
cat(sprintf("  gene ~ %s\n", gene_rhs))
cat(sprintf("  int  ~ %s\n", int_rhs))
cat(sprintf("Parameters: 8 covariates + 9 gene mains + 9 interactions = 26\n\n"))

# ── Define brms priors ──────────────────────────────────
model_prior <- prior(normal(0, 2), class = "b", nlpar = "cov")

# ── stanvar() injection for hierarchical structure ──────
sv_parameters <- stanvar(
  scode = "  real<lower=0> sigma_gene;\n  real<lower=0> sigma_int;",
  block = "parameters"
)

sv_model <- stanvar(
  scode = paste0(
    sprintf("  // Hierarchical hyperpriors (half-Normal(0, %.2f))\n", sigma_scale),
    sprintf("  target += normal_lpdf(sigma_gene | 0, %.2f) + log(2);\n", sigma_scale),
    sprintf("  target += normal_lpdf(sigma_int  | 0, %.2f) + log(2);\n", sigma_scale),
    "  // Hierarchical priors on gene and interaction coefficients\n",
    "  target += normal_lpdf(b_gene | 0, sigma_gene);\n",
    "  target += normal_lpdf(b_int  | 0, sigma_int);"
  ),
  block = "model"
)

all_stanvars <- sv_parameters + sv_model

# ── Verify Stan code ────────────────────────────────────
cat("Verifying Stan code via make_stancode()...\n")
stan_code <- make_stancode(
  formula = model_formula,
  data = imputed_data[[1]],
  family = cumulative("logit"),
  prior = model_prior,
  stanvars = all_stanvars
)

has_b_gene <- grepl("\\bb_gene\\b", stan_code)
has_b_int  <- grepl("\\bb_int\\b", stan_code)
if (!has_b_gene || !has_b_int) {
  stop("FATAL: brms does not use expected coefficient names (b_gene/b_int)")
}
cat("  Confirmed: b_gene and b_int found in Stan code.\n")

stan_code_file <- file.path(out_dir, "stan_code.stan")
writeLines(stan_code, stan_code_file)
cat(sprintf("  Stan code saved to: %s\n\n", stan_code_file))

# ── Fit + LOO loop ──────────────────────────────────────
n_fitted <- 0
n_skipped <- 0
fit_times <- c()
loo_rows <- list()

for (i in start_idx:end_idx) {
  out_file <- file.path(out_dir, sprintf("imp_%03d.rds", i))
  out_file_loo <- file.path(out_dir, sprintf("loo_%03d.rds", i))

  if (file.exists(out_file_loo)) {
    cat(sprintf("[%s] Imp %03d: LOO exists, loading...\n", parenting_var, i))
    loo_result <- readRDS(out_file_loo)
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      parenting = parenting_var, model = "m4_v3", imputation = i,
      elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
      se_elpd = loo_result$estimates["elpd_loo", "SE"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      n_high_pareto_k = sum(loo_result$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
    n_skipped <- n_skipped + 1
    next
  }

  cat(sprintf("[%s] Imp %03d: fitting... ", parenting_var, i))
  t0 <- Sys.time()

  fit <- tryCatch({
    brm(
      formula = model_formula,
      data = imputed_data[[i]],
      family = cumulative("logit"),
      prior = model_prior,
      stanvars = all_stanvars,
      chains = 4, cores = 2,
      iter = 8000, warmup = 4000,
      init = 0.5,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      seed = 20250310 + i,
      silent = 2, refresh = 0
    )
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(fit)) next

  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  fit_times <- c(fit_times, elapsed)

  saveRDS(fit, out_file)
  n_fitted <- n_fitted + 1

  np <- nuts_params(fit)
  n_div <- sum(np$Value[np$Parameter == "divergent__"])
  cat(sprintf("fit %.0fs (%d div), ", elapsed, n_div))

  # Compute LOO inline
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
    n_bad <- sum(loo_result$diagnostics$pareto_k > 0.7)
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      parenting = parenting_var, model = "m4_v3", imputation = i,
      elpd_loo = elpd,
      se_elpd = loo_result$estimates["elpd_loo", "SE"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      n_high_pareto_k = n_bad,
      stringsAsFactors = FALSE
    )
    cat(sprintf("LOO=%.1f (bad_k=%d)\n", elpd, n_bad))
  }

  rm(fit); gc(verbose = FALSE)
}

# ── Save combined LOO results ────────────────────────────
if (length(loo_rows) > 0) {
  loo_df <- do.call(rbind, loo_rows)
  write.csv(loo_df, file.path(out_dir, "loo_results.csv"), row.names = FALSE)
  cat(sprintf("\nLOO: %d results, mean ELPD=%.1f\n",
              nrow(loo_df), mean(loo_df$elpd_loo)))
}

# ── Generate summary CSV (pool posteriors across imputations) ─
cat(sprintf("\n[%s] Generating summary...\n", parenting_var))

all_files <- sort(list.files(out_dir, pattern = "^imp_\\d+\\.rds$", full.names = TRUE))

if (length(all_files) > 0) {
  all_draws <- data.frame()
  for (f in all_files) {
    fit <- readRDS(f)
    draws <- as.data.frame(fit)
    all_draws <- rbind(all_draws, draws)
  }

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

  # Add hyperparameter rows
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

  write.csv(results, file.path(out_dir, "summary.csv"), row.names = FALSE)
  cat(sprintf("[%s] Summary saved (%d effects + hyperparameters)\n",
              parenting_var, sum(results$effect_type != "hyperparameter")))

  # Show top interactions
  interactions <- results[results$effect_type == "interaction", ]
  interactions <- interactions[order(-interactions$post_prob), ]
  cat(sprintf("\n[%s] Top 5 interactions:\n", parenting_var))
  for (j in seq_len(min(5, nrow(interactions)))) {
    r <- interactions[j, ]
    cat(sprintf("  %s: PP=%.3f (est=%.3f, %s)\n",
                r$term, r$post_prob, r$estimate, r$direction))
  }
}

# ── Final summary ────────────────────────────────────────
cat(sprintf("\n=================================================================\n"))
cat(sprintf("[%s] COMPLETE: fitted=%d, skipped=%d, total=%d\n",
            parenting_var, n_fitted, n_skipped, length(all_files)))
if (length(fit_times) > 0) {
  cat(sprintf("  Mean fit time: %.0fs, Total: %.1f min\n",
              mean(fit_times), sum(fit_times) / 60))
}
cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n")

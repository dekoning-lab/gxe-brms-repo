#!/usr/bin/env Rscript
############################################################
# 52_fit_m3_v3.R
# M3 (v3): Demographics + child behavior + 9 genes (hierarchical)
#
# Model:
#   Final7PointLikertScaleFactor ~ cov + gene
#   cov  ~ 0 + PC1 + PC2 + PC3 + Sex + infant_age + pass + diff  (7 terms)
#   gene ~ 0 + BDNF + ... + X5HTTLPR                              (9 terms)
#
# Priors:
#   cov:  Normal(0, 2)
#   gene: Normal(0, sigma_gene), sigma_gene ~ half-Normal(0, 0.5)
#
# Key change from 49_fit_m3_genes_only.R:
#   - Adds pass + diff as child behavior covariates
#   - Uses augmented v2 imputed data
#   - Parenting-independent: 100 fits serve all comparisons
#   - Computes LOO inline
#
# Run from: the repository root (locally, or on a cluster after sourcing config.sh; see README).
# Usage:    Rscript scripts/52_fit_m3_v3.R [start] [end]
#
# Output:   results/v3_m3/imp_NNN.rds
#           results/v3_m3/loo_NNN.rds
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
start_idx <- if (length(args) >= 1) as.integer(args[1]) else 1L
end_idx   <- if (length(args) >= 2) as.integer(args[2]) else 100L

cat("=================================================================\n")
cat(sprintf("M3 (v3): DEMOGRAPHICS + CHILD BEHAVIOR + GENES: imp %d-%d\n",
            start_idx, end_idx))
cat("  cov  ~ PC1 + PC2 + PC3 + Sex + infant_age + pass + diff (7 terms)\n")
cat("  gene ~ BDNF + CNR1.77 + ... + X5HTTLPR (9 terms)\n")
cat("  Priors: Normal(0,2) on cov; Normal(0, sigma_gene) on genes\n")
cat("          sigma_gene ~ half-Normal(0, 0.5)\n")
cat(sprintf("Start time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

# ── Load augmented imputed data ─────────────────────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets (n=%d each)\n",
            length(imputed_data), nrow(imputed_data[[1]])))
stopifnot("pass" %in% names(imputed_data[[1]]))
stopifnot("diff" %in% names(imputed_data[[1]]))
cat("  Confirmed: pass and diff present in data\n\n")

# ── Output directory ─────────────────────────────────────
out_dir <- "results/v3_m3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Build non-linear model formula (2 nlpar blocks) ─────
genes <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
           "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

cov_rhs  <- paste(c("0", "PC1", "PC2", "PC3", "Sex", "infant_age",
                     "pass", "diff"), collapse = " + ")
gene_rhs <- paste(c("0", genes), collapse = " + ")

model_formula <- bf(
  Final7PointLikertScaleFactor ~ cov + gene,
  as.formula(paste("cov ~", cov_rhs)),
  as.formula(paste("gene ~", gene_rhs)),
  nl = TRUE,
  decomp = "QR"
)

cat("Formula:\n")
cat(sprintf("  cov  ~ %s\n", cov_rhs))
cat(sprintf("  gene ~ %s\n", gene_rhs))
cat(sprintf("Parameters: 7 covariates + 9 gene mains = 16\n\n"))

# Prior for covariates only
model_prior <- prior(normal(0, 2), class = "b", nlpar = "cov")

# stanvar: single sigma_gene hyperparameter
sigma_scale <- 0.5

sv_parameters <- stanvar(
  scode = "  real<lower=0> sigma_gene;",
  block = "parameters"
)

sv_model <- stanvar(
  scode = paste0(
    sprintf("  // Hierarchical hyperprior: half-Normal(0, %.2f)\n", sigma_scale),
    sprintf("  target += normal_lpdf(sigma_gene | 0, %.2f) + log(2);\n", sigma_scale),
    "  // Hierarchical prior on gene coefficients\n",
    "  target += normal_lpdf(b_gene | 0, sigma_gene);"
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
if (!has_b_gene) {
  stop("FATAL: brms does not use 'b_gene' in Stan code")
}
cat("  Confirmed: b_gene found in Stan code.\n")

stan_code_file <- file.path(out_dir, "stan_code.stan")
writeLines(stan_code, stan_code_file)
cat(sprintf("  Stan code saved to: %s\n\n", stan_code_file))

# ── Fit + LOO loop ──────────────────────────────────────
loo_rows <- list()

for (i in start_idx:end_idx) {
  out_file <- file.path(out_dir, sprintf("imp_%03d.rds", i))
  out_file_loo <- file.path(out_dir, sprintf("loo_%03d.rds", i))

  if (file.exists(out_file_loo)) {
    cat(sprintf("  Imp %03d: LOO exists, loading...\n", i))
    loo_result <- readRDS(out_file_loo)
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      model = "m3_v3", imputation = i,
      elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
      se_elpd = loo_result$estimates["elpd_loo", "SE"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      n_high_pareto_k = sum(loo_result$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
    next
  }

  cat(sprintf("  Imp %03d: fitting... ", i))
  t0 <- Sys.time()

  fit <- tryCatch({
    brm(
      formula = model_formula,
      data = imputed_data[[i]],
      family = cumulative("logit"),
      prior = model_prior,
      stanvars = all_stanvars,
      chains = 4, cores = 4,
      iter = 8000, warmup = 4000,
      init = 0.5,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      seed = 20250310L + i,
      silent = 2, refresh = 0
    )
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(fit)) next

  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(fit, out_file)

  np <- nuts_params(fit)
  n_div <- sum(np$Value[np$Parameter == "divergent__"])
  cat(sprintf("fit %.0fs (%d div), ", t_fit, n_div))

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
      model = "m3_v3", imputation = i,
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
  cat(sprintf("\n=== SUMMARY: %d LOO results, mean ELPD=%.1f ===\n",
              nrow(loo_df), mean(loo_df$elpd_loo)))
}

cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

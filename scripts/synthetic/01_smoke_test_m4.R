#!/usr/bin/env Rscript
############################################################
# 01_smoke_test_m4.R
#
# Lightweight smoke test of the brms M4 model on the synthetic
# dataset. Fits ONE imputation only (out of 100) at a reduced
# sampling budget (2 chains × 1000 iter, 500 warmup, 500 sampling),
# to verify that:
#   - All R/brms/Stan dependencies are installed correctly
#   - The data + model formula compile
#   - The hierarchical stanvar() injection parses
#   - One full fit completes without errors
#
# This is NOT the production model — see 52_fit_m4_v3.R for the
# real fit (4 chains × 8000 iter, all 100 imputations).
#
# Run from: repository root
# Usage:    Rscript scripts/synthetic/01_smoke_test_m4.R [parenting]
#           parenting defaults to "sens"
############################################################

suppressPackageStartupMessages({
  library(brms); library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
parenting_var <- if (length(args) >= 1) args[1] else "sens"
stopifnot(parenting_var %in% c("sens","cont","unre"))

cat(sprintf("Smoke test: M4 brms fit on synthetic data, parenting=%s\n",
            parenting_var))
cat(sprintf("Start time: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets (n=%d each)\n",
            length(imputed_data), nrow(imputed_data[[1]])))
stopifnot("pass" %in% names(imputed_data[[1]]))
stopifnot("diff" %in% names(imputed_data[[1]]))

genes <- c("BDNF","CNR1.77","CNR1.10","DRD2","DRD4",
           "MAOA","SLC6A3.9R","SLC6A3.10R","X5HTTLPR")

cov_rhs  <- paste(c("0","PC1","PC2","PC3","Sex","infant_age",
                    "pass","diff",parenting_var), collapse=" + ")
gene_rhs <- paste(c("0", genes), collapse=" + ")
int_rhs  <- paste(c("0", paste0(parenting_var,":",genes)), collapse=" + ")

model_formula <- bf(
  as.formula("Final7PointLikertScaleFactor ~ cov + gene + int"),
  as.formula(paste("cov ~", cov_rhs)),
  as.formula(paste("gene ~", gene_rhs)),
  as.formula(paste("int  ~", int_rhs)),
  decomp = "QR",
  nl = TRUE
)

sigma_scale <- 0.5
model_prior <- prior(normal(0, 2), class = "b", nlpar = "cov")
sv_parameters <- stanvar(
  scode = "  real<lower=0> sigma_gene;\n  real<lower=0> sigma_int;",
  block = "parameters"
)
sv_model <- stanvar(
  scode = sprintf(paste(
    "  sigma_gene ~ normal(0, %f);",
    "  sigma_int  ~ normal(0, %f);",
    "  b_gene ~ normal(0, sigma_gene);",
    "  b_int  ~ normal(0, sigma_int);", sep="\n"),
    sigma_scale, sigma_scale),
  block = "model"
)
all_stanvars <- sv_parameters + sv_model

cat("\nFitting imputation 1 only (reduced budget)...\n")
t0 <- Sys.time()
fit <- brm(
  formula = model_formula,
  data = imputed_data[[1]],
  family = cumulative("logit"),
  prior = model_prior,
  stanvars = all_stanvars,
  chains = 2, cores = 2,
  iter = 1000, warmup = 500,
  init = 0.5,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = 20250310,
  silent = 2, refresh = 0
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Fit completed in %.0f seconds\n", elapsed))

np <- nuts_params(fit)
n_div <- sum(np$Value[np$Parameter == "divergent__"])
sm <- summary(fit)
cat(sprintf("Divergent transitions: %d / 1000\n", n_div))
cat(sprintf("Model parameters: %d\n", nrow(sm$fixed)))

cat("\n[Smoke test passed.]  The full production fit lives in\n")
cat("scripts/52_fit_m4_v3.R (4 chains x 8000 iter, all 100 imputations).\n")
cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

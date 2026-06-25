#!/usr/bin/env Rscript
############################################################
# Step 4: Multiple Imputation with MICE
#
# Reads the fully harmonized dataset (with child + maternal PCs)
# and produces M=100 multiply-imputed datasets ready for brms.
#
# Key design decisions:
#   - Child PC1-PC3 are ANALYSIS variables (included in brms model)
#   - Maternal PC1-PC3 are AUXILIARY variables (used in MICE to help
#     impute missing child PCs, then dropped before brms)
#   - Outcome (Final7PointLikertScale) is NOT imputed but IS used
#     as a predictor in the imputation model
#   - Genetic variants use pmm (predictive mean matching)
#   - Continuous ancestry PCs use pmm
#
# Run from: the repository root.
# Usage:    Rscript scripts/04_run_mice_imputation.R
#
# Input:  data/Fulldata_with_PCs_and_maternal_PCs.txt
# Output: data/mice_imputed_m100_maxit30.rds
#         data/imputed_datasets_for_brms_m100.rds
#         results/mice_full_convergence.pdf
#         results/mice_imputation_report.txt
############################################################

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
  library(purrr)
  library(future)
})

# ── Parallelization settings ────────────────────────────────────────
n_cores <- 18
cat("Using", n_cores, "cores for parallel imputation\n\n")

sink_file <- "results/mice_imputation_report.txt"
sink(sink_file, split = TRUE)

cat("=================================================================\n")
cat("STEP 4: MULTIPLE IMPUTATION WITH MICE\n")
cat("=================================================================\n\n")

# ── Load harmonized data ─────────────────────────────────────────────
input_file <- "data/Fulldata_with_PCs_and_maternal_PCs.txt"
dat <- read.table(input_file, header = TRUE, sep = "\t",
                  check.names = FALSE, encoding = "latin-1")

cat("Input file:", input_file, "\n")
cat("Raw data: n =", nrow(dat), ", p =", ncol(dat), "\n\n")

# ── Clean column names for R compatibility ───────────────────────────
names(dat) <- gsub("-", ".", names(dat))
names(dat) <- make.names(names(dat), unique = TRUE)

# ── Select variables for imputation ──────────────────────────────────
# Analysis variables (will be in brms model)
analysis_vars <- c(
  "Final7PointLikertScale",
  "ethnicity", "Sex", "infant_age",
  "sens", "cont", "unre",
  "DRD2", "BDNF", "CNR1.77", "DRD4", "MAOA",
  "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR", "CNR1.10",
  "PC1", "PC2", "PC3"
)

# Auxiliary variables (help imputation, dropped before brms)
auxiliary_vars <- c("maternal_PC1", "maternal_PC2", "maternal_PC3")

all_vars <- c(analysis_vars, auxiliary_vars)

# Verify all variables exist
missing_vars <- setdiff(all_vars, names(dat))
if (length(missing_vars) > 0) {
  stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
}

dat_analysis <- dat[, all_vars] %>%
  mutate(
    ethnicity = factor(ethnicity),
    Sex = factor(Sex)
  )

cat("Variables selected for imputation:\n")
cat("  Analysis variables:", length(analysis_vars), "\n")
cat("    ", paste(analysis_vars, collapse = ", "), "\n")
cat("  Auxiliary variables:", length(auxiliary_vars), "\n")
cat("    ", paste(auxiliary_vars, collapse = ", "), "\n")
cat("  Total:", length(all_vars), "\n\n")

# ── Missingness summary ──────────────────────────────────────────────
cat("Missingness summary:\n")
miss_counts <- colSums(is.na(dat_analysis))
miss_pct <- round(100 * miss_counts / nrow(dat_analysis), 1)
for (v in all_vars) {
  if (miss_counts[v] > 0) {
    cat(sprintf("  %-30s %3d missing (%4.1f%%)\n",
                v, miss_counts[v], miss_pct[v]))
  }
}
cat(sprintf("  %-30s %3d missing\n", "Total rows with missing outcome",
            sum(is.na(dat_analysis$Final7PointLikertScale))))
cat("\n")

# ── Configure imputation ─────────────────────────────────────────────
init <- mice(dat_analysis, maxit = 0, print = FALSE)
imp_methods <- init$method

# Outcome: do NOT impute (but use as predictor)
imp_methods["Final7PointLikertScale"] <- ""

# Genetic variants: use predictive mean matching
genetic_vars <- c("DRD2", "BDNF", "CNR1.77", "DRD4", "MAOA",
                  "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR", "CNR1.10")
for (var in genetic_vars) {
  if (var %in% names(imp_methods) && imp_methods[var] != "") {
    imp_methods[var] <- "pmm"
  }
}

# Ancestry PCs: use predictive mean matching (continuous)
pc_vars <- c("PC1", "PC2", "PC3",
             "maternal_PC1", "maternal_PC2", "maternal_PC3")
for (var in pc_vars) {
  if (var %in% names(imp_methods) && imp_methods[var] != "") {
    imp_methods[var] <- "pmm"
  }
}

# Predictor matrix: outcome predicts others but is not predicted
pred <- init$predictorMatrix
pred[, "Final7PointLikertScale"] <- 1
pred["Final7PointLikertScale", ] <- 0

cat("Imputation configuration:\n")
cat("  M = 100 imputations\n")
cat("  Cores:", n_cores, "\n")
cat("  maxit = 30 iterations\n")
cat("  Outcome: NOT imputed (used as predictor)\n")
cat("  Genetic variants: pmm\n")
cat("  Ancestry PCs: pmm\n")
cat("  Maternal PCs: pmm (auxiliary)\n\n")

cat("Methods:\n")
for (v in names(imp_methods)) {
  if (imp_methods[v] != "") {
    cat(sprintf("  %-25s %s\n", v, imp_methods[v]))
  }
}
cat("\n")

# ── Run imputation (parallel across 18 cores) ────────────────────────
cat("Starting parallel imputation (M=100, maxit=30, cores=", n_cores, ")...\n")
cat("  Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# futuremice distributes m imputations across n.core workers.
# Each worker runs independent mice() calls with unique seeds.
start_time <- Sys.time()

imp_full <- futuremice(
  data = dat_analysis,
  m = 100,
  method = imp_methods,
  predictorMatrix = pred,
  maxit = 30,
  parallelseed = 20250217,
  n.core = n_cores,
  future.plan = "multisession"
)

end_time <- Sys.time()
duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\nImputation completed!\n")
cat("  End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Duration:", round(duration, 2), "minutes\n")
cat("  Cores used:", n_cores, "\n\n")

# ── Save raw MICE object ─────────────────────────────────────────────
mice_file <- "data/mice_imputed_m100_maxit30.rds"
saveRDS(imp_full, mice_file)
cat("Saved raw MICE object:", mice_file, "\n")

# ── Convergence plots ────────────────────────────────────────────────
conv_file <- "results/mice_full_convergence.pdf"
cat("Creating convergence plots...\n")
pdf(conv_file, width = 14, height = 12)
plot(imp_full, layout = c(4, 4))
dev.off()
cat("Saved:", conv_file, "\n\n")

# ── Prepare datasets for brms ────────────────────────────────────────
cat("Preparing datasets for brms...\n")
imputed_list <- complete(imp_full, action = "all")

imputed_list_brms <- map(imputed_list, function(df) {
  # Create ordered factor outcome
  df$Final7PointLikertScaleFactor <- ordered(
    df$Final7PointLikertScale,
    levels = sort(unique(na.omit(df$Final7PointLikertScale)))
  )

  # Ensure factor types
  df$ethnicity <- factor(df$ethnicity)
  df$Sex <- factor(df$Sex)

  # Scale continuous variables
  df$infant_age <- as.numeric(scale(df$infant_age))
  df$sens <- as.numeric(scale(df$sens))
  df$cont <- as.numeric(scale(df$cont))
  df$unre <- as.numeric(scale(df$unre))
  df$PC1 <- as.numeric(scale(df$PC1))
  df$PC2 <- as.numeric(scale(df$PC2))
  df$PC3 <- as.numeric(scale(df$PC3))

  # Drop auxiliary variables (maternal PCs)
  df$maternal_PC1 <- NULL
  df$maternal_PC2 <- NULL
  df$maternal_PC3 <- NULL

  # Drop rows without outcome
  df <- df[!is.na(df$Final7PointLikertScale), ]

  return(df)
})

sample_sizes <- map_int(imputed_list_brms, nrow)

brms_file <- "data/imputed_datasets_for_brms_m100.rds"
saveRDS(imputed_list_brms, brms_file)
cat("Saved:", brms_file, "\n\n")

# ── Summary ──────────────────────────────────────────────────────────
cat("=================================================================\n")
cat("SUMMARY\n")
cat("=================================================================\n\n")
cat("Imputations:", imp_full$m, "\n")
cat("Iterations per imputation:", imp_full$iteration, "\n")
cat("Time:", round(duration, 2), "minutes\n")
cat("Final sample size (after dropping missing outcome):",
    sample_sizes[1], "\n")
cat("\n")
cat("Variables in brms-ready datasets:\n")
cat("  ", paste(names(imputed_list_brms[[1]]), collapse = ", "), "\n")
cat("\n")
cat("Scaled continuous variables: infant_age, sens, cont, unre,",
    "PC1, PC2, PC3\n")
cat("Auxiliary variables (maternal_PC1-3) have been DROPPED from",
    "brms datasets.\n")
cat("\n")
cat("Ready for brms!\n")
cat("=================================================================\n")

sink()
cat("Report saved to:", sink_file, "\n")

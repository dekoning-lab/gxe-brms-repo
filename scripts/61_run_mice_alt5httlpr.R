#!/usr/bin/env Rscript
############################################################
# 61_run_mice_alt5httlpr.R
#
# ALTERNATE-CODING RUN of the production MICE imputation pipeline,
# using the Hu et al. (2006) triallelic-functional coding of 5-HTTLPR
# instead of the original biallelic S-count coding.
#
# This is a near-copy of 04_run_mice_imputation.R with the following
# changes only:
#   - input_file  changed to data_alt_5httlpr/...alt5httlpr.txt
#   - mice_file   changed to data_alt_5httlpr/mice_imputed_m100_maxit30.rds
#   - brms_file   changed to data_alt_5httlpr/imputed_datasets_for_brms_m100.rds
#   - report file changed to results_alt_5httlpr/mice_imputation_report.txt
#   - conv_file   changed to results_alt_5httlpr/mice_full_convergence.pdf
# Same MICE seed (parallelseed = 20250217), same n.core, same maxit,
# same predictor matrix, same methods. Goal: maximise comparability
# with the production run; only the input column values change.
#
# Run from: the repository root (locally, or on a cluster after sourcing config.sh; see README).
# Usage:    Rscript scripts/61_run_mice_alt5httlpr.R
############################################################

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
  library(purrr)
  library(future)
})

n_cores <- 18
cat("Using", n_cores, "cores for parallel imputation\n\n")

dir.create("results_alt_5httlpr", recursive = TRUE, showWarnings = FALSE)
dir.create("data_alt_5httlpr", recursive = TRUE, showWarnings = FALSE)

sink_file <- "results_alt_5httlpr/mice_imputation_report.txt"
sink(sink_file, split = TRUE)

cat("=================================================================\n")
cat("STEP 61: MICE on ALTERNATE 5-HTTLPR CODING (Hu et al. 2006)\n")
cat("=================================================================\n\n")

input_file <- "data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt"
dat <- read.table(input_file, header = TRUE, sep = "\t",
                  check.names = FALSE, encoding = "latin-1")

cat("Input file:", input_file, "\n")
cat("Raw data: n =", nrow(dat), ", p =", ncol(dat), "\n\n")

names(dat) <- gsub("-", ".", names(dat))
names(dat) <- make.names(names(dat), unique = TRUE)

analysis_vars <- c(
  "Final7PointLikertScale",
  "ethnicity", "Sex", "infant_age",
  "sens", "cont", "unre",
  "DRD2", "BDNF", "CNR1.77", "DRD4", "MAOA",
  "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR", "CNR1.10",
  "PC1", "PC2", "PC3"
)
auxiliary_vars <- c("maternal_PC1", "maternal_PC2", "maternal_PC3")
all_vars <- c(analysis_vars, auxiliary_vars)

missing_vars <- setdiff(all_vars, names(dat))
if (length(missing_vars) > 0) {
  stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
}

dat_analysis <- dat[, all_vars] %>%
  mutate(ethnicity = factor(ethnicity), Sex = factor(Sex))

cat("Variables selected for imputation:\n")
cat("  Analysis variables:", length(analysis_vars), "\n")
cat("    ", paste(analysis_vars, collapse = ", "), "\n")
cat("  Auxiliary variables:", length(auxiliary_vars), "\n")
cat("    ", paste(auxiliary_vars, collapse = ", "), "\n")
cat("  Total:", length(all_vars), "\n\n")

cat("Missingness summary:\n")
miss_counts <- colSums(is.na(dat_analysis))
miss_pct <- round(100 * miss_counts / nrow(dat_analysis), 1)
for (v in all_vars) {
  if (miss_counts[v] > 0) {
    cat(sprintf("  %-30s %3d missing (%4.1f%%)\n",
                v, miss_counts[v], miss_pct[v]))
  }
}
cat(sprintf("  %-30s %3d missing\n",
            "Total rows with missing outcome",
            sum(is.na(dat_analysis$Final7PointLikertScale))))
cat("\nX5HTTLPR distribution (post-recoding):\n")
print(table(dat_analysis$X5HTTLPR, useNA = "ifany"))

init <- mice(dat_analysis, maxit = 0, print = FALSE)
imp_methods <- init$method
imp_methods["Final7PointLikertScale"] <- ""

genetic_vars <- c("DRD2", "BDNF", "CNR1.77", "DRD4", "MAOA",
                  "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR", "CNR1.10")
for (var in genetic_vars) {
  if (var %in% names(imp_methods) && imp_methods[var] != "") {
    imp_methods[var] <- "pmm"
  }
}
pc_vars <- c("PC1", "PC2", "PC3",
             "maternal_PC1", "maternal_PC2", "maternal_PC3")
for (var in pc_vars) {
  if (var %in% names(imp_methods) && imp_methods[var] != "") {
    imp_methods[var] <- "pmm"
  }
}

pred <- init$predictorMatrix
pred[, "Final7PointLikertScale"] <- 1
pred["Final7PointLikertScale", ] <- 0

cat("\nImputation configuration:\n")
cat("  M = 100 imputations\n")
cat("  Cores:", n_cores, "\n")
cat("  maxit = 30\n")
cat("  Seed (parallelseed): 20250217  (same as production for comparability)\n\n")

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
cat(sprintf("\nImputation completed in %.2f min\n", duration))

mice_file <- "data_alt_5httlpr/mice_imputed_m100_maxit30.rds"
saveRDS(imp_full, mice_file)
cat("Saved raw MICE object:", mice_file, "\n")

conv_file <- "results_alt_5httlpr/mice_full_convergence.pdf"
pdf(conv_file, width = 14, height = 12)
plot(imp_full, layout = c(4, 4))
dev.off()
cat("Saved:", conv_file, "\n\n")

imputed_list <- complete(imp_full, action = "all")
imputed_list_brms <- map(imputed_list, function(df) {
  df$Final7PointLikertScaleFactor <- ordered(
    df$Final7PointLikertScale,
    levels = sort(unique(na.omit(df$Final7PointLikertScale)))
  )
  df$ethnicity <- factor(df$ethnicity)
  df$Sex <- factor(df$Sex)
  df$infant_age <- as.numeric(scale(df$infant_age))
  df$sens <- as.numeric(scale(df$sens))
  df$cont <- as.numeric(scale(df$cont))
  df$unre <- as.numeric(scale(df$unre))
  df$PC1 <- as.numeric(scale(df$PC1))
  df$PC2 <- as.numeric(scale(df$PC2))
  df$PC3 <- as.numeric(scale(df$PC3))
  df$maternal_PC1 <- NULL
  df$maternal_PC2 <- NULL
  df$maternal_PC3 <- NULL
  df <- df[!is.na(df$Final7PointLikertScale), ]
  return(df)
})

brms_file <- "data_alt_5httlpr/imputed_datasets_for_brms_m100.rds"
saveRDS(imputed_list_brms, brms_file)
cat("Saved:", brms_file, "\n\n")

cat("X5HTTLPR distribution in first imputed dataset:\n")
print(table(imputed_list_brms[[1]]$X5HTTLPR))
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()

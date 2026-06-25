#!/usr/bin/env Rscript
############################################################
# 52_compute_rpss_v3.R
# Compute Ranked Probability Skill Score (RPSS) for v3 models
#
# RPSS = 1 - RPS_model / RPS_reference
# Reference = "climatological" forecast using marginal category
# frequencies (predicting base rates for everyone).
#
# Run from: the repository root.
# Usage: Rscript scripts/52_compute_rpss_v3.R
############################################################

cat("=================================================================\n")
cat("V3 RANKED PROBABILITY SKILL SCORE (RPSS)\n")
cat("=================================================================\n\n")

# ── Load RPS results ──────────────────────────────────────
rps_all <- read.csv("results/v3_loo_rps/all_rps_combined.csv",
                    stringsAsFactors = FALSE)
cat(sprintf("Loaded %d RPS results\n", nrow(rps_all)))

# ── Load imputed data to compute base-rate RPS ───────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
cat(sprintf("Loaded %d imputed datasets\n\n", length(imputed_data)))

K <- 7  # 7-point Likert scale

# ── Compute reference RPS per imputation ──────────────────
# The reference forecast is the marginal category distribution
# (i.e., predict p_k = observed frequency of category k for everyone)
cat("Computing reference (base-rate) RPS per imputation...\n")

ref_rps <- data.frame()

for (imp in 1:100) {
  dat <- imputed_data[[imp]]
  if (!"Final7PointLikertScaleFactor" %in% names(dat)) {
    dat$Final7PointLikertScaleFactor <- factor(dat$Final7PointLikertScale)
  }
  dat <- dat[!is.na(dat$Final7PointLikertScaleFactor), ]
  y_num <- as.numeric(dat$Final7PointLikertScaleFactor)
  n <- length(y_num)

  # Marginal category frequencies
  p_marginal <- tabulate(y_num, nbins = K) / n

  # RPS for each observation under the marginal forecast
  rps_obs <- numeric(n)
  F_marginal <- cumsum(p_marginal)

  for (i in 1:n) {
    obs_cum <- as.numeric(y_num[i] <= seq_len(K))
    rps_obs[i] <- sum((F_marginal - obs_cum)^2) / (K - 1)
  }

  ref_rps <- rbind(ref_rps, data.frame(
    imputation = imp,
    ref_mean_rps = mean(rps_obs),
    ref_total_rps = sum(rps_obs),
    n = n,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("  Reference RPS: mean = %.4f (SD = %.4f) across %d imputations\n",
            mean(ref_rps$ref_mean_rps), sd(ref_rps$ref_mean_rps), nrow(ref_rps)))

# ── Compute RPSS for each model × parenting × imputation ─
cat("\nComputing RPSS...\n\n")

# Merge RPS results with reference
rpss_data <- merge(rps_all, ref_rps, by = "imputation")

# RPSS = 1 - RPS_model / RPS_reference
rpss_data$rpss <- 1 - rpss_data$mean_rps / rpss_data$ref_mean_rps

# ── Summarize ────────────────────────────────────────────
cat("=== RPSS SUMMARY (averaged across imputations) ===\n\n")

parenting_label <- c(sens = "Sensitivity", cont = "Controlling", unre = "Unresponsiveness")

cat(sprintf("%-8s %-8s %8s %8s %8s %10s %10s\n",
            "Parent.", "Model", "RPS", "Ref RPS", "RPSS", "RPSS 2.5%", "RPSS 97.5%"))
cat(paste(rep("-", 75), collapse = ""), "\n")

rpss_summary <- list()

for (pv in c("sens", "cont", "unre")) {
  for (m in c("M1", "M2", "M3", "M4")) {
    sub <- rpss_data[rpss_data$parenting == pv & rpss_data$model == m, ]
    if (nrow(sub) == 0) next

    rpss_summary[[length(rpss_summary) + 1]] <- data.frame(
      parenting = pv,
      parenting_label = parenting_label[pv],
      model = m,
      mean_rps = mean(sub$mean_rps),
      mean_ref_rps = mean(sub$ref_mean_rps),
      mean_rpss = mean(sub$rpss),
      sd_rpss = sd(sub$rpss),
      q025_rpss = quantile(sub$rpss, 0.025),
      q975_rpss = quantile(sub$rpss, 0.975),
      n_imps = nrow(sub),
      stringsAsFactors = FALSE
    )

    cat(sprintf("%-8s %-8s %8.4f %8.4f %8.4f %10.4f %10.4f\n",
                pv, m, mean(sub$mean_rps), mean(sub$ref_mean_rps),
                mean(sub$rpss),
                quantile(sub$rpss, 0.025), quantile(sub$rpss, 0.975)))
  }
  cat("\n")
}

rpss_df <- do.call(rbind, rpss_summary)

# ── Pairwise RPSS differences ─────────────────────────────
cat("=== RPSS DIFFERENCES (does adding predictors improve skill?) ===\n\n")

for (pv in c("sens", "cont", "unre")) {
  cat(sprintf("--- %s ---\n", parenting_label[pv]))

  for (pair in list(c("M1", "M2"), c("M1", "M3"), c("M2", "M4"), c("M1", "M4"))) {
    base <- rpss_data[rpss_data$parenting == pv & rpss_data$model == pair[1], ]
    comp <- rpss_data[rpss_data$parenting == pv & rpss_data$model == pair[2], ]
    shared <- intersect(base$imputation, comp$imputation)

    if (length(shared) > 0) {
      base_rpss <- base$rpss[match(shared, base$imputation)]
      comp_rpss <- comp$rpss[match(shared, comp$imputation)]
      d <- comp_rpss - base_rpss

      cat(sprintf("  %s vs %s: dRPSS = %+.4f (SD %.4f) [%.0f%% complex better]\n",
                  pair[2], pair[1], mean(d), sd(d), mean(d > 0) * 100))
    }
  }
  cat("\n")
}

# ── Interpretation ────────────────────────────────────────
cat("=== INTERPRETATION ===\n\n")
cat("RPSS ranges from -Inf to 1:\n")
cat("  RPSS = 0:  model is no better than predicting base rates\n")
cat("  RPSS > 0:  model has skill beyond base rates\n")
cat("  RPSS < 0:  model is worse than base rates (overfitting)\n")
cat("  RPSS = 1:  perfect prediction\n\n")

overall_m1 <- rpss_df[rpss_df$model == "M1", ]
cat(sprintf("M1 (demographics + child behavior): RPSS = %.4f [%.4f, %.4f]\n",
            mean(overall_m1$mean_rpss),
            mean(overall_m1$q025_rpss),
            mean(overall_m1$q975_rpss)))
cat("This means M1 reduces the RPS by ~X%% relative to predicting base rates.\n\n")

# ── Save ──────────────────────────────────────────────────
dir.create("results/v3_loo_rps", recursive = TRUE, showWarnings = FALSE)
write.csv(rpss_df, "results/v3_loo_rps/rpss_summary.csv", row.names = FALSE)
write.csv(rpss_data[, c("model", "parenting", "imputation", "mean_rps",
                         "ref_mean_rps", "rpss")],
          "results/v3_loo_rps/rpss_all_imputations.csv", row.names = FALSE)
cat("Saved: results/v3_loo_rps/rpss_summary.csv\n")
cat("Saved: results/v3_loo_rps/rpss_all_imputations.csv\n")

cat("\n=== DONE ===\n")

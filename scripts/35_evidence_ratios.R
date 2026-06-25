#!/usr/bin/env Rscript
############################################################
# 35_evidence_ratios.R
# Evidence calibration for hierarchical model results
#
# Computes:
#   1. Directional Bayes factors from posterior probabilities
#   2. ROPE (Region of Practical Equivalence) analysis
#   3. Evidence summary table + figure
#
# Key insight: Under symmetric Normal(0, sigma) prior,
#   prior P(beta > 0) = 0.5, so directional BF = PP/(1-PP)
#
# Note: Savage-Dickey point-null BF (H0: beta = 0) is
#   ill-defined here because the marginal prior density at 0
#   is infinite when sigma ~ half-Normal (density ~ 1/sigma
#   as sigma -> 0). We therefore use directional BFs only.
#
# Run from: the repository root.
# Usage:    Rscript scripts/35_evidence_ratios.R
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

cat("=== EVIDENCE CALIBRATION ANALYSIS ===\n\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ── Load primary results (v2, s=0.50) ────────────────────
all_results <- list()
for (p in c("sens", "cont", "unre")) {
  f <- file.path("results", "hierarchical_v2", p, "summary.csv")
  if (file.exists(f)) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    all_results[[p]] <- d
  }
}

df <- do.call(rbind, all_results)
rownames(df) <- NULL

parenting_label <- c(sens = "Sensitivity", cont = "Controlling", unre = "Unresponsiveness")

cat(sprintf("Loaded %d rows across %d parenting models\n\n",
            nrow(df), length(all_results)))

# ── 1. Directional Bayes Factors ─────────────────────────
# Under symmetric prior centered at 0:
#   Prior P(beta > 0) = P(beta < 0) = 0.5 (prior odds = 1)
#   Posterior odds = PP / (1 - PP)
#   Directional BF = posterior odds / prior odds = PP / (1 - PP)
#
# This is the BF for H1: "beta is in the observed direction"
# vs H0: "beta is in the opposite direction"

cat("=", rep("=", 60), "\n")
cat("1. DIRECTIONAL BAYES FACTORS\n")
cat("=", rep("=", 60), "\n\n")

effects <- df %>%
  filter(effect_type %in% c("gene_main", "interaction"))

effects$bf_directional <- effects$post_prob / (1 - effects$post_prob)

# Kass & Raftery (1995) scale based on 2*ln(BF):
#   0-2:  Not worth more than a bare mention  (BF 1.0 - 2.7)
#   2-6:  Positive evidence                    (BF 2.7 - 20.1)
#   6-10: Strong evidence                      (BF 20.1 - 148.4)
#   >10:  Very strong evidence                 (BF > 148.4)
classify_kr <- function(bf) {
  twoln <- 2 * log(bf)
  ifelse(twoln > 10, "Very strong",
  ifelse(twoln > 6,  "Strong",
  ifelse(twoln > 2,  "Positive",
                      "Bare mention")))
}

effects$kr_evidence <- classify_kr(effects$bf_directional)

# Print summary for gene + interaction effects
cat("Kass & Raftery (1995) evidence categories [2 ln(BF)]:\n")
cat("  0-2:  Not worth more than a bare mention (BF 1.0 - 2.7)\n")
cat("  2-6:  Positive evidence                   (BF 2.7 - 20.1)\n")
cat("  6-10: Strong evidence                     (BF 20.1 - 148)\n")
cat("  >10:  Very strong evidence                (BF > 148)\n\n")

cat("PP thresholds and corresponding directional BFs:\n")
cat("  PP = 60%  -> BF = 1.50, 2ln(BF) = 0.81 [Bare mention]\n")
cat("  PP = 70%  -> BF = 2.33, 2ln(BF) = 1.69 [Bare mention]\n")
cat("  PP = 73%  -> BF = 2.70, 2ln(BF) = 1.99 [~Positive threshold]\n")
cat("  PP = 80%  -> BF = 4.00, 2ln(BF) = 2.77 [Positive]\n")
cat("  PP = 90%  -> BF = 9.00, 2ln(BF) = 4.39 [Positive]\n")
cat("  PP = 95%  -> BF = 19.0, 2ln(BF) = 5.89 [Positive]\n")
cat("  PP = 97%  -> BF = 32.3, 2ln(BF) = 6.95 [Strong]\n\n")

# Show key effects
key_effects <- effects %>%
  filter(post_prob >= 0.65) %>%
  arrange(desc(bf_directional))

if (nrow(key_effects) > 0) {
  cat(sprintf("%-10s %-22s %7s %7s %8s  %-14s\n",
              "Model", "Term", "PP(%)", "BF", "2ln(BF)", "Evidence"))
  cat(paste(rep("-", 78), collapse = ""), "\n")
  for (i in seq_len(nrow(key_effects))) {
    r <- key_effects[i, ]
    cat(sprintf("%-10s %-22s %6.1f %7.2f %8.2f  %-14s\n",
                parenting_label[r$parenting], r$term,
                r$post_prob * 100, r$bf_directional,
                2 * log(r$bf_directional), r$kr_evidence))
  }
}

# ── 2. ROPE Analysis ─────────────────────────────────────
# Region of Practical Equivalence: what proportion of the
# posterior falls outside a region near zero?
#
# Using Normal(estimate, sd) approximation.
# ROPE = [-delta, delta] on log-odds scale.

cat("\n\n")
cat("=", rep("=", 60), "\n")
cat("2. ROPE (REGION OF PRACTICAL EQUIVALENCE) ANALYSIS\n")
cat("=", rep("=", 60), "\n\n")

# For ordinal logit with 7 categories:
#   0.10 log-odds ~ 2-3 percentage point shift in cumulative prob
#   0.05 log-odds ~ 1-1.5 percentage point shift
rope_delta <- 0.10

effects$p_in_rope <- pnorm(rope_delta, effects$estimate, effects$sd) -
                     pnorm(-rope_delta, effects$estimate, effects$sd)
effects$p_outside_rope <- 1 - effects$p_in_rope

cat(sprintf("ROPE = [%.2f, +%.2f] log-odds\n", -rope_delta, rope_delta))
cat("(~2-3 percentage point shift in cumulative probability per category)\n\n")

rope_key <- effects %>%
  filter(post_prob >= 0.65) %>%
  arrange(desc(p_outside_rope))

cat(sprintf("%-10s %-22s %7s %8s %10s\n",
            "Model", "Term", "PP(%)", "BF", "P(|β|>0.1)"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_len(nrow(rope_key))) {
  r <- rope_key[i, ]
  cat(sprintf("%-10s %-22s %6.1f %8.2f %9.1f%%\n",
              parenting_label[r$parenting], r$term,
              r$post_prob * 100, r$bf_directional, r$p_outside_rope * 100))
}

cat("\nNote: ROPE probability tells us something different from PP.\n")
cat("PP asks: which side of zero? ROPE asks: is the effect non-negligible?\n")
cat("An effect can have high PP (consistently one-sided) but low ROPE\n")
cat("probability (the effect is small and much of the posterior is near zero).\n")

# ── 3. Pooled gene main effects ──────────────────────────
# Since gene mains are nearly identical across models, pool them
cat("\n\n")
cat("=", rep("=", 60), "\n")
cat("3. POOLED GENE MAIN EFFECT EVIDENCE\n")
cat("=", rep("=", 60), "\n\n")

gene_df <- effects %>% filter(effect_type == "gene_main")

pooled_genes <- gene_df %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    sd = mean(sd),
    q025 = mean(q025),
    q975 = mean(q975),
    post_prob = mean(post_prob),
    direction = first(direction),
    .groups = "drop"
  ) %>%
  mutate(
    bf_directional = post_prob / (1 - post_prob),
    kr_evidence = classify_kr(bf_directional),
    p_in_rope = pnorm(rope_delta, estimate, sd) - pnorm(-rope_delta, estimate, sd),
    p_outside_rope = 1 - p_in_rope
  ) %>%
  arrange(desc(bf_directional))

cat(sprintf("%-14s %7s %7s %7s %8s  %-14s %10s\n",
            "Gene", "Est", "PP(%)", "BF", "2ln(BF)", "Evidence", "P(|β|>0.1)"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in seq_len(nrow(pooled_genes))) {
  r <- pooled_genes[i, ]
  cat(sprintf("%-14s %7.3f %6.1f %7.2f %8.2f  %-14s %9.1f%%\n",
              r$term, r$estimate, r$post_prob * 100, r$bf_directional,
              2 * log(r$bf_directional), r$kr_evidence, r$p_outside_rope * 100))
}

# ── 4. Full evidence summary table (for export) ──────────
cat("\n\n")
cat("=", rep("=", 60), "\n")
cat("4. COMPLETE EVIDENCE SUMMARY\n")
cat("=", rep("=", 60), "\n\n")

# Build a tidy summary for all effects
evidence_summary <- effects %>%
  mutate(
    model = parenting_label[parenting],
    p_in_rope = pnorm(rope_delta, estimate, sd) - pnorm(-rope_delta, estimate, sd),
    p_outside_rope = 1 - p_in_rope,
    twoln_bf = 2 * log(bf_directional)
  ) %>%
  select(model, term, effect_type, estimate, sd, q025, q975,
         post_prob, direction, bf_directional, twoln_bf, kr_evidence,
         p_outside_rope) %>%
  arrange(effect_type, desc(bf_directional))

write.csv(evidence_summary, "results/figures/evidence_summary.csv", row.names = FALSE)
cat("Saved: results/figures/evidence_summary.csv\n")

# ── 5. Evidence calibration figure ────────────────────────
cat("\nGenerating evidence calibration figure...\n")

# Figure: BF scale with our effects mapped onto it
# Shows where each effect falls on the evidence scale

# All effects with PP > 55% (otherwise it's just noise)
plot_effects <- effects %>%
  filter(post_prob >= 0.55) %>%
  mutate(
    model = parenting_label[parenting],
    gene = sub("^(sens|cont|unre):", "", term),
    twoln_bf = 2 * log(bf_directional),
    eff_label = ifelse(effect_type == "interaction",
                       paste0(model, ": P × ", gene),
                       paste0(gene, " (pooled)")),
    kr_cat = factor(kr_evidence,
                    levels = c("Bare mention", "Positive", "Strong", "Very strong"))
  )

# For gene mains, show pooled version only
gene_pooled_plot <- pooled_genes %>%
  filter(post_prob >= 0.55) %>%
  mutate(
    model = "Pooled",
    gene = term,
    twoln_bf = 2 * log(bf_directional),
    eff_label = paste0(gene, " (gene main)"),
    effect_type = "gene_main",
    kr_cat = factor(kr_evidence,
                    levels = c("Bare mention", "Positive", "Strong", "Very strong"))
  )

int_plot <- plot_effects %>%
  filter(effect_type == "interaction") %>%
  mutate(eff_label = paste0(sub(" \\(pooled\\)", "", model), ": P × ", gene))

combined_plot <- bind_rows(
  gene_pooled_plot %>% select(eff_label, post_prob, bf_directional, twoln_bf, kr_cat, effect_type),
  int_plot %>% select(eff_label, post_prob, bf_directional, twoln_bf, kr_cat, effect_type)
) %>%
  arrange(twoln_bf) %>%
  mutate(eff_label = factor(eff_label, levels = eff_label))

p_bf <- ggplot(combined_plot, aes(x = twoln_bf, y = eff_label, colour = kr_cat)) +
  # Evidence threshold bands
  annotate("rect", xmin = -Inf, xmax = 2, ymin = -Inf, ymax = Inf,
           fill = "#fee0d2", alpha = 0.3) +
  annotate("rect", xmin = 2, xmax = 6, ymin = -Inf, ymax = Inf,
           fill = "#deebf7", alpha = 0.3) +
  annotate("rect", xmin = 6, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#e5f5e0", alpha = 0.3) +
  # Threshold lines
  geom_vline(xintercept = 2, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  # Points
  geom_point(aes(shape = effect_type), size = 3) +
  # Labels at top
  annotate("text", x = 1, y = Inf, label = "Bare mention", vjust = -0.5,
           size = 3, colour = "grey40") +
  annotate("text", x = 4, y = Inf, label = "Positive", vjust = -0.5,
           size = 3, colour = "grey40") +
  annotate("text", x = 7, y = Inf, label = "Strong", vjust = -0.5,
           size = 3, colour = "grey40") +
  # Scale
  scale_colour_manual(
    values = c("Bare mention" = "#d95f02", "Positive" = "#1b9e77",
               "Strong" = "#7570b3", "Very strong" = "#e7298a"),
    name = "Kass & Raftery\nevidence"
  ) +
  scale_shape_manual(
    values = c("gene_main" = 16, "interaction" = 17),
    labels = c("Gene main", "Interaction"),
    name = "Effect type"
  ) +
  labs(
    title = "Evidence calibration: Directional Bayes factors",
    subtitle = "Under symmetric hierarchical prior, BF = PP / (1 - PP)",
    x = expression(2 ~ ln(BF[directional])),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/evidence_calibration.pdf", p_bf,
       width = 9, height = 7)
cat("  Saved: results/figures/evidence_calibration.pdf\n")

# ── 6. PP-BF conversion reference table ──────────────────
cat("\n")
cat("=", rep("=", 60), "\n")
cat("5. CALIBRATION SUMMARY\n")
cat("=", rep("=", 60), "\n\n")

cat("For this hierarchical model with symmetric Normal(0, sigma) priors:\n")
cat("  - Prior PP = 50% (by symmetry)\n")
cat("  - Directional BF = PP / (1 - PP)\n")
cat("  - The Savage-Dickey point-null BF is undefined (marginal prior\n")
cat("    density at beta=0 is infinite due to half-Normal hyperprior on sigma)\n\n")

cat("Our strongest effects:\n\n")

top3 <- effects %>%
  arrange(desc(bf_directional)) %>%
  head(5)

for (i in seq_len(nrow(top3))) {
  r <- top3[i, ]
  cat(sprintf("  %s in %s model:\n", r$term, parenting_label[r$parenting]))
  cat(sprintf("    PP = %.1f%% -> BF = %.2f -> 2ln(BF) = %.2f -> %s\n",
              r$post_prob * 100, r$bf_directional,
              2 * log(r$bf_directional), r$kr_evidence))
  cat(sprintf("    P(|beta| > 0.10) = %.1f%% (ROPE)\n\n",
              r$p_outside_rope * 100))
}

cat("\nConclusion: The strongest effects in this analysis (PP ~ 80%)\n")
cat("correspond to directional BF ~ 4, which Kass & Raftery classify as\n")
cat("'positive' evidence — meaningful but not strong (BF > 20) or\n")
cat("decisive (BF > 150). This is consistent with a weakly informative\n")
cat("sample (n = 168) for detecting gene-environment interactions.\n")
cat("The effects that survived hierarchical shrinkage carry real signal\n")
cat("but would benefit from replication in a larger sample.\n")

cat("\n=== DONE ===\n")

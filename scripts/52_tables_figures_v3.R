#!/usr/bin/env Rscript
############################################################
# 52_tables_figures_v3.R
# Tables and figures for v3 model comparison
# (pass + diff as child behavior covariates in all models)
#
# Generates:
#   1. Table 1: Coefficient estimates (3 parenting models)
#   2. Table 2: LOO pairwise comparisons
#   3. Figure 1: Forest plot (gene mains + interactions)
#   4. Figure 2: Evidence calibration (BF scale)
#   5. Figure 3: Posterior densities for PP >= 70%
#   6. Figure 4: LOO delta ELPD distribution
#
# Run from: the repository root.
# Usage:    Rscript scripts/52_tables_figures_v3.R
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

cat("=================================================================\n")
cat("V3 TABLES AND FIGURES\n")
cat("(pass + diff as child behavior covariates in all models)\n")
cat("=================================================================\n\n")

dir.create("results/latex_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v3_figures", recursive = TRUE, showWarnings = FALSE)

parenting_vars  <- c("sens", "cont", "unre")
parenting_label <- c(sens = "Sensitivity", cont = "Controlling",
                     unre = "Unresponsiveness")

# ── Load M4 summary CSVs ──────────────────────────────────────────
all_results <- list()
for (pv in parenting_vars) {
  f <- file.path("results", "v3_m4", pv, "summary.csv")
  if (file.exists(f)) {
    all_results[[pv]] <- read.csv(f, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded %s: %d rows\n", pv, nrow(all_results[[pv]])))
  } else {
    cat(sprintf("  WARNING: %s not found\n", f))
  }
}

df <- do.call(rbind, all_results)
rownames(df) <- NULL

effects <- df %>% filter(effect_type != "hyperparameter")
hyper   <- df %>% filter(effect_type == "hyperparameter")

cat(sprintf("\nTotal: %d effect rows, %d hyperparameter rows\n\n",
            nrow(effects), nrow(hyper)))

# ── Load LOO data ──────────────────────────────────────────────────
loo_all <- read.csv("results/v3_all_loo.csv", stringsAsFactors = FALSE)
loo_pairwise <- read.csv("results/v3_pairwise.csv", stringsAsFactors = FALSE)
cat(sprintf("LOO: %d observations, %d pairwise comparisons\n\n",
            nrow(loo_all), nrow(loo_pairwise)))

# ── Helpers ────────────────────────────────────────────────────────
fmt_term <- function(t) {
  t <- gsub("_", "\\_", t, fixed = TRUE)
  t <- gsub(":", "$\\times$", t, fixed = TRUE)
  t
}

fmt_cell <- function(est, prob, threshold = 0.70) {
  if (is.na(est) || is.na(prob)) return("---")
  e <- sprintf("%.2f", est)
  p <- round(prob * 100)
  # Compare rounded value so display and bolding are consistent
  if (p >= threshold * 100)
    sprintf("\\textbf{%s} (\\textbf{%d})", e, p)
  else
    sprintf("%s (%d)", e, p)
}

fmt_hyper_cell <- function(est, q025, q975) {
  if (is.na(est)) return("---")
  sprintf("%.3f $[%.3f,\\, %.3f]$", est, q025, q975)
}

# Gene names
genes_unsorted <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
                    "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

# Sort genes by max |estimate| across models (descending)
gene_abs <- effects %>%
  filter(effect_type == "gene_main") %>%
  group_by(term) %>%
  summarise(max_abs = max(abs(estimate)), .groups = "drop") %>%
  arrange(desc(max_abs))
genes <- gene_abs$term

# Sort interactions by max |estimate| across models (descending)
int_abs <- effects %>%
  filter(effect_type == "interaction") %>%
  mutate(gene_name = sub("^(sens|cont|unre):", "", term)) %>%
  group_by(gene_name) %>%
  summarise(max_abs = max(abs(estimate)), .groups = "drop") %>%
  arrange(desc(max_abs))
int_gene_order <- int_abs$gene_name

# Covariates: demographics + child behavior + parenting
# In v3, pass and diff are covariates (not parenting variables)
common_covs_unsorted <- c("PC1", "PC2", "PC3", "Sex0", "Sex1",
                          "infant_age", "pass", "diff")

cov_abs <- effects %>%
  filter(effect_type == "covariate", term %in% common_covs_unsorted) %>%
  group_by(term) %>%
  summarise(max_abs = max(abs(estimate)), .groups = "drop") %>%
  arrange(desc(max_abs))
common_covs <- cov_abs$term

# ══════════════════════════════════════════════════════════════════
# TABLE 1: Combined Main Results
# ══════════════════════════════════════════════════════════════════
cat("Creating Table 1: Combined main results...\n")

make_lookup <- function(pv) {
  d <- effects %>% filter(parenting == pv)
  list(
    cov  = d %>% filter(effect_type == "covariate"),
    gene = d %>% filter(effect_type == "gene_main"),
    int  = d %>% filter(effect_type == "interaction")
  )
}

lookups <- lapply(setNames(parenting_vars, parenting_vars), make_lookup)

hyper_lookup <- function(pv, param) {
  r <- hyper %>% filter(parenting == pv, term == param)
  if (nrow(r) > 0) r[1, ] else NULL
}

table1_rows <- character()

# ─── Panel A: Demographics ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{A. Demographic Covariates}}} \\\\",
  "\\addlinespace[3pt]"
)

demo_covs <- intersect(common_covs, c("PC1", "PC2", "PC3", "Sex0", "Sex1", "infant_age"))
for (cv in demo_covs) {
  cells <- sapply(parenting_vars, function(pv) {
    r <- lookups[[pv]]$cov %>% filter(term == cv)
    if (nrow(r) > 0) fmt_cell(r$estimate[1], r$post_prob[1]) else "---"
  })
  table1_rows <- c(table1_rows,
    sprintf("%s & %s & %s & %s \\\\", fmt_term(cv), cells[1], cells[2], cells[3]))
}

table1_rows <- c(table1_rows, "\\addlinespace[6pt]")

# ─── Panel A2: Child Behavior Covariates ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{B. Child Behavior Covariates (CARE Index)}}} \\\\",
  "\\addlinespace[3pt]"
)

child_covs <- intersect(common_covs, c("pass", "diff"))
child_labels <- c(pass = "Passivity", diff = "Difficultness")
for (cv in child_covs) {
  cells <- sapply(parenting_vars, function(pv) {
    r <- lookups[[pv]]$cov %>% filter(term == cv)
    if (nrow(r) > 0) fmt_cell(r$estimate[1], r$post_prob[1]) else "---"
  })
  lab <- if (cv %in% names(child_labels)) child_labels[cv] else fmt_term(cv)
  table1_rows <- c(table1_rows,
    sprintf("%s & %s & %s & %s \\\\", lab, cells[1], cells[2], cells[3]))
}

table1_rows <- c(table1_rows, "\\addlinespace[6pt]")

# ─── Panel C: Parenting Variable (model-specific) ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{C. Parenting Variable (model-specific)}}} \\\\",
  "\\addlinespace[3pt]"
)

pcells <- sapply(parenting_vars, function(pv) {
  r <- lookups[[pv]]$cov %>% filter(term == pv)
  if (nrow(r) > 0) fmt_cell(r$estimate[1], r$post_prob[1]) else "---"
})
table1_rows <- c(table1_rows,
  sprintf("Parenting & %s & %s & %s \\\\", pcells[1], pcells[2], pcells[3]))

table1_rows <- c(table1_rows, "\\addlinespace[6pt]")

# ─── Panel D: Gene Main Effects ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{D. Gene Main Effects}}} \\\\",
  "\\addlinespace[3pt]"
)

for (g in genes) {
  cells <- sapply(parenting_vars, function(pv) {
    r <- lookups[[pv]]$gene %>% filter(term == g)
    if (nrow(r) > 0) fmt_cell(r$estimate[1], r$post_prob[1]) else "---"
  })
  table1_rows <- c(table1_rows,
    sprintf("%s & %s & %s & %s \\\\", fmt_term(g), cells[1], cells[2], cells[3]))
}

table1_rows <- c(table1_rows, "\\addlinespace[6pt]")

# ─── Panel D: Gene × Parenting Interactions ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{E. Gene $\\times$ Parenting Interactions}}} \\\\",
  "\\addlinespace[3pt]"
)

for (g in int_gene_order) {
  cells <- sapply(parenting_vars, function(pv) {
    int_term <- paste0(pv, ":", g)
    r <- lookups[[pv]]$int %>% filter(term == int_term)
    if (nrow(r) > 0) fmt_cell(r$estimate[1], r$post_prob[1]) else "---"
  })
  table1_rows <- c(table1_rows,
    sprintf("P $\\times$ %s & %s & %s & %s \\\\", fmt_term(g),
            cells[1], cells[2], cells[3]))
}

table1_rows <- c(table1_rows, "\\addlinespace[6pt]")

# ─── Panel E: Hyperparameters ───
table1_rows <- c(table1_rows,
  "\\multicolumn{4}{l}{\\textit{\\textbf{F. Hierarchical Scale Parameters}}} \\\\",
  "\\addlinespace[3pt]"
)

for (param in c("sigma_gene", "sigma_int")) {
  cells <- sapply(parenting_vars, function(pv) {
    r <- hyper_lookup(pv, param)
    if (!is.null(r)) fmt_hyper_cell(r$estimate, r$q025, r$q975) else "---"
  })
  param_label <- if (param == "sigma_gene") "$\\sigma_{\\text{gene}}$" else "$\\sigma_{\\text{int}}$"
  table1_rows <- c(table1_rows,
    sprintf("%s & %s & %s & %s \\\\", param_label, cells[1], cells[2], cells[3]))
}

note_main <- paste0(
  "\\textit{Note:} Panels A--E show posterior mean estimate with posterior probability (\\%) in parentheses. ",
  "Posterior probability = max(P($\\beta>0$), P($\\beta<0$)). ",
  "\\textbf{Bold} indicates PP $\\geq$ 70\\%. ",
  "Panel F shows posterior mean with 95\\% credible interval. ",
  "Estimates are log-odds from Bayesian cumulative logit regression ",
  "($n = 168$, 100 multiply-imputed datasets, draw-level pooling). ",
  "Child behavior covariates (passivity, difficultness) from the CARE Index ",
  "are included in all models as controls. ",
  "Gene main effects share scale $\\sigma_{\\text{gene}}$; ",
  "interactions share $\\sigma_{\\text{int}}$ (hierarchical partial pooling). ",
  "Covariates: Normal(0, 2); hyperpriors: half-Normal(0, 0.5). ",
  "``Parenting'' denotes the model-specific parenting variable ",
  "(sensitivity, controlling, or unresponsiveness). ",
  "P $\\times$ Gene denotes the interaction of that parenting variable with each gene. ",
  "Within each panel, rows are sorted by maximum absolute estimate across models (descending)."
)

table1 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Bayesian Ordinal Regression Results: Hierarchical Partial Pooling with Child Behavior Controls}",
  "\\label{tab:v3_main}",
  "\\begin{tabular}{l ccc}",
  "\\toprule",
  " & \\textbf{Sensitivity} & \\textbf{Controlling} & \\textbf{Unresponsiveness} \\\\",
  " & Est (PP\\%) & Est (PP\\%) & Est (PP\\%) \\\\",
  "\\midrule",
  table1_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\vspace{0.3cm}",
  "",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\footnotesize",
  note_main,
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(table1, "results/latex_tables/v3_main_results.tex")
cat("  Saved: results/latex_tables/v3_main_results.tex\n")

# ══════════════════════════════════════════════════════════════════
# TABLE 2: Three-Scoring-Rule Model Comparison
# ══════════════════════════════════════════════════════════════════
cat("Creating Table 2: Three-scoring-rule model comparison...\n")

# Load RPS data if available
rps_file <- "results/v3_loo_rps/all_rps_combined.csv"
rps_pw_file <- "results/v3_loo_rps/pairwise_summary.csv"

if (file.exists(rps_file) && file.exists(rps_pw_file)) {
  rps_all <- read.csv(rps_file, stringsAsFactors = FALSE)
  rps_pw  <- read.csv(rps_pw_file, stringsAsFactors = FALSE)

  # ── Panel A: Model-level summaries ──
  # Average across parenting models for M1 and M3 (same across all)
  rps_summary <- rps_all %>%
    group_by(model, parenting) %>%
    summarise(
      mean_elpd = mean(loo_elpd),
      mean_rps  = mean(mean_rps),
      mean_bin  = mean(binary_elpd),
      mean_ploo = mean(p_loo),
      mean_highk = mean(n_high_k),
      n = n(),
      .groups = "drop"
    )

  t2_rows <- character()

  # M1 (same across parenting)
  m1_rps <- rps_summary %>% filter(model == "M1") %>%
    summarise(across(c(mean_elpd, mean_rps, mean_bin, mean_ploo), mean))
  t2_rows <- c(t2_rows,
    sprintf("M1 (demog.~+ child behav.) & --- & %.1f & %.4f & %.1f & --- \\\\",
            m1_rps$mean_elpd, m1_rps$mean_rps, m1_rps$mean_bin))

  # M3 (same across parenting)
  m3_rps <- rps_summary %>% filter(model == "M3") %>%
    summarise(across(c(mean_elpd, mean_rps, mean_bin, mean_ploo), mean))
  t2_rows <- c(t2_rows,
    sprintf("M3 (M1 + 9 genes) & --- & %.1f & %.4f & %.1f & --- \\\\",
            m3_rps$mean_elpd, m3_rps$mean_rps, m3_rps$mean_bin))

  t2_rows <- c(t2_rows, "\\midrule")

  # M2 and M4 by parenting
  for (pv in parenting_vars) {
    m2r <- rps_summary %>% filter(model == "M2", parenting == pv)
    m4r <- rps_summary %>% filter(model == "M4", parenting == pv)
    lab <- parenting_label[pv]
    if (nrow(m2r) > 0) {
      t2_rows <- c(t2_rows,
        sprintf("M2 (M1 + parenting) & %s & %.1f & %.4f & %.1f & --- \\\\",
                lab, m2r$mean_elpd, m2r$mean_rps, m2r$mean_bin))
    }
    if (nrow(m4r) > 0) {
      t2_rows <- c(t2_rows,
        sprintf("M4 (full) & %s & %.1f & %.4f & %.1f & --- \\\\",
                lab, m4r$mean_elpd, m4r$mean_rps, m4r$mean_bin))
    }
  }

  # ── Panel B: Pairwise comparisons under all 3 rules ──
  t2_rows <- c(t2_rows, "\\midrule",
    "\\multicolumn{6}{l}{\\textit{\\textbf{Pairwise Comparisons}}} \\\\",
    "\\addlinespace[3pt]",
    " & & $\\Delta$ELPD & $\\Delta$RPS & $\\Delta$Binary & \\% pos. \\\\",
    "\\addlinespace[2pt]"
  )

  for (comp in c("M2 vs M1", "M3 vs M1", "M4 vs M2", "M4 vs M1")) {
    pw_sub <- rps_pw %>% filter(comparison == comp)
    for (i in seq_len(nrow(pw_sub))) {
      r <- pw_sub[i, ]
      lab <- if (r$parenting %in% names(parenting_label)) parenting_label[r$parenting] else "---"
      # For % positive, show the max across all 3 rules (all are 0 here)
      max_pct <- max(r$pct_complex_better_elpd, r$pct_complex_better_rps,
                     r$pct_complex_better_binary)
      t2_rows <- c(t2_rows,
        sprintf("%s & %s & %+.2f & %+.4f & %+.2f & %.0f\\%% \\\\",
                comp, lab, r$mean_d_elpd, r$mean_d_rps, r$mean_d_binary, max_pct))
    }
    if (comp != "M4 vs M1") t2_rows <- c(t2_rows, "\\addlinespace[2pt]")
  }

  note_loo <- paste0(
    "\\textit{Note:} Three proper scoring rules evaluated under LOO-CV across 99 multiply-imputed datasets. ",
    "ELPD = expected log pointwise predictive density (higher = better); evaluates predicted probability of the exact observed category. ",
    "RPS = Ranked Probability Score (lower = better); evaluates the full cumulative distribution, penalizing predictions ",
    "far from the truth more than those that are close---the appropriate metric for ordinal outcomes. ",
    "Binary ELPD = log predictive density under a secure (categories 4--7) vs.\\ insecure (1--3) dichotomization. ",
    "$\\Delta$ = complex $-$ simple for ELPD and binary (negative = simpler better); ",
    "simple $-$ complex for RPS (negative = simpler better). ",
    "\\% pos.\\ = maximum percentage of imputations where the complex model was better under any rule. ",
    "M1: PC1--3 + Sex + infant age + passivity + difficultness. ",
    "M2: M1 + parenting. M3: M1 + 9 genes (hierarchical). ",
    "M4: M2 + 9 genes + 9 gene $\\times$ parenting interactions (hierarchical)."
  )

  table2 <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\small",
    "\\caption{LOO Cross-Validation: Three-Scoring-Rule Model Comparison}",
    "\\label{tab:v3_loo}",
    "\\begin{tabular}{l l rrrr}",
    "\\toprule",
    "Model & Parenting & ELPD & RPS & Binary & \\\\",
    "\\midrule",
    t2_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.3cm}",
    "",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize",
    note_loo,
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(table2, "results/latex_tables/v3_loo_comparison.tex")
  cat("  Saved: results/latex_tables/v3_loo_comparison.tex\n")

} else {
  cat("  RPS data not found, falling back to ELPD-only table.\n")
  # Fallback ELPD-only table omitted for brevity
}

# ══════════════════════════════════════════════════════════════════
# FIGURE 1: Forest Plot
# ══════════════════════════════════════════════════════════════════
cat("Creating Figure 1: Forest plot...\n")

# Homogeneity test for gene main effects
gene_df <- df %>% filter(effect_type == "gene_main")
pooled_gene <- gene_df %>%
  group_by(term) %>%
  summarise(
    mean_est  = mean(estimate),
    mean_sd   = mean(sd),
    mean_q025 = mean(q025),
    mean_q975 = mean(q975),
    mean_pp   = mean(post_prob),
    between_sd = sd(estimate),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(mean_est)))

# Pooled gene main effects for forest plot
pooled_gene_forest <- pooled_gene %>%
  mutate(
    display_term = term,
    estimate = mean_est,
    q025 = mean_q025,
    q975 = mean_q975,
    post_prob = mean_pp,
    pp_cat = case_when(
      post_prob >= 0.70 ~ "PP >= 70%",
      post_prob >= 0.60 ~ "PP 60-70%",
      TRUE              ~ "PP < 60%"
    ),
    pp_cat = factor(pp_cat, levels = c("PP >= 70%", "PP 60-70%", "PP < 60%"))
  )

# Interactions
int_df <- df %>%
  filter(effect_type == "interaction") %>%
  mutate(
    gene_name = sub("^(sens|cont|unre):", "", term),
    display_term = paste0("P \u00d7 ", gene_name),
    model = parenting_label[parenting],
    pp_cat = case_when(
      post_prob >= 0.70 ~ "PP >= 70%",
      post_prob >= 0.60 ~ "PP 60-70%",
      TRUE              ~ "PP < 60%"
    ),
    pp_cat = factor(pp_cat, levels = c("PP >= 70%", "PP 60-70%", "PP < 60%"))
  )

gene_order <- pooled_gene_forest %>%
  arrange(abs(estimate)) %>%
  pull(display_term)

int_order <- int_df %>%
  group_by(display_term) %>%
  summarise(max_abs = max(abs(estimate)), .groups = "drop") %>%
  arrange(max_abs) %>%
  pull(display_term)

pooled_gene_forest$display_term <- factor(pooled_gene_forest$display_term,
                                          levels = gene_order)
int_df$display_term <- factor(int_df$display_term, levels = int_order)
int_df$model <- factor(int_df$model,
                       levels = c("Sensitivity", "Controlling", "Unresponsiveness"))

p_gene <- ggplot(pooled_gene_forest,
       aes(x = estimate, y = display_term, color = pp_cat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.3, linewidth = 0.6) +
  geom_point(aes(size = post_prob), shape = 16) +
  scale_color_manual(
    values = c("PP >= 70%" = "#d62728", "PP 60-70%" = "#1f77b4",
               "PP < 60%" = "grey55"),
    name = NULL, drop = FALSE
  ) +
  scale_size_continuous(range = c(2, 4), guide = "none") +
  labs(x = "Posterior Mean (log-odds)", y = NULL,
       title = "A. Gene Main Effects (pooled across parenting models)") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

p_int <- ggplot(int_df,
       aes(x = estimate, y = display_term, color = pp_cat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.3, linewidth = 0.5) +
  geom_point(aes(size = post_prob), shape = 16) +
  scale_color_manual(
    values = c("PP >= 70%" = "#d62728", "PP 60-70%" = "#1f77b4",
               "PP < 60%" = "grey55"),
    name = NULL, drop = FALSE
  ) +
  scale_size_continuous(range = c(1.5, 3.5), guide = "none") +
  facet_wrap(~model, ncol = 3) +
  labs(x = "Posterior Mean (log-odds)", y = NULL,
       title = "B. Gene \u00d7 Parenting Interactions (model-specific)") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 12)
  )

p_forest <- arrangeGrob(p_gene, p_int, ncol = 1, heights = c(0.8, 1))
ggsave("results/v3_figures/forest_plot_combined.pdf", p_forest, width = 10, height = 9)
ggsave("results/v3_figures/forest_gene_mains_pooled.pdf", p_gene, width = 6, height = 4)
ggsave("results/v3_figures/forest_interactions.pdf", p_int, width = 10, height = 4.5)
cat("  Saved: results/v3_figures/forest_plot_combined.pdf\n")

# ══════════════════════════════════════════════════════════════════
# FIGURE 2: Evidence Calibration
# ══════════════════════════════════════════════════════════════════
cat("Creating Figure 2: Evidence calibration...\n")

all_effects <- effects %>%
  filter(effect_type %in% c("gene_main", "interaction")) %>%
  mutate(
    bf_directional = post_prob / (1 - post_prob),
    twoln_bf = 2 * log(bf_directional)
  )

classify_kr <- function(bf) {
  twoln <- 2 * log(bf)
  ifelse(twoln > 10, "Very strong",
  ifelse(twoln > 6,  "Strong",
  ifelse(twoln > 2,  "Positive",
                      "Bare mention")))
}

all_effects$kr_evidence <- classify_kr(all_effects$bf_directional)

# Pooled gene mains
gene_pooled <- pooled_gene %>%
  mutate(
    post_prob = mean_pp,
    bf_directional = post_prob / (1 - post_prob),
    twoln_bf = 2 * log(bf_directional),
    kr_evidence = classify_kr(bf_directional)
  )

# For plot: pooled genes + all interactions with PP > 55%
gene_plot <- gene_pooled %>%
  filter(post_prob >= 0.55) %>%
  mutate(
    eff_label = paste0(term, " (gene main)"),
    effect_type = "gene_main",
    kr_cat = factor(kr_evidence,
                    levels = c("Bare mention", "Positive", "Strong", "Very strong"))
  )

int_plot <- all_effects %>%
  filter(effect_type == "interaction", post_prob >= 0.55) %>%
  mutate(
    model = parenting_label[parenting],
    gene = sub("^(sens|cont|unre):", "", term),
    eff_label = paste0(model, ": P \u00d7 ", gene),
    kr_cat = factor(kr_evidence,
                    levels = c("Bare mention", "Positive", "Strong", "Very strong"))
  )

combined_bf <- bind_rows(
  gene_plot %>% select(eff_label, post_prob, bf_directional, twoln_bf, kr_cat, effect_type),
  int_plot %>% select(eff_label, post_prob, bf_directional, twoln_bf, kr_cat, effect_type)
) %>%
  arrange(twoln_bf) %>%
  mutate(eff_label = factor(eff_label, levels = eff_label))

p_bf <- ggplot(combined_bf, aes(x = twoln_bf, y = eff_label, colour = kr_cat)) +
  annotate("rect", xmin = -Inf, xmax = 2, ymin = -Inf, ymax = Inf,
           fill = "#fee0d2", alpha = 0.3) +
  annotate("rect", xmin = 2, xmax = 6, ymin = -Inf, ymax = Inf,
           fill = "#deebf7", alpha = 0.3) +
  annotate("rect", xmin = 6, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#e5f5e0", alpha = 0.3) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_point(aes(shape = effect_type), size = 3) +
  annotate("text", x = 1, y = Inf, label = "Bare mention", vjust = 1.5,
           size = 3, colour = "grey40", fontface = "italic") +
  annotate("text", x = 4, y = Inf, label = "Positive", vjust = 1.5,
           size = 3, colour = "grey40", fontface = "italic") +
  annotate("text", x = 7, y = Inf, label = "Strong", vjust = 1.5,
           size = 3, colour = "grey40", fontface = "italic") +
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
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  labs(
    title = "Evidence calibration: Directional Bayes factors",
    subtitle = "Under symmetric hierarchical prior, BF = PP / (1 - PP); model includes child behavior controls",
    x = expression(2 ~ ln(BF[directional])),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("results/v3_figures/evidence_calibration.pdf", p_bf, width = 9, height = 7)
cat("  Saved: results/v3_figures/evidence_calibration.pdf\n")

# Print evidence summary
cat("\n  Key effects (PP >= 65%):\n")
key_effs <- all_effects %>%
  filter(post_prob >= 0.65) %>%
  arrange(desc(bf_directional))
if (nrow(key_effs) > 0) {
  for (i in seq_len(nrow(key_effs))) {
    r <- key_effs[i, ]
    cat(sprintf("    %s %-22s PP=%.0f%% BF=%.2f [%s]\n",
                parenting_label[r$parenting], r$term,
                r$post_prob * 100, r$bf_directional, r$kr_evidence))
  }
}

# ══════════════════════════════════════════════════════════════════
# FIGURE 3: Posterior Densities for PP >= 70%
# ══════════════════════════════════════════════════════════════════
cat("\nCreating Figure 3: Posterior densities for PP >= 70%...\n")

high_pp <- df %>%
  filter(post_prob >= 0.70,
         effect_type %in% c("gene_main", "interaction", "covariate"),
         !(effect_type == "covariate" & !(term %in% c("pass", "diff", parenting_vars)))) %>%
  mutate(
    model = parenting_label[parenting],
    gene_name = case_when(
      effect_type == "interaction" ~ sub("^(sens|cont|unre):", "P \u00d7 ", term),
      effect_type == "covariate" & term == "pass" ~ "Passivity",
      effect_type == "covariate" & term == "diff" ~ "Difficultness",
      effect_type == "covariate" & term %in% parenting_vars ~ paste0("Parenting (", model, ")"),
      TRUE ~ term
    ),
    label = paste0(gene_name, "\n", model, " model")
  ) %>%
  arrange(desc(post_prob))

cat(sprintf("  Found %d effects with PP >= 70%%\n", nrow(high_pp)))

if (nrow(high_pp) > 0) {
  set.seed(42)
  density_data <- do.call(rbind, lapply(seq_len(nrow(high_pp)), function(i) {
    r <- high_pp[i, ]
    x <- seq(r$estimate - 4*r$sd, r$estimate + 4*r$sd, length.out = 500)
    y <- dnorm(x, mean = r$estimate, sd = r$sd)
    data.frame(
      x = x, y = y,
      label = r$label,
      gene_name = r$gene_name,
      model = r$model,
      estimate = r$estimate,
      q025 = r$q025, q975 = r$q975,
      post_prob = r$post_prob,
      direction = r$direction,
      shade = ifelse(r$direction == "positive", x > 0, x < 0),
      stringsAsFactors = FALSE
    )
  }))

  density_data$label <- factor(density_data$label,
    levels = rev(unique(high_pp$label)))

  p_density <- ggplot(density_data, aes(x = x, y = y)) +
    geom_area(data = density_data %>% filter(shade),
              fill = "#3182bd", alpha = 0.3) +
    geom_line(color = "grey20", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(aes(xintercept = estimate), color = "#d62728", linewidth = 0.5) +
    geom_segment(aes(x = q025, xend = q975, y = 0, yend = 0),
                 color = "#d62728", linewidth = 1.5) +
    geom_text(
      data = high_pp %>% mutate(label = factor(label, levels = rev(unique(label)))),
      aes(x = Inf, y = Inf, label = sprintf("PP = %.0f%%", post_prob * 100)),
      hjust = 1.1, vjust = 1.5, size = 3.2, color = "grey30", fontface = "italic"
    ) +
    facet_wrap(~label, scales = "free", ncol = 3) +
    labs(
      x = "Log-odds coefficient",
      y = "Density",
      title = "Approximate Posterior Distributions for Effects with PP \u2265 70%",
      subtitle = "Blue shading = direction of evidence; red line = posterior mean; red bar = 95% CI"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  n_panels <- nrow(high_pp)
  fig_height <- ceiling(n_panels / 3) * 2.5 + 1.2
  ggsave("results/v3_figures/posterior_densities_high_pp.pdf", p_density,
         width = 10, height = fig_height)
  cat(sprintf("  Saved: results/v3_figures/posterior_densities_high_pp.pdf (%d panels)\n",
              n_panels))
}

# ══════════════════════════════════════════════════════════════════
# FIGURE 4: LOO Delta ELPD Distribution
# ══════════════════════════════════════════════════════════════════
cat("Creating Figure 4: LOO delta ELPD distribution...\n")

# Compute per-imputation deltas for M4 vs M2
delta_data <- list()

for (pv in parenting_vars) {
  m1_sub <- loo_all %>% filter(model == "M1")
  m2_sub <- loo_all %>% filter(model == "M2", parenting == pv)
  m3_sub <- loo_all %>% filter(model == "M3")
  m4_sub <- loo_all %>% filter(model == "M4", parenting == pv)

  # M4 vs M2
  shared <- intersect(m2_sub$imputation, m4_sub$imputation)
  if (length(shared) > 0) {
    d42 <- m4_sub$elpd_loo[match(shared, m4_sub$imputation)] -
           m2_sub$elpd_loo[match(shared, m2_sub$imputation)]
    delta_data[[paste0(pv, "_M4vM2")]] <- data.frame(
      parenting = parenting_label[pv],
      comparison = "M4 vs M2",
      delta = d42, stringsAsFactors = FALSE
    )
  }

  # M4 vs M1
  shared <- intersect(m1_sub$imputation, m4_sub$imputation)
  if (length(shared) > 0) {
    d41 <- m4_sub$elpd_loo[match(shared, m4_sub$imputation)] -
           m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
    delta_data[[paste0(pv, "_M4vM1")]] <- data.frame(
      parenting = parenting_label[pv],
      comparison = "M4 vs M1",
      delta = d41, stringsAsFactors = FALSE
    )
  }

  # M2 vs M1
  shared <- intersect(m1_sub$imputation, m2_sub$imputation)
  if (length(shared) > 0) {
    d21 <- m2_sub$elpd_loo[match(shared, m2_sub$imputation)] -
           m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
    delta_data[[paste0(pv, "_M2vM1")]] <- data.frame(
      parenting = parenting_label[pv],
      comparison = "M2 vs M1",
      delta = d21, stringsAsFactors = FALSE
    )
  }
}

# M3 vs M1 (same across all parenting)
m1_sub <- loo_all %>% filter(model == "M1")
m3_sub <- loo_all %>% filter(model == "M3")
shared <- intersect(m1_sub$imputation, m3_sub$imputation)
if (length(shared) > 0) {
  d31 <- m3_sub$elpd_loo[match(shared, m3_sub$imputation)] -
         m1_sub$elpd_loo[match(shared, m1_sub$imputation)]
  delta_data[["M3vM1"]] <- data.frame(
    parenting = "All",
    comparison = "M3 vs M1",
    delta = d31, stringsAsFactors = FALSE
  )
}

delta_df <- do.call(rbind, delta_data)
rownames(delta_df) <- NULL

delta_df$comparison <- factor(delta_df$comparison,
  levels = c("M2 vs M1", "M3 vs M1", "M4 vs M2", "M4 vs M1"))

p_delta <- ggplot(delta_df, aes(x = comparison, y = delta, fill = parenting)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(
    values = c("Sensitivity" = "#1b9e77", "Controlling" = "#d95f02",
               "Unresponsiveness" = "#7570b3", "All" = "grey70"),
    name = "Parenting model"
  ) +
  labs(
    title = expression(paste(Delta, "ELPD distribution across 100 imputations")),
    subtitle = "Negative values = simpler model predicts better; models include child behavior controls",
    x = NULL,
    y = expression(paste(Delta, "ELPD (complex - simple)"))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave("results/v3_figures/loo_delta_distribution.pdf", p_delta, width = 9, height = 5.5)
cat("  Saved: results/v3_figures/loo_delta_distribution.pdf\n")

# ══════════════════════════════════════════════════════════════════
# TABLE 3: Explained Variation (Bayesian R², McKelvey-Zavoina R², RPSS)
# ══════════════════════════════════════════════════════════════════
cat("Creating Table 3: Explained variation measures...\n")

r2_file <- "results/v3_r2/r2_summary.csv"
rpss_file <- "results/v3_loo_rps/rpss_summary.csv"

if (file.exists(r2_file) && file.exists(rpss_file)) {
  r2_sum  <- read.csv(r2_file, stringsAsFactors = FALSE)
  rpss_sum <- read.csv(rpss_file, stringsAsFactors = FALSE)

  t3_rows <- character()

  for (pv in parenting_vars) {
    lab <- parenting_label[pv]
    t3_rows <- c(t3_rows,
      sprintf("\\multicolumn{6}{l}{\\textit{%s}} \\\\", lab),
      "\\addlinespace[2pt]"
    )

    for (m in c("M1", "M2", "M3", "M4")) {
      r2r <- r2_sum[r2_sum$model == m & r2_sum$parenting == pv, ]
      rr  <- rpss_sum[rpss_sum$model == m & rpss_sum$parenting == pv, ]

      if (nrow(r2r) > 0 && nrow(rr) > 0) {
        t3_rows <- c(t3_rows,
          sprintf("\\quad %s & %.3f $[%.3f,\\, %.3f]$ & %.3f $[%.3f,\\, %.3f]$ & %+.3f $[%+.3f,\\, %+.3f]$ \\\\",
                  m,
                  r2r$bayes_r2_mean[1], r2r$bayes_r2_q025[1], r2r$bayes_r2_q975[1],
                  r2r$mz_r2_mean[1], r2r$mz_r2_q025[1], r2r$mz_r2_q975[1],
                  rr$mean_rpss[1], rr$q025_rpss[1], rr$q975_rpss[1]))
      }
    }
    t3_rows <- c(t3_rows, "\\addlinespace[4pt]")
  }

  note_r2 <- paste0(
    "\\textit{Note:} Three measures of explained variation. ",
    "Bayesian $R^2$ (Gelman et al., 2019) = Var(predicted) / [Var(predicted) + Var(residual)], ",
    "estimated from posterior predictive draws. ",
    "McKelvey--Zavoina $R^2$ = Var($X\\beta$) / [Var($X\\beta$) + $\\pi^2/3$], ",
    "the natural $R^2$ for latent-variable ordinal models. ",
    "RPSS (Ranked Probability Skill Score) = $1 - \\text{RPS}_{\\text{model}} / \\text{RPS}_{\\text{reference}}$, ",
    "where the reference is the marginal category-frequency forecast; evaluated out-of-sample via LOO-CV. ",
    "RPSS $> 0$ indicates skill beyond base rates; RPSS $< 0$ indicates worse-than-base-rate prediction. ",
    "Values in brackets are 95\\% credible intervals ($R^2$ measures) or 95\\% quantile range across imputations (RPSS). ",
    "M1: demographics + child behavior. M2: M1 + parenting. M3: M1 + 9 genes. M4: M2 + genes + interactions."
  )

  table3 <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\small",
    "\\caption{Explained Variation: Bayesian $R^2$, McKelvey--Zavoina $R^2$, and RPSS}",
    "\\label{tab:v3_r2}",
    "\\begin{tabular}{l ccc}",
    "\\toprule",
    " & Bayesian $R^2$ & McKelvey--Zavoina $R^2$ & RPSS (out-of-sample) \\\\",
    " & Mean $[95\\%~\\text{CI}]$ & Mean $[95\\%~\\text{CI}]$ & Mean $[95\\%~\\text{range}]$ \\\\",
    "\\midrule",
    t3_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.3cm}",
    "",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize",
    note_r2,
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(table3, "results/latex_tables/v3_explained_variation.tex")
  cat("  Saved: results/latex_tables/v3_explained_variation.tex\n")
} else {
  cat("  R² or RPSS data not found, skipping Table 3.\n")
}

# ══════════════════════════════════════════════════════════════════
# TABLE 4: Prior Sensitivity Analysis
# ══════════════════════════════════════════════════════════════════
cat("Creating Table 4: Prior sensitivity analysis...\n")

sens_file <- "results/v3_sensitivity/all_sensitivity.csv"
pp_stab_file <- "results/v3_sensitivity/pp_stability.csv"

if (file.exists(sens_file)) {
  sens_all <- read.csv(sens_file, stringsAsFactors = FALSE)

  hyper_sens <- sens_all[sens_all$effect_type == "hyperparameter", ]
  effects_sens <- sens_all[sens_all$effect_type %in% c("gene_main", "interaction"), ]

  # ── Panel A: Hyperparameter posteriors ──
  t4_rows <- character()
  t4_rows <- c(t4_rows,
    "\\multicolumn{4}{l}{\\textit{\\textbf{A. Hierarchical Scale Parameter Posteriors}}} \\\\",
    "\\addlinespace[3pt]",
    " & $s = 0.25$ & $s = 0.50$ & $s = 1.00$ \\\\",
    "\\addlinespace[2pt]"
  )

  for (pv in parenting_vars) {
    lab <- parenting_label[pv]
    t4_rows <- c(t4_rows,
      sprintf("\\multicolumn{4}{l}{\\quad\\textit{%s}} \\\\", lab))

    for (param in c("sigma_gene", "sigma_int")) {
      cells <- sapply(c("s025", "s050", "s100"), function(sl) {
        r <- hyper_sens[hyper_sens$parenting == pv &
                        hyper_sens$term == param &
                        hyper_sens$scale_label == sl, ]
        if (nrow(r) > 0) {
          sprintf("%.3f $[%.3f,\\, %.3f]$", r$estimate[1], r$q025[1], r$q975[1])
        } else "---"
      })
      param_label <- if (param == "sigma_gene") "$\\sigma_{\\text{gene}}$" else "$\\sigma_{\\text{int}}$"
      t4_rows <- c(t4_rows,
        sprintf("\\quad\\quad %s & %s & %s & %s \\\\",
                param_label, cells[1], cells[2], cells[3]))
    }
    t4_rows <- c(t4_rows, "\\addlinespace[3pt]")
  }

  # ── Panel B: PP stability for key effects ──
  t4_rows <- c(t4_rows,
    "\\addlinespace[4pt]",
    "\\multicolumn{4}{l}{\\textit{\\textbf{B. Posterior Probability Stability (PP \\%, effects with PP $\\geq$ 70\\% at $s=0.50$)}}} \\\\",
    "\\addlinespace[3pt]",
    " & $s = 0.25$ & $s = 0.50$ & $s = 1.00$ \\\\",
    "\\addlinespace[2pt]"
  )

  # Filter to effects with PP >= 70% at s=0.50
  eff_s050 <- effects_sens[effects_sens$scale_label == "s050" & effects_sens$post_prob >= 0.70, ]

  if (nrow(eff_s050) > 0) {
    # Get PP at other scales for these effects
    for (i in seq_len(nrow(eff_s050))) {
      pv <- eff_s050$parenting[i]
      tm <- eff_s050$term[i]

      cells <- sapply(c("s025", "s050", "s100"), function(sl) {
        r <- effects_sens[effects_sens$parenting == pv &
                          effects_sens$term == tm &
                          effects_sens$scale_label == sl, ]
        if (nrow(r) > 0) sprintf("%.1f", r$post_prob[1] * 100) else "---"
      })

      # Format term name
      display_term <- tm
      if (eff_s050$effect_type[i] == "interaction") {
        gene <- sub("^(sens|cont|unre):", "", tm)
        display_term <- paste0("P $\\times$ ", gsub("_", "\\_", gene, fixed = TRUE))
      } else {
        display_term <- gsub("_", "\\_", tm, fixed = TRUE)
      }

      lab <- substr(parenting_label[pv], 1, 4)
      t4_rows <- c(t4_rows,
        sprintf("%s (%s) & %s & \\textbf{%s} & %s \\\\",
                display_term, lab, cells[1], cells[2], cells[3]))
    }
  }

  # PP range summary
  if (file.exists(pp_stab_file)) {
    pp_stab <- read.csv(pp_stab_file, stringsAsFactors = FALSE)
    max_range <- max(pp_stab$pp_range, na.rm = TRUE) * 100
    mean_range <- mean(pp_stab$pp_range, na.rm = TRUE) * 100
    t4_rows <- c(t4_rows,
      "\\addlinespace[4pt]",
      sprintf("\\multicolumn{4}{l}{\\quad Maximum PP variation across scales: %.1f percentage points} \\\\", max_range),
      sprintf("\\multicolumn{4}{l}{\\quad Mean PP variation across scales: %.1f percentage points} \\\\", mean_range)
    )
  }

  note_sens <- paste0(
    "\\textit{Note:} Prior sensitivity analysis for the M4 (full hierarchical) model ",
    "under three hyperprior scales: half-Normal$(0, s)$ with $s \\in \\{0.25, 0.50, 1.00\\}$. ",
    "Panel A shows posterior mean and 95\\% credible interval for the hierarchical scale parameters. ",
    "Panel B shows posterior probabilities (\\%) for effects that reached PP $\\geq$ 70\\% ",
    "at the primary scale ($s = 0.50$, bolded). ",
    "Maximum PP variation of $< 5$ percentage points across a 4$\\times$ range of prior scales ",
    "indicates that substantive conclusions are robust to the choice of hyperprior."
  )

  table4 <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\small",
    "\\caption{Prior Sensitivity Analysis: M4 under Three Hyperprior Scales}",
    "\\label{tab:v3_sensitivity}",
    "\\begin{tabular}{l ccc}",
    "\\toprule",
    t4_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.3cm}",
    "",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize",
    note_sens,
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(table4, "results/latex_tables/v3_prior_sensitivity.tex")
  cat("  Saved: results/latex_tables/v3_prior_sensitivity.tex\n")
} else {
  cat("  Sensitivity data not found, skipping Table 4.\n")
}

# ══════════════════════════════════════════════════════════════════
# Complete compilable LaTeX document
# ══════════════════════════════════════════════════════════════════
cat("\nCreating complete LaTeX document...\n")

# Copy figures for compilation
for (fig in c("forest_plot_combined.pdf", "evidence_calibration.pdf",
              "posterior_densities_high_pp.pdf", "loo_delta_distribution.pdf")) {
  src <- file.path("results/v3_figures", fig)
  if (file.exists(src)) {
    file.copy(src, file.path("results/latex_tables", fig), overwrite = TRUE)
  }
}

doc <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage{booktabs}",
  "\\usepackage{geometry}",
  "\\usepackage{amsmath}",
  "\\usepackage{makecell}",
  "\\usepackage{multirow}",
  "\\usepackage{graphicx}",
  "\\usepackage{float}",
  "\\geometry{margin=0.75in}",
  "",
  "\\title{Gene--Environment Interaction Analysis (v3):\\\\",
  "       Bayesian Ordinal Regression with Child Behavior Controls}",
  "\\author{Potter-Dickey, Letourneau \\& de Koning}",
  "\\date{\\today}",
  "",
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\begin{abstract}",
  "\\noindent",
  "We applied Bayesian ordinal regression with hierarchical partial pooling to",
  "examine gene--environment interactions between nine candidate genes and three",
  "parenting measures from the CARE Index (sensitivity, controlling, unresponsiveness)",
  "in predicting infant attachment classification ($n = 168$). In this analysis,",
  "two child behavior sub-indices of the CARE Index---passivity and difficultness---are",
  "included as covariates in all models, providing a stronger baseline that controls",
  "for infant behavioral contributions to the dyadic interaction.",
  "Results show that no model including parenting, genes, or gene $\\times$ parenting",
  "interactions improves out-of-sample prediction over the baseline model",
  "(demographics + child behavior). All $\\Delta$ELPD values are negative across",
  "100 multiply-imputed datasets.",
  "\\end{abstract}",
  "",
  "\\clearpage",
  "",
  "% ── TABLE 1 ──",
  readLines("results/latex_tables/v3_main_results.tex"),
  "",
  "\\clearpage",
  "",
  "% ── TABLE 2 ──",
  readLines("results/latex_tables/v3_loo_comparison.tex"),
  "",
  "\\clearpage",
  ""
)

# Add Table 3 if it exists
if (file.exists("results/latex_tables/v3_explained_variation.tex")) {
  doc <- c(doc,
    "% ── TABLE 3 ──",
    readLines("results/latex_tables/v3_explained_variation.tex"),
    "",
    "\\clearpage",
    ""
  )
}

# Add Table 4 if it exists
if (file.exists("results/latex_tables/v3_prior_sensitivity.tex")) {
  doc <- c(doc,
    "% ── TABLE 4 ──",
    readLines("results/latex_tables/v3_prior_sensitivity.tex"),
    "",
    "\\clearpage",
    ""
  )
}

doc <- c(doc,
  "% ── FIGURE 1 ──",
  "\\begin{figure}[H]",
  "\\centering",
  "\\includegraphics[width=\\textwidth]{forest_plot_combined.pdf}",
  "\\caption{Forest plot of gene main effects (pooled across parenting models)",
  "and gene $\\times$ parenting interactions. Point size proportional to posterior",
  "probability. Child behavior covariates (passivity, difficultness) controlled",
  "in all models.}",
  "\\label{fig:v3_forest}",
  "\\end{figure}",
  "",
  "\\clearpage",
  "",
  "% ── FIGURE 2 ──",
  "\\begin{figure}[H]",
  "\\centering",
  "\\includegraphics[width=\\textwidth]{evidence_calibration.pdf}",
  "\\caption{Evidence calibration: directional Bayes factors mapped onto the",
  "Kass \\& Raftery (1995) scale. Under the symmetric hierarchical prior,",
  "$\\text{BF} = \\text{PP} / (1 - \\text{PP})$. All models include child behavior controls.}",
  "\\label{fig:v3_evidence}",
  "\\end{figure}",
  "",
  "\\clearpage",
  ""
)

# Add posterior density figure if it exists
if (file.exists("results/latex_tables/posterior_densities_high_pp.pdf")) {
  doc <- c(doc,
    "% ── FIGURE 3 ──",
    "\\begin{figure}[H]",
    "\\centering",
    "\\includegraphics[width=\\textwidth]{posterior_densities_high_pp.pdf}",
    "\\caption{Approximate posterior distributions for effects with PP $\\geq$ 70\\%.",
    "Blue shading indicates direction of evidence; red line = posterior mean;",
    "red bar = 95\\% credible interval.}",
    "\\label{fig:v3_densities}",
    "\\end{figure}",
    "",
    "\\clearpage",
    ""
  )
}

doc <- c(doc,
  "% ── FIGURE 4 ──",
  "\\begin{figure}[H]",
  "\\centering",
  "\\includegraphics[width=\\textwidth]{loo_delta_distribution.pdf}",
  "\\caption{Distribution of $\\Delta$ELPD across 100 multiply-imputed datasets.",
  "All values are negative, indicating that no model with additional predictors",
  "improves out-of-sample prediction over M1 (demographics + child behavior).}",
  "\\label{fig:v3_loo}",
  "\\end{figure}",
  "",
  "\\end{document}"
)

writeLines(doc, "results/latex_tables/v3_all_tables.tex")
cat("  Saved: results/latex_tables/v3_all_tables.tex\n")

# ── Summary ──────────────────────────────────────────────────────
cat("\n=================================================================\n")
cat("V3 TABLES AND FIGURES COMPLETE\n")
cat("=================================================================\n")
cat("Tables:\n")
cat("  v3_main_results.tex        -- Coefficient estimates (Table 1)\n")
cat("  v3_loo_comparison.tex      -- LOO model comparison (Table 2)\n")
cat("  v3_explained_variation.tex -- R² + RPSS (Table 3)\n")
cat("  v3_prior_sensitivity.tex   -- Prior sensitivity (Table 4)\n")
cat("  v3_all_tables.tex          -- Complete compilable document\n")
cat("Figures:\n")
cat("  forest_plot_combined.pdf   -- Gene mains + interactions\n")
cat("  evidence_calibration.pdf   -- Directional BF scale\n")
cat("  posterior_densities_high_pp.pdf -- Densities for PP >= 70%\n")
cat("  loo_delta_distribution.pdf -- dELPD boxplots\n")
cat("=================================================================\n")

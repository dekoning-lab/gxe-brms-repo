#!/usr/bin/env Rscript
############################################################
# 99_predicted_probs_v3_cluster.R
# Predicted probability plots for three key v3 sensitivity effects:
#   1. SLC6A3.10R gene main effect
#   2. Sensitivity x CNR1.77 interaction
#   3. Sensitivity x 5-HTTLPR interaction
#
# Pools predictions across all 100 imputations.
# Run from: the repository root (locally or on a cluster; see config.sh).
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(dplyr)
})

cat("=================================================================\n")
cat("PREDICTED PROBABILITY PLOTS — V3 SENSITIVITY MODEL\n")
cat(sprintf("Start: %s\n", Sys.time()))
cat("=================================================================\n\n")

out_dir <- "results/v3_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
dat <- imputed_data[[50]]
cat(sprintf("Data loaded: %d obs, %d vars\n", nrow(dat), ncol(dat)))

# ── Gene list and helper ─────────────────────────────────
genes <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
           "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

get_mode <- function(x) as.numeric(names(sort(table(x), decreasing = TRUE))[1])

# ── Imputations to pool ──────────────────────────────────
imp_indices <- 1:100

# ══════════════════════════════════════════════════════════
# Helper: get P(secure) = P(Y>=4) from multiple imputations
# Memory-efficient: computes P(secure) per imputation, stacks only that
# ══════════════════════════════════════════════════════════
get_pooled_p_secure <- function(nd, imp_indices) {
  all_p_secure <- list()
  for (imp_i in imp_indices) {
    fit_file <- sprintf("results/v3_m4/sens/imp_%03d.rds", imp_i)
    if (!file.exists(fit_file)) {
      cat(sprintf("  Skipping imp %d (not found)\n", imp_i))
      next
    }
    cat(sprintf("  Loading imp %d... ", imp_i))
    fit <- readRDS(fit_file)
    pred <- fitted(fit, newdata = nd, summary = FALSE)  # [draws x obs x cats]
    n_cat <- dim(pred)[3]
    # Sum P(Y=4) + ... + P(Y=n_cat) for each draw and observation
    ps <- apply(pred[, , 4:n_cat, drop = FALSE], c(1, 2), sum)  # [draws x obs]
    cat(sprintf("draws=%d, obs=%d\n", nrow(ps), ncol(ps)))
    all_p_secure[[length(all_p_secure) + 1]] <- ps
    rm(fit, pred); gc(verbose = FALSE)
  }
  # Stack: [all_draws x obs] — much smaller than full [draws x obs x cats]
  do.call(rbind, all_p_secure)
}

# ══════════════════════════════════════════════════════════
# PLOT 1: SLC6A3.10R gene main effect
# ══════════════════════════════════════════════════════════
cat("\n── Plot 1: SLC6A3.10R gene main effect ──\n")

gene_levels_10r <- sort(unique(dat$SLC6A3.10R))
cat(sprintf("  SLC6A3.10R levels: %s\n", paste(gene_levels_10r, collapse = ", ")))

nd1 <- data.frame()
for (gl in gene_levels_10r) {
  row <- data.frame(
    PC1 = 0, PC2 = 0, PC3 = 0,
    Sex = factor("1", levels = levels(dat$Sex)),
    infant_age = 0, pass = 0, diff = 0, sens = 0
  )
  row$SLC6A3.10R <- gl
  for (og in setdiff(genes, "SLC6A3.10R")) row[[og]] <- get_mode(dat[[og]])
  row$gene_group <- gl
  nd1 <- rbind(nd1, row)
}

p_secure1 <- get_pooled_p_secure(nd1, imp_indices)
cat(sprintf("  Pooled: %d draws, %d obs\n", nrow(p_secure1), ncol(p_secure1)))

plot_df1 <- data.frame(
  genotype = factor(gene_levels_10r, labels = paste0(gene_levels_10r, " copies")),
  mean = colMeans(p_secure1),
  lo = apply(p_secure1, 2, quantile, 0.025),
  hi = apply(p_secure1, 2, quantile, 0.975)
)
cat("  SLC6A3.10R predicted P(secure):\n")
print(plot_df1)

p1 <- ggplot(plot_df1, aes(x = genotype, y = mean)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  size = 1.2, linewidth = 0.9, colour = "#2166ac") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = expression(italic("SLC6A3") ~ "(DAT1 10R): Predicted probability of secure attachment"),
       subtitle = "Covariates at mean/mode; pooled across 100 imputations",
       x = "Number of 10R allele copies",
       y = "P(Attachment category >= 4)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(out_dir, "predicted_probs_slc6a3_10r.pdf"), p1, width = 6, height = 5)
cat("  Saved: predicted_probs_slc6a3_10r.pdf\n")
rm(pred1, p_secure1); gc(verbose = FALSE)

# ══════════════════════════════════════════════════════════
# PLOTS 2-3: Sensitivity x gene interactions
# ══════════════════════════════════════════════════════════
sens_range <- seq(min(dat$sens), max(dat$sens), length.out = 40)

for (gene_var in c("CNR1.77", "X5HTTLPR")) {
  cat(sprintf("\n── Plot: Sensitivity x %s ──\n", gene_var))

  gene_levels <- sort(unique(dat[[gene_var]]))
  cat(sprintf("  %s levels: %s\n", gene_var, paste(gene_levels, collapse = ", ")))

  nd <- data.frame()
  for (gl in gene_levels) {
    for (s in sens_range) {
      row <- data.frame(
        PC1 = 0, PC2 = 0, PC3 = 0,
        Sex = factor("1", levels = levels(dat$Sex)),
        infant_age = 0, pass = 0, diff = 0, sens = s
      )
      row[[gene_var]] <- gl
      for (og in setdiff(genes, gene_var)) row[[og]] <- get_mode(dat[[og]])
      row$gene_group <- gl
      row$sens_val <- s
      nd <- rbind(nd, row)
    }
  }
  cat(sprintf("  Newdata: %d rows\n", nrow(nd)))

  p_secure <- get_pooled_p_secure(nd, imp_indices)

  plot_df <- data.frame(
    sens = nd$sens_val,
    gene_group = nd$gene_group,
    mean = colMeans(p_secure),
    lo = apply(p_secure, 2, quantile, 0.025),
    hi = apply(p_secure, 2, quantile, 0.975)
  )

  # Gene labels
  if (gene_var == "X5HTTLPR") {
    plot_df$gene_label <- factor(plot_df$gene_group,
                                  levels = sort(unique(plot_df$gene_group)),
                                  labels = c("L/L", "S/L", "S/S")[seq_along(unique(plot_df$gene_group))])
    title_expr <- "Predicted P(Secure): Sensitivity × 5-HTTLPR"
  } else {
    plot_df$gene_label <- factor(plot_df$gene_group,
                                  levels = sort(unique(plot_df$gene_group)),
                                  labels = paste0(sort(unique(plot_df$gene_group)), " copies"))
    title_expr <- "Predicted P(Secure): Sensitivity × CNR1 (rs806377)"
  }

  p <- ggplot(plot_df, aes(x = sens, y = mean,
                            colour = gene_label, fill = gene_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1.1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = title_expr,
         subtitle = "Covariates at mean/mode; pooled across 100 imputations",
         x = "Maternal Sensitivity (z-scored)",
         y = "P(Attachment category >= 4)",
         colour = gene_var, fill = gene_var) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom")

  outf <- sprintf("predicted_probs_sens_%s.pdf",
                   tolower(gsub("\\.", "", gene_var)))
  ggsave(file.path(out_dir, outf), p, width = 7.5, height = 5.5)
  cat(sprintf("  Saved: %s\n", outf))
  rm(p_secure); gc(verbose = FALSE)
}

cat(sprintf("\n=== DONE: %s ===\n", Sys.time()))

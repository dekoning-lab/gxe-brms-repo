#!/usr/bin/env Rscript
############################################################
# 99_predicted_probs_v3.R
# Predicted probability plots for three key effects in the
# v3 sensitivity model (M4, with child behavior controls):
#   1. SLC6A3.10R gene main effect
#   2. Sensitivity x CNR1.77 interaction
#   3. Sensitivity x 5-HTTLPR interaction
#
# Pools predictions across multiple imputations for stability.
# Requires brms model fits in results/v3_m4/sens/
#
# Run from: the repository root (locally or on a cluster, if fits are available).
# Usage: Rscript scripts/99_predicted_probs_v3.R
############################################################

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(dplyr)
})

cat("=================================================================\n")
cat("PREDICTED PROBABILITY PLOTS — V3 SENSITIVITY MODEL\n")
cat("=================================================================\n\n")

out_dir <- "results/v3_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────
imputed_data <- readRDS("data/imputed_datasets_for_brms_m100_v2.rds")
dat <- imputed_data[[50]]  # reference dataset for covariate ranges

# ── Gene list and helper ─────────────────────────────────
genes <- c("BDNF", "CNR1.77", "CNR1.10", "DRD2", "DRD4",
           "MAOA", "SLC6A3.9R", "SLC6A3.10R", "X5HTTLPR")

get_mode <- function(x) as.numeric(names(sort(table(x), decreasing = TRUE))[1])

# ── Select imputations to pool ───────────────────────────
# Use 10 evenly spaced imputations for stability
imp_indices <- seq(5, 100, by = 10)
cat(sprintf("Pooling across %d imputations: %s\n\n",
            length(imp_indices), paste(imp_indices, collapse = ", ")))

# ── Base newdata (covariates at mean/mode) ───────────────
make_base_newdata <- function(n_rows) {
  data.frame(
    PC1 = 0, PC2 = 0, PC3 = 0,
    Sex = factor("1", levels = levels(dat$Sex)),
    infant_age = 0,
    pass = 0, diff = 0  # v3: child behavior controls at mean (z-scored)
  )
}

# ══════════════════════════════════════════════════════════
# PLOT 1: SLC6A3.10R gene main effect
# ══════════════════════════════════════════════════════════
cat("── Plot 1: SLC6A3.10R gene main effect ──\n")

gene_levels <- sort(unique(dat$SLC6A3.10R))
cat(sprintf("  SLC6A3.10R levels: %s\n", paste(gene_levels, collapse = ", ")))

nd <- data.frame()
for (gl in gene_levels) {
  row <- make_base_newdata(1)
  row$sens <- 0  # parenting at mean
  row$SLC6A3.10R <- gl
  for (og in setdiff(genes, "SLC6A3.10R")) row[[og]] <- get_mode(dat[[og]])
  row$gene_group <- gl
  nd <- rbind(nd, row)
}

# Pool predictions across imputations
all_preds <- list()
for (imp_i in imp_indices) {
  fit_file <- sprintf("results/v3_m4/sens/imp_%03d.rds", imp_i)
  if (!file.exists(fit_file)) {
    cat(sprintf("  Skipping imp %d (file not found)\n", imp_i))
    next
  }
  cat(sprintf("  Loading imp %d...\n", imp_i))
  fit <- readRDS(fit_file)
  pred <- fitted(fit, newdata = nd, summary = FALSE)  # [draws x obs x cats]
  all_preds[[length(all_preds) + 1]] <- pred
}

if (length(all_preds) > 0) {
  # Average across imputations
  pred_array <- abind::abind(all_preds, along = 1)  # [all_draws x obs x cats]
  n_cat <- dim(pred_array)[3]

  # Compute P(secure) = P(Y >= 4) for each draw and genotype
  p_secure <- apply(pred_array[, , 4:n_cat, drop = FALSE], c(1, 2), sum)

  plot_df <- data.frame(
    genotype = factor(rep(gene_levels, each = nrow(p_secure)),
                      labels = paste0(gene_levels, " copies")),
    p_secure_mean = colMeans(p_secure),
    p_secure_lo = apply(p_secure, 2, quantile, 0.025),
    p_secure_hi = apply(p_secure, 2, quantile, 0.975)
  )

  p1 <- ggplot(plot_df, aes(x = genotype, y = p_secure_mean)) +
    geom_pointrange(aes(ymin = p_secure_lo, ymax = p_secure_hi),
                    size = 1, linewidth = 0.8, colour = "#2166ac") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = expression(italic("SLC6A3") ~ "(10R): Predicted P(Secure Attachment)"),
         subtitle = "Other variables at mean/mode; pooled across imputations",
         x = "Number of 10R allele copies",
         y = "P(Attachment ≥ 4)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "predicted_probs_slc6a3_10r.pdf"),
         p1, width = 6, height = 5)
  cat("  Saved: predicted_probs_slc6a3_10r.pdf\n")
}

# ══════════════════════════════════════════════════════════
# PLOTS 2-3: Sensitivity x gene interactions
# ══════════════════════════════════════════════════════════
sens_range <- seq(min(dat$sens), max(dat$sens), length.out = 50)

for (gene_var in c("CNR1.77", "X5HTTLPR")) {
  cat(sprintf("\n── Plot: Sensitivity x %s interaction ──\n", gene_var))

  gene_levels <- sort(unique(dat[[gene_var]]))
  cat(sprintf("  %s levels: %s\n", gene_var, paste(gene_levels, collapse = ", ")))

  # Build prediction grid
  nd <- data.frame()
  for (gl in gene_levels) {
    for (s in sens_range) {
      row <- make_base_newdata(1)
      row$sens <- s
      row[[gene_var]] <- gl
      for (og in setdiff(genes, gene_var)) row[[og]] <- get_mode(dat[[og]])
      row$gene_group <- gl
      row$sens_val <- s
      nd <- rbind(nd, row)
    }
  }

  # Pool predictions
  all_preds <- list()
  for (imp_i in imp_indices) {
    fit_file <- sprintf("results/v3_m4/sens/imp_%03d.rds", imp_i)
    if (!file.exists(fit_file)) next
    cat(sprintf("  Loading imp %d...\n", imp_i))
    fit <- readRDS(fit_file)
    pred <- fitted(fit, newdata = nd, summary = FALSE)
    all_preds[[length(all_preds) + 1]] <- pred
  }

  if (length(all_preds) == 0) {
    cat("  No model fits found, skipping\n")
    next
  }

  pred_array <- abind::abind(all_preds, along = 1)
  n_cat <- dim(pred_array)[3]

  # P(secure) = P(Y >= 4)
  p_secure <- apply(pred_array[, , 4:n_cat, drop = FALSE], c(1, 2), sum)

  plot_df <- data.frame(
    sens = nd$sens_val,
    gene_group = nd$gene_group,
    p_mean = colMeans(p_secure),
    p_lo = apply(p_secure, 2, quantile, 0.025),
    p_hi = apply(p_secure, 2, quantile, 0.975)
  )

  # Gene labels
  if (gene_var == "X5HTTLPR") {
    plot_df$gene_label <- factor(plot_df$gene_group,
                                  levels = c(0, 1, 2),
                                  labels = c("L/L", "S/L", "S/S"))
    gene_italic <- expression(italic("5-HTTLPR"))
  } else if (gene_var == "CNR1.77") {
    plot_df$gene_label <- factor(plot_df$gene_group,
                                  levels = sort(unique(plot_df$gene_group)))
    gene_italic <- expression(italic("CNR1") ~ "(rs806377)")
  }

  p <- ggplot(plot_df, aes(x = sens, y = p_mean,
                            colour = gene_label, fill = gene_label)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = bquote("Predicted P(Secure Attachment): Sensitivity" ~ "\u00d7" ~ .(gene_var)),
         subtitle = "Other variables at mean/mode; pooled across imputations",
         x = "Maternal Sensitivity (z-scored)",
         y = "P(Attachment ≥ 4)",
         colour = gene_var, fill = gene_var) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  outf <- sprintf("predicted_probs_sens_%s.pdf",
                   tolower(gsub("\\.", "", gene_var)))
  ggsave(file.path(out_dir, outf), p, width = 8, height = 5.5)
  cat(sprintf("  Saved: %s\n", outf))
}

cat("\n=== DONE ===\n")

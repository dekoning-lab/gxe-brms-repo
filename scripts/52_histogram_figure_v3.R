#!/usr/bin/env Rscript
############################################################
# 52_histogram_figure_v3.R
# Posterior density plots for key effects in the
# sensitivity model (M4, parenting = sens)
#
# Four panels showing prior vs posterior:
#   - SLC6A3.10R (gene main, PP=84%) — strongest signal
#   - Sensitivity x CNR1.77 (interaction, PP=77%)
#   - Sensitivity x 5-HTTLPR (interaction, PP=79%)
#   - CNR1.77 (gene main, PP=54%) — near-null comparison
#
# Uses actual MCMC draws pooled across 100 imputations.
#
# Run from: the repository root.
# Usage: Rscript scripts/52_histogram_figure_v3.R
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

cat("=================================================================\n")
cat("POSTERIOR DENSITY PLOTS FOR KEY SENSITIVITY MODEL EFFECTS\n")
cat("=================================================================\n\n")

# ── Load draws ────────────────────────────────────────────
draws1 <- read.csv("results/v3_figures/posterior_draws_sens_sub.csv",
                   stringsAsFactors = FALSE)
draws2 <- read.csv("results/v3_figures/posterior_draws_cnr177_sub.csv",
                   stringsAsFactors = FALSE)
draws <- rbind(draws1, draws2)
cat(sprintf("Loaded %d total draws\n", nrow(draws)))

# ── Rename parameters for display ────────────────────────
draws$parameter <- recode(draws$parameter,
  "Sens x CNR1.77"  = "Sensitivity x CNR1.77",
  "Sens x 5-HTTLPR" = "Sensitivity x 5-HTTLPR"
)

# ── Summaries ─────────────────────────────────────────────
summ <- draws %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    pp = max(mean(value > 0), mean(value < 0)),
    n = n(),
    .groups = "drop"
  )
cat("\nParameter summaries:\n")
for (i in seq_len(nrow(summ))) {
  r <- summ[i, ]
  cat(sprintf("  %s: mean=%.3f, SD=%.3f, PP=%.1f%%, n=%d draws\n",
              r$parameter, r$mean, r$sd, r$pp * 100, r$n))
}

# Panel order: strongest signal first, null comparison last
panel_order <- c("SLC6A3.10R (gene main)", "Sensitivity x CNR1.77",
                 "Sensitivity x 5-HTTLPR", "CNR1.77 (gene main)")
summ <- summ %>% mutate(parameter = factor(parameter, levels = panel_order))
draws$parameter <- factor(draws$parameter, levels = panel_order)

# ── Hyperparameters for prior curves ──────────────────────
# From results/v3_m4/sens/summary.csv
sigma_gene <- 0.2039
sigma_int  <- 0.1935

prior_sigma <- c(
  "SLC6A3.10R (gene main)"    = sigma_gene,
  "Sensitivity x CNR1.77"     = sigma_int,
  "Sensitivity x 5-HTTLPR"    = sigma_int,
  "CNR1.77 (gene main)"       = sigma_gene
)

# ── Colour palette: blue for signals, grey for null ───────
fill_colors <- c(
  "SLC6A3.10R (gene main)"    = "#2166ac",
  "Sensitivity x CNR1.77"     = "#2166ac",
  "Sensitivity x 5-HTTLPR"    = "#2166ac",
  "CNR1.77 (gene main)"       = "#999999"
)

outline_colors <- c(
  "SLC6A3.10R (gene main)"    = "black",
  "Sensitivity x CNR1.77"     = "black",
  "Sensitivity x 5-HTTLPR"    = "black",
  "CNR1.77 (gene main)"       = "#777777"
)

mean_line_colors <- c(
  "SLC6A3.10R (gene main)"    = "#d62728",
  "Sensitivity x CNR1.77"     = "#d62728",
  "Sensitivity x 5-HTTLPR"    = "#d62728",
  "CNR1.77 (gene main)"       = "#666666"
)

# ── Build prior density curves ────────────────────────────
x_range <- range(draws$value)
x_pad <- diff(x_range) * 0.08
x_lims <- c(x_range[1] - x_pad, x_range[2] + x_pad)
x_seq <- seq(x_lims[1], x_lims[2], length.out = 500)

prior_df <- do.call(rbind, lapply(panel_order, function(p) {
  data.frame(
    parameter = factor(p, levels = panel_order),
    x = x_seq,
    y = dnorm(x_seq, mean = 0, sd = prior_sigma[p]),
    stringsAsFactors = FALSE
  )
}))

# ── Create figure ─────────────────────────────────────────
cat("\nCreating density panel figure...\n")

# Annotation data for posterior mean
# Nudge the null panel text further left to avoid overlap
ann <- summ %>%
  mutate(
    label = sprintf("Posterior mean = %.3f\nPP = %.0f%%", mean, pp * 100),
    nudge_x = ifelse(parameter == "CNR1.77 (gene main)", -0.15, 0)
  )

p <- ggplot(draws, aes(x = value)) +
  # Posterior density (filled)
  geom_density(aes(fill = parameter), alpha = 0.4, color = NA,
               adjust = 2) +
  # Posterior density outline (use outline_colors, not the color aes)
  geom_density(aes(group = parameter, colour = parameter), linewidth = 0.6,
               adjust = 2, fill = NA) +
  # Prior curve
  geom_line(data = prior_df, aes(x = x, y = y),
            color = "grey30", linewidth = 0.65, linetype = "dashed") +
  # Prior mean (zero)
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5,
             linetype = "dotted") +
  # Posterior mean line (per-panel colour via loop-generated geom layers)
  lapply(seq_len(nrow(summ)), function(i) {
    r <- summ[i, ]
    geom_vline(data = r, aes(xintercept = mean),
               colour = mean_line_colors[as.character(r$parameter)],
               linewidth = 0.7, linetype = "solid")
  }) +
  # Annotations (to the LEFT of the posterior mean line)
  geom_text(data = ann,
            aes(x = mean + nudge_x, y = Inf, label = label),
            vjust = 1.3, hjust = 1.08, size = 3, color = "grey20",
            lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~parameter, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_colour_manual(values = outline_colors, guide = "none") +
  coord_cartesian(xlim = x_lims) +
  labs(
    x = "Coefficient (log-odds)",
    y = "Posterior probability",
    title = "Prior vs. posterior distributions: Sensitivity model effects",
    subtitle = "Shaded = posterior (pooled across 100 imputations); dashed grey = prior N(0, sigma)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

ggsave("results/v3_figures/posterior_density_panels.pdf",
       p, width = 9, height = 6)
cat("  Saved: results/v3_figures/posterior_density_panels.pdf\n")

cat("\n=== DONE ===\n")

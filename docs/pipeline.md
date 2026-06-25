# Pipeline runbook

Sequential commands to reproduce the paper's main results. All commands
are run from the repository root.

The numbered prefixes on script names trace the chronological order of
the analysis in `scripts/`; the manuscript Supplementary Extended
Methods refers to specific scripts by these names. The pipeline is
organised into six phases plus a dedicated diagnostics + sensitivity
phase.

A single-shot end-to-end driver is provided as
`scripts/run_pipeline.sh` (locally) and `scripts/slurm/run_pipeline_slurm.sh`
(HPC). The phase-by-phase breakdown below documents what each step
does and where its output lands.

---

## Phase 0 — One-time setup

```bash
# Install R packages (see docs/INSTALL.md)
Rscript -e 'install.packages(c("brms","rstan","posterior","loo","mice",
                               "furrr","future","dplyr","purrr","tidyr",
                               "ggplot2","grid","gridExtra","ordinal"))'

# Install Python packages
pip install openpyxl pandas numpy
```

**Decide which input dataset to use.**

- **Real cohort data** (researchers with REB approval): place
  `Fulldata_with_PCs_and_maternal_PCs.txt` at `data/`.
- **Synthetic data** (everyone else): generate and copy in place:
  ```bash
  python3 scripts/synthetic/00_generate_synthetic_data.py
  cp data/synthetic/Fulldata_synthetic.txt data/Fulldata_with_PCs_and_maternal_PCs.txt
  ```

---

## Phase 1 — Data harmonization and imputation

Joins child + maternal ancestry PCs (already harmonized in the shipped
input file), runs MICE multiple imputation (M=100, 30 iterations, pmm),
and adds the two child-behaviour covariates (`pass`, `diff`) to the
imputed datasets.

```bash
Rscript scripts/04_run_mice_imputation.R
Rscript scripts/45_augment_imputed_data.R
```

Outputs:
- `data/mice_imputed_m100_maxit30.rds` — raw MICE object
- `data/imputed_datasets_for_brms_m100.rds` — 100 imputed datasets
- `data/imputed_datasets_for_brms_m100_v2.rds` — augmented with `pass`, `diff`
- `results/mice_imputation_report.txt`
- `results/mice_full_convergence.pdf`

---

## Phase 2 — Bayesian ordinal regression (HPC recommended)

Fits the four model levels for each parenting variable and prior scale.
Each script accepts the parenting variable as its first argument.

```bash
# On a single machine (slow — ~ tens of hours total):
for pv in sens cont unre; do
  Rscript scripts/52_fit_m1_v3.R $pv
  Rscript scripts/52_fit_m2_v3.R $pv
  Rscript scripts/52_fit_m3_v3.R $pv
  Rscript scripts/52_fit_m4_v3.R $pv
done
# Prior sensitivity for M4 at the two off-primary scales:
for s in 0.25 1.00; do
  for pv in sens cont unre; do
    Rscript scripts/52_fit_m4_v3_sensitivity.R $pv $s 1 100
  done
done

# Or on SLURM (submit.sh injects --partition / --account from config.sh):
scripts/slurm/submit.sh scripts/slurm/52_slurm_m1_v3.sh
scripts/slurm/submit.sh scripts/slurm/52_slurm_m2_v3.sh
scripts/slurm/submit.sh scripts/slurm/52_slurm_m3_v3.sh
scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3.sh
scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3_s025.sh
scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3_s100.sh
```

Outputs:
- `results/v3_m{1,2,3,4}/{sens,cont,unre}/imp_NNN.rds` — 100 saved brms fits per cell
- `results/v3_m{1,2,3,4}/{sens,cont,unre}/loo_NNN.rds`
- `results/v3_m{1,2,3,4}/{sens,cont,unre}/summary.csv` — pooled posterior summary
- `results/v3_m4_s{025,100}/{sens,cont,unre}/` — prior-sensitivity variants

---

## Phase 3 — Aggregation and model comparison

```bash
Rscript scripts/52_aggregate_v3.R              # cross-model LOO table
Rscript scripts/52_aggregate_sensitivity_v3.R  # prior-sensitivity aggregation
Rscript scripts/52_loo_rps_v3.R                # LOO + RPS + binary-threshold scoring
Rscript scripts/52_compute_r2_v3.R sens
Rscript scripts/52_compute_r2_v3.R cont
Rscript scripts/52_compute_r2_v3.R unre
Rscript scripts/52_aggregate_r2_v3.R
```

Outputs:
- `results/v3_pairwise.csv` — pairwise model comparison
- `results/v3_loo_rps/` — predictive performance results
- `results/v3_r2/` — Bayesian and McKelvey–Zavoina R²

---

## Phase 4 — Posterior extraction and figures

```bash
Rscript scripts/52_extract_draws_v3.R       # pooled posterior draws for figures
Rscript scripts/52_histogram_figure_v3.R    # Figure 2: prior + posterior densities
Rscript scripts/35_evidence_ratios.R        # Figure 3: directional Bayes factors
Rscript scripts/99_predicted_probs_v3.R     # Figure 4: predicted-prob plots
Rscript scripts/52_tables_figures_v3.R      # Tables 1, 2 + forest plots
```

Outputs:
- `results/v3_figures/` — all figure PDFs
- `results/latex_tables/` — LaTeX-ready tables

---

## Phase 5 — Proportional-odds and binary-cutpoint checks (optional)

```bash
Rscript scripts/99_test_proportional_odds.R
```

For the binary-cutpoint sensitivity analysis (Stan model integrating
out the secure/insecure cutpoint), see
`scripts/51_fit_binary_cutpoint.R` and the supplementary methods.

---

## Phase 6 — MCMC diagnostics

Extracts comprehensive convergence diagnostics from every saved fit:

```bash
for pv in sens cont unre; do
  Rscript scripts/53_diagnostics_v3_m4.R $pv
done
# Or via SLURM (3-task array):
scripts/slurm/submit.sh scripts/slurm/53_slurm_diagnostics_v3.sh
```

Outputs:
- `results/diagnostics_v3/{sens,cont,unre}/diagnostics_per_fit.csv`
- `results/diagnostics_v3/{sens,cont,unre}/diagnostics_per_param.csv`
- `results/diagnostics_v3/{sens,cont,unre}/diagnostics_param_summary.csv`
- `results/diagnostics_v3/{sens,cont,unre}/diagnostics_summary.txt`

---

## Phase 7 — Outlier-exclusion sensitivity

Re-pools the posterior with the single most-pathological fit per
parenting model removed (identified from Phase 6).

```bash
for pv in sens cont unre; do
  Rscript scripts/54_sensitivity_excl_outlier.R $pv
done
# Or via SLURM:
scripts/slurm/submit.sh scripts/slurm/54_slurm_sensitivity_excl.sh

python3 scripts/55_compare_excl_outlier.py
```

Outputs:
- `results/v3_m4_excl_outlier/{sens,cont,unre}/summary.csv`
- `results/v3_m4_excl_outlier/comparison_original_vs_excluded.csv`
- `results/v3_m4_excl_outlier/comparison_report.txt`

---

## Phase 8 — Alternative 5-HTTLPR coding (optional)

Re-runs the entire pipeline with the Hu et al. (2006)
triallelic-functional coding instead of the biallelic S-count coding.
Requires the original Kobor lab triallelic genotype spreadsheet
(`data/Kobor_lab_genotype_summary_5HTTLPR_triallelic.xlsx`).

```bash
python3 scripts/60_recode_5httlpr_functional.py
scripts/slurm/submit.sh scripts/slurm/61_slurm_mice_alt5httlpr.sh
scripts/slurm/submit.sh scripts/slurm/63_slurm_m4_alt5httlpr.sh    # depends on 61
python3 scripts/64_compare_alt5httlpr.py
```

Outputs:
- `data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt`
- `results_alt_5httlpr/v3_m4/{sens,cont,unre}/summary.csv`
- `results_alt_5httlpr/comparison_original_vs_alt5httlpr.csv`

---

## Mapping to manuscript objects

| Manuscript element | Produced by |
|---|---|
| Table 1 (main coefficients)             | `52_tables_figures_v3.R` |
| Table 2 (prior sensitivity)             | `52_aggregate_sensitivity_v3.R` + `52_tables_figures_v3.R` |
| Figure 1 (model schematic)              | manually drawn (TikZ in manuscript source) |
| Figure 2 (posterior densities)          | `52_histogram_figure_v3.R` |
| Figure 3 (evidence calibration)         | `35_evidence_ratios.R` |
| Figure 4 (predicted probabilities)      | `99_predicted_probs_v3.R` |
| Table S5 + S6 (MCMC diagnostics)        | `53_diagnostics_v3_m4.R` |
| Table S7 (outlier-exclusion sensitivity)| `55_compare_excl_outlier.py` |
| Table 4 (variant codebook)              | static (in manuscript source) |

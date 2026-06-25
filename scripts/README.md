# Scripts

All analysis code for the paper. Scripts are organised by the
chronological phase of the pipeline; their numeric prefixes match the
labels used in the manuscript's Supplementary Extended Methods.

## Phase 1: Data harmonization and imputation
| Script | Purpose |
|---|---|
| `01_harmonize_child_pcs.py`        | Join infant ancestry PCs (PC1–PC50) to phenotype dataset |
| `02_harmonize_maternal_pcs.py`     | Join maternal ancestry PCs as auxiliary imputation variables |
| `03_validate_harmonization.py`     | Cross-check join integrity, row counts, duplicate IDs |
| `04_run_mice_imputation.R`         | MICE: m=100 imputations, maxit=30, pmm, parallel |
| `45_augment_imputed_data.R`        | Add child-behaviour covariates (`pass`, `diff`) to each imputed dataset |

## Phase 2: Bayesian ordinal regression (brms)
| Script | Purpose |
|---|---|
| `52_fit_m1_v3.R`                   | M1: demographics + child-behaviour covariates only |
| `52_fit_m2_v3.R`                   | M2: adds the parenting variable |
| `52_fit_m3_v3.R`                   | M3: adds hierarchical gene main effects |
| `52_fit_m4_v3.R`                   | M4: full GxE model (gene mains + gene × parenting interactions) |
| `52_fit_m4_v3_sensitivity.R`       | M4 under non-primary hyperprior scales (s = 0.25, 1.00) |

Each fit script takes the parenting variable (`sens`, `cont`, or
`unre`) as its first argument; the sensitivity variant additionally
takes the hyperprior scale and optionally a start/end imputation range.

## Phase 3: Aggregation and model comparison
| Script | Purpose |
|---|---|
| `52_aggregate_v3.R`                | LOO pairwise model comparison (M1 vs M2, M3 vs M4, etc.) |
| `52_aggregate_sensitivity_v3.R`    | Pool prior-sensitivity variants |
| `52_loo_rps_v3.R`                  | LOO-CV + ranked-probability score + binary-threshold log-score |
| `52_aggregate_rps_v3.R`            | Aggregate the RPS computations |
| `52_compute_r2_v3.R`               | Per-imputation Bayesian and McKelvey–Zavoina R² |
| `52_compute_rpss_v3.R`             | Per-imputation ranked-probability skill scores |
| `52_aggregate_r2_v3.R`             | Pool R² across imputations |

## Phase 4: Posterior extraction and figures
| Script | Purpose |
|---|---|
| `52_extract_draws_v3.R`            | Concatenate pooled posterior draws for key parameters |
| `52_histogram_figure_v3.R`         | Figure 2: prior + posterior density panels |
| `52_tables_figures_v3.R`           | Tables 1, 2; forest plots |
| `35_evidence_ratios.R`             | Figure 3: directional Bayes-factor evidence calibration |
| `99_predicted_probs_v3.R`          | Figure 4: predicted-probability plots |
| `99_predicted_probs_v3_cluster.R`  | HPC variant of the predicted-probability extraction (computationally heavier) |

## Phase 5: Robustness checks
| Script | Purpose |
|---|---|
| `99_test_proportional_odds.R`      | LR tests for proportional-odds assumption on sampled imputations |

## Phase 6: MCMC diagnostics
| Script | Purpose |
|---|---|
| `53_diagnostics_v3_m4.R`           | Per-fit and per-parameter R-hat, ESS, divergences, BFMI from saved fits |

## Phase 7: Outlier-exclusion sensitivity
| Script | Purpose |
|---|---|
| `54_sensitivity_excl_outlier.R`    | Re-pool with the catastrophic-pathology fit excluded |
| `55_compare_excl_outlier.py`       | Original vs excluded posterior summary comparison |

## Phase 8: Alternate 5-HTTLPR coding
| Script | Purpose |
|---|---|
| `60_recode_5httlpr_functional.py`  | Build the Hu et al. (2006) triallelic-functional dataset from Kobor lab Excel |
| `61_run_mice_alt5httlpr.R`         | MICE on alt-coded data |
| `62_augment_imputed_alt5httlpr.R`  | Augment with pass/diff (alt) |
| `63_fit_m4_alt5httlpr.R`           | M4 fit on alt-coded data |
| `64_compare_alt5httlpr.py`         | Primary vs alt-coding posterior comparison |

## SLURM drivers — `slurm/`

Cluster job-submission scripts for the compute-heavy phases. Each
mirrors the underlying R/Python script. Cluster-specific settings
(partition, account, module loads, project root) are read from the
top-level `config.sh` — edit that one file to match your site. Submit
jobs with `slurm/submit.sh`, which injects `--partition`/`--account`
from `config.sh`, e.g.:

    scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3.sh

## Synthetic data generator — `synthetic/`

`synthetic/00_generate_synthetic_data.py` produces a
`Fulldata_synthetic.txt` file matching the input schema. See the
top-level `docs/DATA.md` for details.

## Single-shot driver — `run_pipeline.sh`

`run_pipeline.sh` is a thin shell wrapper that invokes the local
(non-SLURM) pipeline in order. It is intended for smoke-testing on
the synthetic dataset; for production use, submit each phase as
a SLURM job (or, on the cluster, source `slurm/run_pipeline_slurm.sh`).

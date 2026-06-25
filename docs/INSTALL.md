# Installation and environment

The pipeline runs on R 4.4.x and Python 3.11 (the versions used for the
paper). Model fitting on a workstation is feasible but slow; we
recommend an HPC environment with SLURM for the full M1–M4 production
runs.

## Software versions

The paper used:

| Tool / package | Version |
|---|---|
| R              | 4.4.1   |
| Stan (via rstan) | 2.32.6 |
| brms           | 2.21.0  |
| posterior      | 1.6.0   |
| loo            | 2.8.0   |
| mice           | 3.16.0  |
| futuremice (via mice/furrr) | bundled |
| Python         | 3.11    |
| openpyxl       | 3.x     |
| pandas         | 2.x     |
| numpy          | 1.26.x  |
| pdflatex (TeX Live) | 2025 |

## R packages

Install everything used by the production pipeline:

```r
install.packages(c(
  "brms",        # Bayesian regression in Stan
  "rstan",       # Stan backend (compile-toolchain dependent)
  "posterior",   # draws manipulation + diagnostics
  "loo",         # PSIS-LOO
  "mice",        # multiple imputation
  "furrr",       # parallel mice (used by futuremice)
  "future",
  "dplyr",
  "purrr",
  "tidyr",
  "ggplot2",
  "grid",
  "gridExtra",
  "ordinal"      # ordinal::clm for proportional-odds check (Phase 6)
))
```

`rstan` requires a working C++ toolchain. On Linux, GCC ≥ 13 is
recommended (some pre-built Stan binaries need `GLIBCXX_3.4.32`); on
macOS, install the Xcode command line tools.

## Python packages

```bash
pip install -r requirements.txt
```

or directly:

```bash
pip install openpyxl pandas numpy
```

## HPC environment (optional, recommended)

For the full M1–M4 production runs we recommend a SLURM-managed
cluster. The SLURM driver scripts in `scripts/slurm/` are written to be
site-agnostic: instead of hard-coding any partition, account, or
module names, they read those settings from `config.sh` at the
repository root. To adapt the pipeline to your cluster, set the
relevant variables there:

- `PARTITION` — the SLURM partition (queue) to submit to (leave empty
  for your cluster's default).
- `ACCOUNT` — your SLURM account / allocation, if your site requires
  one (leave empty if not).
- `MODULE_LOADS` — the environment-module command(s) that put R (with
  a C++17 toolchain) on `PATH`. This is highly site-specific; leave it
  empty if R is already available. As an **example only**, on one
  particular Lmod-based cluster this was:
  `MODULE_LOADS="module load R/4.4.1 gcc/13.3.0"` — substitute the
  exact module names for your own system (or omit entirely).

Submit jobs through the wrapper, which injects `--partition` /
`--account` from `config.sh`:

```bash
scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3.sh
```

Any extra arguments (e.g. `--dependency=...`) are passed straight
through to `sbatch`. R **≥ 4.4** with a working **C++17** toolchain is
required for `rstan`/`brms`.

Compute requirements per M4 fit (16 000 post-warmup draws, 4 chains
NUTS, hierarchical funnel with `adapt_delta = 0.99`,
`max_treedepth = 15`): roughly 3–5 minutes on 4 CPU cores, ≤8 GB RAM.

The complete pipeline (M1–M4 across three parenting variables and 100
imputed datasets, three prior scales for M4, plus LOO and aggregation)
took approximately 4 300 core-hours over ~2 000 SLURM array tasks.

## Local-only execution

The smaller pieces of the pipeline — MICE imputation, post-fit
aggregation, comparison scripts, figure generation — run comfortably on
a laptop. Full Bayesian fitting locally is feasible if you have
patience: M4 fits at the production sampler budget are 3–5 minutes
each, 100 imputations × 3 parenting models = 300 fits ≈ 15–25 hours
serially. Use `scripts/04_run_mice_imputation.R`'s `n_cores` setting
(or its equivalents) to parallelise where possible.

#!/usr/bin/env bash
###############################################################################
# run_pipeline.sh
#
# Sequential local execution of the full analysis pipeline. Intended for
# smoke-testing on the synthetic dataset; production runs of the brms
# fits (Phase 2) should use the SLURM drivers in scripts/slurm/.
#
# Usage (from the repository root):
#
#   # smoke test on synthetic data:
#   bash scripts/run_pipeline.sh smoke
#
#   # real pipeline (Phase 2 will take many hours locally):
#   bash scripts/run_pipeline.sh full
###############################################################################
set -euo pipefail

MODE="${1:-smoke}"
[[ "$MODE" == "smoke" || "$MODE" == "full" ]] || {
  echo "Usage: $0 [smoke|full]" >&2; exit 1
}

cd "$(dirname "$0")/.."

# Load local settings (data filename, and cluster settings if present).
source ./config.sh

echo "================================================================="
echo "PIPELINE: mode=$MODE  ($(date))"
echo "================================================================="

# ---- Phase 0: synthetic-data staging (smoke mode only) ----------------
if [[ "$MODE" == "smoke" ]]; then
  python3 scripts/synthetic/00_generate_synthetic_data.py
  cp data/synthetic/Fulldata_synthetic.txt data/Fulldata_with_PCs_and_maternal_PCs.txt
  echo
fi

# ---- Phase 1: MICE imputation + augmentation -------------------------
echo "## Phase 1: MICE + augmentation"
Rscript scripts/04_run_mice_imputation.R
Rscript scripts/45_augment_imputed_data.R
echo

# ---- Phase 2: brms fits ----------------------------------------------
if [[ "$MODE" == "smoke" ]]; then
  echo "## Phase 2 (smoke): single M4 fit at reduced budget"
  Rscript scripts/synthetic/01_smoke_test_m4.R sens
  echo
  echo "[Smoke complete.]  See docs/pipeline.md for the production pipeline."
  exit 0
fi

echo "## Phase 2 (full): M1-M4 and prior sensitivity, all parenting vars"
echo "   This will take many hours on a single machine. Consider"
echo "   submitting scripts/slurm/52_slurm_m4_v3.sh on a cluster instead."

for pv in sens cont unre; do
  Rscript scripts/52_fit_m1_v3.R "$pv"
  Rscript scripts/52_fit_m2_v3.R "$pv"
  Rscript scripts/52_fit_m3_v3.R "$pv"
  Rscript scripts/52_fit_m4_v3.R "$pv"
done

for s in 0.25 1.00; do
  for pv in sens cont unre; do
    Rscript scripts/52_fit_m4_v3_sensitivity.R "$pv" "$s" 1 100
  done
done

# ---- Phase 3: aggregation + model comparison --------------------------
echo "## Phase 3: aggregation, LOO, RPS, R^2"
Rscript scripts/52_aggregate_v3.R
Rscript scripts/52_aggregate_sensitivity_v3.R
Rscript scripts/52_loo_rps_v3.R
for pv in sens cont unre; do
  Rscript scripts/52_compute_r2_v3.R "$pv"
done
Rscript scripts/52_aggregate_r2_v3.R

# ---- Phase 4: figures and tables --------------------------------------
echo "## Phase 4: extract draws, build figures and tables"
Rscript scripts/52_extract_draws_v3.R
Rscript scripts/52_histogram_figure_v3.R
Rscript scripts/35_evidence_ratios.R
Rscript scripts/99_predicted_probs_v3.R
Rscript scripts/52_tables_figures_v3.R

# ---- Phase 5: proportional-odds check ---------------------------------
echo "## Phase 5: proportional-odds assumption check"
Rscript scripts/99_test_proportional_odds.R

# ---- Phase 6: MCMC diagnostics ----------------------------------------
echo "## Phase 6: MCMC diagnostics"
for pv in sens cont unre; do
  Rscript scripts/53_diagnostics_v3_m4.R "$pv"
done

# ---- Phase 7: outlier-exclusion sensitivity ---------------------------
echo "## Phase 7: outlier-exclusion sensitivity"
for pv in sens cont unre; do
  Rscript scripts/54_sensitivity_excl_outlier.R "$pv"
done
python3 scripts/55_compare_excl_outlier.py

echo
echo "================================================================="
echo "PIPELINE complete.  Results under results/."
echo "================================================================="

#!/bin/bash
#SBATCH --job-name=mice_alt5httlpr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=02:00:00
#SBATCH --mem=64G
# Partition/account are injected at submit time from config.sh via
# scripts/slurm/submit.sh. To hard-code a partition instead, uncomment:
# #SBATCH --partition=<your_partition>
#SBATCH --output=logs/mice_alt5httlpr_%j.out
#SBATCH --error=logs/mice_alt5httlpr_%j.err

# Run from: the repository root (the job cd's to $PROJECT_ROOT from config.sh).
# Runs MICE (m=100, maxit=30, parallel over 18 cores) then augments
# with pass/diff. Produces data_alt_5httlpr/imputed_datasets_for_brms_m100_v2.rds.

set -euo pipefail

# Load local/cluster settings (project root + module loads from config.sh;
# partition/account were applied at submit time by scripts/slurm/submit.sh).
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
cd "$PROJECT_ROOT"

source /etc/profile.d/modules.sh 2>/dev/null || true
eval "${MODULE_LOADS:-:}"   # site-specific module load(s), from config.sh

mkdir -p logs data_alt_5httlpr results_alt_5httlpr

echo "============================================================"
echo "ALT 5-HTTLPR: MICE + augment"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

Rscript scripts/61_run_mice_alt5httlpr.R
Rscript scripts/62_augment_imputed_alt5httlpr.R

echo "Finished: $(date)"

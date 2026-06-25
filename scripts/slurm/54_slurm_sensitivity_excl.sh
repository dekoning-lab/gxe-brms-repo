#!/bin/bash
#SBATCH --job-name=sens_excl
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=24G
# Partition/account are injected at submit time from config.sh via
# scripts/slurm/submit.sh. To hard-code a partition instead, uncomment:
# #SBATCH --partition=<your_partition>
#SBATCH --output=logs/sens_excl_%A_%a.out
#SBATCH --error=logs/sens_excl_%A_%a.err

# Run from: the repository root (the job cd's to $PROJECT_ROOT from config.sh).
# Submits 3-task array: 1=sens (excl imp 3), 2=cont (excl imp 1), 3=unre (excl imp 38)

set -euo pipefail

# Load local/cluster settings (project root + module loads from config.sh;
# partition/account were applied at submit time by scripts/slurm/submit.sh).
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
cd "$PROJECT_ROOT"

source /etc/profile.d/modules.sh 2>/dev/null || true
eval "${MODULE_LOADS:-:}"   # site-specific module load(s), from config.sh

mkdir -p logs results/v3_m4_excl_outlier

case "$SLURM_ARRAY_TASK_ID" in
  1) PARENTING=sens ;;
  2) PARENTING=cont ;;
  3) PARENTING=unre ;;
  *) echo "Bad task id $SLURM_ARRAY_TASK_ID"; exit 1 ;;
esac

echo "============================================================"
echo "Task $SLURM_ARRAY_TASK_ID: outlier-exclusion sensitivity for $PARENTING"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "============================================================"

Rscript scripts/54_sensitivity_excl_outlier.R "$PARENTING"

echo "Finished: $(date)"

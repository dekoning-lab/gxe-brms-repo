#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=16G
# Partition/account are injected at submit time from config.sh via
# scripts/slurm/submit.sh. To hard-code a partition instead, uncomment:
# #SBATCH --partition=<your_partition>
#SBATCH --job-name=v3_rps
#SBATCH --array=1-300
#SBATCH --output=logs/v3_rps_%A_%a.out
#SBATCH --error=logs/v3_rps_%A_%a.err

# 300 tasks: 3 parenting × 100 imputations (1 per task)
# Loads pre-fitted M1, M2, M3, M4 and computes RPS + binary ELPD
# ~10-20 min per task (load 4 models + posterior_epred + LOO scoring)
#
# Tasks   1-100: sens
# Tasks 101-200: cont
# Tasks 201-300: unre

# Load local/cluster settings (project root + module loads from config.sh;
# partition/account were applied at submit time by scripts/slurm/submit.sh).
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

source /etc/profile.d/modules.sh 2>/dev/null || true
eval "${MODULE_LOADS:-:}"   # site-specific module load(s), from config.sh

cd "$PROJECT_ROOT"

IDX=$((SLURM_ARRAY_TASK_ID - 1))
P_IDX=$((IDX / 100))
IMP=$((IDX % 100 + 1))

case $P_IDX in
  0) PV="sens" ;;
  1) PV="cont" ;;
  2) PV="unre" ;;
  *) echo "ERROR: invalid P_IDX=$P_IDX"; exit 1 ;;
esac

echo "RPS v3: parenting=${PV}, imp=${IMP} (task ${SLURM_ARRAY_TASK_ID})"
Rscript scripts/52_loo_rps_v3.R ${PV} ${IMP} ${IMP}

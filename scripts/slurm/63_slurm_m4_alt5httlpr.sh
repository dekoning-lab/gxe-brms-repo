#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=16G
# Partition/account are injected at submit time from config.sh via
# scripts/slurm/submit.sh. To hard-code a partition instead, uncomment:
# #SBATCH --partition=<your_partition>
#SBATCH --job-name=m4alt5
#SBATCH --array=1-300
#SBATCH --output=logs/m4_alt5httlpr_%A_%a.out
#SBATCH --error=logs/m4_alt5httlpr_%A_%a.err

# 300 tasks: 3 parenting × 100 imputations (1 per task), mirroring the
# production v3 layout in 52_slurm_m4_v3.sh.
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

echo "ALT 5-HTTLPR M4: parenting=${PV}, imp=${IMP} (task ${SLURM_ARRAY_TASK_ID})"
Rscript scripts/63_fit_m4_alt5httlpr.R ${PV} ${IMP} ${IMP}

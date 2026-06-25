#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
# Partition/account are injected at submit time from config.sh via
# scripts/slurm/submit.sh. To hard-code a partition instead, uncomment:
# #SBATCH --partition=<your_partition>
#SBATCH --job-name=v3_m3
#SBATCH --array=1-100
#SBATCH --output=logs/v3_m3_%A_%a.out
#SBATCH --error=logs/v3_m3_%A_%a.err

# 100 tasks: 1 imputation each (maximum parallelism)
# M3 (v3): PC1-3 + Sex + infant_age + pass + diff + 9 genes (hierarchical)
# ~20-30 min per fit + LOO

# Load local/cluster settings (project root + module loads from config.sh;
# partition/account were applied at submit time by scripts/slurm/submit.sh).
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

source /etc/profile.d/modules.sh 2>/dev/null || true
eval "${MODULE_LOADS:-:}"   # site-specific module load(s), from config.sh

cd "$PROJECT_ROOT"

IMP=$SLURM_ARRAY_TASK_ID
echo "M3 v3: imputation ${IMP} (task ${SLURM_ARRAY_TASK_ID})"
Rscript scripts/52_fit_m3_v3.R ${IMP} ${IMP}

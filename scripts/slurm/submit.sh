#!/usr/bin/env bash
###############################################################################
# submit.sh — thin wrapper around `sbatch` that injects cluster settings
#             (partition, account) from config.sh, so the individual job
#             scripts stay free of any site-specific names.
#
# Usage (from the repository root):
#   scripts/slurm/submit.sh scripts/slurm/52_slurm_m4_v3.sh
#   scripts/slurm/submit.sh --dependency=afterok:123 scripts/slurm/63_slurm_m4_alt5httlpr.sh
#
# Any extra arguments are passed straight through to sbatch, so you can still
# add --dependency, --time, etc. Command-line --partition/--account override
# the values from config.sh.
###############################################################################
set -euo pipefail

# Repository root = two levels up from this script (scripts/slurm/ -> repo).
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Load settings (prefer an un-tracked config.local.sh if present).
if [[ -f "$REPO_ROOT/config.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "$REPO_ROOT/config.local.sh"
else
  # shellcheck source=/dev/null
  source "$REPO_ROOT/config.sh"
fi

# Build optional flags only when the corresponding setting is non-empty.
extra_flags=()
[[ -n "${PARTITION:-}" ]] && extra_flags+=(--partition="$PARTITION")
[[ -n "${ACCOUNT:-}"   ]] && extra_flags+=(--account="$ACCOUNT")

echo "submit.sh: sbatch ${extra_flags[*]:-(cluster defaults)} $*" >&2
exec sbatch "${extra_flags[@]}" "$@"

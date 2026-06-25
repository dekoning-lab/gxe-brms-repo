#!/usr/bin/env bash
###############################################################################
# config.sh — local / cluster settings for the analysis pipeline.
#
# Edit the values below to match YOUR environment. This file is sourced by:
#   * scripts/run_pipeline.sh        (local end-to-end runs)
#   * scripts/slurm/submit.sh        (injects --partition / --account at submit)
#   * scripts/slurm/*.sh             (module loads + project root, at run time)
#   * the Python harmonization scripts (via the exported RAW_DATA_FILE)
#
# It contains NO secrets and is safe to commit. If you prefer to keep your
# personal cluster settings out of version control, copy it to config.local.sh
# and source that instead.
#
# Every value uses the "${VAR:-default}" idiom, so you can also override any
# setting from the environment without editing this file, e.g.:
#   PARTITION=cpu2023 ACCOUNT=def-yourname sbatch ...
###############################################################################

# --- HPC / SLURM ------------------------------------------------------------

# SLURM partition (queue) to submit jobs to. Leave EMPTY to use your cluster's
# default partition. scripts/slurm/submit.sh passes this as --partition=<...>.
# Example: PARTITION="compute"
export PARTITION="${PARTITION:-}"

# SLURM account / allocation, if your cluster requires one (many Slurm sites do,
# e.g. "def-yourname"). Leave EMPTY if not needed. submit.sh passes it as
# --account=<...>.
export ACCOUNT="${ACCOUNT:-}"

# Directory the batch jobs change into before running. Defaults to the directory
# you ran `sbatch` from (SLURM_SUBMIT_DIR), which is correct if you always submit
# from the repository root. Override only if your jobs must run elsewhere.
export PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"

# Environment-module command(s) that make R (with a C++17 toolchain for
# rstan/brms) available on your cluster. This is highly site-specific. Set it to
# the exact `module load ...` line(s) for your system, or leave EMPTY to skip
# (e.g. if R is already on your PATH). Multiple commands can be separated by ';'.
# Example (Lmod):  MODULE_LOADS="module load R/4.4.1 gcc/13.3.0"
export MODULE_LOADS="${MODULE_LOADS:-}"

# --- Data -------------------------------------------------------------------

# Path (relative to the repository root) of the raw input dataset consumed by
# the harmonization scripts (01–03). The synthetic smoke test does NOT use this
# file. Place the real, private dataset here, or point this at its location.
export RAW_DATA_FILE="${RAW_DATA_FILE:-data/Fulldata_raw.txt}"

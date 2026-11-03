#!/bin/bash
# Purpose: SLURM array wrapper for 55_calibrate-ccs-gate.R; task i fits one regime x seed of <run_dir>/jobs.tsv.
# Run:     cd <repo root>; sbatch --array=1-<n_jobs> scripts/slurm_ccs_gate.sh <run_dir>
# Inputs:  <run_dir>/jobs.tsv (regime_id, seed), <run_dir>/regimes.tsv [55 setup]
# Outputs: <run_dir>/reps/<regime_id>/seed_NN.tsv; logs in <run_dir>/logs
#SBATCH --job-name=cj_ccsgate
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail

RUN_DIR=${1:?usage: sbatch --array=1-N slurm_ccs_gate.sh <run_dir>}
TASK=${SLURM_ARRAY_TASK_ID:?run through sbatch --array}
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/job${SLURM_ARRAY_JOB_ID}_task${TASK}.out" 2>&1

R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
read -r REGIME SEED < <(awk -F'\t' -v t="$TASK" 'NR == t + 1 { print $1, $2 }' "$RUN_DIR/jobs.tsv")
echo "task $TASK: $REGIME seed $SEED on $(hostname)"
CJ_PKG_DIR=${CJ_PKG_DIR:-cytoJointR-wt/ccs-gate} NOT_CRAN=true PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/55_calibrate-ccs-gate.R "$RUN_DIR" fit "$REGIME" "$SEED"

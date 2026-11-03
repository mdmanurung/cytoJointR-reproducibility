#!/bin/bash
# Purpose: RV144 posterior predictive check (60_diagnose-estimand.R fit then stats) for one candidate.
# Run:     CJ_PKG_DIR=<worktree> CJ_FIT_ARGS=<cj_fit args> sbatch scripts/slurm_ppc_c9.sh <run_dir>
# Outputs: <run_dir>/rv144/ (see 60_diagnose-estimand.R); log in <run_dir>/logs
#SBATCH --job-name=cj_ppc_c9
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/dev/null
set -uo pipefail
RUN_DIR=${1:?run_dir}
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/ppc_${SLURM_JOB_ID}.out" 2>&1
R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
export PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib
Rscript scripts/60_diagnose-estimand.R fit "$RUN_DIR" rv144 && Rscript scripts/60_diagnose-estimand.R stats "$RUN_DIR" rv144

#!/bin/bash
# Purpose: SLURM array wrapper for the fit stage of 60_diagnose-estimand.R; task i fits data set i.
# Run:     cd <repo root>; sbatch --array=1-4 scripts/slurm_estimand.sh runs/2026-11-03_estimand
# Inputs:  data/processed/{rv144,tb_s,tb_u,tb_simplecompass}/
# Outputs: <run_dir>/<dataset>/ (see 60_diagnose-estimand.R); logs in <run_dir>/logs
#SBATCH --job-name=cj_estimand
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail
RUN_DIR=${1:?run_dir}
DATASETS=(rv144 tb_s tb_u tb_simplecompass)
DS=${DATASETS[$((${SLURM_ARRAY_TASK_ID:?run through sbatch --array} - 1))]}
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/fit_${DS}_${SLURM_JOB_ID}.out" 2>&1
R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/60_diagnose-estimand.R fit "$RUN_DIR" "$DS"

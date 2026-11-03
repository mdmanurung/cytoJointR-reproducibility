#!/bin/bash
# Purpose: SLURM array wrapper for the SARS-CoV-2 ICS analysis; task 1 = CD4, task 2 = CD8. Each task runs stage
#          cytojointr (environment A) then stage polyics (environment B) of 32_analyse-sarscov2.R.
# Run:     cd <repo root>; sbatch --array=1-2 scripts/slurm_realdata_sarscov2.sh runs/2026-11-03_realdata/sarscov2
#          Compile Stan once first (see slurm_benchmark.sh). Do not edit scripts or the package while tasks run.
# Inputs:  data/processed/sarscov2/counts.tsv from 24_prepare-sarscov2.R
# Outputs: <run_dir>/<pop>/ (see 32); logs in <run_dir>/logs
#SBATCH --job-name=cj_sarscov2
#SBATCH --partition=medium   # 4 h limit; the QOS caps a user at 4 CPUs and 8G, so the two array tasks run one after the other (partition all queued for hours)
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/dev/null
set -uo pipefail

RUN_DIR=${1:?usage: sbatch --array=1-2 slurm_realdata_sarscov2.sh <run_dir>}
POPS=(CD4 CD8)
POP=${POPS[$((${SLURM_ARRAY_TASK_ID:?run through sbatch --array} - 1))]}

cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/job${SLURM_ARRAY_JOB_ID}_${POP}.out" 2>&1

R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51   # cytoJointR, COMPASS
R_B=/exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript      # R 4.4.3, polyICSFlow

echo "$POP on $(hostname)"
PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/32_analyse-sarscov2.R "$RUN_DIR" cytojointr "$POP" || { echo "cytojointr stage failed"; exit 1; }
R_LIBS=$PWD/renv/lib-r44 "$R_B" scripts/32_analyse-sarscov2.R "$RUN_DIR" polyics "$POP" || echo "polyics stage failed"

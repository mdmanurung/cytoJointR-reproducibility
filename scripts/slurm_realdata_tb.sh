#!/bin/bash
# Purpose: SLURM array wrapper for the TB real-data benchmark; task i runs the three steps of 31_analyse-tb.R
#          (fit in environment A, polyics in environment B, summary in environment A) on data set i.
# Run:     cd <repo root>; sbatch --array=1-3 scripts/slurm_realdata_tb.sh runs/2026-11-03_realdata/tb fit
#          sbatch --array=1-3 --dependency=afterok:<fit job id> scripts/slurm_realdata_tb.sh runs/2026-11-03_realdata/tb post
#          (fit = step fit; post = steps polyics and summary; "all" runs the three in one job). One step set per
#          job keeps each under the 1 h limit of the short partition.
#          Compile Stan once first (see slurm_benchmark.sh); do not edit scripts or the package while tasks run.
# Inputs:  data/processed/{tb_s,tb_u,tb_simplecompass}/ from 22_prepare-tb-compass.R and 23_prepare-tb-simplecompass.R
# Outputs: <run_dir>/<dataset>/ (see 31_analyse-tb.R); logs in <run_dir>/logs
#SBATCH --job-name=cj_tb
#SBATCH --partition=short   # 1 h limit; the all partition queue was 700 jobs deep
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G           # polyICSFlow expands about 4 million cells
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail

RUN_DIR=${1:?usage: sbatch --array=1-3 slurm_realdata_tb.sh <run_dir> <fit|post|all>}
WHAT=${2:?fit, post or all}
DATASETS=(tb_s tb_u tb_simplecompass)
DATASET=${DATASETS[$((${SLURM_ARRAY_TASK_ID:?run through sbatch --array} - 1))]}

cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/${DATASET}_${WHAT}_job${SLURM_ARRAY_JOB_ID}.out" 2>&1

R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51   # cytoJointR, COMPASS
R_B=/exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript      # R 4.4.3, polyICSFlow

echo "$DATASET on $(hostname)"
if [[ $WHAT == fit || $WHAT == all ]]; then
  PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/31_analyse-tb.R "$RUN_DIR" "$DATASET" fit || { echo "fit failed"; exit 1; }
fi
if [[ $WHAT == post || $WHAT == all ]]; then
  R_LIBS=$PWD/renv/lib-r44 "$R_B" scripts/31_analyse-tb.R "$RUN_DIR" "$DATASET" polyics || { echo "polyics failed"; exit 1; }
  PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/31_analyse-tb.R "$RUN_DIR" "$DATASET" summary || { echo "summary failed"; exit 1; }
fi

#!/bin/bash
# Purpose: SLURM wrapper for the vaccine-trial stages of 30_analyse-vaccine.R; one task per dataset and stage.
# Run:     cd <repo root>; sbatch scripts/slurm_vaccine.sh <fit|compass|polyics|stats> <run_dir> <rv144|hvtn078>
#          Compile Stan once first (not inside tasks, so they do not race to build):
#            PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib \
#              Rscript -e 'pkgload::load_all("cytoJointR-wt/redesign-run"); invisible(cytoJointR:::m1_compile())'
#          Override the limits on the command line, e.g. sbatch --partition=all --time=08:00:00 --mem=48G ...
# Inputs:  data/processed/<dataset>/ from 20_prepare-rv144.R, 21_prepare-hvtn078.R
# Outputs: <run_dir>/<dataset>/ (see 30_analyse-vaccine.R); log in <run_dir>/logs
#SBATCH --job-name=cj_vaccine
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail
STAGE=${1:?stage}; RUN_DIR=${2:?run_dir}; DS=${3:?dataset}
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/${STAGE}_${DS}_${SLURM_JOB_ID}.out" 2>&1
R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
R_B=/exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript
if [ "$STAGE" = polyics ]; then
  R_LIBS=$PWD/renv/lib-r44 "$R_B" scripts/30_analyse-vaccine.R "$STAGE" "$RUN_DIR" "$DS"
else
  PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/30_analyse-vaccine.R "$STAGE" "$RUN_DIR" "$DS"
fi

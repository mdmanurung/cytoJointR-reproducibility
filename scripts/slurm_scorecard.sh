#!/bin/bash
# Purpose: SLURM array wrapper for the fit stages of 70_run-scorecard.R; task i runs one stage of one candidate.
#          Tasks: 1 rv144, 2 tb_s, 3 tb_u, 4-10 sim (the seven scenarios of benchmark/sim_scenarios.tsv, in file order).
# Run:     cd <repo root>; compile Stan once first (tasks must not race to build it):
#            PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript -e 'pkgload::load_all("<pkg>"); invisible(cytoJointR:::m1_compile())'
#          sbatch --array=1-10 scripts/slurm_scorecard.sh <run_dir> <label> <pkg_worktree> <fit_args>
#          then, after the array has finished (--dependency=afterok:<job id>):
#          sbatch scripts/slurm_scorecard.sh <run_dir> <label> <pkg_worktree> <fit_args> collect
#          <fit_args> is a cj_fit argument set, e.g. model=m1dm,level=sample. Do not edit scripts or the package while tasks run.
# Inputs:  see 70_run-scorecard.R
# Outputs: <run_dir>/<label>/ (see 70_run-scorecard.R); logs in <run_dir>/<label>/logs
#SBATCH --job-name=cj_scorecard
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail
RUN_DIR=${1:?run_dir}; LABEL=${2:?label}; PKG=${3:?pkg worktree}; FIT=${4:?fit args}; MODE=${5:-fit}
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/$LABEL/logs"
R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
if [ "$MODE" = collect ]; then
  STAGE=collect; EXTRA=()
else
  SCENARIOS=($(awk -F'\t' 'NR > 1 { print $6 }' benchmark/sim_scenarios.tsv))
  TASK=${SLURM_ARRAY_TASK_ID:?run through sbatch --array}
  case $TASK in
    1) STAGE=rv144; EXTRA=() ;;
    2) STAGE=tb_s; EXTRA=() ;;
    3) STAGE=tb_u; EXTRA=() ;;
    *) STAGE=sim; EXTRA=("--scenario=${SCENARIOS[$((TASK - 4))]}") ;;
  esac
fi
exec > "$RUN_DIR/$LABEL/logs/${STAGE}${EXTRA[*]:+_${EXTRA[0]#--scenario=}}_${SLURM_JOB_ID}.out" 2>&1
PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/70_run-scorecard.R "$STAGE" "$RUN_DIR" --pkg="$PKG" --fit="$FIT" --label="$LABEL" "${EXTRA[@]}"

#!/bin/bash
# Purpose: SLURM array wrapper; task i runs every replicate that jobs.tsv assigns to task i (several per task,
#          so the array stays near 100 tasks): 06 (environment A: simulate, candidate models, COMPASS) then 07
#          (environment B: polyICSFlow) for each scenario x seed.
# Run:     cd <repo root>; sbatch --array=1-<n_tasks> scripts/slurm_benchmark.sh <run_dir>
#          n_tasks = the largest task in <run_dir>/jobs.tsv. Compile Stan once first, and do not edit scripts or the
#          package while tasks run:
#            PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib \
#              Rscript -e 'pkgload::load_all("cytoJointR-wt/redesign"); invisible(cytoJointR:::m1_compile())'
# Inputs:  <run_dir>/scenarios.tsv, jobs.tsv from 05_define-benchmark.R. A finished step is skipped, so a failed or
#          timed-out task can be resubmitted: sbatch --array=<task> ...
# Override on the command line (sbatch --partition=all --time=...) for tasks longer than 1 h.
# Outputs: <run_dir>/reps/<scenario_id>/seed_<NN>/ (see 06 and 07); logs in <run_dir>/logs
#SBATCH --job-name=cj_redesign
#SBATCH --partition=short   # 1 h limit; medium caps a user at 4 CPUs and 8G (QOS restrictmedium), so no array
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -uo pipefail

RUN_DIR=${1:?usage: sbatch --array=1-N slurm_benchmark.sh <run_dir>}
TASK=${SLURM_ARRAY_TASK_ID:?run through sbatch --array}
BUDGET=2700   # seconds: start no new replicate after 45 min, so the 1 h limit never kills a fit mid-write

cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/job${SLURM_ARRAY_JOB_ID}_task${TASK}.out" 2>&1

R_A=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51   # cytoJointR, COMPASS
R_B=/exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript      # R 4.4.3, polyICSFlow

awk -F'\t' -v t="$TASK" 'NR > 1 && $1 == t { print $2 "\t" $3 }' "$RUN_DIR/jobs.tsv" | while IFS=$'\t' read -r SCENARIO SEED; do
  REP=$RUN_DIR/reps/$SCENARIO/$(printf 'seed_%02d' "$SEED")
  if (( SECONDS > BUDGET )); then echo "skipped $SCENARIO seed $SEED: time budget"; continue; fi
  echo "task $TASK: $SCENARIO seed $SEED on $(hostname), ${SECONDS}s in"
  if [ ! -f "$REP/diagnostics.tsv" ]; then
    PATH=$R_A/bin:$PATH R_LIBS=$PWD/renv/lib "$R_A/bin/Rscript" scripts/06_run-benchmark-replicate.R "$RUN_DIR" "$SCENARIO" "$SEED" || { echo "06 failed: $SCENARIO seed $SEED"; continue; }
  fi
  if [ ! -f "$REP/polyicsflow/runtime.tsv" ]; then
    R_LIBS=$PWD/renv/lib-r44 "$R_B" scripts/07_run-polyicsflow.R "$RUN_DIR" "$SCENARIO" "$SEED" || echo "07 failed: $SCENARIO seed $SEED"
  fi
done

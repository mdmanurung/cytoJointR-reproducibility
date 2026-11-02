#!/bin/bash
# Purpose: SLURM array wrapper; task i runs 03_run-replicate.R for one scenario x seed.
# Run:     cd <repo root>; sbatch --array=1-<=100 scripts/slurm_pilot.sh <run_dir> <n_seeds> <in_smoke|in_pilot> [offset]
# Inputs:  <run_dir>/scenarios.tsv from 02_define-scenarios.R; the third argument names the
#          TRUE/FALSE column that selects scenarios. Task i + offset -> scenario (i-1)/n_seeds, seed (i-1)%n_seeds+1.
# Override on the command line (sbatch --partition=all --time=...) for tasks longer than 1 h.
# Outputs: <run_dir>/reps/<scenario_id>/seed_<NN>/ (see 03_run-replicate.R); logs in <run_dir>/logs
#SBATCH --job-name=cj_pilot
#SBATCH --partition=short   # medium caps a user at 4 CPUs and 8G (QOS restrictmedium), so no array
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
set -euo pipefail

RUN_DIR=${1:?usage: sbatch --array=1-N slurm_pilot.sh <run_dir> <n_seeds> <in_smoke|in_pilot>}
N_SEEDS=${2:?usage: sbatch --array=1-N slurm_pilot.sh <run_dir> <n_seeds> <in_smoke|in_pilot>}
COLUMN=${3:?usage: sbatch --array=1-N slurm_pilot.sh <run_dir> <n_seeds> <in_smoke|in_pilot>}
OFFSET=${4:-0}  # the cluster caps array indices at 124, so larger sets are submitted in chunks of 100
TASK=$(( ${SLURM_ARRAY_TASK_ID:?run through sbatch --array} + OFFSET ))

cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p "$RUN_DIR/logs"
exec > "$RUN_DIR/logs/job${SLURM_ARRAY_JOB_ID}_task${TASK}.out" 2>&1

R_ENV=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51
export PATH=$R_ENV/bin:$PATH
export R_LIBS=$PWD/renv/lib

# Scenario ids whose COLUMN is TRUE, in file order.
mapfile -t IDS < <(awk -F'\t' -v col="$COLUMN" '
  NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "scenario_id") s = i; if ($i == col) f = i }; next }
  $f == "TRUE" { print $s }' "$RUN_DIR/scenarios.tsv")
SCENARIO=${IDS[$(( (TASK - 1) / N_SEEDS ))]:?task index beyond the scenario list}
SEED=$(( (TASK - 1) % N_SEEDS + 1 ))

echo "task $TASK: $SCENARIO seed $SEED on $(hostname)"
"$R_ENV/bin/Rscript" scripts/03_run-replicate.R "$RUN_DIR" "$SCENARIO" "$SEED"

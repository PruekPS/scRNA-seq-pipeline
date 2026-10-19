#!/bin/bash
#SBATCH --job-name=split-pipe-comb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#
# split-pipe --mode comb: merge per-sublibrary outputs (each run through at
# least the dge step) into one combined dataset. Chemistry/kit/samples/genome
# are read from the sublibrary runs. Usually chained after the array job:
#
#   sbatch --account=<GROUP> --qos=<GROUP>-b --dependency=afterok:<ARRAY_JOBID> \
#          --output=<OUT_ROOT>/logs/split-pipe-comb-%j.out \
#          pipeline/split-pipe/run_split_pipe_comb.slurm <SUBLIBS_TSV> <OUT_ROOT>
#
# Combined output goes to <OUT_ROOT>/combined/.

set -euo pipefail

SUBLIBS_TSV="${1:?usage: run_split_pipe_comb.slurm <SUBLIBS_TSV> <OUT_ROOT>}"
OUT_ROOT="${2:?usage: run_split_pipe_comb.slurm <SUBLIBS_TSV> <OUT_ROOT>}"
ENV_NAME="${SPLITPIPE_ENV:-splitpipe}"

mapfile -t SUBLIB_DIRS < <(cut -f1 "${SUBLIBS_TSV}" | sed "s|^|${OUT_ROOT}/|")
for d in "${SUBLIB_DIRS[@]}"; do
  [[ -d "${d}" ]] || { echo "ERROR: missing sublibrary output ${d}" >&2; exit 1; }
done

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
echo "[$(date)] combining ${#SUBLIB_DIRS[@]} sublibraries"

conda run --no-capture-output -n "${ENV_NAME}" split-pipe \
  --mode comb \
  --sublibraries "${SUBLIB_DIRS[@]}" \
  --output_dir "${OUT_ROOT}/combined" \
  --nthreads "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] combine done"

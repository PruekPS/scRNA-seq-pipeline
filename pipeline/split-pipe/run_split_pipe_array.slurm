#!/bin/bash
#SBATCH --job-name=split-pipe-all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL
#
# split-pipe --mode all, one SLURM array task per sublibrary. Sized for
# Parse's 500M-1B+ reads/sublibrary tier (32 CPUs, 256 GB). Non-GPU.
# Account/QoS/log path/array range are passed at submit time; throttle the
# array (%N) so N x 32 CPUs fits your QoS (e.g. %3 = 96 CPUs):
#
#   sbatch --account=<GROUP> --qos=<GROUP>-b --array=1-<N>%3 \
#          --output=<OUT_ROOT>/logs/split-pipe-%A_%a.out \
#          pipeline/split-pipe/run_split_pipe_array.slurm \
#          <SUBLIBS_TSV> <SAMP_LIST> <GENOME_DIR> <OUT_ROOT> <CHEMISTRY> <KIT>
#
# SUBLIBS_TSV: one sublibrary per line, tab-separated: <name> <fq1> <fq2>
#   (fq1 = mRNA read R1, fq2 = barcode read R2; array task i uses line i).
# SAMP_LIST:   split-pipe --samp_list file, "<sample_name> <wells>" per line,
#   wells in split-pipe notation (see `split-pipe --explain`), e.g. A11-B3.
# Output for each sublibrary goes to <OUT_ROOT>/<name>/; combine them
# afterwards with run_split_pipe_comb.slurm.

set -euo pipefail

usage="usage: run_split_pipe_array.slurm <SUBLIBS_TSV> <SAMP_LIST> <GENOME_DIR> <OUT_ROOT> <CHEMISTRY> <KIT>"
SUBLIBS_TSV="${1:?${usage}}"
SAMP_LIST="${2:?${usage}}"
GENOME_DIR="${3:?${usage}}"
OUT_ROOT="${4:?${usage}}"
CHEMISTRY="${5:?${usage}}"
KIT="${6:?${usage}}"
ENV_NAME="${SPLITPIPE_ENV:-splitpipe}"
TASK="${SLURM_ARRAY_TASK_ID:?submit as an array job (--array=1-N)}"

IFS=$'\t' read -r NAME FQ1 FQ2 < <(sed -n "${TASK}p" "${SUBLIBS_TSV}")
if [[ -z "${NAME:-}" || ! -s "${FQ1}" || ! -s "${FQ2}" ]]; then
  echo "ERROR: line ${TASK} of ${SUBLIBS_TSV} is empty or its FASTQs are missing" >&2
  exit 1
fi

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
mkdir -p "${OUT_ROOT}/${NAME}"
echo "[$(date)] task=${TASK} sublib=${NAME} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK}"

conda run --no-capture-output -n "${ENV_NAME}" split-pipe \
  --mode all \
  --chemistry "${CHEMISTRY}" \
  --kit "${KIT}" \
  --genome_dir "${GENOME_DIR}" \
  --fq1 "${FQ1}" \
  --fq2 "${FQ2}" \
  --samp_list "${SAMP_LIST}" \
  --output_dir "${OUT_ROOT}/${NAME}" \
  --nthreads "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] ${NAME} done"

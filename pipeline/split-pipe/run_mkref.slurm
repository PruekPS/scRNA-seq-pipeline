#!/bin/bash
#SBATCH --job-name=split-pipe-mkref
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb          # STAR genomeGenerate on human needs ~32 GB; 64 GB leaves headroom
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#
# Build a split-pipe genome index (--mode mkref). Non-GPU, CPU + RAM bound.
# Account/QoS/log path are passed at submit time so this file stays generic:
#
#   sbatch --account=<GROUP> --qos=<GROUP> \
#          --output=<OUT_PARENT>/logs/mkref-%j.out \
#          pipeline/split-pipe/run_mkref.slurm <GENOME_NAME> <FASTA> <GTF> <OUT_DIR>
#
# FASTA and GTF must share seqnames (e.g. both "chr1", as in a GENCODE
# primary_assembly pair) and the GTF must not be UCSC-formatted -- see
# README.md "Building a genome reference". Outputs built by split-pipe
# versions older than the installed one are not compatible; rebuild.

set -euo pipefail

GENOME_NAME="${1:?usage: run_mkref.slurm <GENOME_NAME> <FASTA> <GTF> <OUT_DIR>}"
FASTA="${2:?usage: run_mkref.slurm <GENOME_NAME> <FASTA> <GTF> <OUT_DIR>}"
GTF="${3:?usage: run_mkref.slurm <GENOME_NAME> <FASTA> <GTF> <OUT_DIR>}"
OUT_DIR="${4:?usage: run_mkref.slurm <GENOME_NAME> <FASTA> <GTF> <OUT_DIR>}"
ENV_NAME="${SPLITPIPE_ENV:-splitpipe}"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true

# Fail fast on the most common mkref error: FASTA/GTF seqname mismatch.
fasta_chr1=$(grep -m1 '^>' "${FASTA}" | sed 's/^>//; s/[[:space:]].*//')
if ! awk -F'\t' -v c="${fasta_chr1}" '!/^#/ && $1 == c {found = 1; exit} END {exit !found}' "${GTF}"; then
  echo "ERROR: first FASTA seqname '${fasta_chr1}' not found in GTF column 1." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
echo "[$(date)] split-pipe $(conda run -n "${ENV_NAME}" split-pipe --version 2>&1 | tail -1)"

conda run -n "${ENV_NAME}" split-pipe \
  --mode mkref \
  --genome_name "${GENOME_NAME}" \
  --fasta "${FASTA}" \
  --genes "${GTF}" \
  --output_dir "${OUT_DIR}" \
  --nthreads "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] mkref finished"
ls -la "${OUT_DIR}"

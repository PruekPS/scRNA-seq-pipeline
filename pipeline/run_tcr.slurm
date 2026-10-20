#!/bin/bash
#SBATCH --job-name=scrna-tcr
#SBATCH --account=<GROUP>   # set via `sbatch --account=... --qos=...` at submission
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb          # UNVALIDATED on real data; ~150k cells x one TCR table is modest
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#
# TCR stage, run AFTER the WT chain (annotated .h5ad exists) and AFTER
# split-pipe --tcr_analysis (+ --mode comb) has produced the AIRR tables:
#   12 attach TCR to WT cells, chain QC, clonotypes (TRA+TRB CDR3 aa by default)
#   13 diversity (raw + rarefied), sharing across samples, phenotype of expanded clones
#
# Usage (from the repo root):
#   sbatch --account=<g> --qos=<q> --output=<logs>/scrna-tcr-%j.out pipeline/run_tcr.slurm \
#     OUT_DIR WT_H5AD AIRR_ROOT [sample_col] [group_col] [annotation_col]
# Extra flags for step 12 (e.g. "--tcr-suffix-format __{label}") go in the
# TCR_ATTACH_EXTRA environment variable; for step 13 in TCR_ANALYSIS_EXTRA:
#   TCR_ATTACH_EXTRA="--strip-tcr-suffix" sbatch ... pipeline/run_tcr.slurm ...
#
# NOT YET RUN ON REAL DATA -- scripts 12/13 are smoke-tested on
# tests/synthetic_tcr.py only. Read 12's match-rate line and the
# tcr_fraction_by_celltype plot before trusting anything downstream.

set -euo pipefail

OUT_DIR="${1:?Usage: run_tcr.slurm OUT_DIR WT_H5AD AIRR_ROOT [sample_col] [group_col] [annotation_col]}"
WT_H5AD="${2:?need the annotated WT .h5ad}"
AIRR_ROOT="${3:?need the split-pipe TCR output root (searched for */TCR_filtered/*tcr_annotation_airr.tsv)}"
SAMPLE_COL="${4:-sample}"
GROUP_COL="${5:-}"
ANNOTATION_COL="${6:-cell_type_v2}"
ENV_NAME="scrna"

cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p "${OUT_DIR}"
module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

python scripts/12_tcr_attach.py \
  --wt-input "${WT_H5AD}" \
  --airr-root "${AIRR_ROOT}" \
  --output "${OUT_DIR}/12_tcr.h5ad" \
  --plot-dir "${OUT_DIR}/12_tcr_plots" \
  --sample-col "${SAMPLE_COL}" \
  --annotation-col "${ANNOTATION_COL}" \
  ${TCR_ATTACH_EXTRA:-}

GROUP_ARG=()
[ -n "${GROUP_COL}" ] && GROUP_ARG=(--group-col "${GROUP_COL}")
python scripts/13_tcr_analysis.py \
  --input "${OUT_DIR}/12_tcr.h5ad" \
  --output-dir "${OUT_DIR}/13_tcr_results" \
  --sample-col "${SAMPLE_COL}" \
  --annotation-col "${ANNOTATION_COL}" \
  "${GROUP_ARG[@]}" \
  ${TCR_ANALYSIS_EXTRA:-}

echo "TCR stage done: ${OUT_DIR}"

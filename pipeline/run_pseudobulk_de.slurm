#!/bin/bash
#SBATCH --job-name=scrna-pseudobulk-de
#SBATCH --account=<GROUP>   # TODO: set via `sbatch --account=... --qos=...` at submission, not by editing this file
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#
# Optional pseudobulk differential-expression chain, run AFTER an annotated
# .h5ad already exists (run_pipeline.slurm + optionally run_downstream.slurm):
#   10 aggregate cells into pseudobulk pseudo-samples (Python)
#   -> 11 DESeq2 DE per cell type + volcano plots + top-DE-gene heatmap (R)
#
# Dataset-agnostic: --condition-col/--condition-a/--condition-b are supplied
# at submission time, not hardcoded, since this isn't tied to any one
# comparison. NEITHER step has been test-run as of this commit -- there is
# no bundled real condition comparison to validate against by design (see
# scripts/10_pseudobulk_prep.py docstring). Smoke-test against tests/data/
# with a synthetic condition column, or a small real subset, before
# trusting output.
#
# Usage (run sbatch from inside the repo):
#   cd /path/to/scRNA-pipeline
#   sbatch --account=<your-account> --qos=<your-qos> \
#     --output=/path/to/logs/scrna-pseudobulk-de-%j.out \
#     pipeline/run_pseudobulk_de.slurm \
#     /path/to/OUT_DIR CONDITION_COL CONDITION_A CONDITION_B [annotation_col] [sample_col]
#
# OUT_DIR must already contain an annotated .h5ad (04b_annotated.h5ad
# preferred, falls back to 04_annotated.h5ad) with --condition-col actually
# present in its .obs (i.e. it needs to have survived from the raw Parse
# cell_metadata.csv through every step -- same mechanism as the existing
# 'sample' column, see scripts/00_load_parse_output.py).
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna"

OUT_DIR="${1:?Usage: sbatch run_pseudobulk_de.slurm /path/to/OUT_DIR CONDITION_COL CONDITION_A CONDITION_B [annotation_col] [sample_col]}"
CONDITION_COL="${2:?A --condition-col value is required (e.g. disease_status) -- pass as the 2nd positional arg}"
CONDITION_A="${3:?A first condition level is required (e.g. disease) -- pass as the 3rd positional arg}"
CONDITION_B="${4:?A second condition level is required (e.g. healthy) -- pass as the 4th positional arg}"
ANNOTATION_COL="${5:-cell_type_v2}"
SAMPLE_COL="${6:-sample}"

ANNOTATED_INPUT="${OUT_DIR}/04b_annotated.h5ad"
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  ANNOTATED_INPUT="${OUT_DIR}/04_annotated.h5ad"
fi
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  echo "ERROR: no annotated .h5ad found in ${OUT_DIR} (looked for 04b_annotated.h5ad, 04_annotated.h5ad)." >&2
  exit 1
fi

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"

echo "=== [1/2] Pseudobulk aggregation (condition-col=${CONDITION_COL}) ==="
python scripts/10_pseudobulk_prep.py \
  --input "${ANNOTATED_INPUT}" \
  --output-dir "${OUT_DIR}/pseudobulk" \
  --condition-col "${CONDITION_COL}" \
  --sample-col "${SAMPLE_COL}" \
  --annotation-col "${ANNOTATION_COL}"

echo "=== [2/2] DESeq2 (${CONDITION_A} vs ${CONDITION_B}) + volcano + heatmap ==="
Rscript scripts/11_pseudobulk_deseq2.R \
  --input-dir="${OUT_DIR}/pseudobulk" \
  --output-dir="${OUT_DIR}/pseudobulk_de" \
  --condition-a="${CONDITION_A}" \
  --condition-b="${CONDITION_B}"

echo "=== Done. Final outputs: ${OUT_DIR}/pseudobulk_de/ ==="

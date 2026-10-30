#!/bin/bash
#SBATCH --job-name=scrna-pseudobulk-de
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
# at submission time, not hardcoded, since this isn't tied to any one comparison.
#
# Usage (run sbatch from inside the repo):
#   cd /path/to/scRNA-pipeline
#   sbatch --account=<your-account> --qos=<your-qos> \
#     --output=/path/to/logs/scrna-pseudobulk-de-%j.out \
#     pipeline/run_pseudobulk_de.slurm \
#     /path/to/OUT_DIR CONDITION_COL CONDITION_A CONDITION_B [annotation_col] [sample_col]
#
# OUT_DIR must already contain an annotated .h5ad: 04c_celltypist.h5ad if present,
# else 04b_annotated.h5ad, or name another with INPUT_H5AD=/path/file.h5ad (for
# example 05b's final labels). --condition-col must be present in
# its .obs -- normally added at load time from the sample sheet
# (scripts/00_load_parse_output.py --sample-sheet, or SAMPLE_SHEET= for
# pipeline/run_pipeline.slurm).
#
# Design is per experiment: put sample-level covariates in COVARIATES
# (comma-separated) to fit ~ covariates + condition, e.g. a paired
# before/after study:
#   COVARIATES=patient sbatch ... pipeline/run_pseudobulk_de.slurm OUT_DIR timepoint 6h pre
# or, when some samples have no partner (uses only pairs seen in both conditions):
#   PAIR_COL=patient sbatch ... pipeline/run_pseudobulk_de.slurm OUT_DIR timepoint 6h pre
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna"

OUT_DIR="${1:?Usage: sbatch run_pseudobulk_de.slurm /path/to/OUT_DIR CONDITION_COL CONDITION_A CONDITION_B [annotation_col] [sample_col]}"
CONDITION_COL="${2:?A --condition-col value is required (e.g. disease_status) -- pass as the 2nd positional arg}"
CONDITION_A="${3:?A first condition level is required (e.g. disease) -- pass as the 3rd positional arg}"
CONDITION_B="${4:?A second condition level is required (e.g. healthy) -- pass as the 4th positional arg}"
ANNOTATION_COL="${5:-cell_type_v2}"
SAMPLE_COL="${6:-sample}"

ANNOTATED_INPUT="${INPUT_H5AD:-}"
if [ -z "${ANNOTATED_INPUT}" ]; then
  for candidate in 04c_celltypist.h5ad 04b_annotated.h5ad 04_annotated.h5ad; do
    if [ -f "${OUT_DIR}/${candidate}" ]; then ANNOTATED_INPUT="${OUT_DIR}/${candidate}"; break; fi
  done
fi
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  echo "ERROR: no annotated .h5ad found in ${OUT_DIR} (or INPUT_H5AD does not exist)." >&2
  exit 1
fi
PAIR_COL="${PAIR_COL:-}"
PREP_PAIR=(); DE_PAIR=()
if [ -n "${PAIR_COL}" ]; then PREP_PAIR=(--pair-col "${PAIR_COL}"); DE_PAIR=(--pair-col=pair); fi

cd "${REPO_ROOT}"
source pipeline/activate_env.sh "${ENV_NAME}"

echo "=== [1/2] Pseudobulk aggregation of ${ANNOTATED_INPUT} (condition-col=${CONDITION_COL}) ==="
python scripts/10_pseudobulk_prep.py \
  --input "${ANNOTATED_INPUT}" \
  --output-dir "${OUT_DIR}/pseudobulk" \
  --condition-col "${CONDITION_COL}" \
  --sample-col "${SAMPLE_COL}" \
  --annotation-col "${ANNOTATION_COL}" \
  --covariate-cols "${COVARIATES:-}" \
  ${PREP_PAIR[@]+"${PREP_PAIR[@]}"}

echo "=== [2/2] DESeq2 (${CONDITION_A} vs ${CONDITION_B}) + volcano + heatmap ==="
Rscript scripts/11_pseudobulk_deseq2.R \
  --input-dir="${OUT_DIR}/pseudobulk" \
  --output-dir="${OUT_DIR}/pseudobulk_de" \
  --condition-a="${CONDITION_A}" \
  --condition-b="${CONDITION_B}" \
  --covariates="${COVARIATES:-}" \
  ${DE_PAIR[@]+"${DE_PAIR[@]}"}

echo "=== Done. Final outputs: ${OUT_DIR}/pseudobulk_de/ ==="

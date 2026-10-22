#!/bin/bash
#SBATCH --job-name=scrna-pipeline
#SBATCH --account=<GROUP>   # TODO: set via `sbatch --account=... --qos=...` at submission, not by editing this file
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#
# Chains scripts/00_load_parse_output.py -> 01_qc_filter.py -> 02_normalize_pca.py
# -> 03_cluster.py -> 04b_better_annotation.py as one batch job.
#
# Default resources (16 CPU / 128GB) suit a ~100-200k-cell run. For a ~1M-cell
# run override at submit time: a 987,585-cell x 60,664-gene matrix is ~1.18B
# non-zero entries, ~17GB just to parse, and at 128GB it OOM-killed inside
# scripts/01_qc_filter.py (Scrublet densifies working arrays; peak spiked
# well past the last-sampled MaxRSS of ~67.6GB) -- use --mem=400gb there.
# 16 CPU exceeds some group investment QoS caps; submit to the burst QoS
# (<group>-b) if the regular one rejects it with QOSGrpCpuLimit.
#
# Usage (run sbatch from inside the repo -- REPO_ROOT is derived from
# $SLURM_SUBMIT_DIR, not hardcoded; --output also needs a real path
# since #SBATCH lines are parsed by SLURM itself, not this script, so
# a placeholder there can't be filled in from a shell variable below):
#   cd /path/to/scRNA-pipeline
#   sbatch --account=<your-account> --qos=<your-qos> \
#     --output=/path/to/logs/scrna-pipeline-%j.out \
#     pipeline/run_pipeline.slurm \
#     /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier] [sample_sheet.csv]
#
# For split-pipe combined output (.../combined/all-sample/DGE_filtered), the
# per-cell 'sample' column in cell_metadata.csv is kept and sample_name is
# ignored. The sample sheet (study design per sample; docs/03) is the 5th
# argument or SAMPLE_SHEET=...; its sample-name column is 'sample' unless
# SAMPLE_SHEET_KEY=... (older name: SAMPLE_METADATA_KEY). SPECIES (default
# human; set mouse for mouse data -- QC's mitochondrial prefix and 04b's
# markers depend on it) and RESOLUTION (default 0.6, the Leiden resolution
# 04b annotates) can also be set as environment variables at submit time.
#
# If split-pipe's per-sample cell calls are unreliable (check with
# scripts/00a_check_cell_calls.py), give DGE_unfiltered instead and one
# cutoff for every sample via MIN_COUNTS / MIN_GENES (docs/03-your-data.md):
#   SAMPLE_SHEET=samples.csv MIN_COUNTS=500 MIN_GENES=300 sbatch ... \
#     pipeline/run_pipeline.slurm <comb_out>/all-sample/DGE_unfiltered out/ all-sample
#
# Optional steps, switched on by environment variables at submit time:
#   EXCLUDE_SUBLIBS="s1 ..."  drop split-pipe sublibraries in step 0 (library artifact)
#   HARMONY_KEY=sample        run 02b Harmony on that obs column and cluster on
#                             X_pca_harmony (embedding only; counts untouched)
#   CELLTYPIST_MODEL=/path/Immune_All_Low.pkl
#                             run 04c CellTypist after 04b (majority vote over
#                             leiden_r${RESOLUTION}); CELLTYPIST_MODEL2 adds a
#                             second model (e.g. Immune_All_High) as celltypist_high_*
#   e.g. sbatch --export=ALL,HARMONY_KEY=sample,CELLTYPIST_MODEL=... pipeline/run_pipeline.slurm ...
#
# QC (step 1) runs Scrublet per-'sample' batch by default (01_qc_filter.py's
# --batch-key defaults to "sample") -- this is both more correct (Scrublet
# assumes a single library per run) and far lower peak memory than the
# original pooled call, which OOM-killed a 987k-cell run even at 400GB.
# Single-sample data falls back to pooled mode automatically.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna"

RAW_DGE_DIR="${1:?Usage: sbatch run_pipeline.slurm /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier]}"
OUT_DIR="${2:?Usage: sbatch run_pipeline.slurm /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier]}"
SAMPLE_NAME="${3:-parse_run}"
MAD="${4:-5}"
SAMPLE_SHEET="${5:-${SAMPLE_SHEET:-}}"
SAMPLE_SHEET_KEY="${SAMPLE_SHEET_KEY:-${SAMPLE_METADATA_KEY:-sample}}"
MIN_COUNTS="${MIN_COUNTS:-0}"
MIN_GENES="${MIN_GENES:-0}"
SPECIES="${SPECIES:-human}"
RESOLUTION="${RESOLUTION:-0.6}"
EXCLUDE_SUBLIBS="${EXCLUDE_SUBLIBS:-}"
HARMONY_KEY="${HARMONY_KEY:-}"
CELLTYPIST_MODEL="${CELLTYPIST_MODEL:-}"
CELLTYPIST_MODEL2="${CELLTYPIST_MODEL2:-}"

mkdir -p "${OUT_DIR}"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"
echo "repo: ${REPO_ROOT} @ $(cat REVISION 2>/dev/null || git rev-parse --short HEAD 2>/dev/null || echo unknown)"

META_ARGS=()
if [[ -n "${SAMPLE_SHEET}" ]]; then
  META_ARGS=(--sample-sheet "${SAMPLE_SHEET}" --sample-sheet-key "${SAMPLE_SHEET_KEY}")
fi
if [[ -n "${EXCLUDE_SUBLIBS}" ]]; then
  read -r -a _sublibs <<< "${EXCLUDE_SUBLIBS}"
  META_ARGS+=(--exclude-sublibraries "${_sublibs[@]}")
fi

echo "=== [0/4] Load raw Parse output (${RAW_DGE_DIR}) ==="
python scripts/00_load_parse_output.py \
  --sample-dir "${RAW_DGE_DIR}" "${SAMPLE_NAME}" \
  ${META_ARGS[@]+"${META_ARGS[@]}"} \
  --min-counts "${MIN_COUNTS}" --min-genes "${MIN_GENES}" \
  --output "${OUT_DIR}/00_raw.h5ad"

echo "=== [1/4] QC filtering (MAD=${MAD}) ==="
python scripts/01_qc_filter.py \
  --input "${OUT_DIR}/00_raw.h5ad" \
  --output "${OUT_DIR}/01_filtered.h5ad" \
  --plot-dir "${OUT_DIR}/qc_plots" \
  --species "${SPECIES}" \
  --mad "${MAD}"

echo "=== [2/4] Normalize / HVG / PCA ==="
python scripts/02_normalize_pca.py \
  --input "${OUT_DIR}/01_filtered.h5ad" \
  --output "${OUT_DIR}/02_pca.h5ad" \
  --plot-dir "${OUT_DIR}/pca_plots"

CLUSTER_INPUT="${OUT_DIR}/02_pca.h5ad"
USE_REP="X_pca"
if [[ -n "${HARMONY_KEY}" ]]; then
  echo "=== [2b] Harmony integration on '${HARMONY_KEY}' ==="
  python scripts/02b_harmony_integration.py \
    --input "${OUT_DIR}/02_pca.h5ad" \
    --output "${OUT_DIR}/02b_harmony.h5ad" \
    --plot-dir "${OUT_DIR}/harmony_plots" \
    --batch-key "${HARMONY_KEY}"
  CLUSTER_INPUT="${OUT_DIR}/02b_harmony.h5ad"
  USE_REP="X_pca_harmony"
fi

echo "=== [3/4] Clustering (Leiden resolution scan + UMAP on ${USE_REP}) ==="
python scripts/03_cluster.py \
  --input "${CLUSTER_INPUT}" \
  --output "${OUT_DIR}/03_clustered.h5ad" \
  --plot-dir "${OUT_DIR}/cluster_plots" \
  --use-rep "${USE_REP}"

echo "=== [4/4] Confidence-gated annotation (04b, ${SPECIES}, leiden_r${RESOLUTION}) ==="
python scripts/04b_better_annotation.py \
  --input "${OUT_DIR}/03_clustered.h5ad" \
  --output "${OUT_DIR}/04b_annotated.h5ad" \
  --plot-dir "${OUT_DIR}/marker_plots_v2" \
  --resolution "${RESOLUTION}" \
  --species "${SPECIES}" \
  --markers-csv "${OUT_DIR}/top_markers.csv"

FINAL="${OUT_DIR}/04b_annotated.h5ad"
if [[ -n "${CELLTYPIST_MODEL}" ]]; then
  echo "=== [4c] CellTypist ($(basename "${CELLTYPIST_MODEL}"), majority vote over leiden_r${RESOLUTION}) ==="
  python scripts/04c_celltypist_annotation.py \
    --input "${FINAL}" \
    --output "${OUT_DIR}/04c_celltypist.h5ad" \
    --model "${CELLTYPIST_MODEL}" \
    --over-clustering "leiden_r${RESOLUTION}" \
    --plot-dir "${OUT_DIR}/celltypist_plots"
  FINAL="${OUT_DIR}/04c_celltypist.h5ad"
  if [[ -n "${CELLTYPIST_MODEL2}" ]]; then
    echo "=== [4c] CellTypist second model ($(basename "${CELLTYPIST_MODEL2}")) ==="
    python scripts/04c_celltypist_annotation.py \
      --input "${FINAL}" \
      --output "${FINAL}" \
      --model "${CELLTYPIST_MODEL2}" \
      --over-clustering "leiden_r${RESOLUTION}" \
      --prefix celltypist_high \
      --plot-dir "${OUT_DIR}/celltypist_plots"
  fi
fi

echo "=== Done. Final object: ${FINAL} ==="

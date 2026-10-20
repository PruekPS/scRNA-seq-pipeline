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
#     /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier] [sample_metadata.csv]
#
# For split-pipe combined output (.../combined/all-sample/DGE_filtered), the
# per-cell 'sample' column in cell_metadata.csv is kept and sample_name is
# ignored. sample_metadata.csv (optional) is joined on by sample name via its
# 'sample_id' column (override with SAMPLE_METADATA_KEY=...). SPECIES
# (default human) and RESOLUTION (default 0.6, the Leiden resolution 04b
# annotates) can also be set as environment variables at submit time.
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
SAMPLE_METADATA="${5:-}"
SAMPLE_METADATA_KEY="${SAMPLE_METADATA_KEY:-sample_id}"
SPECIES="${SPECIES:-human}"
RESOLUTION="${RESOLUTION:-0.6}"

mkdir -p "${OUT_DIR}"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"
echo "repo: ${REPO_ROOT} @ $(cat REVISION 2>/dev/null || git rev-parse --short HEAD 2>/dev/null || echo unknown)"

META_ARGS=()
if [[ -n "${SAMPLE_METADATA}" ]]; then
  META_ARGS=(--sample-metadata "${SAMPLE_METADATA}" --sample-metadata-key "${SAMPLE_METADATA_KEY}")
fi

echo "=== [0/4] Load raw Parse output (${RAW_DGE_DIR}) ==="
python scripts/00_load_parse_output.py \
  --sample-dir "${RAW_DGE_DIR}" "${SAMPLE_NAME}" \
  ${META_ARGS[@]+"${META_ARGS[@]}"} \
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

echo "=== [3/4] Clustering (Leiden resolution scan + UMAP) ==="
python scripts/03_cluster.py \
  --input "${OUT_DIR}/02_pca.h5ad" \
  --output "${OUT_DIR}/03_clustered.h5ad" \
  --plot-dir "${OUT_DIR}/cluster_plots"

echo "=== [4/4] Confidence-gated annotation (04b, ${SPECIES}, leiden_r${RESOLUTION}) ==="
python scripts/04b_better_annotation.py \
  --input "${OUT_DIR}/03_clustered.h5ad" \
  --output "${OUT_DIR}/04b_annotated.h5ad" \
  --plot-dir "${OUT_DIR}/marker_plots_v2" \
  --resolution "${RESOLUTION}" \
  --species "${SPECIES}" \
  --markers-csv "${OUT_DIR}/top_markers.csv"

echo "=== Done. Final object: ${OUT_DIR}/04b_annotated.h5ad ==="

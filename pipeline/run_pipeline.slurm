#!/bin/bash
#SBATCH --job-name=scrna-pipeline
#SBATCH --account=<GROUP>   # TODO: set via `sbatch --account=... --qos=...` at submission, not by editing this file
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=<PROJECT_ROOT>/logs/scrna-pipeline-%j.out   # TODO: set to your /blue/<group>/<user> path
#SBATCH --mail-type=END,FAIL
#
# Chains scripts/00_load_parse_output.py -> 01_qc_filter.py -> 02_normalize_pca.py
# -> 03_cluster.py -> 04_markers_annotation.py as one batch job.
#
# Resources sized for a real, full-scale Parse run (~1M cells): a
# 987,585-cell x 60,664-gene matrix is ~1.18B non-zero entries, ~17GB just
# to parse -- 128GB/16 CPU/24h is a real allocation for real data, not the
# 64GB/8CPU/12h used for local smoke-testing at 50k cells. Resize down if
# your actual dataset is smaller than this.
#
# Usage:
#   sbatch --account=<your-account> --qos=<your-qos> run_pipeline.slurm \
#     /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier]
set -euo pipefail

PROJECT_ROOT="<PROJECT_ROOT>"   # TODO: set to your /blue/<group>/<user> path
REPO_ROOT="${PROJECT_ROOT}/scRNA-pipeline"   # TODO: adjust to wherever this repo lives on HPC
ENV_NAME="scrna"

RAW_DGE_DIR="${1:?Usage: sbatch run_pipeline.slurm /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier]}"
OUT_DIR="${2:?Usage: sbatch run_pipeline.slurm /path/to/DGE_filtered_dir /path/to/output_dir [sample_name] [mad_multiplier]}"
SAMPLE_NAME="${3:-parse_run}"
MAD="${4:-5}"

mkdir -p "${OUT_DIR}" "${PROJECT_ROOT}/logs"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"

echo "=== [0/4] Load raw Parse output (${RAW_DGE_DIR}) ==="
python scripts/00_load_parse_output.py \
  --sample-dir "${RAW_DGE_DIR}" "${SAMPLE_NAME}" \
  --output "${OUT_DIR}/00_raw.h5ad"

echo "=== [1/4] QC filtering (MAD=${MAD}) ==="
python scripts/01_qc_filter.py \
  --input "${OUT_DIR}/00_raw.h5ad" \
  --output "${OUT_DIR}/01_filtered.h5ad" \
  --plot-dir "${OUT_DIR}/qc_plots" \
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

echo "=== [4/4] Marker detection (demo annotation is PBMC-specific -- see script docstring) ==="
python scripts/04_markers_annotation.py \
  --input "${OUT_DIR}/03_clustered.h5ad" \
  --output "${OUT_DIR}/04_annotated.h5ad" \
  --plot-dir "${OUT_DIR}/marker_plots" \
  --markers-csv "${OUT_DIR}/top_markers.csv"

echo "=== Done. Final object: ${OUT_DIR}/04_annotated.h5ad ==="

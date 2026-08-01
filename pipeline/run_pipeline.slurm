#!/bin/bash
#SBATCH --job-name=scrna-pipeline
#SBATCH --account=<GROUP>   # TODO: set to your HiPerGator group; confirm with `sacctmgr show assoc user=$USER format=account,qos`
#SBATCH --qos=<GROUP>       # TODO: set to your HiPerGator QoS -- often <group> or <group>-b (burst)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=<PROJECT_ROOT>/logs/scrna-pipeline-%j.out   # TODO: set to your /blue/<group>/<user> path
#SBATCH --mail-type=END,FAIL
#
# Chains scripts/01_qc_filter.py -> 02_normalize_pca.py -> 03_cluster.py ->
# 04_markers_annotation.py as one batch job. Sized modestly (8 CPU, 64GB,
# no GPU -- nothing in this default path needs one); resize --mem upward
# if the input dataset is large enough to need it, rather than requesting
# more by default. Re-run any single step manually (see each script's
# --help) if you just want to iterate on one stage without redoing earlier
# ones -- this driver is for a full, unattended run only.
#
# Usage:
#   sbatch run_pipeline.slurm /path/to/raw.h5ad /path/to/output_dir [mad_multiplier]
set -euo pipefail

PROJECT_ROOT="<PROJECT_ROOT>"   # TODO: set to your /blue/<group>/<user> path
REPO_ROOT="${PROJECT_ROOT}/scRNA-pipeline"   # TODO: adjust to wherever this repo lives on HPC
ENV_NAME="scrna"

RAW_H5AD="${1:?Usage: sbatch run_pipeline.slurm /path/to/raw.h5ad /path/to/output_dir [mad_multiplier]}"
OUT_DIR="${2:?Usage: sbatch run_pipeline.slurm /path/to/raw.h5ad /path/to/output_dir [mad_multiplier]}"
MAD="${3:-5}"

mkdir -p "${OUT_DIR}" "${PROJECT_ROOT}/logs"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"

echo "=== [1/4] QC filtering (MAD=${MAD}) ==="
python scripts/01_qc_filter.py \
  --input "${RAW_H5AD}" \
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

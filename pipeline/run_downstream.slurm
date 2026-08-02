#!/bin/bash
#SBATCH --job-name=scrna-downstream
#SBATCH --account=<GROUP>   # TODO: set via `sbatch --account=... --qos=...` at submission, not by editing this file
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#
# Optional downstream chain, run AFTER pipeline/run_pipeline.slurm has
# already produced 02_pca.h5ad in an OUT_DIR:
#   02b harmony (optional, skip for single-sample data)
#   -> 03 cluster (reclustered on the Harmony-corrected embedding if used)
#   -> 04b better annotation (finer marker set + confidence gating)
#   -> 06 export for CellChat (flat MTX/CSV, no live R<->Python bridge)
#   -> 07 CellChat analysis (R)
#
# Validated end-to-end on the real full-scale run (2026-08-02): 907,696
# cells post-QC -> Harmony -> 694,803 confidently-annotated cells (after
# excluding Ambiguous clusters) -> CellChat found 291 significant L-R
# interactions across 19 pathways. 400gb/16cpu is proven for steps 1-3;
# CellChat itself (step 4) only needed ~300gb and completed in ~1.5h with
# --workers=1 -- see the comment above the Rscript call below for why it's
# forced to 1 rather than following cpus-per-task.
#
# Resources sized for a real, full-scale Parse run (~1M cells) -- resize
# down for smaller/local-scale data.
#
# Usage (run sbatch from inside the repo -- REPO_ROOT is derived from
# $SLURM_SUBMIT_DIR, not hardcoded; --output also needs a real path
# since #SBATCH lines are parsed by SLURM itself, not this script):
#   cd /path/to/scRNA-pipeline
#   sbatch --account=<your-account> --qos=<your-qos> \
#     --output=/path/to/logs/scrna-downstream-%j.out \
#     pipeline/run_downstream.slurm \
#     /path/to/OUT_DIR [batch_key] [resolution] [skip_harmony: 0|1] [species: human|mouse|canine]
#
# OUT_DIR must already contain 02_pca.h5ad from run_pipeline.slurm. All
# downstream outputs are written into the same OUT_DIR.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna"

OUT_DIR="${1:?Usage: sbatch run_downstream.slurm /path/to/OUT_DIR [batch_key] [resolution] [skip_harmony: 0|1] [species]}"
BATCH_KEY="${2:-sample}"
RESOLUTION="${3:-0.6}"
SKIP_HARMONY="${4:-0}"
SPECIES="${5:-human}"

PCA_H5AD="${OUT_DIR}/02_pca.h5ad"
if [ ! -f "${PCA_H5AD}" ]; then
  echo "ERROR: ${PCA_H5AD} not found -- run pipeline/run_pipeline.slurm first." >&2
  exit 1
fi

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"

CLUSTER_INPUT="${PCA_H5AD}"
USE_REP="X_pca"

if [ "${SKIP_HARMONY}" != "1" ]; then
  echo "=== [1/4] Harmony integration (batch key: ${BATCH_KEY}) ==="
  python scripts/02b_harmony_integration.py \
    --input "${PCA_H5AD}" \
    --output "${OUT_DIR}/02b_harmony.h5ad" \
    --plot-dir "${OUT_DIR}/harmony_plots" \
    --batch-key "${BATCH_KEY}"
  CLUSTER_INPUT="${OUT_DIR}/02b_harmony.h5ad"
  USE_REP="X_pca_harmony"
else
  echo "=== [1/4] Skipping Harmony integration (single-sample data) ==="
fi

echo "=== [2/4] Clustering on ${USE_REP} ==="
python scripts/03_cluster.py \
  --input "${CLUSTER_INPUT}" \
  --output "${OUT_DIR}/03b_clustered.h5ad" \
  --plot-dir "${OUT_DIR}/cluster_plots_v2" \
  --use-rep "${USE_REP}"

echo "=== [3/4] Better annotation (resolution=${RESOLUTION}, species=${SPECIES}) ==="
python scripts/04b_better_annotation.py \
  --input "${OUT_DIR}/03b_clustered.h5ad" \
  --output "${OUT_DIR}/04b_annotated.h5ad" \
  --plot-dir "${OUT_DIR}/marker_plots_v2" \
  --resolution "${RESOLUTION}" \
  --species "${SPECIES}" \
  --markers-csv "${OUT_DIR}/top_markers_v2.csv"

echo "=== [4/4] Export for CellChat + run CellChat ==="
python scripts/06_export_for_cellchat.py \
  --input "${OUT_DIR}/04b_annotated.h5ad" \
  --output-dir "${OUT_DIR}/cellchat_input" \
  --group-col cell_type_v2 \
  --exclude-ambiguous

# --workers=1 deliberately, not $SLURM_CPUS_PER_TASK: future::multisession
# spawns separate R processes that each get their own copy of CellChat's
# large intermediate data (a ~7.2GiB dense matrix here) rather than sharing
# it -- 8 workers OOM-killed a 128gb job at this cell count where 1 worker
# completed fine in ~300gb. Confirmed on the real run, not a guess.
Rscript scripts/07_cellchat_analysis.R \
  --input-dir="${OUT_DIR}/cellchat_input" \
  --output-dir="${OUT_DIR}/cellchat_results" \
  --species="${SPECIES}" \
  --workers=1

echo "=== Done. Final objects: ${OUT_DIR}/04b_annotated.h5ad, ${OUT_DIR}/cellchat_results/ ==="

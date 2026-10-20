#!/bin/bash
#SBATCH --job-name=scrna-trajectory
#SBATCH --account=<GROUP>   # TODO: set via `sbatch --account=... --qos=...` at submission, not by editing this file
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#
# Optional trajectory-analysis chain, run AFTER pipeline/run_downstream.slurm
# (or run_pipeline.slurm alone) has already produced an annotated,
# clustered .h5ad:
#   08 PAGA + diffusion pseudotime (cheap, no new data)
#   -> 09 CellRank fate mapping (CytoTRACE + connectivity [+ DPT] kernels,
#      GPCCA estimator) -- NO RNA velocity, see
#      docs/trajectory_cellrank_overview.html for why and what that gives up.
#
# NEITHER step has been test-run as of this commit -- no real bone-marrow+TME
# dataset exists yet to validate against. See
# docs/trajectory_cellrank_overview.html for the experimental-design
# prerequisites (paired compartment sampling, replicate counts) this
# analysis actually needs before its output means anything. Resource sizing
# below is a starting guess, not proven -- CellRank's Schur decomposition
# and macrostate computation have a different profile (dense eigenvalue
# problems on a graph Laplacian) than anything else in this pipeline.
#
# Usage (run sbatch from inside the repo):
#   cd /path/to/scRNA-pipeline
#   sbatch --account=<your-account> --qos=<your-qos> \
#     --output=/path/to/logs/scrna-trajectory-%j.out \
#     pipeline/run_trajectory.slurm \
#     /path/to/OUT_DIR [annotation_col] [root_cluster] [n_states]
#
# OUT_DIR must already contain an annotated .h5ad from run_downstream.slurm
# (04b_annotated.h5ad) or run_pipeline.slurm (04_annotated.h5ad) -- this
# script prefers 04b_annotated.h5ad if present.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna-traj"   # separate python-3.12 env, see envs/environment_trajectory.yml

OUT_DIR="${1:?Usage: sbatch run_trajectory.slurm /path/to/OUT_DIR [annotation_col] [root_cluster] [n_states]}"
ANNOTATION_COL="${2:-cell_type_v2}"
ROOT_CLUSTER="${3:?A --root-cluster value is required (e.g. HSPC) -- pass it as the 3rd positional arg}"
N_STATES="${4:-8}"

ANNOTATED_INPUT="${OUT_DIR}/04b_annotated.h5ad"
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  ANNOTATED_INPUT="${OUT_DIR}/04_annotated.h5ad"
fi
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  echo "ERROR: no annotated .h5ad found in ${OUT_DIR} (looked for 04b_annotated.h5ad, 04_annotated.h5ad)." >&2
  echo "Run pipeline/run_pipeline.slurm and/or pipeline/run_downstream.slurm first." >&2
  exit 1
fi

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
conda activate "${ENV_NAME}"

cd "${REPO_ROOT}"

echo "=== [1/2] PAGA + diffusion pseudotime (root=${ROOT_CLUSTER}) ==="
python scripts/08_trajectory_paga.py \
  --input "${ANNOTATED_INPUT}" \
  --output "${OUT_DIR}/08_paga.h5ad" \
  --plot-dir "${OUT_DIR}/paga_plots" \
  --annotation-col "${ANNOTATION_COL}" \
  --root-cluster "${ROOT_CLUSTER}"

echo "=== [2/2] CellRank fate mapping (n_states=${N_STATES}, CytoTRACE+connectivity+DPT) ==="
python scripts/09_cellrank_fate_mapping.py \
  --input "${OUT_DIR}/08_paga.h5ad" \
  --output "${OUT_DIR}/09_cellrank.h5ad" \
  --plot-dir "${OUT_DIR}/cellrank_plots" \
  --drivers-csv "${OUT_DIR}/lineage_drivers.csv" \
  --annotation-col "${ANNOTATION_COL}" \
  --n-states "${N_STATES}" \
  --use-dpt-pseudotime

echo "=== Done. Final objects: ${OUT_DIR}/09_cellrank.h5ad, ${OUT_DIR}/lineage_drivers.csv ==="

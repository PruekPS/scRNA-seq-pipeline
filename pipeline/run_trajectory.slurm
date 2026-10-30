#!/bin/bash
#SBATCH --job-name=scrna-trajectory
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
#      docs/reference/trajectory-design.md for why and what that gives up.
#
# Tested on synthetic branching data only (tests/synthetic_trajectory.py). See
# docs/reference/trajectory-design.md for the experimental-design
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
#     /path/to/OUT_DIR ANNOTATION_COL ROOT_CLUSTER [n_states]
#
# ANNOTATION_COL is usually cell_type_v2; ROOT_CLUSTER is a value of it (e.g. HSPC).
# OUT_DIR must already contain an annotated .h5ad: 04c_celltypist.h5ad if present,
# else 04b_annotated.h5ad, or name another with INPUT_H5AD=/path/file.h5ad.
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
ENV_NAME="scrna-traj"   # separate python-3.12 env, see envs/environment_trajectory.yml

OUT_DIR="${1:?Usage: sbatch run_trajectory.slurm /path/to/OUT_DIR ANNOTATION_COL ROOT_CLUSTER [n_states]}"
ANNOTATION_COL="${2:?An annotation column is required (usually cell_type_v2) -- pass it as the 2nd argument}"
ROOT_CLUSTER="${3:?A root cluster is required (e.g. HSPC) -- pass it as the 3rd argument}"
N_STATES="${4:-8}"

ANNOTATED_INPUT="${INPUT_H5AD:-}"
if [ -z "${ANNOTATED_INPUT}" ]; then
  for candidate in 04c_celltypist.h5ad 04b_annotated.h5ad 04_annotated.h5ad; do
    if [ -f "${OUT_DIR}/${candidate}" ]; then ANNOTATED_INPUT="${OUT_DIR}/${candidate}"; break; fi
  done
fi
if [ ! -f "${ANNOTATED_INPUT}" ]; then
  echo "ERROR: no annotated .h5ad found in ${OUT_DIR} (or INPUT_H5AD does not exist)." >&2
  echo "Run pipeline/run_pipeline.slurm and/or pipeline/run_downstream.slurm first." >&2
  exit 1
fi

cd "${REPO_ROOT}"
[[ -f pipeline/activate_env.sh ]] || { echo "ERROR: submit this job from the pipeline folder: cd /path/to/scRNA-seq-pipeline" >&2; exit 1; }
source pipeline/activate_env.sh "${ENV_NAME}"

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

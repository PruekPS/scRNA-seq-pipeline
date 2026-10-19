#!/bin/bash
# Install Parse Biosciences split-pipe into a dedicated conda env, entirely
# under your HiPerGator group/user project root. Run this ON HiPerGator
# (ssh into your cluster alias) via `sbatch` to a compute node -- never on
# the login node, even though this is "just" an install: it compiles Cython
# extensions and runs a conda solve, both real compute.
#
# Follows Parse's own v1.9.0 install procedure (README.md bundled in the
# pipeline zip): create an env with python>=3.12.8, activate it, run their
# install_dependencies_conda.sh into that active env (there is no
# env_conda.yml in this version -- that was true of an older split-pipe
# release, not this one), then `pip install .` (build isolation pulls in
# Cython/numpy itself via pyproject.toml, so no --no-build-isolation flag).
#
# Usage: SPLITPIPE_PROJECT_ROOT=/blue/<group>/<user> bash install_hipergator.sh /path/to/ParseBiosciences-Pipeline.<version>.zip [immune]
#   Pass "immune" as the second argument for Evercode TCR/BCR data.
set -euo pipefail

PROJECT_ROOT="${SPLITPIPE_PROJECT_ROOT:?Set SPLITPIPE_PROJECT_ROOT to your /blue/<group>/<user> path}"
SOFTWARE_DIR="${PROJECT_ROOT}/software/split-pipe"
ENV_NAME="splitpipe"
PY_VERSION="3.12.8"   # exact version Parse's own install docs test against (setup.py requires >=3.12.8)

PIPELINE_ZIP="${1:?Usage: $0 /path/to/ParseBiosciences-Pipeline.<version>.zip [immune]}"
WITH_IMMUNE="${2:-}"

if [[ "${PIPELINE_ZIP}" != "${PROJECT_ROOT}"* ]]; then
  echo "Refusing to run: pipeline zip must live under ${PROJECT_ROOT}" >&2
  exit 1
fi

module load conda 2>/dev/null || module load mamba 2>/dev/null || true

mkdir -p "${SOFTWARE_DIR}"
cd "${SOFTWARE_DIR}"
unzip -o "${PIPELINE_ZIP}"

PIPELINE_DIR=$(find "${SOFTWARE_DIR}" -maxdepth 1 -type d -name "ParseBiosciences-Pipeline*" | sort | tail -n1)
echo "Installing from: ${PIPELINE_DIR}"

conda create -y -n "${ENV_NAME}" "python=${PY_VERSION}"

# install_dependencies_conda.sh installs into whatever conda env is
# *currently active* (it has no env-creation logic of its own -- it calls
# `conda install` directly), so the env must exist and be activated first.
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

cd "${PIPELINE_DIR}"
# This is Parse's "Option 1" (conda-scoped, no sudo) install route -- the
# right one for a HiPerGator account. Never also run install_dependencies.sh
# (the sudo/system-wide "Option 2") in the same env: Parse's own docs warn
# this can cause non-obvious version conflicts.
# -i = install, -y = actually execute (default is a dry run). If this fails
# with a "CondaToSNonInteractiveError", accept the conda ToS once first:
#   conda tos accept --override-channels --channel defaults
#   conda tos accept --override-channels --channel conda-forge
bash ./install_dependencies_conda.sh -i -y

# BCR/TCR (Evercode TCR/BCR kits): the IgBlast databases land in
# ./splitpipe/immune_dir/, which setup.py packages as package_data -- so per
# Parse's README they must be installed BEFORE `pip install .`, not after.
if [[ "${WITH_IMMUNE}" == "immune" ]]; then
  bash ./install_immune_dbs.sh -i -y
  bash ./install_immune_dbs.sh -c -y
fi

pip install . --no-cache-dir

echo "Done. Activate with: conda activate ${ENV_NAME}"
echo "Verify with: conda run -n ${ENV_NAME} split-pipe --version"

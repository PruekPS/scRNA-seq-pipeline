#!/bin/bash
# Install Parse Biosciences split-pipe into a dedicated mamba env,
# entirely under /blue/esayour/puttawatsuphapru/. Run this ON HiPerGator
# (ssh hipergator), not on your laptop.
#
# Usage: bash install_hipergator.sh /path/to/ParseBiosciences-Pipeline.<version>.zip
set -euo pipefail

PROJECT_ROOT="/blue/esayour/puttawatsuphapru"
SOFTWARE_DIR="${PROJECT_ROOT}/software/split-pipe"
ENV_NAME="splitpipe"

PIPELINE_ZIP="${1:?Usage: $0 /path/to/ParseBiosciences-Pipeline.<version>.zip}"

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

conda env create -n "${ENV_NAME}" -f "${PIPELINE_DIR}/env_conda.yml" \
  || conda env update -n "${ENV_NAME}" -f "${PIPELINE_DIR}/env_conda.yml"

conda run -n "${ENV_NAME}" pip install --no-cache-dir --no-build-isolation "${PIPELINE_DIR}"

echo "Done. Activate with: conda activate ${ENV_NAME}"
echo "Verify with: conda run -n ${ENV_NAME} split-pipe --version"

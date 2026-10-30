#!/bin/bash
# Sourced by the job scripts: `source pipeline/activate_env.sh ENV_NAME`. Not run on its own.
#
# Checks the job was submitted from the top of this repository, then makes `conda activate`
# work inside a batch job, where a plain `conda activate` often fails because a job's shell
# has not run conda's setup.
if [[ ! -f scripts/01_qc_filter.py ]]; then
  echo "ERROR: this job was submitted from $(pwd), which is not the pipeline folder." >&2
  echo "Run: cd /path/to/scRNA-seq-pipeline   and then submit the job again." >&2
  exit 1
fi
module load conda 2>/dev/null || module load mamba 2>/dev/null || true
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Ask your cluster's support how to load it (docs/02-installation.md)." >&2
  exit 1
fi
eval "$(conda shell.bash hook)"
# Some packages' activation scripts read unset variables, so `set -u` is off while activating.
set +u
if ! conda activate "$1"; then
  echo "ERROR: could not activate the conda environment '$1'. See docs/02-installation.md." >&2
  echo "Environments found:" >&2
  conda env list >&2
  exit 1
fi
set -u
echo "repo: $(pwd) @ $(cat REVISION 2>/dev/null || git rev-parse --short HEAD 2>/dev/null || echo unknown); env: $1 on $(hostname)"

#!/bin/bash
# Post-`mamba env create -f envs/environment.yml` setup for CellChat.
#
# CellChat and two of its dependencies aren't installable purely from
# environment.yml -- they need this separate step, run ONCE after the conda
# env exists. Run as a SLURM job on a compute node, never the login node
# (login-node compute is against policy on most clusters).
#
# Every command below reflects an actual failure hit and fixed on a real
# HiPerGator run (2026-08-01), not a hypothetical:
#   - ComplexHeatmap/BiocNeighbors are Bioconductor-only -- not resolvable
#     by remotes::install_github's CRAN-based dependency resolver.
#   - svglite failed compiling from source on HiPerGator's toolchain
#     (linker error) -- the conda-forge precompiled build works fine.
#   - NMF needs >=0.23.0 for CellChat; conda's r-nmf caps at 0.21.0 on
#     every channel/R-version combination (confirmed via `conda search`) --
#     has to come from CRAN.
#   - presto (CellChat's fast-Wilcoxon backend) isn't on CRAN or conda at
#     all -- GitHub only. Without it, identifyOverExpressedGenes() errors
#     outright rather than falling back to a slow implementation.
#
# Usage (from the top of this repository; the job's output goes to r-extras-<jobid>.out):
#   sbatch --account=<GROUP> --qos=<QOS> --cpus-per-task=4 --mem=16gb --time=01:30:00 \
#     --output=r-extras-%j.out --wrap="bash envs/install_r_extras.sh"
set -euo pipefail

ENV_NAME="${1:-scrna}"
module load conda 2>/dev/null || module load mamba 2>/dev/null || true

echo "=== conda-installable deps (precompiled, avoids source-build issues) ==="
conda install -n "${ENV_NAME}" -c bioconda -c conda-forge \
  r-svglite bioconductor-complexheatmap bioconductor-biocneighbors r-circlize -y

echo "=== NMF from CRAN (conda's r-nmf caps at 0.21.0, CellChat needs >=0.23.0) ==="
conda run -n "${ENV_NAME}" Rscript -e \
  'remotes::install_version("NMF", version="0.28", repos="https://cloud.r-project.org", upgrade="never")'

echo "=== presto (CellChat's fast-Wilcoxon backend, GitHub-only) ==="
conda run -n "${ENV_NAME}" Rscript -e \
  'remotes::install_github("immunogenomics/presto", upgrade="never")'

echo "=== CellChat itself ==="
conda run -n "${ENV_NAME}" Rscript -e \
  'remotes::install_github("jinworks/CellChat", upgrade="never", dependencies=TRUE)'

echo "=== verifying ==="
conda run -n "${ENV_NAME}" Rscript -e \
  'library(CellChat); library(presto); cat("CELLCHAT_AND_PRESTO_OK", as.character(packageVersion("CellChat")), "\n")'

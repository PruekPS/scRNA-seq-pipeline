#!/bin/bash
#SBATCH --job-name=scrna-env-extras
#SBATCH --account=<GROUP>   # set at submission: sbatch --account=... --qos=... envs/install_trajectory_de.slurm <ENV_PREFIX>
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#
# Part 1: adds bioconductor-deseq2 (script 11) to the EXISTING scrna env.
# Part 2: creates a SEPARATE python-3.12 env `scrna-traj` for cellrank +
#         pygpcca (scripts 08/09), because cellrank>=2.3 needs python>=3.12
#         and the validated scrna env is python 3.11. (A combined install
#         into scrna was tried first and failed to solve; nothing changed.)
# Runs as a compute-node job (never run installs on the login node).
#
# Part 1 is safe by construction:
#   1. Records the conda revision first, so `conda install --revision N`
#      can roll the env back if anything goes wrong.
#   2. Solves with --dry-run --json and ABORTS, installing nothing, if the
#      plan would replace/remove any package the validated pipeline depends
#      on (python, r-base, harmonypy, scanpy, ...). harmonypy must stay
#      0.0.10 -- 2.0.0 corrupts Harmony output (see envs/environment.yml).
#   3. After installing, runs a tiny real DESeq2 test and library(CellChat)
#      to prove the R stack was not disturbed.
#
# Usage (from the repo root, so SLURM_SUBMIT_DIR is the repo):
#   sbatch --account=<g> --qos=<q> envs/install_trajectory_de.slurm /path/to/envs/scrna

set -euo pipefail

ENV_PREFIX="${1:?Usage: sbatch ... envs/install_trajectory_de.slurm /path/to/env/prefix}"
PROTECTED="python r-base harmonypy scanpy anndata numpy pandas scipy leidenalg r-seurat r-matrix scikit-learn umap-learn"
TRAJ_PREFIX="$(dirname "${ENV_PREFIX}")/scrna-traj"
REPO_DIR="${SLURM_SUBMIT_DIR:?run sbatch from the repo root}"

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
echo "=== node: $(hostname)  env: ${ENV_PREFIX} ==="

echo "=== baseline revision (roll back with: conda install -p ${ENV_PREFIX} --revision <N>) ==="
conda list -p "${ENV_PREFIX}" --revisions | tail -3

echo "=== PART 1: dry-run solve for DESeq2 in scrna ==="
PLAN=$(mktemp)
if ! conda install -p "${ENV_PREFIX}" -c conda-forge -c bioconda \
  "bioconductor-deseq2>=1.42" \
  --dry-run --json > "${PLAN}" 2> "${PLAN}.err"; then
  echo "SOLVE FAILED -- nothing was installed. Solver output follows:"
  head -c 6000 "${PLAN}"; echo; head -c 6000 "${PLAN}.err"
  exit 4
fi

python - "${PLAN}" "${PROTECTED}" <<'PYEOF'
import json, sys
plan = json.load(open(sys.argv[1]))
protected = set(sys.argv[2].split())
acts = plan.get("actions", {})
unlink = {p["name"] if isinstance(p, dict) else str(p).split("::")[-1].rsplit("-", 2)[0]
          for p in acts.get("UNLINK", [])}
link = acts.get("LINK", [])
print(f"plan: {len(link)} packages to link, {len(unlink)} to unlink/replace")
hit = sorted(unlink & protected)
if hit:
    print("ABORT: plan would replace/remove protected packages:", ", ".join(hit))
    print("Nothing was installed. Unlinked set:", ", ".join(sorted(unlink)))
    sys.exit(3)
print("protected packages untouched; unlinked set:", ", ".join(sorted(unlink)) or "(none)")
PYEOF

echo "=== PART 1: installing DESeq2 ==="
conda install -y -p "${ENV_PREFIX}" -c conda-forge -c bioconda "bioconductor-deseq2>=1.42"

echo "=== PART 1 verify: harmonypy pin, DESeq2 tiny run, CellChat still loads ==="
"${ENV_PREFIX}/bin/python" - <<'PYEOF'
import importlib.metadata as md
for pkg in ("scanpy", "harmonypy", "numpy", "pandas", "anndata"):
    print(f"{pkg:10s} {md.version(pkg)}")
assert md.version("harmonypy").startswith("0.0."), "harmonypy pin broken"
PYEOF
"${ENV_PREFIX}/bin/Rscript" -e '
suppressMessages({library(DESeq2); library(CellChat); library(ComplexHeatmap)})
cat("DESeq2", as.character(packageVersion("DESeq2")), " CellChat", as.character(packageVersion("CellChat")), "\n")
set.seed(1)
dds <- DESeq(makeExampleDESeqDataSet(n = 300, m = 8), quiet = TRUE)
stopifnot(nrow(results(dds)) == 300)
cat("DESeq2 tiny run OK\n")
'

echo "=== PART 2: separate python-3.12 env for cellrank at ${TRAJ_PREFIX} ==="
if [ -d "${TRAJ_PREFIX}" ]; then
  echo "exists already, updating"
  conda env update -p "${TRAJ_PREFIX}" -f "${REPO_DIR}/envs/environment_trajectory.yml"
else
  conda env create -y -p "${TRAJ_PREFIX}" -f "${REPO_DIR}/envs/environment_trajectory.yml"
fi

echo "=== PART 2 verify: imports + the cellrank API script 09 uses ==="
"${TRAJ_PREFIX}/bin/python" - <<'PYEOF'
import importlib.metadata as md
import cellrank as cr, pygpcca, scanpy
for pkg in ("cellrank", "pygpcca", "scanpy", "numpy", "pandas", "anndata"):
    print(f"{pkg:10s} {md.version(pkg)}")
assert hasattr(cr.kernels, "CytoTRACEKernel") and hasattr(cr.kernels, "ConnectivityKernel")
assert hasattr(cr.estimators, "GPCCA")
print("cellrank API surface used by script 09: OK")
PYEOF
echo "=== DONE: all checks passed ==="

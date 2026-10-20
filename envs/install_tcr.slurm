#!/bin/bash
#SBATCH --job-name=scrna-env-tcr
#SBATCH --account=<GROUP>   # set at submission: sbatch --account=... --qos=...
#SBATCH --qos=<GROUP>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#
# Adds scirpy (scripts 12/13, TCR) to the EXISTING scrna env on a compute node
# (never install on the login node), with the same safety rules as
# envs/install_trajectory_de.slurm:
#   1. record the conda revision (roll back: conda install -p <env> --revision N)
#   2. dry-run solve; ABORT, installing nothing, if any protected package would
#      be replaced (python, r-base, harmonypy, scanpy, anndata, numpy, ...)
#   3. after installing, read a tiny AIRR table through scirpy and define a
#      clonotype, and confirm the R side (CellChat, DESeq2) still loads.
# setuptools<81 is required: scirpy's `airr` dependency imports pkg_resources.
#
# Usage: sbatch --account=<g> --qos=<q> envs/install_tcr.slurm /path/to/envs/scrna

set -euo pipefail

ENV_PREFIX="${1:?Usage: sbatch ... envs/install_tcr.slurm /path/to/env/prefix}"
PROTECTED="python r-base harmonypy scanpy anndata numpy pandas scipy leidenalg r-seurat r-matrix scikit-learn umap-learn"
PKGS=("scirpy=0.22" "setuptools<81")

module load conda 2>/dev/null || module load mamba 2>/dev/null || true
echo "=== node: $(hostname)  env: ${ENV_PREFIX} ==="
conda list -p "${ENV_PREFIX}" --revisions | tail -3

echo "=== dry-run solve ==="
PLAN=$(mktemp)
if ! conda install -p "${ENV_PREFIX}" -c conda-forge -c bioconda "${PKGS[@]}" \
     --dry-run --json > "${PLAN}" 2> "${PLAN}.err"; then
  echo "SOLVE FAILED -- nothing was installed:"; head -c 6000 "${PLAN}"; echo; head -c 6000 "${PLAN}.err"
  exit 4
fi
python - "${PLAN}" "${PROTECTED}" <<'PYEOF'
import json, sys
plan = json.load(open(sys.argv[1]))
protected = set(sys.argv[2].split())
acts = plan.get("actions", {})
unlink = {p["name"] if isinstance(p, dict) else str(p).split("::")[-1].rsplit("-", 2)[0]
          for p in acts.get("UNLINK", [])}
print(f"plan: {len(acts.get('LINK', []))} to link, {len(unlink)} to unlink/replace: {', '.join(sorted(unlink)) or '(none)'}")
hit = sorted(unlink & protected)
if hit:
    print("ABORT: plan would replace protected packages:", ", ".join(hit), "-- nothing installed")
    sys.exit(3)
PYEOF

echo "=== installing ==="
conda install -y -p "${ENV_PREFIX}" -c conda-forge -c bioconda "${PKGS[@]}"

echo "=== verify ==="
"${ENV_PREFIX}/bin/python" -W ignore - <<'PYEOF'
import importlib.metadata as md
import pandas as pd, scirpy as ir
for pkg in ("scirpy", "airr", "setuptools", "scanpy", "anndata", "numpy", "harmonypy"):
    print(f"{pkg:10s} {md.version(pkg)}")
assert md.version("harmonypy").startswith("0.0."), "harmonypy pin broken"
rows = []
for cell, (a, b) in {"c1": ("CAVAAF", "CASSAF"), "c2": ("CAVAAF", "CASSAF"), "c3": ("CAVKKF", "CASSKF")}.items():
    for locus, cdr3 in (("TRA", a), ("TRB", b)):
        rows.append(dict(cell_id=cell, sequence_id=f"{cell}_{locus}", locus=locus, productive=True,
                         junction_aa=cdr3, junction="TGT" * len(cdr3), v_call=locus + "V1", j_call=locus + "J1",
                         duplicate_count=5))
a = ir.io.read_airr(pd.DataFrame(rows))
ir.pp.index_chains(a); ir.tl.chain_qc(a)
ir.pp.ir_dist(a, metric="identity", sequence="aa")
ir.tl.define_clonotype_clusters(a, sequence="aa", metric="identity", receptor_arms="all",
                                dual_ir="primary_only", key_added="clone_id")
assert a.obs["clone_id"].nunique() == 2, a.obs
print("scirpy AIRR read + aa clonotype definition: OK (3 cells -> 2 clones)")
PYEOF
"${ENV_PREFIX}/bin/Rscript" -e 'suppressMessages({library(CellChat); library(DESeq2)}); cat("R stack still loads: CellChat", as.character(packageVersion("CellChat")), "DESeq2", as.character(packageVersion("DESeq2")), "\n")'
echo "=== DONE: all checks passed ==="

#!/usr/bin/env python3
"""Runs scripts 08 and 09 on synthetic branching data, with a planted order. Use the scrna-traj
environment. This proves the code and CellRank's plotting calls run, and that pseudotime follows the
planted order. It says nothing about whether a real trajectory is meaningful.

Usage (from the top of the repository, on a compute node): python tests/synthetic_trajectory.py OUT_DIR
"""
import pathlib
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)

# One stem state splitting into two branches; progress t runs 0 -> 1. More genes are on early in
# the stem, as CytoTRACE expects (less differentiated cells express more distinct genes).
n, n_genes = 1800, 600
t = rng.uniform(0, 1, n)
branch = rng.integers(0, 2, n)
means = np.full((n, n_genes), 0.2)
means[:, :200] += 3 * (1 - t)[:, None]                                     # stem program, fades
means[:, 200:300] += 6 * (t * (branch == 0))[:, None]                      # branch A program
means[:, 300:400] += 6 * (t * (branch == 1))[:, None]                      # branch B program
X = sp.csr_matrix(rng.poisson(means), dtype=np.float32)
labels = np.where(t < 0.3, "Stem", np.where(branch == 0, "Branch A", "Branch B"))
a = ad.AnnData(X, obs=pd.DataFrame({"cell_type_v2": labels, "truth_t": t}, index=[f"c{i}" for i in range(n)]))
a.var_names = [f"G{i}" for i in range(n_genes)]
a.layers["counts"] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
a.raw = a
sc.pp.pca(a, n_comps=20)
sc.pp.neighbors(a)
sc.tl.umap(a)
a.obs["cell_type_v2"] = a.obs["cell_type_v2"].astype("category")
a.write_h5ad(out / "input.h5ad")

subprocess.run([sys.executable, "scripts/08_trajectory_paga.py", "--input", str(out / "input.h5ad"), "--output", str(out / "08.h5ad"),
                "--plot-dir", str(out / "paga"), "--root-cluster", "Stem"], check=True)
subprocess.run([sys.executable, "scripts/09_cellrank_fate_mapping.py", "--input", str(out / "08.h5ad"),
                "--output", str(out / "new_folder" / "09.h5ad"), "--plot-dir", str(out / "cellrank"),
                "--drivers-csv", str(out / "drivers.csv"), "--n-states", "4", "--use-dpt-pseudotime"], check=True)

r = ad.read_h5ad(out / "new_folder" / "09.h5ad")
rho = stats.spearmanr(r.obs["dpt_pseudotime"], r.obs["truth_t"]).statistic
results = [("08 pseudotime follows the planted order", rho > 0.6, f"Spearman {rho:.2f}"),
           ("08 figures", (out / "paga" / "umap_dpt_pseudotime.pdf").exists(), "pdf written"),
           ("09 wrote its result into a folder that did not exist", (out / "new_folder" / "09.gpcca.pickle").exists(), ""),
           ("09 driver genes", (out / "drivers.csv").exists(), ""),
           ("09 fate figure", (out / "cellrank" / "fate_probabilities.pdf").exists(), "")]
for name, ok, detail in results:
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
print(f"\n{sum(r[1] for r in results)}/{len(results)} checks passed")
sys.exit(0 if all(r[1] for r in results) else 1)

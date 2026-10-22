#!/usr/bin/env python3
"""Test scripts/05c_score_gene_programs.py on a toy dataset with planted programs.

Two samples x 2 cell groups; 'treg' cells carry the Treg-effector genes and 'th1'
cells the Th1 genes. Checks that each planted program scores highest in its own
group (built-in programs in mouse and human symbol styles), that a custom program
CSV works, and that the per-sample subset summary is written correctly.

Usage: python tests/synthetic_programs.py OUT_DIR   (seconds, no internet)
"""
import pathlib
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
results = []


def check(name, ok, detail):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def make(species):
    rng = np.random.default_rng(3)
    treg = ["Foxp3", "Ikzf2", "Klrg1", "Prdm1", "Il10", "Gzmb", "Icos", "Tnfrsf18"]
    th1 = ["Tbx21", "Ifng", "Cxcr3", "Il12rb2", "Ccl5", "Nkg7", "Fasl"]
    genes = treg + th1 + [f"Fill{i}" for i in range(400)]
    if species == "human":
        genes = [g.upper() for g in genes]
    n = 400
    group = np.array(["treg"] * (n // 2) + ["th1"] * (n // 2))
    sample = np.tile(["S1", "S2"], n // 2)
    rate = np.full((n, len(genes)), 0.5)
    rate[:, len(treg) + len(th1):] = rng.lognormal(-0.5, 0.8, len(genes) - len(treg) - len(th1))
    rate[group == "treg", :len(treg)] = 6
    rate[group == "th1", len(treg):len(treg) + len(th1)] = 6
    X = sp.csr_matrix(rng.poisson(rate).astype(np.float32))
    a = ad.AnnData(X=X.copy(), obs=pd.DataFrame({"sample": sample, "group": group,
                                                  "condition": np.where(sample == "S1", "A", "B")},
                                                 index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    a.layers["counts"] = X
    return a


for species in ("mouse", "human"):
    f = out / f"toy_{species}.h5ad"
    make(species).write_h5ad(f)
    o = out / f"scored_{species}.h5ad"
    csv = out / f"treg_summary_{species}.csv"
    subprocess.run([sys.executable, str(ROOT / "scripts/05c_score_gene_programs.py"), "--input", str(f),
                    "--output", str(o), "--species", species, "--subset-col", "group", "--subset-values", "treg",
                    "--summary-csv", str(csv), "--group-col", "condition"], check=True)
    s = ad.read_h5ad(o).obs
    m = s.groupby("group", observed=True)[["prog_Treg effector", "prog_Th1"]].mean()
    check(f"{species}: Treg-effector program highest in planted Treg cells",
          m.loc["treg", "prog_Treg effector"] > m.loc["th1", "prog_Treg effector"] + 1, m["prog_Treg effector"].round(2).to_dict())
    check(f"{species}: Th1 program highest in planted Th1 cells",
          m.loc["th1", "prog_Th1"] > m.loc["treg", "prog_Th1"] + 1, m["prog_Th1"].round(2).to_dict())
    t = pd.read_csv(csv, index_col=0)
    check(f"{species}: subset summary counts only Treg cells, carries condition",
          int(t["n_cells"].sum()) == 200 and set(t["condition"]) == {"A", "B"}, t[["condition", "n_cells"]].to_dict())

custom = out / "custom_programs.csv"
pd.DataFrame({"program": ["MyTreg"] * 4, "gene": ["Foxp3", "Ikzf2", "Il10", "NotAGene"]}).to_csv(custom, index=False)
o = out / "scored_custom.h5ad"
subprocess.run([sys.executable, str(ROOT / "scripts/05c_score_gene_programs.py"), "--input",
                str(out / "toy_mouse.h5ad"), "--output", str(o), "--programs", str(custom)], check=True)
s = ad.read_h5ad(o)
check("custom program CSV scored (missing gene ignored)", "prog_MyTreg" in s.obs and
      list(s.uns["gene_programs"]["MyTreg"]) == ["Foxp3", "Ikzf2", "Il10"], list(s.uns["gene_programs"]["MyTreg"]))

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)

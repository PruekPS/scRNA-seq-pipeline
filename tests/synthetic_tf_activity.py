"""Synthetic check for scripts/14_tf_activity_ulm.py.

Writes a toy signed network and inputs in the formats of scripts/10 and 11,
with planted TF signal, runs script 14 on them, and checks the scores:
  DE mode    -- 'UP' targets raised, 'DOWN' targets lowered, 'MIXED' has its
                repressed targets lowered (so its activity is positive),
                'NULL' untouched, in contrast A; no signal in contrast B
  pseudobulk -- 'UP' targets raised only in pseudo-samples of condition 'treated'
No network access and no real data needed.

Usage:
    python tests/synthetic_tf_activity.py tests/data/tf_activity
"""
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "tests/data/tf_activity")
script = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "14_tf_activity_ulm.py"
rng = np.random.default_rng(0)
genes = [f"G{i}" for i in range(3000)]

net = pd.DataFrame({
    "source": ["UP"] * 40 + ["DOWN"] * 40 + ["NULL"] * 40 + ["MIXED"] * 40 + ["TINY"] * 3,
    "target": genes[:40] + genes[40:80] + genes[80:120] + genes[120:160] + genes[160:163],
    "weight": [1] * 120 + [1] * 20 + [-1] * 20 + [1] * 3,
})
(out / "de").mkdir(parents=True, exist_ok=True)
(out / "pb").mkdir(parents=True, exist_ok=True)
net.to_csv(out / "net.csv", index=False)

# --- DE mode: two contrasts written like scripts/11 output
for ct, planted in (("Cell A", True), ("Cell B", False)):
    stat = rng.normal(0, 1, len(genes))
    if planted:
        stat[:40] += 3; stat[40:80] -= 3
        stat[120:140] += 3; stat[140:160] -= 3   # activated targets up, repressed down -> MIXED active
    pd.DataFrame({"gene": genes, "stat": stat, "log2FoldChange": stat / 3, "cell_type": ct}).to_csv(
        out / "de" / f"de_{ct.replace(' ', '_')}.csv", index=False)
# the combined file must be ignored
pd.DataFrame({"gene": genes, "stat": 0.0, "cell_type": "combined"}).to_csv(
    out / "de" / "de_all_celltypes_combined.csv", index=False)

# --- pseudobulk mode: 6 pseudo-samples written like scripts/10 output (genes x pseudo-samples)
meta = pd.DataFrame({"pseudo_sample_id": [f"s{i}_Mono" for i in range(6)],
                     "sample": [f"s{i}" for i in range(6)], "cell_type": "Mono",
                     "condition": ["control"] * 3 + ["treated"] * 3})
mu = rng.gamma(2.0, 50.0, len(genes))
boost = np.r_[np.full(40, 4.0), np.ones(len(genes) - 40)]
counts = np.column_stack([rng.poisson(mu * boost if c == "treated" else mu) for c in meta["condition"]])
scipy.io.mmwrite(str(out / "pb" / "pseudobulk_counts.mtx"), sp.csr_matrix(counts))
pd.Series(genes, name="gene").to_csv(out / "pb" / "genes.csv", index=False)
meta.to_csv(out / "pb" / "pseudobulk_metadata.csv", index=False)

subprocess.run([sys.executable, str(script), "--net", str(out / "net.csv"), "--de-dir", str(out / "de"),
                "--pseudobulk-dir", str(out / "pb"), "--output-dir", str(out / "result")], check=True)

de = pd.read_csv(out / "result" / "tf_activity_de.csv", keep_default_na=False).pivot(index="tf", columns="sample_or_contrast", values="score")
pb = pd.read_csv(out / "result" / "tf_activity_pseudobulk.csv", keep_default_na=False)
up_pb = pb[pb.tf == "UP"].set_index("sample_or_contrast")

checks = [
    ("combined DE file ignored", set(de.columns) == {"Cell A", "Cell B"}),
    ("TF with < --min-targets targets dropped", "TINY" not in de.index),
    ("UP active in contrast A (t > 8)", de.loc["UP", "Cell A"] > 8),
    ("DOWN inactive in contrast A (t < -8)", de.loc["DOWN", "Cell A"] < -8),
    ("MIXED (signed targets) active in contrast A (t > 8)", de.loc["MIXED", "Cell A"] > 8),
    ("NULL ~0 in contrast A (|t| < 3)", abs(de.loc["NULL", "Cell A"]) < 3),
    ("no TF called in unplanted contrast B (|t| < 3)", (de["Cell B"].abs() < 3).all()),
    ("pseudobulk metadata joined", "condition" in pb.columns),
    ("pseudobulk UP higher in every treated than every control pseudo-sample",
     up_pb.loc[up_pb.condition == "treated", "score"].min() > up_pb.loc[up_pb.condition == "control", "score"].max()),
]
print(de.round(2).to_string())
failed = 0
for name, ok in checks:
    print(("PASS " if ok else "FAIL ") + name)
    failed += not ok
print(f"{len(checks) - failed}/{len(checks)} checks passed")
sys.exit(1 if failed else 0)

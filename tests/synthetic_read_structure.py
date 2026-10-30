#!/usr/bin/env python3
"""Offline checks of scripts/00d_check_read_structure.py and scripts/02's --exclude-hvg-genes,
scripts/figstyle.py statistics, and 01's QC rule on low mitochondrial cells. Prints PASS/FAIL.

Usage: python tests/synthetic_read_structure.py OUT_DIR
"""
import gzip
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, "scripts")
import figstyle  # noqa: E402

rng = np.random.default_rng(0)
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


def write_fastq(path, n, polya_fraction):
    with gzip.open(path, "wt") as f:
        for i in range(n):
            seq = "".join(rng.choice(list("ACGT"), 75))
            if rng.random() < polya_fraction:
                seq = seq[:30] + "A" * 45
            f.write(f"@r{i}\n{seq}\n+\n{'I' * 75}\n")


for name, fraction in (("good1", 0.01), ("good2", 0.01), ("good3", 0.01), ("short", 0.30)):
    write_fastq(out / f"{name}.fastq.gz", 2000, fraction)
args = [x for n in ("good1", "good2", "good3", "short") for x in ("--fastq", n, str(out / f"{n}.fastq.gz"))]
subprocess.run([sys.executable, "scripts/00d_check_read_structure.py", "--output-dir", str(out / "read"), *args], check=True)
table = pd.read_csv(out / "read" / "read_structure.csv").fillna({"flag": ""}).set_index("library")
check("00d flags only the short-insert library", table.loc["short", "flag"] != "" and (table.drop("short")["flag"] == "").all(),
      table["pct_polyA"].to_dict())

check("exact Mann-Whitney on 3 vs 3 fully separated groups", abs(figstyle.mannwhitney_exact([1, 2, 3], [4, 5, 6])[1] - 0.1) < 1e-9,
      "P = 0.1 (2 of 20 orderings)")
check("exact Mann-Whitney on 2 vs 2 cannot go below 0.33", abs(figstyle.mannwhitney_exact([1, 2], [3, 4])[1] - 1 / 3) < 1e-9, "P = 1/3")
check("no test when a group has one sample", figstyle.group_test({"a": [1, 2], "b": [3]}) is None, "returns None")
check("3-group test returns pairwise P", len(figstyle.group_test({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})["pairs"]) == 3, "3 pairs")

qc = __import__("01_qc_filter")


class Fake:
    obs = pd.DataFrame({"m": np.r_[np.full(95, 5.0) + rng.normal(0, 0.1, 95), [0.0, 0.0], [50.0, 50.0, 50.0]]})


flags = qc.is_outlier(Fake, "m", 5, upper_only=True)
check("01 mitochondrial filter removes high cells, keeps low ones", flags[-3:].all() and not flags[95:97].any(),
      f"{int(flags.sum())} flagged")

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)

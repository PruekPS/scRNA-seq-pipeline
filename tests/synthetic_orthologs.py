#!/usr/bin/env python3
"""Test the human -> mouse network conversion in scripts/14a_fetch_collectri.py offline.

A toy MGI homology table (same columns as HOM_MouseHumanSequence.rpt) with one-to-one,
one-to-many and renamed genes, and a toy human network. Checks that only one-to-one
orthologs are kept, renamed genes map to the MGI mouse symbol (not the case-converted
human one), one-to-many families are dropped, and the AP1/NFKB pseudo-TFs survive.

Usage: python tests/synthetic_orthologs.py OUT_DIR   (seconds, no internet)
"""
import importlib.util
import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fetch14a", ROOT / "scripts/14a_fetch_collectri.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
results = []


def check(name, ok, detail):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


rows = []  # cls, organism, symbol
for cls, h, ms in [("1", "FOXP3", ["Foxp3"]), ("2", "TBX21", ["Tbx21"]), ("3", "IFNG", ["Ifng"]),
                   ("4", "ZNF683", ["Zfp683"]),              # renamed in mouse
                   ("5", "CCL3", ["Ccl3", "Ccl3l1"]),        # one-to-many: dropped
                   ("6", "IL10", ["Il10"])]:
    rows.append({"DB Class Key": cls, "Common Organism Name": "human", "Symbol": h})
    for s in ms:
        rows.append({"DB Class Key": cls, "Common Organism Name": "mouse, laboratory", "Symbol": s})
hom = out / "toy_HOM_MouseHumanSequence.rpt"
pd.DataFrame(rows).to_csv(hom, sep="\t", index=False)
net = pd.DataFrame({"source": ["FOXP3", "FOXP3", "TBX21", "TBX21", "AP1", "NFKB", "ZNF683"],
                    "target": ["IL10", "CCL3", "IFNG", "ZNF683", "IL10", "IFNG", "IFNG"],
                    "weight": [1, 1, 1, -1, 1, 1, 1]})

# OmniPath-style rows: complexes flagged in the UniProt column, members joined in the symbol column
omni = pd.DataFrame({
    "source": ["Q01543", "COMPLEX:P05412_P17535", "COMPLEX:P01100_P17275", "COMPLEX:P19838_Q04206", "COMPLEX:P1_P2", "P0"],
    "source_genesymbol": ["FOXP3", "JUN_JUND", "FOS_JUNB", "NFKB1_RELA", "SMAD2_SMAD3", "TBX21"],
    "target": ["P22301", "P22301", "P01579", "P01579", "P01579", "COMPLEX:P9_P8"],
    "target_genesymbol": ["IL10", "IL10", "IFNG", "IFNG", "IFNG", "A_B"],
    "consensus_inhibition": [0, 0, 0, 0, 0, 0], "consensus_stimulation": [1, 1, 1, 1, 1, 1],
    "is_inhibition": [0, 0, 0, 0, 0, 0], "is_stimulation": [1, 1, 1, 1, 1, 1]})
pn = m.parse_collectri(omni)
pe = set(zip(pn["source"], pn["target"]))
check("OmniPath complexes collapse to AP1 / NFKB; other complexes and complex targets dropped",
      pe == {("FOXP3", "IL10"), ("AP1", "IL10"), ("AP1", "IFNG"), ("NFKB", "IFNG")}, sorted(pe))

o2o = m.mgi_one_to_one(str(hom))
check("one-to-one map excludes one-to-many families", "CCL3" not in o2o.index and len(o2o) == 5, o2o.to_dict())
mm = m.to_mouse(net, o2o)
edges = set(zip(mm["source"], mm["target"]))
check("renamed gene uses the MGI mouse symbol", ("Tbx21", "Zfp683") in edges and ("Zfp683", "Ifng") in edges, sorted(edges))
check("edges to one-to-many targets dropped", not any(t.lower().startswith("ccl3") for _, t in edges), sorted(edges))
check("AP1 / NFKB pseudo-TFs kept", ("AP1", "Il10") in edges and ("NFKB", "Ifng") in edges, sorted(edges))
check("signs preserved", int(mm.set_index(["source", "target"]).loc[("Tbx21", "Zfp683"), "weight"]) == -1, "Tbx21 -| Zfp683")

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)

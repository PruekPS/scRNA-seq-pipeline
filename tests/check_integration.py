#!/usr/bin/env python3
"""Check tests/integration_smoke.sh outputs against the planted truth from
tests/synthetic_splitpipe.py. Prints one PASS/FAIL line per check, exits 1 on any FAIL.

Usage: python tests/check_integration.py OUT_DIR
"""
import json
import pathlib
import sys

import anndata as ad
import pandas as pd

out = pathlib.Path(sys.argv[1])
F, R = out / "fixture", out / "run"
truth = pd.read_csv(F / "truth_cells.csv").set_index("bc_wells")
planted = json.load(open(F / "truth.json"))
results = []


def check(name, ok, detail):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")


# 00: sample sheet + sublibrary + barcode names survive to the end
a = ad.read_h5ad(R / "04b_annotated.h5ad")
obs = a.obs
check("barcode names kept (NN_NN_NN__sN)", obs.index.str.match(r"^\d{2}_\d{2}_\d{2}__s\d+$").all(),
      f"e.g. {obs.index[0]}")
check("sample sheet joined", {"patient", "timepoint"} <= set(obs.columns) and obs["patient"].notna().all(),
      f"patients {sorted(obs['patient'].astype(str).unique())}")
check("sublibrary tagged from barcode", "sublibrary" in obs and obs["sublibrary"].nunique() == 3,
      f"{sorted(obs['sublibrary'].astype(str).unique())}")

# 00a + optional unfiltered path: the planted bad cutoff is flagged, and one cutoff on
# DGE_unfiltered recovers exactly the real cells (no empty barcodes)
cc = pd.read_csv(R / "cell_call_check" / "cell_call_check.csv").fillna({"flag": ""})
flagged = {(r["run"], r["sample"]) for _, r in cc[cc["flag"] != ""].iterrows()}
want = planted["unstable_cutoff"]
check("00a flags the planted unstable cutoff", (want["run"], want["sample"]) in flagged
      and {s for _, s in flagged} == {want["sample"]}, f"flagged {sorted(flagged)}")
u = ad.read_h5ad(R / "00_raw_unfiltered.h5ad")
check("00 --min-counts on DGE_unfiltered keeps exactly the real cells",
      set(u.obs_names) == set(truth.index) and u.uns.get("load_cell_cutoff", {}).get("min_counts") == 200,
      f"{u.n_obs:,} kept of {planted['n_cells'] + planted['n_empty_barcodes']:,} barcodes "
      f"({planted['n_cells']:,} real cells)")

# 04b: annotation agrees with the planted cell types (confident labels only)
lab = obs["cell_type_v2"].astype(str)
conf = ~lab.str.startswith("Ambiguous")
acc = (lab[conf] == truth.loc[lab[conf].index, "true_type"]).mean()
check("04b annotation matches planted types", acc >= 0.9,
      f"{acc:.1%} of {conf.sum():,} confidently labeled cells ({(~conf).sum():,} Ambiguous)")

# 12: TCR matching, both routes
s12 = json.load(open(R / "tcr_plots" / "tcr_attach_summary.json"))
# not 100%: TCR cells whose WT cell was removed by QC/doublet filtering have nothing to match
n_qc_removed = planted["n_cells"] - a.n_obs
check("12 combined: match rate", s12["match_rate"] >= 0.95,
      f"{s12['match_rate']:.1%} of {s12['n_tcr_cells']:,} TCR cells (QC removed {n_qc_removed:,} WT cells)")
check("12 combined: TCR cells sit in T-lineage clusters", s12["paired_tcr_cells_in_t_lineage_fraction"] >= 0.9,
      f"{s12['paired_tcr_cells_in_t_lineage_fraction']:.1%}")
s12b = json.load(open(R / "tcr_plots_persub" / "tcr_attach_summary.json"))
check("12 per-sublibrary route gives the same cells", s12b["n_matched_to_wt"] == s12["n_matched_to_wt"],
      f"{s12b['n_matched_to_wt']:,} vs {s12['n_matched_to_wt']:,}")
check("12 all-sample folder not double-loaded", s12["n_tcr_cells"] == planted["n_tcr_cells"],
      f"{s12['n_tcr_cells']:,} TCR cells read (planted {planted['n_tcr_cells']:,})")
check("12 mis-ordered sublibraries refused", not (R / "12_should_not_exist.h5ad").exists(), "no output written")

# 13: planted clones recovered
t = ad.read_h5ad(R / "13_tcr" / "13_tcr_analyzed.h5ad").obs
t = t.join(truth[["true_clone"]], how="left")
shared = pd.read_csv(R / "13_tcr" / "shared_clones.csv", dtype={"clone_id": str})
ok_exp = []
for k in ("EXP0", "EXP1", "EXP2"):
    ids = t.loc[t["true_clone"] == k, "clone_id"].dropna().astype(str).unique()
    row = shared[shared["clone_id"].isin(ids)]
    ok_exp.append(len(ids) == 1 and len(row) == 1 and set(row.iloc[0]["samples"].split(";")) == {"P1_pre", "P1_6h"})
check("13 three expanded P1 clones found once each and shared pre/6h", all(ok_exp), f"{ok_exp}")
pub_ids = t.loc[t["true_clone"] == "PUBLIC", "clone_id"].dropna().astype(str).unique()
pub_row = shared[shared["clone_id"].isin(pub_ids)]
check("13 public clone flagged across patients", len(pub_row) == 1 and bool(pub_row.iloc[0]["public"]),
      f"{pub_row[['clone_id', 'groups']].to_dict('records')}")
exp = pd.read_csv(R / "13_tcr" / "expanded_by_celltype_pooled.csv", index_col=0)["fraction_expanded"]
top = exp.idxmax()
check("13 expansion concentrated in CD8 T cells", "CD8" in top, f"top: {top} ({exp.max():.1%})")

# 10/11: paired DE recovers planted genes without false hits
de = pd.read_csv(R / "de" / "de_all_celltypes_combined.csv")
sig = de[(de["padj"] < 0.05) & (de["log2FoldChange"].abs() > 1)]
false = sorted(set(sig["gene"]) - set(planted["planted_de_genes"]))
found = sorted(set(sig["gene"]) & set(planted["planted_de_genes"]))
check("11 significant genes are all planted 6h genes", not false, f"false hits: {false[:10]}")
check("11 planted genes recovered", len(found) >= 8, f"{len(found)}/10 planted genes significant in >=1 cell type")
check("11 design recorded", (de["design"] == "~patient + condition").all(), de["design"].iloc[0])

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)

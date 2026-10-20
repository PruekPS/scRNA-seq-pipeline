"""Step 5b: apply a curated subcluster -> final label map to 05's output.

05_refine_annotation.py proposes labels (majority CellTypist call per
lineage subcluster, doublet flags, rescued rare types); a person then
reviews marker dotplots and writes a map (configs/label_maps/*.csv) with
the final name and broad lineage for every subcluster. This script applies
it, so curated labels live in a reviewable file rather than in code.

Two maps:
  --subcluster-map  subcluster,final_label,broad_lineage[,evidence]
                    must cover every subcluster in the object (fails
                    otherwise, so a map written for a different run can't
                    be applied silently).
  --rescued-map     celltypist_label,final_label,broad_lineage (optional)
                    for cells 05 rescued as a confident rare type inside a
                    subcluster (their cell_type_refined differs from the
                    subcluster's majority); unlisted rescued labels fall
                    back to the subcluster's final label.

Doublet subclusters keep a "Doublet ..." final label and broad_lineage
"Doublet" so downstream steps can drop them explicitly.

Usage:
    python 05b_apply_label_map.py --input refined/refined.h5ad --output final_labels.h5ad \
        --subcluster-map configs/label_maps/x_subcluster_labels.csv \
        --rescued-map configs/label_maps/x_rescued_labels.csv
"""
import argparse

import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--subcluster-map", required=True, dest="subcluster_map")
    ap.add_argument("--rescued-map", default=None, dest="rescued_map")
    ap.add_argument("--sample-col", default="sample", dest="sample_col")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    for col in ("subcluster", "cell_type_refined"):
        if col not in adata.obs:
            raise SystemExit(f"'{col}' missing -- run 05_refine_annotation.py first")

    smap = pd.read_csv(args.subcluster_map).set_index("subcluster")
    sub = adata.obs["subcluster"].astype(str)
    missing = sorted(set(sub) - set(smap.index))
    if missing:
        raise SystemExit(f"subclusters not in {args.subcluster_map}: {missing}")
    extra = sorted(set(smap.index) - set(sub))
    if extra:
        print(f"Note: map rows with no cells (ignored): {extra}")

    final = sub.map(smap["final_label"]).astype(str)
    broad = sub.map(smap["broad_lineage"]).astype(str)

    refined = adata.obs["cell_type_refined"].astype(str)
    majority = refined.groupby(sub).agg(lambda x: x.value_counts().index[0])
    rescued = (refined != sub.map(majority)) & ~broad.eq("Doublet")
    n_res = int(rescued.sum())
    if args.rescued_map and n_res:
        rmap = pd.read_csv(args.rescued_map).set_index("celltypist_label")
        hit = rescued & refined.isin(rmap.index)
        final[hit] = refined[hit].map(rmap["final_label"])
        broad[hit] = refined[hit].map(rmap["broad_lineage"])
        unmapped = refined[rescued & ~hit].value_counts()
        print(f"Rescued rare-type cells: {n_res:,}; relabeled via rescued map: {int(hit.sum()):,}")
        if len(unmapped):
            print("  rescued labels without a map entry (kept subcluster label):\n" + unmapped.to_string())

    adata.obs["final_label"] = pd.Categorical(final)
    adata.obs["broad_lineage"] = pd.Categorical(broad)
    print("\nFinal labels:\n" + adata.obs["final_label"].value_counts().to_string())
    print("\nBroad lineage:\n" + adata.obs["broad_lineage"].value_counts().to_string())
    adata.write_h5ad(args.output)
    pd.crosstab(adata.obs["final_label"], adata.obs[args.sample_col]).to_csv(
        args.output.replace(".h5ad", "_counts_by_sample.csv"))
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()

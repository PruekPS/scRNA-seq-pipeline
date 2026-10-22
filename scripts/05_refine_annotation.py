"""Step 5 (optional): refine a CellTypist-annotated object by subclustering
each lineage, flagging cross-lineage doublet subclusters, and rescuing rare
cell types that a whole-dataset majority vote absorbs.

Why: 04c's majority vote runs over global Leiden clusters (~25 for ~140k
PBMC). That is robust for major populations but (a) merges real subtypes
inside one big cluster, (b) assigns rare types smaller than a cluster
(plasma cells, platelets, progenitors, basophils) to their neighbors'
label, and (c) gives doublet pockets (e.g. T-monocyte) the label of
whichever parent dominates. Re-clustering within each lineage on its own
highly variable genes resolves (a); per-cell CellTypist labels with high
confidence handle (b); a cross-lineage marker score combined with the
Scrublet score from 01 handles (c).

Per lineage group (default groups: T/NK, B/plasma, myeloid, other):
  1. log1p(CP10k) from the counts layer, HVGs (seurat_v3 on counts,
     batch-aware), scale, PCA, Harmony on --batch-key, Leiden, UMAP.
  2. Each subcluster gets the majority per-cell CellTypist label
     (--fine-col) and its agreement fraction.
  3. Doublet flag, relative to the lineage's own baseline (absolute
     thresholds don't transfer: Scrublet scores and ambient marker levels
     differ by lineage and dataset): a subcluster is flagged when its mean
     score for a *different* lineage's markers is >= --cross-fold x that
     score's median across the lineage's subclusters AND its median
     doublet_score is >= --doublet-fold x the lineage median (and >=
     --doublet-floor). Flagged cells get cell_type_refined =
     "Doublet (A+B)" so they can be dropped downstream rather than
     silently deleted here.
  4. Rare-type rescue: cells whose per-cell label is one of --rare-types
     with confidence >= --rare-conf keep that per-cell label regardless of
     their subcluster's majority -- but only if the label belongs to the
     lineage being processed (or to no lineage, e.g. HSC/MPP); a confident
     "MAIT cells" call inside the myeloid group is a doublet signal, not a
     rare myeloid type.

Outputs (in --outdir):
  refined.h5ad                     full object + lineage_group, subcluster,
                                   cell_type_refined, doublet_flag columns
  sub_<group>.h5ad                 per-lineage object with its own UMAP
  subcluster_summary.csv           per subcluster: n, majority label,
                                   agreement, doublet score, cross-lineage
                                   scores, top markers, top sample share
  umap_<group>.png, refined_counts_by_sample.csv

Usage:
    python 05_refine_annotation.py --input 04c_celltypist.h5ad --outdir refined/ \
        --broad-col celltypist_high_majority --fine-col celltypist_label --batch-key sample
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.sparse as sp

# Broad CellTypist (Immune_All_High) label -> lineage group. Matching is by
# lowercase substring, first hit wins; unmatched labels go to "other".
GROUP_RULES = [
    ("myeloid", ["monocyte", "macrophage", "dc", "dendritic", "mast", "granulocyte", "neutrophil"]),
    ("B_plasma", ["b cell", "plasma"]),
    ("T_NK", ["t cell", "nk", "ilc", "mait", "thymocyte"]),
]
# Marker sets used to score cross-lineage contamination of a subcluster.
LINEAGE_MARKERS = {
    "T_NK": ["CD3E", "CD3D", "CD2", "IL7R", "NKG7", "GNLY"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1"],
    "myeloid": ["LYZ", "CD14", "S100A8", "S100A9", "FCN1", "CST3"],
}
DEFAULT_RARE = ["Plasma cells", "Plasmablasts", "Megakaryocytes/platelets", "HSC/MPP",
                "Mast cells", "pDC", "DC1", "DC2", "MAIT cells", "Regulatory T cells",
                "gamma-delta T cells", "ILC3", "Cycling T cells", "Cycling NK cells"]


def group_of(label):
    l = str(label).lower()
    for g, keys in GROUP_RULES:
        if any(k in l for k in keys):
            return g
    return "other"


def lognorm(counts):
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    X = sp.diags(1e4 / lib) @ counts
    X.data = np.log1p(X.data)
    return X.tocsr()


def subcluster(sub, batch_key, resolution, n_hvg, n_pcs, seed):
    sub.X = lognorm(sub.layers["counts"])
    try:
        sc.pp.highly_variable_genes(sub, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts",
                                    batch_key=batch_key)
    except Exception as e:  # seurat_v3 can fail on tiny batches; fall back to per-dataset HVGs
        print(f"    batch-aware HVG failed ({e}); using dataset-wide HVGs")
        sc.pp.highly_variable_genes(sub, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts")
    hv = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    n_pcs = min(n_pcs, hv.n_obs - 1, hv.n_vars - 1)
    sc.tl.pca(hv, n_comps=n_pcs, random_state=seed)
    use_rep = "X_pca"
    if batch_key and hv.obs[batch_key].nunique() > 1:
        sce.pp.harmony_integrate(hv, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony",
                                 random_state=seed)
        use_rep = "X_pca_harmony"
    sc.pp.neighbors(hv, use_rep=use_rep, random_state=seed)
    sc.tl.leiden(hv, resolution=resolution, key_added="subcluster", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(hv, random_state=seed)
    sub.obs["subcluster"] = hv.obs["subcluster"].values
    sub.obsm["X_umap"] = hv.obsm["X_umap"]
    sub.obsm[use_rep] = hv.obsm[use_rep]
    return sub


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--broad-col", default="celltypist_high_majority", dest="broad_col")
    ap.add_argument("--fine-col", default="celltypist_label", dest="fine_col",
                    help="per-cell fine label (default: celltypist_label from 04c)")
    ap.add_argument("--conf-col", default="celltypist_conf", dest="conf_col")
    ap.add_argument("--batch-key", default="sample", dest="batch_key")
    ap.add_argument("--resolution", type=float, default=0.8)
    ap.add_argument("--n-hvg", type=int, default=2000, dest="n_hvg")
    ap.add_argument("--n-pcs", type=int, default=30, dest="n_pcs")
    ap.add_argument("--cross-fold", type=float, default=3.0, dest="cross_fold",
                    help="other-lineage marker score fold over the lineage's subcluster median to be suspect")
    ap.add_argument("--doublet-fold", type=float, default=2.5, dest="doublet_fold",
                    help="median doublet_score fold over the lineage's subcluster median to be flagged")
    ap.add_argument("--doublet-floor", type=float, default=0.1, dest="doublet_floor",
                    help="minimum median doublet_score for a flag, whatever the fold")
    ap.add_argument("--rare-types", nargs="*", default=DEFAULT_RARE, dest="rare_types")
    ap.add_argument("--rare-conf", type=float, default=0.8, dest="rare_conf")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    for col in (args.broad_col, args.fine_col, args.conf_col, args.batch_key):
        if col not in adata.obs:
            raise SystemExit(f"'{col}' not in obs")
    if "counts" not in adata.layers:
        raise SystemExit("needs raw counts in layers['counts']")

    adata.obs["lineage_group"] = adata.obs[args.broad_col].map(group_of).astype("category")
    print("Lineage groups:\n" + adata.obs["lineage_group"].value_counts().to_string())

    # per-cell cross-lineage marker scores (log1p CP10k), computed once on the full object
    all_markers = sorted({g for v in LINEAGE_MARKERS.values() for g in v if g in adata.var_names})
    mi = [adata.var_names.get_loc(g) for g in all_markers]
    lib = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    cm = adata.layers["counts"][:, mi]
    cm = cm.toarray() if sp.issparse(cm) else np.asarray(cm)
    expr = pd.DataFrame(np.log1p(cm / np.maximum(lib, 1)[:, None] * 1e4), columns=all_markers,
                        index=adata.obs_names)
    for lg, genes in LINEAGE_MARKERS.items():
        genes = [g for g in genes if g in expr.columns]
        adata.obs[f"score_{lg}"] = expr[genes].mean(axis=1).values

    adata.obs["subcluster"] = ""
    adata.obs["cell_type_refined"] = adata.obs[args.fine_col].astype(str)
    adata.obs["doublet_flag"] = False
    has_dbl = "doublet_score" in adata.obs
    summaries = []

    for grp in [g for g in ("T_NK", "B_plasma", "myeloid", "other") if (adata.obs["lineage_group"] == g).sum() > 0]:
        idx = np.where(adata.obs["lineage_group"].values == grp)[0]
        print(f"\n=== {grp}: {len(idx):,} cells ===")
        if len(idx) < 50:
            adata.obs.iloc[idx, adata.obs.columns.get_loc("subcluster")] = f"{grp}_0"
            continue
        sub = adata[idx].copy()
        sub = subcluster(sub, args.batch_key, args.resolution, args.n_hvg, args.n_pcs, args.seed)
        sc.tl.rank_genes_groups(sub, "subcluster", method="wilcoxon", n_genes=8)
        top = {c: ", ".join(sub.uns["rank_genes_groups"]["names"][c][:8])
               for c in sub.uns["rank_genes_groups"]["names"].dtype.names}

        refined = sub.obs[args.fine_col].astype(str).copy()
        others = [lg for lg in LINEAGE_MARKERS if lg != grp]
        stats = {}
        for c, g in sub.obs.groupby("subcluster", observed=True):
            stats[c] = ({lg: g[f"score_{lg}"].mean() for lg in others},
                        g["doublet_score"].median() if has_dbl else np.nan)
        base_cross = {lg: max(np.median([s[0][lg] for s in stats.values()]), 1e-3) for lg in others}
        base_dbl = max(np.median([s[1] for s in stats.values()]), 1e-3) if has_dbl else np.nan
        for c, g in sub.obs.groupby("subcluster", observed=True):
            vc = g[args.fine_col].astype(str).value_counts()
            maj, agree = vc.index[0], vc.iloc[0] / len(g)
            cross, dmed = stats[c]
            fold = {lg: cross[lg] / base_cross[lg] for lg in others}
            other_lg = max(fold, key=fold.get) if fold else None
            dfold = dmed / base_dbl if has_dbl else np.nan
            flag = bool(other_lg and has_dbl and fold[other_lg] >= args.cross_fold
                        and dfold >= args.doublet_fold and dmed >= args.doublet_floor)
            label = f"Doublet ({grp}+{other_lg})" if flag else maj
            refined.loc[g.index] = label
            s = g[args.batch_key].value_counts(normalize=True)
            summaries.append({
                "lineage_group": grp, "subcluster": f"{grp}_{c}", "n_cells": len(g),
                "majority_label": maj, "agreement": round(agree, 3),
                "second_label": vc.index[1] if len(vc) > 1 else "",
                "median_conf": round(g[args.conf_col].median(), 3),
                "median_doublet_score": round(dmed, 3) if has_dbl else None,
                **{f"mean_score_{lg}": round(v, 3) for lg, v in cross.items()},
                **{f"fold_{lg}": round(v, 2) for lg, v in fold.items()},
                "doublet_fold": round(dfold, 2) if has_dbl else None,
                "doublet_flag": flag, "refined_label": label,
                "top_sample": s.index[0], "top_sample_pct": round(100 * s.iloc[0], 1),
                "n_samples_ge1pct": int((s >= 0.01).sum()),
                "top_markers": top.get(c, ""),
            })
            print(f"  {grp}_{c:>3}: n={len(g):6,}  {maj:32s} agree={agree:.2f}  "
                  f"dbl={dmed:.2f} (x{dfold:.1f})  {other_lg} x{fold.get(other_lg, 0):.1f}  {'FLAG' if flag else ''}")

        # rare-type rescue on confident per-cell calls (not inside doublet-flagged subclusters)
        lineage_ok = sub.obs[args.fine_col].astype(str).map(lambda l: group_of(l) in (grp, "other"))
        rare = (sub.obs[args.fine_col].astype(str).isin(args.rare_types)
                & (sub.obs[args.conf_col] >= args.rare_conf)
                & lineage_ok
                & ~refined.str.startswith("Doublet"))
        refined.loc[rare] = sub.obs.loc[rare, args.fine_col].astype(str)
        print(f"  rare-type rescue: {int(rare.sum())} cells kept their confident per-cell label")

        sub.obs["cell_type_refined"] = refined.values
        sub.obs["subcluster"] = [f"{grp}_{c}" for c in sub.obs["subcluster"]]
        adata.obs.iloc[idx, adata.obs.columns.get_loc("subcluster")] = sub.obs["subcluster"].values
        adata.obs.iloc[idx, adata.obs.columns.get_loc("cell_type_refined")] = refined.values
        adata.obs.iloc[idx, adata.obs.columns.get_loc("doublet_flag")] = refined.str.startswith("Doublet").values

        fig, axes = plt.subplots(1, 3, figsize=(21, 6.5))
        sc.pl.umap(sub, color="subcluster", ax=axes[0], show=False, legend_loc="on data", legend_fontsize=6,
                   frameon=False, size=2, title=f"{grp}: subclusters")
        sc.pl.umap(sub, color="cell_type_refined", ax=axes[1], show=False, frameon=False, size=2,
                   legend_fontsize=6, title="refined label")
        sc.pl.umap(sub, color=args.batch_key, ax=axes[2], show=False, frameon=False, size=2,
                   legend_fontsize=6, title=args.batch_key)
        fig.savefig(out / f"umap_{grp}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        sub.uns.pop("rank_genes_groups", None)
        sub.write_h5ad(out / f"sub_{grp}.h5ad")

    adata.obs["subcluster"] = adata.obs["subcluster"].astype("category")
    adata.obs["cell_type_refined"] = adata.obs["cell_type_refined"].astype("category")
    pd.DataFrame(summaries).to_csv(out / "subcluster_summary.csv", index=False)
    pd.crosstab(adata.obs["cell_type_refined"], adata.obs[args.batch_key]).to_csv(
        out / "refined_counts_by_sample.csv")
    print("\nRefined labels:\n" + adata.obs["cell_type_refined"].value_counts().to_string())
    print(f"Doublet-flagged cells: {int(adata.obs['doublet_flag'].sum()):,}")
    adata.write_h5ad(out / "refined.h5ad")
    print(f"Wrote {out}/refined.h5ad")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Step 13: clonal-structure analysis of the WT + TCR object from step 12.

Three questions, all dataset-agnostic (sample / group / annotation columns are
runtime arguments, nothing about any one study is hard-coded):

1. Repertoire diversity per sample -- number of clonotypes, fraction of cells in
   expanded clones, largest-clone fraction, and normalized Shannon entropy.
   Diversity depends strongly on how many cells were captured, so it is also
   reported RAREFIED: every sample is repeatedly subsampled to the same number
   of cells and the metric averaged. Compare samples on the rarefied columns.
2. Clones shared across samples -- which clonotypes occur in more than one
   sample (e.g. pre vs post, blood vs tissue), pairwise shared-clone counts and
   Jaccard overlap, and "public" clones seen in more than one --group-col level
   (e.g. more than one patient).
3. Phenotype of expanded clones -- for each cell type, the fraction of TCR cells
   that belong to an expanded clone, per sample and pooled, plus the cell-type
   make-up of the largest clones.

Statistics note: pooled per-cell counts are DESCRIPTIVE. Cells from one sample
are not independent replicates, so any test of "condition A vs B" should treat
per-sample fractions (expanded_by_celltype_by_sample.csv) as the data points,
with samples or patients as replicates -- the same reasoning as pseudobulk DE
(scripts 10/11). No p-values are produced here on purpose.

By default only cells with a full alpha+beta pair enter the clonotype analyses
(`--include-orphans` to keep single-chain cells); cells with one chain cannot be
assigned to a clone reliably.

Usage:
  python 13_tcr_analysis.py --input 12_tcr.h5ad --output-dir tcr_results \\
      --sample-col sample --group-col patient --annotation-col cell_type_v2
"""
import argparse
import itertools
import json
import pathlib
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scirpy as ir

import figstyle

warnings.filterwarnings("ignore", category=FutureWarning)

PAIRED = {"single pair", "extra VJ", "extra VDJ", "two full chains"}


def norm_shannon(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    if len(counts) < 2:
        return 0.0
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum() / np.log(len(counts)))


def diversity_table(tcr: pd.DataFrame, sample_col: str, expanded_min: int, rarefy_to: int,
                    n_iter: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for s, g in tcr.groupby(sample_col, observed=True):
        sizes = g["clone_id"].value_counts().to_numpy()
        row = {sample_col: s, "n_tcr_cells": len(g), "n_clonotypes": len(sizes),
               "frac_cells_in_expanded_clones": float(sizes[sizes >= expanded_min].sum() / len(g)),
               "largest_clone_fraction": float(sizes.max() / len(g)),
               "normalized_shannon": norm_shannon(sizes)}
        if rarefy_to and len(g) >= rarefy_to:
            ids = g["clone_id"].to_numpy()
            ent, exp = [], []
            for _ in range(n_iter):
                sub = pd.Series(rng.choice(ids, rarefy_to, replace=False)).value_counts().to_numpy()
                ent.append(norm_shannon(sub))
                exp.append(sub[sub >= expanded_min].sum() / rarefy_to)
            row.update({"rarefied_to": rarefy_to,
                        "rarefied_normalized_shannon_mean": float(np.mean(ent)),
                        "rarefied_normalized_shannon_sd": float(np.std(ent)),
                        "rarefied_clonality_mean": float(1 - np.mean(ent)),
                        "rarefied_frac_expanded_mean": float(np.mean(exp))})
        else:
            row.update({"rarefied_to": rarefy_to, "rarefied_normalized_shannon_mean": np.nan,
                        "rarefied_normalized_shannon_sd": np.nan, "rarefied_clonality_mean": np.nan,
                        "rarefied_frac_expanded_mean": np.nan})
        rows.append(row)
    return pd.DataFrame(rows)


def sharing_tables(tcr: pd.DataFrame, sample_col: str, group_col: str | None, ann_col: str | None):
    counts = tcr.groupby(["clone_id", sample_col], observed=True).size().unstack(fill_value=0)
    present = counts > 0
    shared = counts[present.sum(axis=1) >= 2]
    rows = []
    for clone, c in shared.iterrows():
        samples = c[c > 0]
        cells = tcr[tcr["clone_id"] == clone]
        row = {"clone_id": clone, "n_samples": len(samples), "samples": ";".join(map(str, samples.index)),
               "total_cells": int(samples.sum()),
               "cells_per_sample": ";".join(f"{k}:{v}" for k, v in samples.items())}
        if group_col:
            groups = sorted(cells[group_col].astype(str).unique())
            row.update({"groups": ";".join(groups), "n_groups": len(groups), "public": len(groups) >= 2})
        if ann_col:
            row["dominant_cell_type"] = cells[ann_col].astype(str).value_counts().index[0]
        rows.append(row)
    long = pd.DataFrame(rows).sort_values(["n_samples", "total_cells"], ascending=False) if rows else pd.DataFrame()

    samples = list(present.columns)
    n_shared = pd.DataFrame(0, index=samples, columns=samples)
    jacc = pd.DataFrame(np.nan, index=samples, columns=samples)
    for a, b in itertools.product(samples, samples):
        inter = int((present[a] & present[b]).sum())
        union = int((present[a] | present[b]).sum())
        n_shared.loc[a, b] = inter
        jacc.loc[a, b] = inter / union if union else np.nan
    return long, n_shared, jacc


def expansion_by_celltype(tcr: pd.DataFrame, sample_col: str, ann_col: str):
    pooled = tcr.groupby(ann_col, observed=True)["expanded"].agg(n_tcr_cells="size", n_in_expanded_clones="sum")
    pooled["fraction_expanded"] = pooled["n_in_expanded_clones"] / pooled["n_tcr_cells"]
    per_sample = (tcr.groupby([sample_col, ann_col], observed=True)["expanded"]
                  .agg(n_tcr_cells="size", n_in_expanded_clones="sum").reset_index())
    per_sample["fraction_expanded"] = per_sample["n_in_expanded_clones"] / per_sample["n_tcr_cells"]
    return pooled.sort_values("fraction_expanded", ascending=False), per_sample


def top_clones(adata: ad.AnnData, tcr: pd.DataFrame, sample_col: str, ann_col: str | None, n: int):
    sizes = tcr["clone_id"].value_counts().head(n)
    chains = ir.get.airr(adata, ["junction_aa", "v_call", "j_call"], ["VJ_1", "VDJ_1"])
    rows, comp = [], {}
    for clone, size in sizes.items():
        cells = tcr.index[tcr["clone_id"] == clone]
        rep = chains.loc[cells[0]]
        row = {"clone_id": clone, "n_cells": int(size),
               "TRA_cdr3_aa": rep.get("VJ_1_junction_aa"), "TRB_cdr3_aa": rep.get("VDJ_1_junction_aa"),
               "TRAV": rep.get("VJ_1_v_call"), "TRBV": rep.get("VDJ_1_v_call"),
               "samples": ";".join(f"{k}:{v}" for k, v in tcr.loc[cells, sample_col].value_counts().items() if v > 0)}
        if ann_col:
            frac = tcr.loc[cells, ann_col].astype(str).value_counts(normalize=True)
            comp[clone] = frac
            row["dominant_cell_type"] = frac.index[0]
            row["dominant_fraction"] = float(frac.iloc[0])
        rows.append(row)
    return pd.DataFrame(rows), (pd.DataFrame(comp).T.fillna(0) if comp else pd.DataFrame())


def heatmap(df: pd.DataFrame, label: str, path: pathlib.Path, fmt: str = "{:.2f}", cmap: str = figstyle.SEQUENTIAL):
    if df.empty:
        return
    fig, axes = figstyle.subplots(1, 1, w_mm=min(183, 25 + 6.5 * df.shape[1]), h_mm=min(240, 20 + 5 * df.shape[0]))
    ax = axes[0, 0]
    values = df.to_numpy(dtype=float)
    im = ax.imshow(values, cmap=cmap, aspect="auto")
    ax.set_xticks(range(df.shape[1]), df.columns.astype(str), rotation=90)
    ax.set_yticks(range(df.shape[0]), df.index.astype(str))
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if df.size <= 400:
        middle = np.nanmin(values) + 0.5 * (np.nanmax(values) - np.nanmin(values))
        for (i, j), v in np.ndenumerate(values):
            if not np.isnan(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=5,
                        color="white" if v < middle else figstyle.INK)
    colorbar = fig.colorbar(im, ax=ax, shrink=0.6, label=label)
    colorbar.outline.set_linewidth(0.4)
    figstyle.save_figure(fig, path, legend=f"{label}.", source_data=df, check=False)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help=".h5ad from scripts/12_tcr_attach.py")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--group-col", default=None, help="e.g. patient/donor; flags clones seen in >1 group as public")
    p.add_argument("--annotation-col", default="cell_type_v2")
    p.add_argument("--expanded-min-size", type=int, default=2,
                   help="a clone counts as expanded if it has at least this many cells within the same sample")
    p.add_argument("--rarefy-to", type=int, default=None,
                   help="cells per sample for rarefied diversity (default: smallest sample with >= --min-cells)")
    p.add_argument("--min-cells", type=int, default=50, help="samples with fewer TCR cells are not rarefied")
    p.add_argument("--n-rarefy", type=int, default=50)
    p.add_argument("--top-n-clones", type=int, default=20)
    p.add_argument("--embedding", default="X_umap")
    p.add_argument("--include-orphans", action="store_true")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    figstyle.apply_style()
    rng = np.random.default_rng(args.seed)

    adata = ad.read_h5ad(args.input)
    obs = adata.obs
    for col in [args.sample_col] + ([args.group_col] if args.group_col else []):
        if col not in obs.columns:
            raise SystemExit(f"--{'group' if col == args.group_col else 'sample'}-col '{col}' not in obs: {list(obs.columns)}")
    ann_col = args.annotation_col if args.annotation_col in obs.columns else None
    if ann_col is None:
        print(f"  note: annotation column '{args.annotation_col}' not found -- skipping phenotype outputs")

    keep = obs["clone_id"].notna()
    if not args.include_orphans:
        keep &= obs["chain_pairing"].isin(PAIRED)
    tcr = obs.loc[keep].copy()
    if tcr.empty:
        raise SystemExit("No cell has a usable TCR for clonotype analysis. Check step 12's match rate and "
                         "chain-pairing report; --include-orphans also keeps single-chain cells.")
    tcr["clone_id"] = tcr["clone_id"].astype(str)
    print(f"Cells used for clonotype analyses: {len(tcr):,} of {adata.n_obs:,} "
          f"({'paired alpha+beta only' if not args.include_orphans else 'including single-chain cells'})")

    # clone size counted WITHIN each sample: a clone with 1 cell in each of 3 samples is not expanded
    tcr["clone_size_in_sample"] = tcr.groupby([args.sample_col, "clone_id"], observed=True)["clone_id"].transform("size")
    tcr["expanded"] = tcr["clone_size_in_sample"] >= args.expanded_min_size

    # 1. diversity
    per_sample_n = tcr.groupby(args.sample_col, observed=True).size()
    rarefy_to = args.rarefy_to or int(per_sample_n[per_sample_n >= args.min_cells].min() if (per_sample_n >= args.min_cells).any() else 0)
    div = diversity_table(tcr, args.sample_col, args.expanded_min_size, rarefy_to, args.n_rarefy, rng)
    if args.group_col:
        div = div.merge(tcr.groupby(args.sample_col, observed=True)[args.group_col].first().reset_index(),
                        on=args.sample_col, how="left")
    div.to_csv(out / "repertoire_by_sample.csv", index=False)
    print(f"  diversity: {len(div)} samples (rarefied to {rarefy_to} cells, {args.n_rarefy} draws)")

    fig, axes = figstyle.subplots(1, 1, w_mm=min(183, 25 + 6 * len(div)), h_mm=50)
    ax = axes[0, 0]
    col = "rarefied_clonality_mean" if div["rarefied_clonality_mean"].notna().any() else None
    y = div[col] if col else 1 - div["normalized_shannon"]
    err = div["rarefied_normalized_shannon_sd"] if col else None
    ax.bar(div[args.sample_col].astype(str), y, yerr=err, color=figstyle.PALETTE[0], capsize=1, linewidth=0,
           error_kw={"elinewidth": 0.5, "capthick": 0.5})
    ax.set_ylabel("Clonality")
    ax.tick_params(axis="x", rotation=90)
    figstyle.save_figure(fig, out / "clonality_by_sample",
                         legend="Clonality = 1 - normalized Shannon entropy of clone sizes, "
                                + (f"mean of {args.n_rarefy} draws of {rarefy_to} cells per sample; bars, s.d."
                                   if col else "NOT rarefied (no sample had enough TCR cells)."),
                         source_data=div.set_index(args.sample_col))

    # 2. sharing
    long, n_shared, jacc = sharing_tables(tcr, args.sample_col, args.group_col, ann_col)
    long.to_csv(out / "shared_clones.csv", index=False)
    n_shared.to_csv(out / "sample_overlap_counts.csv")
    jacc.to_csv(out / "sample_overlap_jaccard.csv")
    heatmap(n_shared, "Clonotypes shared between samples", out / "sample_overlap_counts", fmt="{:.0f}")
    heatmap(jacc, "Clonotype overlap between samples (Jaccard index)", out / "sample_overlap_jaccard")
    n_public = int(long["public"].sum()) if "public" in long.columns and len(long) else 0
    print(f"  sharing: {len(long):,} clonotypes in >= 2 samples"
          + (f"; {n_public} seen in >= 2 '{args.group_col}' levels (public)" if args.group_col else ""))

    # 3. phenotype of expanded clones
    summary = {"n_cells_used": int(len(tcr)), "rarefied_to": rarefy_to,
               "n_clonotypes": int(tcr["clone_id"].nunique()),
               "n_shared_clonotypes": int(len(long)), "n_public_clonotypes": n_public,
               "expanded_min_size": args.expanded_min_size}
    if ann_col:
        pooled, per_sample = expansion_by_celltype(tcr, args.sample_col, ann_col)
        pooled.to_csv(out / "expanded_by_celltype_pooled.csv")
        per_sample.to_csv(out / "expanded_by_celltype_by_sample.csv", index=False)
        fig, axes = figstyle.subplots(1, 1, w_mm=figstyle.SINGLE, h_mm=4 * len(pooled) + 15)
        ax = axes[0, 0]
        order = list(pooled.index.astype(str))
        ax.barh(order, pooled["fraction_expanded"], color="#B3B3B3", label="Pooled", linewidth=0)
        for i, ct in enumerate(order):
            vals = per_sample.loc[per_sample[ann_col].astype(str) == ct, "fraction_expanded"]
            ax.scatter(vals, np.full(len(vals), i) + rng.uniform(-0.15, 0.15, len(vals)), s=4, color=figstyle.PALETTE[0], linewidths=0,
                       label="Per sample" if i == 0 else None, zorder=3)
        ax.invert_yaxis()
        ax.set_xlabel("TCR cells in expanded clones (fraction)")
        ax.legend(loc="best")
        figstyle.save_figure(fig, out / "expanded_by_celltype",
                             legend=f"Expanded clone: at least {args.expanded_min_size} cells within the same sample. "
                                    "Bars, all samples pooled; dots, individual samples.", source_data=pooled)
        summary["fraction_expanded_by_celltype_pooled"] = pooled["fraction_expanded"].round(4).to_dict()

    top, comp = top_clones(adata, tcr, args.sample_col, ann_col, args.top_n_clones)
    top.to_csv(out / "top_clones.csv", index=False)
    if not comp.empty:
        heatmap(comp, f"Fraction of each of the {len(comp)} largest clonotypes per cell type",
                out / "top_clones_phenotype")

    if args.embedding in adata.obsm:
        emb = adata.obsm[args.embedding]
        cat = pd.Series("no paired TCR", index=adata.obs_names)
        cat.loc[tcr.index] = np.where(tcr["expanded"], "expanded", "singleton")
        fig, axes = figstyle.subplots(1, 1, w_mm=figstyle.SINGLE, h_mm=70)
        ax = axes[0, 0]
        base = figstyle.point_size(adata.n_obs)
        for label, color, size in (("no paired TCR", "#E6E6E6", base), ("singleton", figstyle.PALETTE[5], base * 1.5),
                                   ("expanded", figstyle.PALETTE[3], base * 2.5)):
            m = (cat == label).to_numpy()
            ax.scatter(emb[m, 0], emb[m, 1], s=size, c=color, label=f"{label[0].upper() + label[1:]} ({m.sum():,})", linewidths=0,
                       rasterized=True)
        figstyle.clean_embedding(ax, arrows=True, x=emb[:, 0], y=emb[:, 1])
        ax.legend(markerscale=3, loc="upper left", bbox_to_anchor=(1.0, 1.0), handletextpad=0.2)
        figstyle.save_figure(fig, out / "embedding_expansion", check=False,
                             legend="UMAP coloured by clone status; numbers are cells in each category.")

    adata.obs["clone_size_in_sample"] = tcr["clone_size_in_sample"].reindex(adata.obs_names)
    adata.obs["clone_expanded"] = tcr["expanded"].reindex(adata.obs_names).astype("object").astype(str)
    (out / "tcr_analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    adata.write_h5ad(out / "13_tcr_analyzed.h5ad")
    print(f"Wrote outputs to {out}")


if __name__ == "__main__":
    main()

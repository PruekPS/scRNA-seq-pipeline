"""Step 4b (replaces 04's demo annotation): finer marker set + confidence-gated
scoring to address the scale-breakdown documented in
docs/markers_annotation_deep_dive.html and docs/real_data_deep_dive.html.

Two concrete, documented failures in the original 6-category
scripts/04_markers_annotation.py demo dictionary, both from the real 50k-cell
Parse subset run at 20 real clusters (resolution 0.6):
  1. "Least-negative-wins" artifact: argmax over raw mean marker scores always
     picks a label even when every score is weak/negative -- several clusters
     got confidently labeled despite no real marker signal.
  2. Too coarse: 6 categories can't cover real PBMC heterogeneity at ~20
     clusters -- e.g. several genuinely distinct T-cell states all collapsed
     into "T cell", and some got mislabeled "Platelet" off weak PPBP/PF4
     ambient signal.

This script does NOT change scripts/04_markers_annotation.py -- that stays as
the original, already-validated PBMC3k smoke-test demo. This is a separate,
additive step for real-scale data:
  - A finer-grained marker dictionary (T-cell subsets, monocyte subsets,
    pDC vs cDC, plasma cells, HSPCs) instead of 6 broad buckets.
  - Cluster scores are z-scored ACROSS CLUSTERS per label before comparison,
    so labels with naturally higher/lower raw score ranges don't dominate.
  - A cluster is only assigned a label if (a) the top z-score is positive AND
    (b) it beats the runner-up by --margin (default 0.5 SD). Otherwise it's
    marked "Ambiguous (top1 / top2)" for manual review instead of forced.

This is still a marker-dictionary heuristic, not a reference-based classifier
(e.g. CellTypist, Azimuth) -- see docs/pipeline_overview.html section 06 for
that flagged-but-not-yet-adopted alternative. Still human-PBMC-specific;
swap FINE_PBMC_MARKERS for a tissue-appropriate reference for other data.

NOT YET TEST-RUN as of this commit (see project memory) -- validate against
tests/data/ (or a small real subset) before trusting output on a full-scale run.

Input: the clustered .h5ad written by scripts/03_cluster.py (optionally after
02b_harmony_integration.py).

Usage:
    python 04b_better_annotation.py --input clustered.h5ad --output annotated_v2.h5ad \
        --resolution 0.6 --plot-dir marker_plots_v2/ --margin 0.5
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Finer-grained than scripts/04_markers_annotation.py's PBMC_MARKERS -- still
# human PBMC, still canonical (Satija lab / Azimuth PBMC reference markers),
# just more categories so real clusters have somewhere specific to land
# instead of collapsing into one of 6 broad buckets.
FINE_PBMC_MARKERS = {
    "CD4 T cell": ["IL7R", "CD3D", "CD3E", "CD4"],
    "CD8 T cell": ["CD8A", "CD8B", "CD3D", "CD3E"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Plasma cell": ["MZB1", "JCHAIN", "XBP1"],
    "CD14+ Monocyte": ["CD14", "LYZ", "S100A8", "S100A9"],
    "CD16+ Monocyte": ["FCGR3A", "MS4A7", "CDKN1C"],
    "cDC": ["FCER1A", "CLEC10A", "CD1C"],
    "pDC": ["LILRA4", "IL3RA", "GZMB"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
    "HSPC": ["CD34", "SOX4", "SPINK2"],
}


def detect_markers(adata, groupby, n_genes=25):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=True)
    return sc.get.rank_genes_groups_df(adata, group=None).groupby("group").head(n_genes)


def annotate_with_confidence(adata, groupby, marker_dict, margin=0.5):
    present = {
        label: [g for g in genes if g in adata.raw.var_names]
        for label, genes in marker_dict.items()
    }
    missing = {label: genes for label, genes in marker_dict.items() if not present.get(label)}
    if missing:
        print(f"Warning: no marker genes found in var_names for: {list(missing)} -- skipping these labels")

    for label, genes in present.items():
        if genes:
            sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{label}", use_raw=True)

    score_cols = [f"score_{label}" for label in present if present[label]]
    labels = [c.replace("score_", "", 1) for c in score_cols]
    cluster_scores = adata.obs.groupby(groupby)[score_cols].mean()
    cluster_scores.columns = labels

    # z-score each label's column across clusters -- puts every marker set on
    # a comparable scale before we compare them to each other within a row.
    std = cluster_scores.std(axis=0, ddof=0).replace(0, 1.0)
    z = (cluster_scores - cluster_scores.mean(axis=0)) / std

    results = {}
    for cluster, row in z.iterrows():
        ranked = row.sort_values(ascending=False)
        top_label, top_score = ranked.index[0], ranked.iloc[0]
        second_label, second_score = ranked.index[1], ranked.iloc[1]
        confident = (top_score > 0) and (top_score - second_score >= margin)
        if confident:
            results[cluster] = top_label
        else:
            results[cluster] = f"Ambiguous ({top_label} / {second_label})"

    adata.obs["cell_type_v2"] = adata.obs[groupby].map(results).astype(str)
    return adata, cluster_scores, z, results


def save_marker_plots(adata, groupby, marker_dict, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    all_markers = [g for genes in marker_dict.values() for g in genes if g in adata.raw.var_names]
    fig = sc.pl.dotplot(adata, all_markers, groupby=groupby, use_raw=True, show=False, return_fig=True)
    fig.savefig(plot_dir / "marker_dotplot_fine.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color="cell_type_v2", ax=ax, show=False, legend_loc="right margin", legend_fontsize=7)
    fig.tight_layout()
    fig.savefig(plot_dir / "umap_annotated_v2.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to clustered .h5ad from 03_cluster.py")
    parser.add_argument("--output", required=True, help="Path to write the annotated .h5ad")
    parser.add_argument("--plot-dir", default="marker_plots_v2", help="Directory for marker/annotation plots")
    parser.add_argument("--resolution", default="0.6", help="Which leiden_r<resolution> column to annotate")
    parser.add_argument("--margin", type=float, default=0.5,
                         help="Minimum z-score gap (in SDs) between top and runner-up label "
                              "to accept a confident assignment; below this, cluster is 'Ambiguous'")
    parser.add_argument("--markers-csv", default=None, help="Optional path to write the top-markers-per-cluster table")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    groupby = f"leiden_r{args.resolution}"
    if groupby not in adata.obs:
        raise SystemExit(f"{groupby} not found -- available: {[c for c in adata.obs if c.startswith('leiden_r')]}")

    markers_df = detect_markers(adata, groupby)
    print(markers_df.groupby("group").head(5).to_string(index=False))

    adata, cluster_scores, z, results = annotate_with_confidence(adata, groupby, FINE_PBMC_MARKERS, margin=args.margin)
    print()
    print("Per-cluster z-scored marker-set scores:")
    print(z.round(2).to_string())
    print()
    n_ambiguous = sum(1 for v in results.values() if v.startswith("Ambiguous"))
    print(f"Assigned {len(results) - n_ambiguous}/{len(results)} clusters confidently "
          f"(margin >= {args.margin} SD); {n_ambiguous} flagged Ambiguous for manual review.")
    print(pd.Series(results).to_string())

    save_marker_plots(adata, groupby, FINE_PBMC_MARKERS, args.plot_dir)

    if args.markers_csv:
        markers_df.to_csv(args.markers_csv, index=False)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

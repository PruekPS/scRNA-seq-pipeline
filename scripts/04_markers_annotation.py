"""Step 4: Marker detection and (demo) marker-based annotation.

Per docs/pipeline_overview.html section 06 default path, item 5: pairwise
marker detection (Wilcoxon rank-sum, the scanpy equivalent of Seurat's
FindMarkers), then annotation.

The canonical-marker dictionary below is PBMC-specific and exists only to
exercise the annotation step against the PBMC3k smoke-test data (see
tests/README.md) -- it is NOT a general-purpose annotation reference.
Real Parse Evercode data will need tissue-appropriate markers or a
reference-based method (see docs/pipeline_overview.html section 06 for
the flagged, not-yet-adopted foundation-model alternative and its caveats).

Input: the clustered .h5ad written by scripts/03_cluster.py.

Usage:
    python 04_markers_annotation.py --input clustered.h5ad --output annotated.h5ad \
        --resolution 0.6 --plot-dir marker_plots/
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

# Canonical PBMC markers (Satija lab / scanpy PBMC3k tutorials) -- smoke-test only.
# Kept deliberately specific per category -- see docs/markers_annotation_deep_dive.html
# section 04 for why: an earlier version included CST3/HLA-DPB1/HLA-DPA1 under
# "Dendritic cell", which are also broadly expressed in monocytes and ended up
# stealing two correctly-labeled monocyte clusters. FCER1A/CLEC10A are more DC-specific.
PBMC_MARKERS = {
    "T cell": ["CD3D", "CD3E", "CD3G", "IL7R"],
    "B cell": ["MS4A1", "CD79A", "CD79B"],
    "Monocyte": ["CD14", "LYZ", "FCGR3A", "MS4A7", "S100A8", "S100A9"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
    "Platelet": ["PPBP", "PF4"],
    "Dendritic cell": ["FCER1A", "CLEC10A"],
}


def detect_markers(adata, groupby, n_genes=25):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=True)
    return sc.get.rank_genes_groups_df(adata, group=None).groupby("group").head(n_genes)


def annotate_by_markers(adata, groupby, marker_dict):
    present = {
        label: [g for g in genes if g in adata.raw.var_names]
        for label, genes in marker_dict.items()
    }
    for label, genes in present.items():
        if genes:
            sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{label}", use_raw=True)

    score_cols = [f"score_{label}" for label in present if present[label]]
    cluster_scores = adata.obs.groupby(groupby)[score_cols].mean()
    cluster_label = cluster_scores.idxmax(axis=1).str.replace("score_", "", regex=False)

    adata.obs["cell_type_demo"] = adata.obs[groupby].map(cluster_label).astype(str)
    return adata, cluster_scores


def save_marker_plots(adata, groupby, marker_dict, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False, show=False)
    plt.gcf().savefig(plot_dir / "rank_genes_groups.png", dpi=150, bbox_inches="tight")
    plt.close(plt.gcf())

    all_markers = [g for genes in marker_dict.values() for g in genes if g in adata.raw.var_names]
    fig = sc.pl.dotplot(adata, all_markers, groupby=groupby, use_raw=True, show=False, return_fig=True)
    fig.savefig(plot_dir / "marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(adata, color="cell_type_demo", ax=ax, show=False, legend_loc="right margin")
    fig.tight_layout()
    fig.savefig(plot_dir / "umap_annotated_demo.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to clustered .h5ad from 03_cluster.py")
    parser.add_argument("--output", required=True, help="Path to write the annotated .h5ad")
    parser.add_argument("--plot-dir", default="marker_plots", help="Directory for marker/annotation plots")
    parser.add_argument("--resolution", default="0.6", help="Which leiden_r<resolution> column to annotate")
    parser.add_argument("--markers-csv", default=None, help="Optional path to write the top-markers-per-cluster table")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    groupby = f"leiden_r{args.resolution}"
    if groupby not in adata.obs:
        raise SystemExit(f"{groupby} not found -- available: {[c for c in adata.obs if c.startswith('leiden_r')]}")

    markers_df = detect_markers(adata, groupby)
    print(markers_df.groupby("group").head(5).to_string(index=False))

    adata, cluster_scores = annotate_by_markers(adata, groupby, PBMC_MARKERS)
    print()
    print("Demo annotation (mean marker-set score per cluster, argmax label):")
    print(cluster_scores.round(3).to_string())
    print()
    print(adata.obs.groupby(groupby)["cell_type_demo"].first().to_string())

    save_marker_plots(adata, groupby, PBMC_MARKERS, args.plot_dir)

    if args.markers_csv:
        markers_df.to_csv(args.markers_csv, index=False)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

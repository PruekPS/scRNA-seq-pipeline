"""Step 3: Neighbor graph, Leiden clustering (resolution scan), UMAP.

Per docs/pipeline_overview.html section 06 default path, item 4: graph-based
Leiden with a resolution scan rather than one fixed resolution. UMAP is
computed purely for visualization, not as an input to any later
quantitative step.

Input: the PCA-ready .h5ad written by scripts/02_normalize_pca.py.

Usage:
    python 03_cluster.py --input pca.h5ad --output clustered.h5ad --plot-dir cluster_plots/
    python 03_cluster.py --input pca.h5ad --output clustered.h5ad --resolutions 0.2,0.5,1.0
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def save_cluster_plots(adata, resolutions, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    n = len(resolutions)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4))
    axes = [axes] if n == 1 else list(axes)
    for ax, res in zip(axes, resolutions):
        sc.pl.umap(adata, color=f"leiden_r{res}", ax=ax, show=False, legend_loc="on data", legend_fontsize=6)
    fig.tight_layout()
    fig.savefig(plot_dir / "umap_by_resolution.png", dpi=150)
    plt.close(fig)

    # Sanity check: clusters shouldn't just be a proxy for sequencing depth.
    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(adata, color="total_counts", ax=ax, show=False)
    fig.tight_layout()
    fig.savefig(plot_dir / "umap_by_total_counts.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to PCA-ready .h5ad from 02_normalize_pca.py")
    parser.add_argument("--output", required=True, help="Path to write the clustered .h5ad")
    parser.add_argument("--plot-dir", default="cluster_plots", help="Directory for cluster/UMAP plots")
    parser.add_argument("--n-pcs", type=int, default=50, dest="n_pcs", help="Number of PCs to use for the neighbor graph")
    parser.add_argument(
        "--resolutions", default="0.2,0.4,0.6,0.8,1.0",
        help="Comma-separated Leiden resolutions to scan",
    )
    args = parser.parse_args()
    resolutions = [float(r) for r in args.resolutions.split(",")]

    adata = sc.read_h5ad(args.input)

    sc.pp.neighbors(adata, n_pcs=args.n_pcs, use_rep="X_pca")
    sc.tl.umap(adata)

    for res in resolutions:
        sc.tl.leiden(
            adata, resolution=res, key_added=f"leiden_r{res}",
            flavor="igraph", n_iterations=2,
        )
        n_clusters = adata.obs[f"leiden_r{res}"].nunique()
        print(f"resolution={res}: {n_clusters} clusters")

    save_cluster_plots(adata, resolutions, args.plot_dir)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

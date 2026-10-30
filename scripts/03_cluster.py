"""Step 3: Neighbor graph, Leiden clustering (resolution scan), UMAP.

Per docs/reference/design-decisions.md roadmap, item 4: graph-based
Leiden with a resolution scan rather than one fixed resolution. UMAP is
computed purely for visualization, not as an input to any later
quantitative step.

Input: the PCA-ready .h5ad written by scripts/02_normalize_pca.py, or the
Harmony-corrected .h5ad written by scripts/02b_harmony_integration.py (pass
--use-rep X_pca_harmony in that case -- see pipeline/run_downstream.slurm).

Usage:
    python 03_cluster.py --input pca.h5ad --output clustered.h5ad --plot-dir cluster_plots/
    python 03_cluster.py --input pca.h5ad --output clustered.h5ad --resolutions 0.2,0.5,1.0
    python 03_cluster.py --input harmony.h5ad --output clustered.h5ad --use-rep X_pca_harmony
"""
import argparse
import pathlib

import numpy as np
import scanpy as sc

import figstyle


def save_cluster_plots(adata, resolutions, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    ncols = min(3, len(resolutions))
    nrows = int(np.ceil(len(resolutions) / ncols))
    fig, axes = figstyle.subplots(nrows, ncols, w_mm=figstyle.DOUBLE, h_mm=55 * nrows)
    for ax, res in zip(axes.flat, resolutions):
        figstyle.umap(adata, f"leiden_r{res}", ax, title=f"Resolution {res}", legend_loc="on data",
                      legend_fontoutline=1)
    for ax in list(axes.flat)[len(resolutions):]:
        ax.axis("off")
    figstyle.save_figure(fig, plot_dir / "umap_by_resolution", check=False,
                         legend="UMAP coloured by Leiden cluster at each resolution scanned.")

    # Sanity check: clusters shouldn't just be a proxy for sequencing depth.
    fig, axes = figstyle.subplots(1, 1, w_mm=70, h_mm=55)
    figstyle.umap(adata, "total_counts", axes[0, 0])
    figstyle.save_figure(fig, plot_dir / "umap_by_total_counts", legend="UMAP coloured by transcripts per cell.")


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
    parser.add_argument(
        "--use-rep", default="X_pca", dest="use_rep",
        help="obsm key to build the neighbor graph on (default: X_pca). "
             "Pass X_pca_harmony after running 02b_harmony_integration.py.",
    )
    args = parser.parse_args()
    resolutions = [float(r) for r in args.resolutions.split(",")]

    adata = sc.read_h5ad(args.input)

    if args.use_rep not in adata.obsm:
        raise SystemExit(f"--use-rep '{args.use_rep}' not found in adata.obsm (available: {list(adata.obsm.keys())})")

    # A representation can have fewer dimensions than --n-pcs (e.g. 02 was run with --n-pcs 30).
    n_pcs = min(args.n_pcs, adata.obsm[args.use_rep].shape[1])
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep=args.use_rep)
    sc.tl.umap(adata)

    for res in resolutions:
        sc.tl.leiden(
            adata, resolution=res, key_added=f"leiden_r{res}",
            flavor="igraph", n_iterations=2,
        )
        n_clusters = adata.obs[f"leiden_r{res}"].nunique()
        print(f"resolution={res}: {n_clusters} clusters")

    figstyle.apply_style()
    save_cluster_plots(adata, resolutions, args.plot_dir)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

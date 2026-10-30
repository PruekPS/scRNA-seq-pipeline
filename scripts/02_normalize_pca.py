"""Step 2: Normalization, HVG selection, and PCA.

Implements the feature-selection/dimensionality-reduction stage documented
in docs/reference/design-decisions.md roadmap (item 2):
log-normalization -> highly-variable-gene selection on raw counts
(Seurat v3 flavor) -> scaling -> PCA. Batch integration and clustering are
separate downstream scripts (03+).

Input: the filtered .h5ad written by scripts/01_qc_filter.py.

Usage:
    python 02_normalize_pca.py --input filtered.h5ad --output pca.h5ad --plot-dir pca_plots/
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc

import figstyle


def save_pca_plots(adata, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    ratio = adata.uns["pca"]["variance_ratio"]
    fig, axes = figstyle.subplots(1, 1, w_mm=60, h_mm=45)
    ax = axes[0, 0]
    ax.plot(np.arange(1, len(ratio) + 1), ratio, "o-", ms=1.5, lw=0.5, color=figstyle.PALETTE[0])
    ax.set_yscale("log")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (fraction)")
    figstyle.save_figure(fig, plot_dir / "pca_elbow", legend="Variance explained by each principal component.")

    # Sanity check: PC1/PC2 shouldn't be dominated by sequencing depth alone.
    fig, axes = figstyle.subplots(1, 1, w_mm=70, h_mm=55)
    ax = axes[0, 0]
    sc.pl.pca(adata, color="total_counts", ax=ax, show=False, size=figstyle.point_size(adata.n_obs), title="",
              frameon=True)
    figstyle.rasterize_points(ax)
    figstyle.save_figure(fig, plot_dir / "pca_by_total_counts",
                         legend="First two principal components, coloured by transcripts per cell.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to filtered .h5ad from 01_qc_filter.py")
    parser.add_argument("--output", required=True, help="Path to write the PCA-ready .h5ad")
    parser.add_argument("--plot-dir", default="pca_plots", help="Directory for PCA diagnostic plots")
    parser.add_argument("--target-sum", type=float, default=1e4, help="Normalization target sum per cell")
    parser.add_argument("--n-hvg", type=int, default=2000, dest="n_hvg", help="Number of highly variable genes")
    parser.add_argument("--n-pcs", type=int, default=50, dest="n_pcs", help="Number of principal components")
    parser.add_argument("--exclude-hvg-genes", default=None, dest="exclude_hvg",
                        help="Optional CSV with a 'gene' column (and optionally an 'excluded' True/False column, as "
                             "written by 05d_flag_ambient_stress.py): these genes cannot be highly variable genes, so "
                             "they do not shape PCA or clustering. They stay in the data. Default: none.")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    # Keep raw counts around -- HVG selection (seurat_v3) needs them, and
    # downstream marker detection should run on log-normalized, unscaled
    # values rather than the scaled matrix used for PCA.
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata  # log-normalized, unscaled -- for marker detection later

    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=args.n_hvg, layer="counts", subset=False
    )

    if args.exclude_hvg:
        table = pd.read_csv(args.exclude_hvg, keep_default_na=False)
        if "gene" not in table.columns:
            raise SystemExit(f"--exclude-hvg-genes file needs a 'gene' column (has {list(table.columns)})")
        if "excluded" in table.columns:
            table = table[table["excluded"].astype(str).str.lower() == "true"]
        drop = adata.var_names.isin(table["gene"])
        print(f"Excluded from highly variable genes: {int((drop & adata.var['highly_variable']).sum())} "
              f"of {int(drop.sum())} listed genes found in the data")
        adata.var.loc[drop, "highly_variable"] = False

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, svd_solver="arpack")

    # Carry the PCA embedding back onto the full (non-subset) object so
    # downstream steps keep access to all genes, not just the HVG subset.
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]
    # Per-gene PC loadings are only defined on the HVG subset (adata_hvg.varm["PCs"]);
    # not carried onto the full object since downstream clustering only needs obsm["X_pca"].

    figstyle.apply_style()
    save_pca_plots(adata_hvg, args.plot_dir)

    print(f"Input:  {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")
    print(f"PCA: {args.n_pcs} components computed on the HVG subset")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

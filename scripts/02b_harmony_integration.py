"""Step 2b (optional): Harmony batch integration on the PCA embedding.

Motivated by a real finding, not a hypothetical: docs/reference/real-data-validation.md
documents one Leiden cluster (out of 20, resolution 0.6) that was 95%
single-donor in the real 50k-cell Parse subset run WITHOUT batch correction.
This step sits between 02_normalize_pca.py and 03_cluster.py and corrects the
PCA embedding for a batch covariate (default: the 'sample' column, i.e. donor,
already preserved through 00_load_parse_output.py for multi-donor Parse runs)
before clustering.

This is optional and additive: it does not overwrite obsm['X_pca'], it adds
obsm['X_pca_harmony'] alongside it. Skip this step entirely for single-sample
data -- Harmony has nothing to correct for in that case and running it anyway
just adds ~minutes of harmonypy iteration for no effect. Point 03_cluster.py
at the corrected embedding via --use-rep X_pca_harmony.

Validated on the real 907,696-cell full-scale run (2026-08-02): 24 donors,
converged in 6 iterations. Pin harmonypy==0.0.10 in environment.yml --
harmonypy 2.0.0 (conda's newest) returned a corrupted/ragged Z_corr array
at this scale and crashed the obsm assignment below; see environment.yml's
comment for the full story.

Input: the PCA-ready .h5ad written by scripts/02_normalize_pca.py.

Usage (local):
    python 02b_harmony_integration.py --input pca.h5ad --output harmony.h5ad \
        --batch-key sample --plot-dir harmony_plots/

Usage (HPC): see pipeline/run_downstream.slurm, which chains this step into
03_cluster.py --use-rep X_pca_harmony automatically.
"""
import argparse
import pathlib

import scanpy as sc
import scanpy.external as sce

import figstyle


def save_integration_plots(adata, batch_key, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Cheap before/after check on PC1 vs PC2 -- a full neighbors/UMAP re-run
    # belongs in 03_cluster.py (which this step feeds), not duplicated here.
    fig, axes = figstyle.subplots(1, 2, w_mm=figstyle.ONE_AND_HALF, h_mm=50)
    size = figstyle.point_size(adata.n_obs)
    sc.pl.embedding(adata, basis="pca", color=batch_key, ax=axes[0, 0], show=False, size=size,
                    title="Before Harmony", legend_loc=None)
    sc.pl.embedding(adata, basis="pca_harmony", color=batch_key, ax=axes[0, 1], show=False, size=size,
                    title="After Harmony", legend_loc="right margin", legend_fontsize=5)
    figstyle.rasterize_points(fig)
    figstyle.save_figure(fig, plot_dir / "pca_before_after_harmony",
                         legend=f"First two components before and after Harmony, coloured by {batch_key}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to PCA-ready .h5ad from 02_normalize_pca.py")
    parser.add_argument("--output", required=True, help="Path to write the Harmony-corrected .h5ad")
    parser.add_argument("--plot-dir", default="harmony_plots", help="Directory for before/after PCA plots")
    parser.add_argument("--batch-key", default="sample", dest="batch_key",
                         help="adata.obs column to correct for (default: 'sample', i.e. donor)")
    parser.add_argument("--theta", type=float, default=2.0,
                         help="Harmony diversity-clustering penalty (harmonypy default: 2.0). "
                              "Higher = more aggressive mixing across batches.")
    parser.add_argument("--max-iter-harmony", type=int, default=10, dest="max_iter_harmony")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.batch_key not in adata.obs:
        raise SystemExit(
            f"--batch-key '{args.batch_key}' not found in adata.obs "
            f"(available: {list(adata.obs.columns)}). For single-sample data, "
            f"skip this step and feed 02_pca.h5ad straight into 03_cluster.py."
        )
    n_batches = adata.obs[args.batch_key].nunique()
    if n_batches < 2:
        raise SystemExit(
            f"--batch-key '{args.batch_key}' has only {n_batches} unique value(s) -- "
            f"nothing for Harmony to correct. Skip this step for single-sample data."
        )
    print(f"Integrating {adata.n_obs} cells across {n_batches} '{args.batch_key}' batches...")

    try:
        import importlib.metadata
        hp_version = importlib.metadata.version("harmonypy")
        if not hp_version.startswith("0.0."):
            print(f"WARNING: harmonypy {hp_version} is installed, not the pinned 0.0.10 -- "
                  f"harmonypy 2.0.0 was found on a real 907k-cell run to return a ragged/"
                  f"object-dtype Z_corr array (shape (n_pcs,) instead of (n_cells, n_pcs)), "
                  f"which crashes the obsm assignment below with a confusing anndata shape "
                  f"error. If that happens, `conda install harmonypy=0.0.10` and retry.")
    except Exception:
        pass

    sce.pp.harmony_integrate(
        adata, key=args.batch_key,
        basis="X_pca", adjusted_basis="X_pca_harmony",
        theta=args.theta, max_iter_harmony=args.max_iter_harmony,
    )

    figstyle.apply_style()
    save_integration_plots(adata, args.batch_key, args.plot_dir)

    print(f"Wrote obsm['X_pca_harmony'] ({adata.obsm['X_pca_harmony'].shape[1]} dims).")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

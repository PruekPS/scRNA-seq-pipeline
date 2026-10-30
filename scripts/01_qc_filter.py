"""Step 1: QC filtering for a raw scRNA-seq count matrix.

Implements the QC stage documented in docs/reference/design-decisions.md (QC metrics table):
per-cell/per-gene metrics, adaptive median-absolute-deviation (MAD) outlier
detection (Amezquita et al., OSCA; Hong et al., Nat Commun 2022), and
doublet flagging via Scrublet.

Input is the raw-count .h5ad written by scripts/00_load_parse_output.py.

Usage:
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --plot-dir qc_plots/
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --mad 3 --species mouse
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --species canine --mt-prefix MT-
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --batch-key sample
"""
import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

import figstyle

# Mitochondrial gene-symbol prefix by species convention. Human/canine
# reference annotations (Ensembl/NCBI) both typically use an "MT-" prefix
# (HGNC-style); mouse uses lowercase "mt-" (MGI convention). Canine
# annotation conventions vary more by genome build/GTF source than
# human/mouse do -- override with --mt-prefix if this guess is wrong for
# your specific reference.
SPECIES_MT_PREFIX = {"human": "MT-", "mouse": "mt-", "canine": "MT-"}


def is_outlier(adata, metric, n_mads, upper_only=False):
    values = adata.obs[metric].to_numpy()
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    if mad == 0:
        return np.zeros(len(values), dtype=bool)
    high = values > med + n_mads * mad
    return high if upper_only else (values < med - n_mads * mad) | high


def compute_qc_metrics(adata, mt_prefix):
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    n_mt = int(adata.var["mt"].sum())
    if n_mt == 0:
        print(f"WARNING: 0 genes matched mt-prefix '{mt_prefix}' -- pct_counts_mt will be 0 for "
              f"every cell and the mito QC filter will do nothing. Check --species/--mt-prefix "
              f"against the gene symbols actually in this reference (adata.var_names).")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )
    return adata


def flag_outliers(adata, n_mads, hard_mt_ceiling):
    adata.obs["outlier_counts"] = is_outlier(adata, "log1p_total_counts", n_mads)
    adata.obs["outlier_genes"] = is_outlier(adata, "log1p_n_genes_by_counts", n_mads)
    # Upper tail only: a low mitochondrial share is not a sign of a damaged cell.
    adata.obs["outlier_mt"] = is_outlier(adata, "pct_counts_mt", n_mads, upper_only=True) | (
        adata.obs["pct_counts_mt"] > hard_mt_ceiling
    )
    adata.obs["qc_outlier"] = (
        adata.obs["outlier_counts"] | adata.obs["outlier_genes"] | adata.obs["outlier_mt"]
    )
    return adata


def flag_doublets(adata, batch_key=None):
    # Running Scrublet pooled across many donors/samples in one call both (a)
    # violates Scrublet's own methodology -- it simulates doublets from and
    # fits a kNN within the SAME library, so mixing 24 donors treats
    # cross-donor cell pairs as plausible same-library doublets, and its
    # per-library doublet-rate assumption breaks (Wolock, Lopez & Klein,
    # Cell Systems 2019 -- the Scrublet paper -- recommends running per
    # 10x lane/sample) -- and (b) is what OOM-killed the real full-scale
    # HPC run even at 400GB: Scrublet's internal simulated-doublet PCA/kNN
    # step holds effectively 2x the full 987,585-cell matrix in memory at
    # once. scanpy's batch_key loops Scrublet per batch internally, so peak
    # memory is bounded by the LARGEST single sample, not the whole dataset.
    if batch_key and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        print(f"Running Scrublet per-batch on '{batch_key}' ({adata.obs[batch_key].nunique()} batches)")
        sc.pp.scrublet(adata, batch_key=batch_key)
    else:
        sc.pp.scrublet(adata)
    # Scrublet's automatic threshold can fail silently: on a real 165k-cell run it called 2 doublets
    # while clusters of obvious multiplets remained. Say so rather than imply the doublets are gone.
    rate = float(adata.obs["predicted_doublet"].mean())
    if rate < 0.001 and adata.n_obs >= 5000:
        print(f"WARNING: Scrublet called only {rate:.2%} of cells doublets. Its automatic threshold may have "
              f"failed: do not treat doublets as removed. Check for clusters with very high counts or mixed "
              f"markers at step 04 (doublet_score is kept in adata.obs).")
    return adata


def save_qc_plots(adata, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # One violin per sample when there are several, so a bad sample is visible. Individual dots
    # are drawn only for small data: they add nothing but file size at 100k+ cells.
    labels = {"total_counts": "Transcripts per cell", "n_genes_by_counts": "Genes per cell",
              "pct_counts_mt": "Mitochondrial transcripts (%)"}
    by_sample = "sample" in adata.obs and 1 < adata.obs["sample"].nunique() <= 40
    fig, axes = figstyle.subplots(1, 3, w_mm=figstyle.DOUBLE, h_mm=50)
    for ax, (metric, label) in zip(axes[0], labels.items()):
        sc.pl.violin(adata, metric, groupby="sample" if by_sample else None, stripplot=adata.n_obs < 20000,
                     jitter=0.4, size=0.3, ax=ax, show=False, rotation=90 if by_sample else None, linewidth=0.3)
        ax.set_ylabel(label)
        ax.set_xlabel("")
        if not by_sample:
            ax.set_xticks([])
        figstyle.rasterize_points(ax)
    fig.tight_layout()
    figstyle.save_figure(fig, plot_dir / "qc_violin",
                         legend="Quality metrics after the hard floors, before outlier and doublet removal.")

    fig, ax = plt.subplots(figsize=(figstyle.SINGLE * figstyle.MM * 0.7, figstyle.SINGLE * figstyle.MM * 0.6))
    order = np.argsort(adata.obs["pct_counts_mt"].to_numpy())
    pts = ax.scatter(adata.obs["total_counts"].to_numpy()[order], adata.obs["n_genes_by_counts"].to_numpy()[order],
                     c=adata.obs["pct_counts_mt"].to_numpy()[order], s=figstyle.point_size(adata.n_obs), linewidths=0,
                     cmap=figstyle.SEQUENTIAL, rasterized=True)
    ax.set_xlabel("Transcripts per cell")
    ax.set_ylabel("Genes per cell")
    fig.colorbar(pts, ax=ax, shrink=0.6, label="Mitochondrial (%)").outline.set_linewidth(0.4)
    figstyle.save_figure(fig, plot_dir / "qc_scatter_counts_vs_genes",
                         legend="Each dot is one cell, before outlier and doublet removal.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to raw .h5ad count matrix")
    parser.add_argument("--output", required=True, help="Path to write filtered .h5ad")
    parser.add_argument("--plot-dir", default="qc_plots", help="Directory for QC plots")
    parser.add_argument(
        "--mad", type=float, default=5.0,
        help="MAD multiplier for adaptive outlier detection (OSCA examples typically use 3-5)",
    )
    parser.add_argument(
        "--mt-ceiling", type=float, default=8.0, dest="mt_ceiling",
        help="Hard mitochondrial-percent ceiling applied in addition to the adaptive threshold",
    )
    parser.add_argument("--species", default="human", choices=["human", "mouse", "canine"])
    parser.add_argument(
        "--mt-prefix", default=None, dest="mt_prefix",
        help="Override the species-inferred mitochondrial gene-symbol prefix "
             "(e.g. if your reference annotation's convention differs from the default)",
    )
    parser.add_argument("--min-genes", type=int, default=200, help="Hard floor: min genes detected per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Hard floor: min cells a gene must appear in")
    parser.add_argument(
        "--batch-key", default="sample", dest="batch_key",
        help="adata.obs column to run Scrublet separately per batch on (default: 'sample'). "
             "Falls back to a single pooled run automatically if this column is missing or "
             "has only one value (e.g. single-sample data) -- pass '' to force pooled mode.",
    )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()
    n_cells_raw, n_genes_raw = adata.shape

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    mt_prefix = args.mt_prefix or SPECIES_MT_PREFIX[args.species]
    adata = compute_qc_metrics(adata, mt_prefix=mt_prefix)
    adata = flag_outliers(adata, n_mads=args.mad, hard_mt_ceiling=args.mt_ceiling)
    adata = flag_doublets(adata, batch_key=args.batch_key or None)

    figstyle.apply_style()
    save_qc_plots(adata, args.plot_dir)

    n_before = adata.n_obs
    keep = ~adata.obs["qc_outlier"] & ~adata.obs["predicted_doublet"]
    adata_filtered = adata[keep].copy()
    n_after = adata_filtered.n_obs

    print(f"Loaded:            {n_cells_raw} cells x {n_genes_raw} genes")
    print(f"After hard floors: {n_before} cells x {adata.n_vars} genes")
    print(f"QC outliers:       {int(adata.obs['qc_outlier'].sum())} cells "
          f"(counts={int(adata.obs['outlier_counts'].sum())}, "
          f"genes={int(adata.obs['outlier_genes'].sum())}, "
          f"mito={int(adata.obs['outlier_mt'].sum())})")
    print(f"Predicted doublets: {int(adata.obs['predicted_doublet'].sum())} cells")
    print(f"Final:             {n_after} cells x {adata_filtered.n_vars} genes "
          f"({n_before - n_after} removed)")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata_filtered.write_h5ad(args.output)


if __name__ == "__main__":
    main()

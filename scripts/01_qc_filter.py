"""Step 1: QC filtering for a raw scRNA-seq count matrix.

Implements the QC stage documented in docs/pipeline_overview.html section 05:
per-cell/per-gene metrics, adaptive median-absolute-deviation (MAD) outlier
detection (Amezquita et al., OSCA; Hong et al., Nat Commun 2022), and
doublet flagging via Scrublet.

Input is an AnnData .h5ad with raw counts. Parse Trailmaker exports a
count matrix (h5ad/Seurat object download); split-pipe's DGE_filtered
output can be read into AnnData with `sc.read_mtx` + barcodes/genes and
saved as .h5ad before running this script.

Usage:
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --plot-dir qc_plots/
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --mad 3 --species mouse
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --species canine --mt-prefix MT-
    python 01_qc_filter.py --input raw.h5ad --output filtered.h5ad --batch-key sample
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

# Mitochondrial gene-symbol prefix by species convention. Human/canine
# reference annotations (Ensembl/NCBI) both typically use an "MT-" prefix
# (HGNC-style); mouse uses lowercase "mt-" (MGI convention). Canine
# annotation conventions vary more by genome build/GTF source than
# human/mouse do -- override with --mt-prefix if this guess is wrong for
# your specific reference.
SPECIES_MT_PREFIX = {"human": "MT-", "mouse": "mt-", "canine": "MT-"}


def is_outlier(adata, metric, n_mads):
    values = adata.obs[metric]
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    if mad == 0:
        return np.zeros(len(values), dtype=bool)
    return (values < med - n_mads * mad) | (values > med + n_mads * mad)


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
    adata.obs["outlier_mt"] = is_outlier(adata, "pct_counts_mt", n_mads) | (
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
    return adata


def save_qc_plots(adata, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, metric in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
        sc.pl.violin(adata, metric, jitter=0.4, ax=ax, show=False)
    fig.tight_layout()
    fig.savefig(plot_dir / "qc_violin.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", ax=ax, show=False)
    fig.tight_layout()
    fig.savefig(plot_dir / "qc_scatter_counts_vs_genes.png", dpi=150)
    plt.close(fig)


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

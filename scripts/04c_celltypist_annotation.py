"""Step 4c (optional): reference-based cell-type annotation with CellTypist.

Why alongside 04b: 04b scores clusters against a hand-picked marker
dictionary, which breaks down when cell *states* shift whole populations
away from textbook markers -- e.g. a strong, fast vaccine/interferon
response, or disease-skewed myeloid cells. On the first real Parse WT+TCR
PBMC run most 04b labels came out "Ambiguous (CD14+ Monocyte / HSPC)" for
exactly that reason. CellTypist (Dominguez Conde et al., Science 2022,
doi:10.1126/science.abl5197) instead applies logistic-regression models
trained on large curated immune atlases, scoring every gene per cell.

Two labels are written, both from the same model:
  celltypist_label    -- per-cell prediction
  celltypist_majority -- majority vote within each --over-clustering cluster
                         (smooths per-cell noise; the label to use by default)
plus celltypist_conf (per-cell confidence score of the predicted label).

Input expectations, enforced here rather than assumed: CellTypist models
are trained on log1p(counts normalized to 10,000 per cell) over the full
gene set, so this script rebuilds that from adata.layers['counts'] instead
of trusting adata.X (which by step 03 may be scaled / HVG-restricted).

Models are not downloaded here (no network calls inside analysis steps):
fetch them once beforehand, e.g.
    CELLTYPIST_FOLDER=/path/to/models python -c \
      "import celltypist; celltypist.models.download_models(model=['Immune_All_Low.pkl','Immune_All_High.pkl'])"
and pass --model /path/to/models/Immune_All_Low.pkl.

Usage:
    python 04c_celltypist_annotation.py --input 04b_annotated.h5ad --output 04c_celltypist.h5ad \
        --model /path/Immune_All_Low.pkl --over-clustering leiden_r0.6 --plot-dir celltypist_plots
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Clustered (or 04b-annotated) .h5ad with raw counts in a layer")
    parser.add_argument("--output", required=True, help="Path to write the CellTypist-annotated .h5ad")
    parser.add_argument("--model", required=True, help="Path to a downloaded CellTypist model .pkl")
    parser.add_argument("--over-clustering", default="leiden_r0.6", dest="over_clustering",
                        help="obs column used for majority voting (default: leiden_r0.6)")
    parser.add_argument("--counts-layer", default="counts", dest="counts_layer",
                        help="adata.layers key with raw counts (default: counts)")
    parser.add_argument("--prefix", default="celltypist",
                        help="Prefix for output obs columns (default: celltypist) -- change it to "
                             "store a second model's labels alongside the first")
    parser.add_argument("--plot-dir", default="celltypist_plots", dest="plot_dir")
    args = parser.parse_args()

    import celltypist  # imported here so --help works without it installed

    adata = sc.read_h5ad(args.input)
    if args.counts_layer not in adata.layers:
        raise SystemExit(f"--counts-layer '{args.counts_layer}' not in adata.layers ({list(adata.layers.keys())})")
    if args.over_clustering not in adata.obs:
        raise SystemExit(f"--over-clustering '{args.over_clustering}' not in adata.obs")

    # Rebuild exactly what CellTypist models expect: all genes, log1p(CP10k).
    query = sc.AnnData(X=adata.layers[args.counts_layer].copy(), obs=adata.obs[[args.over_clustering]].copy(),
                       var=adata.var[[]].copy())
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)

    model = celltypist.models.Model.load(args.model)
    overlap = np.isin(model.features, query.var_names).sum()
    print(f"Model {pathlib.Path(args.model).name}: {len(model.cell_types)} cell types, "
          f"{overlap}/{len(model.features)} model genes present in data")
    if overlap < 0.5 * len(model.features):
        raise SystemExit("Fewer than half the model's genes are in the data -- wrong species or gene naming?")

    pred = celltypist.annotate(query, model=model, majority_voting=True, over_clustering=args.over_clustering)
    labels = pred.predicted_labels
    p = args.prefix
    adata.obs[f"{p}_label"] = labels["predicted_labels"].astype("category").values
    adata.obs[f"{p}_majority"] = labels["majority_voting"].astype("category").values
    adata.obs[f"{p}_conf"] = pred.probability_matrix.max(axis=1).values

    print(f"\n{p}_majority counts:")
    print(adata.obs[f"{p}_majority"].value_counts().to_string())
    print(f"\nMedian per-cell confidence: {adata.obs[f'{p}_conf'].median():.3f}")

    plot_dir = pathlib.Path(args.plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    if "X_umap" in adata.obsm:
        for col in (f"{p}_majority", f"{p}_label"):
            sc.pl.umap(adata, color=col, show=False, frameon=False, size=1,
                       legend_loc="on data", legend_fontsize=5)
            plt.savefig(plot_dir / f"umap_{col}.png", dpi=150, bbox_inches="tight")
            plt.close("all")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

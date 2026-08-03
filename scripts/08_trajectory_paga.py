"""Step 8 (optional): PAGA coarse-grained trajectory graph + diffusion pseudotime.

The cheap, no-new-data first step toward trajectory analysis (see
docs/trajectory_cellrank_overview.html for the full design rationale). Runs
directly on an already-clustered, already-annotated object -- no raw FASTQ
reprocessing, no RNA velocity, nothing this pipeline doesn't already produce.

What it computes:
  - PAGA (Wolf et al., Genome Biology 2019): a coarse graph where each node
    is a cluster/cell-type and each edge weight reflects how connected two
    populations are in the neighbor graph -- a summary of which cell types
    plausibly transition into which others, not a proof of it.
  - Diffusion pseudotime (DPT, built into scanpy): an ordering of cells
    along inferred developmental distance from a user-specified root
    population (e.g. the most progenitor-like annotated cluster, such as
    "HSPC" from scripts/04b_better_annotation.py's marker dictionary).

This step's pseudotime feeds scripts/09_cellrank_fate_mapping.py as an
OPTIONAL second signal alongside CytoTRACE (--use-dpt-pseudotime there) --
CellRank's CytoTRACE kernel does not require this step to have been run.

NOT YET TEST-RUN as of this commit -- no real bone-marrow+TME dataset
exists yet to validate against (see docs/trajectory_cellrank_overview.html
for the experimental-design prerequisites). Validate against tests/data/
or a small real subset before trusting output.

Input: an annotated, clustered .h5ad (scripts/04b_better_annotation.py
output, ideally after Harmony integration for multi-sample data).

Usage:
    python 08_trajectory_paga.py --input 04b_annotated.h5ad --output 08_paga.h5ad \
        --plot-dir paga_plots/ --annotation-col cell_type_v2 --root-cluster HSPC
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def compute_root_cell(adata, annotation_col, root_cluster):
    mask = (adata.obs[annotation_col] == root_cluster).to_numpy()
    if not mask.any():
        raise SystemExit(
            f"--root-cluster '{root_cluster}' not found in adata.obs['{annotation_col}'] "
            f"(available: {sorted(adata.obs[annotation_col].unique())})"
        )
    # Heuristic, not ground truth: among cells in the root cluster, pick the
    # one furthest out along the first diffusion component as the single
    # representative "root" DPT needs -- typically a purer/more extreme
    # member of that population, less likely a boundary/doublet-like cell.
    # Override with --root-cell-index if a specific cell is known to be better.
    dc1 = adata.obsm["X_diffmap"][:, 1]  # component 0 is the trivial/constant one
    candidates = np.flatnonzero(mask)
    root_idx = candidates[np.argmax(np.abs(dc1[candidates]))]
    return int(root_idx)


def save_plots(adata, annotation_col, plot_dir):
    plot_dir = pathlib.Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.paga(adata, color=annotation_col, ax=ax, show=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(plot_dir / "paga_graph.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False, title="Diffusion pseudotime")
    sc.pl.umap(adata, color=annotation_col, ax=axes[1], show=False, title="Annotation", legend_fontsize=6)
    fig.tight_layout()
    fig.savefig(plot_dir / "umap_dpt_pseudotime.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to annotated, clustered .h5ad")
    parser.add_argument("--output", required=True, help="Path to write the .h5ad with PAGA + DPT results")
    parser.add_argument("--plot-dir", default="paga_plots", help="Directory for PAGA/pseudotime plots")
    parser.add_argument("--annotation-col", default="cell_type_v2", dest="annotation_col",
                         help="adata.obs column PAGA groups cells by (default: cell_type_v2, from 04b)")
    parser.add_argument("--root-cluster", required=True, dest="root_cluster",
                         help="Annotation value identifying the least-differentiated/progenitor "
                              "population (e.g. 'HSPC') -- diffusion pseudotime is measured "
                              "as distance FROM this population, so it materially determines "
                              "the direction of every downstream trajectory result.")
    parser.add_argument("--root-cell-index", type=int, default=None, dest="root_cell_index",
                         help="Override the auto-selected root cell with a specific adata row index")
    parser.add_argument("--n-dcs", type=int, default=15, dest="n_dcs", help="Number of diffusion components")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.annotation_col not in adata.obs:
        raise SystemExit(f"--annotation-col '{args.annotation_col}' not in adata.obs "
                          f"(available: {list(adata.obs.columns)})")
    if "neighbors" not in adata.uns:
        raise SystemExit("adata.uns['neighbors'] not found -- run scripts/03_cluster.py first "
                          "(PAGA and diffusion pseudotime both reuse its neighbor graph).")

    print(f"Computing PAGA on '{args.annotation_col}' ({adata.obs[args.annotation_col].nunique()} groups)...")
    sc.tl.paga(adata, groups=args.annotation_col)

    sc.tl.diffmap(adata, n_comps=args.n_dcs)

    root_idx = args.root_cell_index
    if root_idx is None:
        root_idx = compute_root_cell(adata, args.annotation_col, args.root_cluster)
    adata.uns["iroot"] = root_idx
    print(f"Root cell: index {root_idx} (cluster '{args.root_cluster}', "
          f"{'auto-selected' if args.root_cell_index is None else 'user-specified'})")

    sc.tl.dpt(adata, n_dcs=args.n_dcs)

    print(f"dpt_pseudotime range: [{adata.obs['dpt_pseudotime'].min():.4f}, "
          f"{adata.obs['dpt_pseudotime'].max():.4f}]")

    save_plots(adata, args.annotation_col, args.plot_dir)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()

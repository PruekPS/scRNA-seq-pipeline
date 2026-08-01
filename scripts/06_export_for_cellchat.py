"""Step 6: Export an annotated .h5ad to flat MTX/CSV for CellChat (R).

CellChat is an R package with no Python equivalent, so this is a one-way
bridge out of the Python side of the pipeline. Deliberately does NOT use
zellkonverter/reticulate to hand an AnnData object to R directly -- that
already failed once in this project (see 00b_seurat_to_dge_format.R docstring
and project memory): Bioconductor's `basilisk` tried to compile its own
Python 3.14 from source instead of reusing the `scrna` conda env, and the
build failed. Flat-file handoff sidesteps any live R<->Python bridge, same
pattern used for the Trailmaker Seurat->h5ad direction (00b -> 00).

Exports the LOG-NORMALIZED (not scaled, not raw-count) matrix from
adata.raw -- CellChat's createCellChat() expects normalized data, and its
tutorials work directly from Seurat's "data" slot (log-normalized), not
"counts" or "scale.data". Writes genes-as-rows/cells-as-columns (MTX
transposed from AnnData's cells x genes convention), matching what
Matrix::readMM() + CellChat expect on the R side -- see 07_cellchat_analysis.R.

Input: an annotated .h5ad with a cell-type/group column in .obs -- either
scripts/04_markers_annotation.py's 'cell_type_demo' or
scripts/04b_better_annotation.py's 'cell_type_v2'.

NOT YET TEST-RUN as of this commit (see project memory) -- validate against
tests/data/ (or a small real subset) before trusting output on a full-scale run.

Usage:
    python 06_export_for_cellchat.py --input annotated_v2.h5ad \
        --output-dir cellchat_input/ --group-col cell_type_v2
    python 06_export_for_cellchat.py --input annotated.h5ad \
        --output-dir cellchat_input/ --group-col cell_type_demo --exclude-ambiguous
"""
import argparse
import pathlib

import scanpy as sc
import scipy.io as sio


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to annotated .h5ad (from 04_ or 04b_)")
    parser.add_argument("--output-dir", required=True, dest="output_dir",
                         help="Directory to write normalized_matrix.mtx, genes.csv, cell_metadata.csv")
    parser.add_argument("--group-col", default="cell_type_v2", dest="group_col",
                         help="adata.obs column CellChat will group cells by (default: cell_type_v2, "
                              "from 04b_better_annotation.py; use cell_type_demo for 04_ output)")
    parser.add_argument("--exclude-ambiguous", action="store_true", dest="exclude_ambiguous",
                         help="Drop cells whose --group-col value starts with 'Ambiguous' "
                              "(only meaningful for cell_type_v2 -- see 04b_better_annotation.py)")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.group_col not in adata.obs:
        raise SystemExit(
            f"--group-col '{args.group_col}' not in adata.obs "
            f"(available: {list(adata.obs.columns)})"
        )

    if args.exclude_ambiguous:
        before = adata.n_obs
        adata = adata[~adata.obs[args.group_col].astype(str).str.startswith("Ambiguous")].copy()
        print(f"Excluded {before - adata.n_obs} ambiguous-labeled cells ({adata.n_obs} remain)")

    if adata.raw is None:
        raise SystemExit(
            "adata.raw is empty -- expected the log-normalized layer set by "
            "02_normalize_pca.py (adata.raw = adata before scaling). "
            "Re-run from 02_normalize_pca.py if this object skipped that step."
        )

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = adata.raw.to_adata()
    print(f"Exporting {raw.n_obs} cells x {raw.n_vars} genes (log-normalized, from adata.raw)")

    # Transpose: AnnData is cells x genes, CellChat/Matrix::readMM convention is genes x cells.
    sio.mmwrite(str(out_dir / "normalized_matrix.mtx"), raw.X.T.tocsr())

    raw.var_names.to_series().to_csv(out_dir / "genes.csv", index=False, header=["gene"])

    meta = adata.obs[[args.group_col]].copy()
    meta.columns = ["cell_type"]
    if "sample" in adata.obs.columns:
        meta["sample"] = adata.obs["sample"].values
    meta.index.name = "cell_barcode"
    meta.to_csv(out_dir / "cell_metadata.csv")

    print(f"Wrote {out_dir}/normalized_matrix.mtx, genes.csv, cell_metadata.csv")
    print(f"Group sizes ({args.group_col}):")
    print(meta["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()

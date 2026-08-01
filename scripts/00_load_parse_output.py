"""Step 0: Convert raw Parse split-pipe/Trailmaker output to .h5ad.

Every other script in this pipeline (01-04) expects an AnnData .h5ad with
raw counts. Parse's own output isn't that -- both Trailmaker's "Pipeline"
module and the split-pipe CLI produce a DGE_filtered/ folder per sample
with three files:

    count_matrix.mtx (or DGE.mtx in split-pipe < 1.1.0) -- sparse matrix,
        cells as ROWS, genes as COLUMNS (this is already AnnData's
        obs x var orientation -- no transpose needed, unlike 10x's mtx
        convention which is the other way around)
    all_genes.csv     -- gene_name, gene_id, genome per matrix column
    cell_metadata.csv -- barcode + per-cell metadata per matrix row

This script reads one or more such sample folders and writes a single
concatenated .h5ad, with a "sample" column in .obs so multi-sample runs
can be told apart downstream (e.g. for batch integration later).

Usage:
    python 00_load_parse_output.py --sample-dir DGE_filtered/ --output raw.h5ad
    python 00_load_parse_output.py \
        --sample-dir sample1/DGE_filtered/ sample1 \
        --sample-dir sample2/DGE_filtered/ sample2 \
        --output raw.h5ad
"""
import argparse
import pathlib

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse

BARCODE_COLUMN_CANDIDATES = ["bc_wells", "cell_barcode", "barcode", "cell_barcode_wells"]
GENE_NAME_COLUMN_CANDIDATES = ["gene_name", "gene_symbol", "name"]


def find_matrix_file(sample_dir: pathlib.Path) -> pathlib.Path:
    for fname in ("count_matrix.mtx", "DGE.mtx"):
        candidate = sample_dir / fname
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No count_matrix.mtx or DGE.mtx found in {sample_dir} -- "
        f"is this a split-pipe/Trailmaker DGE_filtered directory?"
    )


def pick_column(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"None of {candidates} found for {what}; columns present: {list(df.columns)}")


def load_one_sample(sample_dir: pathlib.Path, sample_name: str) -> ad.AnnData:
    matrix_path = find_matrix_file(sample_dir)
    genes_path = sample_dir / "all_genes.csv"
    cells_path = sample_dir / "cell_metadata.csv"
    if not genes_path.exists():
        raise FileNotFoundError(f"Missing {genes_path}")
    if not cells_path.exists():
        raise FileNotFoundError(f"Missing {cells_path}")

    X = scipy.io.mmread(matrix_path).tocsr()

    genes = pd.read_csv(genes_path)
    gene_name_col = pick_column(genes, GENE_NAME_COLUMN_CANDIDATES, "gene name")

    cells = pd.read_csv(cells_path)
    barcode_col = pick_column(cells, BARCODE_COLUMN_CANDIDATES, "cell barcode")
    cells = cells.set_index(barcode_col)
    cells.index = cells.index.astype(str)
    cells.index.name = None

    if X.shape[0] != len(cells) or X.shape[1] != len(genes):
        raise ValueError(
            f"Matrix shape {X.shape} doesn't match cell_metadata rows ({len(cells)}) "
            f"/ all_genes rows ({len(genes)}) in {sample_dir} -- wrong files paired together?"
        )

    adata = ad.AnnData(X=X, obs=cells, var=genes.set_index(gene_name_col))
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_name
    return adata


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--sample-dir", nargs="+", action="append", metavar=("PATH", "SAMPLE_NAME"),
        required=True,
        help="One DGE_filtered directory + a sample name, repeatable for multi-sample runs. "
             "If SAMPLE_NAME is omitted, the parent directory name is used.",
    )
    parser.add_argument("--output", required=True, help="Path to write the combined raw .h5ad")
    args = parser.parse_args()

    samples = []
    for entry in args.sample_dir:
        sample_dir = pathlib.Path(entry[0])
        sample_name = entry[1] if len(entry) > 1 else sample_dir.parent.name
        print(f"Loading {sample_dir} as sample '{sample_name}'...")
        samples.append(load_one_sample(sample_dir, sample_name))

    if len(samples) == 1:
        combined = samples[0]
    else:
        combined = ad.concat(samples, join="outer", index_unique="-")

    print(f"Combined: {combined.n_obs} cells x {combined.n_vars} genes "
          f"across {combined.obs['sample'].nunique()} sample(s)")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(args.output)


if __name__ == "__main__":
    main()

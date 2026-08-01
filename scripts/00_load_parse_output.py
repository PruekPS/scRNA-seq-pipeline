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


def find_one(sample_dir: pathlib.Path, patterns: list[str], what: str) -> pathlib.Path:
    """Real Parse output filenames aren't perfectly consistent across
    versions/datasets -- e.g. published datasets use dataset-suffixed names
    like DGE_1M_PBMC.mtx rather than the generic DGE.mtx the docs describe.
    Glob-match instead of assuming an exact filename."""
    for pattern in patterns:
        matches = sorted(sample_dir.glob(pattern))
        if matches:
            if len(matches) > 1:
                print(f"Warning: multiple {what} candidates matched {pattern} in {sample_dir}, using {matches[0].name}")
            return matches[0]
    raise FileNotFoundError(
        f"No {what} found in {sample_dir} (tried patterns {patterns}) -- "
        f"is this a split-pipe/Trailmaker DGE_filtered directory?"
    )


def find_matrix_file(sample_dir: pathlib.Path) -> pathlib.Path:
    return find_one(sample_dir, ["count_matrix*.mtx", "DGE*.mtx", "*.mtx"], "count matrix")


def find_genes_file(sample_dir: pathlib.Path) -> pathlib.Path:
    return find_one(sample_dir, ["all_genes*.csv", "*genes*.csv"], "gene metadata")


def find_cells_file(sample_dir: pathlib.Path) -> pathlib.Path:
    return find_one(sample_dir, ["cell_metadata*.csv", "*cell*metadata*.csv"], "cell metadata")


def pick_column(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"None of {candidates} found for {what}; columns present: {list(df.columns)}")


def load_one_sample(sample_dir: pathlib.Path, sample_name: str) -> ad.AnnData:
    matrix_path = find_matrix_file(sample_dir)
    genes_path = find_genes_file(sample_dir)
    cells_path = find_cells_file(sample_dir)
    print(f"  matrix: {matrix_path.name}, genes: {genes_path.name}, cells: {cells_path.name}")

    X = scipy.io.mmread(matrix_path).tocsr()

    genes = pd.read_csv(genes_path)
    gene_name_col = pick_column(genes, GENE_NAME_COLUMN_CANDIDATES, "gene name")
    if "gene_id" in genes.columns:
        n_missing = genes[gene_name_col].isna().sum()
        if n_missing:
            # Real gene annotation files routinely have no symbol for a chunk of
            # Ensembl IDs (lncRNAs, pseudogenes, etc.) -- e.g. ~1/3 of genes in
            # the Parse 1M PBMC dataset. Falling back to gene_id keeps every
            # gene instead of silently dropping or crashing on the NaN index.
            print(f"  {n_missing}/{len(genes)} genes have no {gene_name_col} -- falling back to gene_id for those")
            genes[gene_name_col] = genes[gene_name_col].fillna(genes["gene_id"])

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
    if "sample" in adata.obs.columns:
        # cell_metadata.csv can already carry real per-cell sample/donor IDs
        # (e.g. multiple donors pooled into one DGE_filtered folder) -- don't
        # clobber that with the folder-level --sample-dir label.
        print(f"  cell_metadata already has a 'sample' column "
              f"({adata.obs['sample'].nunique()} unique values) -- keeping it as-is.")
    else:
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

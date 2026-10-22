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

split-pipe `--mode comb` output: point --sample-dir at either
`<comb>/all-sample/DGE_filtered` (one folder, cell_metadata already carries
each cell's `sample`) OR at the per-sample `<comb>/<sample>/DGE_filtered`
folders -- never both, or every cell is loaded twice. Cell names are the
`bc_wells` barcodes (`NN_NN_NN__s<N>`, where s<N> is the sublibrary's
position in the comb list); they are kept as-is so TCR tables (script 12)
can be matched to them, and the s<N> tag is copied into obs["sublibrary"]
-- the real technical batch in split-pool data (every sublibrary contains
every sample), usable as Harmony's --batch-key.

--sample-sheet CSV: split-pipe only knows sample names, so study design
columns (patient, timepoint, tissue, condition, ...) come from a small CSV
with a `sample` column plus any others; they are joined onto every cell.
Sample names must match the split-pipe sample names exactly; samples with
cells but no sheet row stop the run. The name column is `sample` unless
--sample-sheet-key says otherwise (e.g. `sample_id`).

--exclude-sublibraries s1 ...: drop whole split-pipe sublibraries (by their
__sN tag) when one carries a library-level artifact.

--min-counts / --min-genes (optional): one cell cutoff for every sample.
split-pipe calls cells with a separate knee cutoff per sample and per
sublibrary. When a sample's count distribution has no clean knee (e.g. a
high-UMI tail of activated cells), that cutoff can land far above the main
population and DGE_filtered keeps only a few hundred cells. The fix is to
load DGE_unfiltered (same three files, every barcode) with these options,
so every sample gets the same threshold. Check whether you need it with
00a_check_cell_calls.py; see docs/03-your-data.md. Off (0) by default, so
DGE_filtered input is loaded unchanged.

Usage:
    python 00_load_parse_output.py --sample-dir DGE_filtered/ --output raw.h5ad
    python 00_load_parse_output.py \
        --sample-dir sample1/DGE_filtered/ sample1 \
        --sample-dir sample2/DGE_filtered/ sample2 \
        --output raw.h5ad
    python 00_load_parse_output.py --sample-dir comb/all-sample/DGE_filtered \
        --sample-sheet samples.csv --output raw.h5ad
    # split-pipe's per-sample cell calls unreliable: one cutoff for all samples
    python 00_load_parse_output.py --sample-dir comb/all-sample/DGE_unfiltered \
        --min-counts 500 --min-genes 300 --sample-sheet samples.csv --output raw.h5ad
"""
import argparse
import pathlib
import sys

import anndata as ad
import numpy as np
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


def load_one_sample(sample_dir: pathlib.Path, sample_name: str,
                    min_counts: int = 0, min_genes: int = 0) -> ad.AnnData:
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

    if min_counts or min_genes:
        # Computed from the matrix itself, not cell_metadata's tscp_count/gene_count,
        # so the cutoff means the same thing whichever pipeline wrote the folder.
        counts = np.asarray(X.sum(axis=1)).ravel()
        genes_per_cell = np.diff(X.indptr)
        keep = (counts >= min_counts) & (genes_per_cell >= min_genes)
        print(f"  cell cutoff (>= {min_counts} counts, >= {min_genes} genes): "
              f"kept {int(keep.sum()):,} of {len(keep):,} barcodes")
        X, cells = X[keep], cells.iloc[keep]
    elif "unfiltered" in sample_dir.name.lower():
        print(f"  WARNING: {sample_dir.name} holds every barcode, including empty ones, and no "
              f"--min-counts/--min-genes was given -- set a cutoff (see 00a_check_cell_calls.py).")

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
    parser.add_argument("--sample-sheet", default=None,
                        help="Optional CSV with a sample-name column + design columns (patient, timepoint, ...) "
                             "joined onto every cell by sample name")
    parser.add_argument("--sample-sheet-key", default="sample",
                        help="Column of --sample-sheet holding the split-pipe sample name (default 'sample')")
    parser.add_argument("--sample-metadata", default=None, help=argparse.SUPPRESS)  # older alias of --sample-sheet
    parser.add_argument("--sample-metadata-key", default="sample_id", help=argparse.SUPPRESS)  # alias, old default
    parser.add_argument("--exclude-sublibraries", nargs="+", default=None, metavar="sN",
                        help="Drop cells from these split-pipe sublibraries (e.g. s1), numbered by the __sN "
                             "suffix of --mode comb barcodes -- for a sublibrary with a library-level artifact. "
                             "Order follows the sublibrary list given to split-pipe comb.")
    parser.add_argument("--min-counts", type=int, default=0,
                        help="Keep only barcodes with at least this many counts (transcripts/UMIs). "
                             "Use with DGE_unfiltered to apply one cutoff to every sample instead of "
                             "split-pipe's per-sample cell calls. Default 0 = no cutoff.")
    parser.add_argument("--min-genes", type=int, default=0,
                        help="Keep only barcodes with at least this many detected genes. Default 0.")
    args = parser.parse_args()

    samples = []
    for entry in args.sample_dir:
        sample_dir = pathlib.Path(entry[0])
        sample_name = entry[1] if len(entry) > 1 else sample_dir.parent.name
        print(f"Loading {sample_dir} as sample '{sample_name}'...")
        samples.append(load_one_sample(sample_dir, sample_name, args.min_counts, args.min_genes))

    if len(samples) == 1:
        combined = samples[0]
    else:
        # split-pipe barcodes (with their __s<N> tag) are already unique across samples; keep them
        # untouched so TCR/other per-cell tables can still be matched. Only disambiguate on a clash.
        names = pd.Index([n for a in samples for n in a.obs_names])
        if names.is_unique:
            combined = ad.concat(samples, join="outer")
        else:
            print("  WARNING: cell names repeat across --sample-dir folders (same folder twice, or "
                  "all-sample plus per-sample folders?) -- appending '-<i>' to keep them apart; "
                  "TCR tables will no longer match these names.")
            combined = ad.concat(samples, join="outer", index_unique="-")

    sub = combined.obs_names.to_series().str.extract(r"__(s\d+)$")[0]
    if sub.notna().all():
        combined.obs["sublibrary"] = pd.Categorical(sub.values)
        print(f"  sublibrary tag found in cell names: {combined.obs['sublibrary'].nunique()} sublibraries")

    if args.exclude_sublibraries:
        if "sublibrary" not in combined.obs:
            sys.exit("--exclude-sublibraries given, but cell names carry no __sN sublibrary tag "
                     "(expects split-pipe --mode comb output)")
        unknown = sorted(set(args.exclude_sublibraries) - set(combined.obs["sublibrary"].astype(str)))
        if unknown:
            sys.exit(f"--exclude-sublibraries {unknown} not present; have "
                     f"{sorted(combined.obs['sublibrary'].astype(str).unique())}")
        drop = combined.obs["sublibrary"].astype(str).isin(args.exclude_sublibraries).values
        print(f"  excluding sublibraries {args.exclude_sublibraries}: {int(drop.sum()):,} cells removed")
        combined = combined[~drop].copy()
        combined.obs["sublibrary"] = combined.obs["sublibrary"].cat.remove_unused_categories()

    sheet_path, key = args.sample_sheet, args.sample_sheet_key
    if args.sample_metadata:
        if sheet_path:
            sys.exit("Give --sample-sheet or its older alias --sample-metadata, not both")
        sheet_path, key = args.sample_metadata, args.sample_metadata_key
    if sheet_path:
        sheet = pd.read_csv(sheet_path, dtype=str)
        if key not in sheet.columns:
            sys.exit(f"--sample-sheet needs a '{key}' column (set another with --sample-sheet-key); "
                     f"has {list(sheet.columns)}")
        if sheet[key].duplicated().any():
            sys.exit(f"--sample-sheet lists a sample twice: {sorted(sheet.loc[sheet[key].duplicated(), key])}")
        have = set(combined.obs["sample"].astype(str))
        missing = sorted(have - set(sheet[key]))
        if missing:
            sys.exit(f"Samples with cells but no --sample-sheet row: {missing}. Names must match split-pipe's "
                     f"sample names exactly (sheet has {sorted(sheet[key])}).")
        unused = sorted(set(sheet[key]) - have)
        if unused:
            print(f"  note: sheet rows with no cells (ignored): {unused}")
        clash = [c for c in sheet.columns if c != key and c in combined.obs.columns]
        if clash:
            sys.exit(f"--sample-sheet columns already exist in the cell metadata: {clash}; rename them")
        lookup = sheet.set_index(key)
        for col in lookup.columns:
            combined.obs[col] = pd.Categorical(combined.obs["sample"].astype(str).map(lookup[col]).values)
        print(f"  joined sample sheet columns: {list(lookup.columns)}")

    if args.min_counts or args.min_genes:
        combined.uns["load_cell_cutoff"] = {"min_counts": args.min_counts, "min_genes": args.min_genes}
        per_sample = combined.obs["sample"].astype(str).value_counts().sort_index()
        print("  cells per sample after the cutoff: "
              + ", ".join(f"{s} {n:,}" for s, n in per_sample.items()))

    print(f"Combined: {combined.n_obs} cells x {combined.n_vars} genes "
          f"across {combined.obs['sample'].nunique()} sample(s)")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(args.output)


if __name__ == "__main__":
    main()

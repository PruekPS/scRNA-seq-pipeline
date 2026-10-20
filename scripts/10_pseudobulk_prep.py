"""Step 10 (optional): Aggregate cells into pseudobulk pseudo-samples for DE.

Per-cell differential expression tests (e.g. naively running a t-test or
Wilcoxon test treating every cell as an independent sample) inflate false
discovery rate, because cells from the same biological replicate are not
independent samples -- pseudobulk aggregation (sum counts per
sample x cell-type x condition group, then run a bulk DE method) controls
this properly (Murphy & Skene, Nat Commun 2022, doi:10.1038/s41467-022-35519-4;
see docs/reference/design-decisions.md roadmap, item 6). This script does the
aggregation; scripts/11_pseudobulk_deseq2.R runs DESeq2 on the result.

Dataset-agnostic by design: --condition-col takes any adata.obs column
name at run time (e.g. "disease_status", "treatment", "timepoint") --
nothing here is specific to any one dataset or comparison.

Different experiments need different DE designs (unpaired groups, paired
before/after within a patient, a known processing batch, ...). Any
SAMPLE-LEVEL column named in --covariate-cols (e.g. "patient", "batch") is
carried into pseudobulk_metadata.csv so scripts/11 can put it in the DESeq2
design (--covariates there). A covariate must be constant within each
sample -- per-cell quantities (e.g. a sublibrary tag that mixes cells of
one sample) can't be a pseudobulk covariate and are rejected.

Each pseudo-sample = one (sample, cell_type, condition) combination,
summed from raw counts (not log-normalized data -- DESeq2 models raw
counts directly, that's the whole point of using it over a per-cell test).
Pseudo-samples aggregated from too few cells are dropped
(--min-cells-per-group) since a "pseudobulk" of a handful of cells isn't
meaningfully different from single-cell noise.

NOT YET TEST-RUN as of this commit -- dataset-agnostic by design, so
there's nothing to validate against until a real condition comparison
exists. Validate against tests/data/ (add a synthetic condition column
if needed) or a small real subset before trusting output.

Input: an annotated .h5ad with raw counts in a layer (scripts/02_normalize_pca.py
onward all preserve adata.layers["counts"]) and a cell-type annotation column
(scripts/04_markers_annotation.py's cell_type_demo or 04b's cell_type_v2).

Usage:
    python 10_pseudobulk_prep.py --input 04b_annotated.h5ad --output-dir pseudobulk/ \
        --condition-col disease_status --sample-col sample --annotation-col cell_type_v2
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp


def aggregate_pseudobulk(adata, sample_col, annotation_col, condition_col, counts_layer, min_cells,
                         covariate_cols=()):
    counts = adata.layers[counts_layer]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)

    group_cols = [sample_col, annotation_col, condition_col]
    groups = adata.obs.groupby(group_cols, observed=True).indices

    sample_covs = (adata.obs.groupby(sample_col, observed=True)[list(covariate_cols)].first()
                   if covariate_cols else None)

    pseudo_counts = []
    meta_rows = []
    for (sample_val, celltype_val, condition_val), cell_idx in groups.items():
        n_cells = len(cell_idx)
        if n_cells < min_cells:
            continue
        summed = np.asarray(counts[cell_idx, :].sum(axis=0)).ravel()
        pseudo_id = f"{sample_val}__{celltype_val}__{condition_val}".replace(" ", "_")
        pseudo_counts.append(summed)
        meta_rows.append({
            "pseudo_sample_id": pseudo_id,
            "sample": sample_val,
            "cell_type": celltype_val,
            "condition": condition_val,
            "n_cells_aggregated": n_cells,
            **({c: sample_covs.loc[sample_val, c] for c in covariate_cols} if covariate_cols else {}),
        })

    if not meta_rows:
        raise SystemExit(
            f"No (sample, cell_type, condition) group had >= {min_cells} cells -- "
            f"nothing to aggregate. Lower --min-cells-per-group or check that "
            f"--condition-col/--sample-col/--annotation-col are actually populated."
        )

    pseudo_matrix = np.vstack(pseudo_counts)  # pseudo-samples x genes
    meta_df = pd.DataFrame(meta_rows)
    return pseudo_matrix, meta_df


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to annotated .h5ad")
    parser.add_argument("--output-dir", required=True, dest="output_dir",
                         help="Directory to write pseudobulk_counts.mtx, genes.csv, pseudobulk_metadata.csv")
    parser.add_argument("--condition-col", required=True, dest="condition_col",
                         help="adata.obs column defining the comparison groups (e.g. 'disease_status', "
                              "'treatment') -- not tied to any specific dataset, supply whatever column "
                              "your real data actually has")
    parser.add_argument("--sample-col", default="sample", dest="sample_col",
                         help="adata.obs column identifying biological replicates (default: 'sample')")
    parser.add_argument("--annotation-col", default="cell_type_v2", dest="annotation_col",
                         help="adata.obs column with cell-type labels (default: cell_type_v2, from 04b)")
    parser.add_argument("--counts-layer", default="counts", dest="counts_layer",
                         help="adata.layers key holding raw (unnormalized) counts (default: 'counts')")
    parser.add_argument("--min-cells-per-group", type=int, default=10, dest="min_cells",
                         help="Minimum cells required to keep a (sample, cell_type, condition) "
                              "pseudo-sample (default: 10 -- a pseudobulk of fewer cells is unreliable)")
    parser.add_argument("--covariate-cols", default="", dest="covariate_cols",
                         help="Comma-separated SAMPLE-LEVEL adata.obs columns to carry into the metadata for "
                              "use as DESeq2 covariates in script 11 (e.g. 'patient' for a paired design, "
                              "'batch'). Default: none.")
    args = parser.parse_args()
    covariate_cols = [c.strip() for c in args.covariate_cols.split(",") if c.strip()]

    adata = sc.read_h5ad(args.input)

    for col in (args.sample_col, args.annotation_col, args.condition_col):
        if col not in adata.obs:
            raise SystemExit(f"'{col}' not in adata.obs (available: {list(adata.obs.columns)})")
    for col in covariate_cols:
        if col not in adata.obs:
            raise SystemExit(f"--covariate-cols '{col}' not in adata.obs (available: {list(adata.obs.columns)})")
        if col in (args.sample_col, args.annotation_col, args.condition_col):
            raise SystemExit(f"--covariate-cols '{col}' is already the sample/annotation/condition column")
        varies = adata.obs.groupby(args.sample_col, observed=True)[col].nunique()
        if (varies > 1).any():
            raise SystemExit(f"--covariate-cols '{col}' is not constant within samples "
                              f"(e.g. {list(varies[varies > 1].index[:3])}) -- a pseudobulk covariate must "
                              "describe the whole sample, not individual cells.")
    if args.counts_layer not in adata.layers:
        raise SystemExit(f"--counts-layer '{args.counts_layer}' not in adata.layers "
                          f"(available: {list(adata.layers.keys())})")

    n_conditions = adata.obs[args.condition_col].nunique()
    if n_conditions < 2:
        raise SystemExit(f"--condition-col '{args.condition_col}' has only {n_conditions} unique value(s) -- "
                          f"need at least 2 to compare (found: {adata.obs[args.condition_col].unique().tolist()})")

    print(f"Aggregating by ({args.sample_col}, {args.annotation_col}, {args.condition_col}), "
          f"min {args.min_cells} cells/group...")
    pseudo_matrix, meta_df = aggregate_pseudobulk(
        adata, args.sample_col, args.annotation_col, args.condition_col, args.counts_layer, args.min_cells,
        covariate_cols,
    )

    print(f"Kept {len(meta_df)} pseudo-samples (dropped groups with < {args.min_cells} cells)")
    print("Pseudo-samples per (cell_type, condition):")
    print(meta_df.groupby(["cell_type", "condition"], observed=True).size().to_string())

    out_dir = pathlib.Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # genes x pseudo-samples, matching the R/DESeq2-side convention (same
    # transpose convention as 06_export_for_cellchat.py for the same reason).
    sio.mmwrite(str(out_dir / "pseudobulk_counts.mtx"), sp.csr_matrix(pseudo_matrix.T))
    adata.var_names.to_series().to_csv(out_dir / "genes.csv", index=False, header=["gene"])
    meta_df.to_csv(out_dir / "pseudobulk_metadata.csv", index=False)

    print(f"Wrote {out_dir}/pseudobulk_counts.mtx, genes.csv, pseudobulk_metadata.csv")


if __name__ == "__main__":
    main()

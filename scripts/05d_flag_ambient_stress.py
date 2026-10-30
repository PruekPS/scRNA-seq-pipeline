"""Step 5d (optional): find likely ambient-RNA genes, and flag ambient and stress genes in DE tables.

Two separate problems can look like biology:

Ambient RNA. Damaged cells release RNA into the suspension, and it ends up in other
droplets or wells, so a T cell can appear to express a macrophage or lung gene. This is
likely a bigger issue for tissue, nuclei and FFPE samples. If you have a cluster of
low-quality cells (low counts, no clear identity), its top genes are a good estimate of
the ambient profile. A gene is called ambient when it is among the top --top-n genes of
that cluster and detected in under --max-detection of a clean reference population whose
biology should not include it.

Stress. Heat-shock, immediate-early and hypoxia genes rise with tissue handling and time
on ice (van den Brink et al., Nat Methods 2017, doi:10.1038/nmeth.4437). If samples were
collected on different days, they can dominate a DE result. Such genes are flagged by
name, for human or mouse symbols.

Nothing is removed. Part 1 writes a list you can give to 02_normalize_pca.py
--exclude-hvg-genes, so ambient genes do not shape clustering. Part 2 adds True/False
columns to a DE table, and counts the hits left without them.

Part 1 needs --input, --low-quality-col, --low-quality-values, --reference-col and
--reference-values. Part 2 needs --de-table. Run either or both.

Usage:
    python 05d_flag_ambient_stress.py --output-dir ambient/ --input 04b_annotated.h5ad \
        --low-quality-col leiden_r0.6 --low-quality-values 2 \
        --reference-col leiden_r0.6 --reference-values 4 8 9
    python 05d_flag_ambient_stress.py --output-dir ambient/ \
        --de-table pseudobulk_de/de_all_celltypes_combined.csv --ambient-genes ambient/ambient_genes.csv
"""
import argparse
import pathlib
import re

import numpy as np
import pandas as pd

# Matched without regard to case, so the same pattern covers human and mouse symbols.
STRESS_PATTERN = re.compile(
    r"^(HSP[A-Z0-9]+|DNAJ[ABC]\d+|BAG3|FOS[BL]?\d?|JUN[BD]?|EGR\d|IER\d|ATF3|NR4A\d|DUSP1|ZFP36(L\d)?|KLF[246]|SOCS3|"
    r"PPP1R15A|CEBPD|MT[12]|HIF1A|TXNIP|DDIT4|GADD45[ABG]|BTG2|CYR61|CCN1|UBC)$", re.IGNORECASE)


def is_stress_gene(gene) -> bool:
    return bool(STRESS_PATTERN.match(str(gene)))


def ambient_table(adata, low_quality, reference, top_n, max_detection, counts_layer="counts"):
    if counts_layer not in adata.layers:
        raise SystemExit(f"'{counts_layer}' not in adata.layers ({list(adata.layers.keys())})")
    if not low_quality.any() or not reference.any():
        raise SystemExit(f"Selected {int(low_quality.sum())} low-quality and {int(reference.sum())} reference cells; "
                         "both must be above 0. Check the column names and values.")
    counts = adata.layers[counts_layer].tocsr()
    library = np.maximum(np.asarray(counts[low_quality].sum(axis=1)).ravel(), 1)
    mean_cp10k = np.asarray(counts[low_quality].multiply(1e4 / library[:, None]).mean(axis=0)).ravel()
    detected = np.asarray((counts[reference] > 0).mean(axis=0)).ravel()
    top = np.argsort(mean_cp10k)[::-1][:top_n]
    return pd.DataFrame({"gene": adata.var_names[top], "mean_cp10k_low_quality": mean_cp10k[top].round(2),
                         "fraction_detected_in_reference": detected[top].round(4),
                         "excluded": detected[top] < max_detection})


def flag_de_table(de, ambient_genes):
    de = de.copy()
    de["stress_gene"] = de["gene"].map(is_stress_gene)
    de["possible_ambient"] = de["gene"].astype(str).isin(ambient_genes)
    return de


def summarise_hits(de, padj, lfc):
    sig = de[(de["padj"] < padj) & (de["log2FoldChange"].abs() > lfc)]
    rows = []
    for cell_type, s in sig.groupby("cell_type", observed=True):
        clean = s[~s["possible_ambient"] & ~s["stress_gene"]]
        rows.append({"cell_type": cell_type, "n_significant": len(s), "stress": int(s["stress_gene"].sum()),
                     "possible_ambient": int(s["possible_ambient"].sum()),
                     "clean_up": int((clean["log2FoldChange"] > 0).sum()),
                     "clean_down": int((clean["log2FoldChange"] < 0).sum())})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True, help="Where to write results (created if needed)")
    parser.add_argument("--input", default=None, help="Part 1: annotated .h5ad with raw counts in a layer")
    parser.add_argument("--low-quality-col", default=None, help="Part 1: obs column marking the low-quality cells")
    parser.add_argument("--low-quality-values", nargs="+", default=None, help="Part 1: its value(s), e.g. a cluster")
    parser.add_argument("--reference-col", default=None, help="Part 1: obs column marking the clean reference cells")
    parser.add_argument("--reference-values", nargs="+", default=None, help="Part 1: its value(s)")
    parser.add_argument("--top-n", type=int, default=150, help="Part 1: top low-quality genes considered (default 150)")
    parser.add_argument("--max-detection", type=float, default=0.10,
                        help="Part 1: called ambient if detected in under this fraction of reference cells (default 0.10)")
    parser.add_argument("--de-table", default=None, help="Part 2: de_all_celltypes_combined.csv from step 11")
    parser.add_argument("--ambient-genes", default=None,
                        help="Part 2: ambient_genes.csv from part 1 (default: the one just written, or none)")
    parser.add_argument("--padj", type=float, default=0.05, help="Part 2: significance cutoff (default 0.05)")
    parser.add_argument("--lfc", type=float, default=1.0, help="Part 2: |log2 fold change| cutoff (default 1)")
    args = parser.parse_args()
    if not (args.input or args.de_table):
        raise SystemExit("Nothing to do: give --input (part 1) and/or --de-table (part 2)")

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ambient_csv = args.ambient_genes

    if args.input:
        needed = (args.low_quality_col, args.low_quality_values, args.reference_col, args.reference_values)
        if not all(needed):
            raise SystemExit("Part 1 needs --low-quality-col, --low-quality-values, --reference-col, --reference-values")
        import anndata as ad
        adata = ad.read_h5ad(args.input)
        for col in (args.low_quality_col, args.reference_col):
            if col not in adata.obs:
                raise SystemExit(f"'{col}' not in adata.obs (available: {list(adata.obs.columns)})")
        low = adata.obs[args.low_quality_col].astype(str).isin(args.low_quality_values).to_numpy()
        reference = adata.obs[args.reference_col].astype(str).isin(args.reference_values).to_numpy()
        table = ambient_table(adata, low, reference, args.top_n, args.max_detection)
        ambient_csv = out / "ambient_genes.csv"
        table.to_csv(ambient_csv, index=False)
        print(f"{int(table['excluded'].sum())} likely ambient genes of the top {len(table)}, from {int(low.sum()):,} "
              f"low-quality and {int(reference.sum()):,} reference cells -> {ambient_csv}")
        print("  e.g. " + ", ".join(table.loc[table["excluded"], "gene"].head(15)))

    if args.de_table:
        de = pd.read_csv(args.de_table, keep_default_na=False,
                         na_values={"padj": ["NA", ""], "log2FoldChange": ["NA", ""]})
        ambient_genes = set()
        if ambient_csv:
            listed = pd.read_csv(ambient_csv, keep_default_na=False)
            ambient_genes = set(listed.loc[listed["excluded"].astype(str).str.lower() == "true", "gene"])
        else:
            print("  note: no ambient gene list given; only stress genes are flagged")
        flagged = flag_de_table(de, ambient_genes)
        flagged.to_csv(out / "de_flagged.csv", index=False)
        summary = summarise_hits(flagged, args.padj, args.lfc)
        summary.to_csv(out / "de_summary_counts.csv", index=False)
        print(f"Wrote {out / 'de_flagged.csv'} and {out / 'de_summary_counts.csv'}")
        with pd.option_context("display.width", 200, "display.max_rows", 100):
            print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

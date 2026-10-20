#!/usr/bin/env python3
"""Step 12: attach Parse split-pipe TCR (AIRR) results to the annotated WT AnnData.

Reads the AIRR rearrangement tables that `split-pipe --tcr_analysis` writes
(`<sample>/TCR_filtered/tcr_annotation_airr.tsv`: one row per chain, cells
identified by `cell_id`), matches them to WT cells, and uses scirpy
(Sturm et al., Bioinformatics 2020) to index chains, run chain QC, and define
clonotypes. Adds to `adata.obsm["airr"]` and `adata.obs`:
  has_ir, receptor_type, chain_pairing, clone_id, clone_id_size,
  tcr_source (which AIRR file the cell came from).

Clonotype definition (the user's choice for this project): a clonotype is a
set of cells sharing identical TRA + TRB CDR3 amino-acid sequences
(`--clone-sequence aa`, `receptor_arms="all"`, primary chains only for dual-TCR
cells). `--same-v-gene` and `--clone-sequence nt` make it stricter.

THE ONE THING TO WATCH: the join between TCR cell ids and WT cell names.
split-pipe names cells by their barcode wells, `NN_NN_NN` (the WT `bc_wells`
column, which script 00 uses as the cell name). After `--mode comb` both WT and
TCR append `__s<N>`, N = the sublibrary's position in the --sublibraries list,
so combined TCR matches combined WT directly IF both were combined in the same
order (checked per sublibrary below). This
script matches exactly, reports the match rate, and STOPS if it is below
`--min-match-rate` instead of silently attaching nothing. If the two sides
differ by the sublibrary suffix, prefer --tcr-suffix-format '__s{label}' (rebuild
the suffix on the TCR side from each --airr LABEL); --strip-wt-suffix is
refused whenever the same barcode occurs in more than one sublibrary.

Sanity gate: paired-TCR cells should sit in T-lineage clusters. If most of
them land in B-cell or myeloid clusters, the match is wrong (or the data are
doublets) and the script warns loudly -- don't trust such output.

NOT YET RUN ON REAL SPLIT-PIPE TCR OUTPUT. Validated only on a synthetic
fixture (tests/synthetic_tcr.py) that mimics split-pipe's documented AIRR
columns, which proves the plumbing, not the real file layout.

Usage:
  python 12_tcr_attach.py --wt-input 04b_annotated.h5ad \\
      --airr-root /path/to/split-pipe/tcr_comb --output 12_tcr.h5ad
  # per-sublibrary TCR outputs; LABEL = that sublibrary's position in the WT combine list
  python 12_tcr_attach.py --wt-input 04b_annotated.h5ad --tcr-suffix-format '__s{label}' \\
      --airr tcr_sub1/all-sample/TCR_filtered/tcr_annotation_airr.tsv 1 \\
      --airr tcr_sub2/all-sample/TCR_filtered/tcr_annotation_airr.tsv 2 --output 12_tcr.h5ad
"""
import argparse
import json
import pathlib
import re
import sys
import warnings

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scirpy as ir

warnings.filterwarnings("ignore", category=FutureWarning)

AIRR_GLOB = "**/*tcr_annotation_airr.tsv*"  # split-pipe's clean step may gzip outputs (.tsv.gz)
COMBINED_SAMPLE_DIRS = {"all-sample", "all-well"}
TRUE_STRINGS = {"T", "TRUE", "True", "true", "1", "Y", "yes"}
VALID_TCR_LOCI = {"TRA", "TRB", "TRG", "TRD"}
DEFAULT_TCELL_PATTERN = r"T cell|Treg|CD4|CD8|\bT\b|MAIT|gamma|delta|NKT"


def discover_airr(args) -> list[tuple[pathlib.Path, str]]:
    found = []
    for entry in args.airr or []:
        path = pathlib.Path(entry[0])
        label = entry[1] if len(entry) > 1 else path.parent.parent.name
        found.append((path, label))
    if args.airr_root:
        root = pathlib.Path(args.airr_root)
        paths = [p for p in sorted(root.glob(AIRR_GLOB)) if p.name.endswith((".tsv", ".tsv.gz"))]
        plain = {p for p in paths if p.name.endswith(".tsv")}
        paths = [p for p in paths if not (p.name.endswith(".gz") and p.with_suffix("") in plain)]
        hits = [(path, path.parent.parent.name) for path in paths]
        # <root>/<sample>/TCR_filtered/tcr_annotation_airr.tsv. split-pipe also writes combined
        # "all-sample" (and optionally "all-well") folders holding the SAME cells again; mixing them
        # with per-sample folders would load every cell twice (and fake dual-receptor cells).
        combined = [h for h in hits if h[1] in COMBINED_SAMPLE_DIRS]
        per_sample = [h for h in hits if h[1] not in COMBINED_SAMPLE_DIRS]
        if per_sample:
            if combined:
                print(f"  skipping {len(combined)} combined folder(s) "
                      f"({', '.join(sorted({h[1] for h in combined}))}) -- per-sample folders are used instead")
            found.extend(per_sample)
        else:
            found.extend(combined)
    if not found:
        sys.exit("No AIRR files: pass --airr FILE [LABEL] and/or --airr-root DIR "
                 f"(searched for {AIRR_GLOB}, .tsv or .tsv.gz).")
    missing = [str(p) for p, _ in found if not p.exists()]
    if missing:
        sys.exit(f"AIRR file(s) not found: {missing}")
    return found


def load_airr_tables(files: list[tuple[pathlib.Path, str]]) -> pd.DataFrame:
    frames = []
    for path, label in files:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        df["tcr_source"] = label
        print(f"  {path} -> {len(df):,} chains, {df['cell_id'].nunique():,} cells (label '{label}')")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    if "sequence_id" in df.columns:
        dup = df.duplicated(["cell_id", "sequence_id"])
        if dup.any():
            print(f"  !! dropping {int(dup.sum()):,} duplicated chain rows (same cell_id + sequence_id in more "
                  "than one file) -- check that overlapping folders were not both passed")
            df = df[~dup]
    for col in ("cell_id", "locus", "junction_aa"):
        if col not in df.columns:
            sys.exit(f"AIRR table has no '{col}' column; columns present: {list(df.columns)[:30]}")
    return df


def n_cells(df: pd.DataFrame) -> int:
    # the same bare barcode can occur in two sublibraries, so count (barcode, source) pairs
    return len(df[["cell_id", "tcr_source"]].drop_duplicates())


def normalise_airr(df: pd.DataFrame, keep_multiplets: bool) -> pd.DataFrame:
    n0 = n_cells(df)
    df = df[df["locus"].isin(VALID_TCR_LOCI)].copy()
    df = df[df["junction_aa"].str.len() > 0]
    if "isMultiplet" in df.columns and not keep_multiplets:
        mult = df["isMultiplet"].isin(TRUE_STRINGS)
        print(f"  dropping {n_cells(df[mult]):,} cells flagged isMultiplet by split-pipe")
        df = df[~mult]
    # scirpy's AIRR reader wants real booleans / numbers, split-pipe writes T/F and text
    for col in ("productive", "rev_comp", "stop_codon", "vj_in_frame"):
        if col in df.columns:
            df[col] = df[col].isin(TRUE_STRINGS)
    if "transcript_count" in df.columns:
        df["duplicate_count"] = pd.to_numeric(df["transcript_count"], errors="coerce").fillna(0).astype(int)
    if "read_count" in df.columns:
        df["consensus_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)
    print(f"  kept {n_cells(df):,} of {n0:,} cells with usable TCR chains")
    return df


def wt_keys(adata: ad.AnnData, key_col: str | None, strip: bool) -> pd.Series:
    keys = pd.Series(adata.obs_names if key_col is None else adata.obs[key_col].astype(str).values,
                     index=adata.obs_names)
    if strip:
        keys = keys.str.replace(r"__[^_]+$", "", regex=True)
        if keys.duplicated().any():
            sys.exit("--strip-wt-suffix would make WT cell names collide (the same barcode occurs in "
                     "more than one sublibrary). Match on the full name instead.")
    return keys


def attach(adata: ad.AnnData, airr_df: pd.DataFrame, args) -> tuple[ad.AnnData, dict]:
    keys = wt_keys(adata, args.wt_key_col, args.strip_wt_suffix)
    tcr_id = airr_df["cell_id"].str.replace(r"__[^_]+$", "", regex=True) if args.strip_tcr_suffix else airr_df["cell_id"]
    if args.tcr_suffix_format:
        # per-sublibrary TCR runs give bare barcodes; the same barcode recurs in other sublibraries,
        # so rebuild the combined-style name from the --airr LABEL instead of stripping the WT side
        tcr_id = tcr_id + airr_df["tcr_source"].map(lambda lab: args.tcr_suffix_format.format(label=lab))
    airr_df = airr_df.assign(cell_id=tcr_id)

    tcr_cells = set(airr_df["cell_id"])
    matched = tcr_cells & set(keys)
    rate = len(matched) / max(len(tcr_cells), 1)
    stats = {"n_tcr_cells": len(tcr_cells), "n_matched_to_wt": len(matched), "match_rate": round(rate, 4)}
    print(f"  TCR cells: {len(tcr_cells):,}; matched to a WT cell: {len(matched):,} ({rate:.1%})")
    if rate < args.min_match_rate:
        sys.exit(
            f"Only {rate:.1%} of TCR cell ids match a WT cell name (need >= {args.min_match_rate:.0%}). "
            f"Examples TCR: {sorted(tcr_cells)[:3]}  WT: {list(keys.iloc[:3])}. "
            "Likely a naming difference (sublibrary suffix? different barcode field?). "
            "Try --strip-tcr-suffix / --strip-wt-suffix / --wt-key-col; nothing was written.")

    # split-pipe's --mode comb names cells <bc>__s<N>, where N is the sublibrary's POSITION in the
    # --sublibraries list. If TCR and WT were combined with lists in different orders (or a TCR
    # sublibrary was skipped as empty, which shifts the numbering), one sublibrary's cells look up
    # the wrong WT sublibrary: most fail to match, but by barcode coincidence a few % match the
    # WRONG cell. Overall match rate can still look fine, so check each sublibrary separately.
    group = airr_df["cell_id"].str.extract(r"__(s\d+)$")[0].fillna(airr_df["tcr_source"])
    per = (pd.DataFrame({"cell": airr_df["cell_id"], "group": group}).drop_duplicates()
           .assign(hit=lambda d: d["cell"].isin(matched))
           .groupby("group")["hit"].agg(["mean", "size"]))
    stats["match_rate_by_group"] = per["mean"].round(4).to_dict()
    if len(per) > 1:
        print("  match rate per sublibrary/source: " +
              ", ".join(f"{g} {r:.0%} (n={int(n):,})" for g, (r, n) in per.iterrows()))
    bad = per[(per["mean"] < args.min_match_rate) & (per["size"] >= 20)]
    if len(bad) and not args.allow_partial_match:
        sys.exit(f"Sublibrary/source {list(bad.index)} matched only "
                 f"{', '.join(f'{r:.0%}' for r in bad['mean'])} of its TCR cells while others match well. "
                 "Most likely the TCR and WT `--mode comb --sublibraries` lists were in different orders "
                 "(split-pipe numbers sublibraries __s1, __s2, ... by list position) or a TCR sublibrary was "
                 "skipped. Re-combine in the same order, or attach per-sublibrary outputs with "
                 "--tcr-suffix-format '__s{label}' and the WT sublibrary number as each --airr LABEL. "
                 "(--allow-partial-match overrides; nothing was written.)")

    # point AIRR rows at the WT obs_names so scirpy can align them
    name_for_key = pd.Series(keys.index, index=keys.values)
    airr_df = airr_df[airr_df["cell_id"].isin(matched)].copy()
    airr_df["cell_id"] = airr_df["cell_id"].map(name_for_key)

    # scirpy's merge_airr needs AIRR data on BOTH sides; the WT object has none, so rebuild a
    # cell list in WT order (empty AirrCell for cells with no TCR) and store it in obsm["airr"].
    adata_tcr = ir.io.read_airr(airr_df)
    tcr_cells = {c.cell_id: c for c in ir.io.to_airr_cells(adata_tcr)}
    cells = [tcr_cells[n] if n in tcr_cells else ir.io.AirrCell(n) for n in adata.obs_names]
    adata.obsm["airr"] = ir.io.from_airr_cells(cells).obsm["airr"]
    src = airr_df.drop_duplicates("cell_id").set_index("cell_id")["tcr_source"]
    adata.obs["tcr_source"] = src.reindex(adata.obs_names).values
    return adata, stats


def tcr_qc_and_clonotypes(adata: ad.AnnData, args) -> None:
    ir.pp.index_chains(adata)
    ir.tl.chain_qc(adata)
    adata.obs["has_ir"] = (adata.obs["receptor_type"] != "no IR").astype(str)
    ir.pp.ir_dist(adata, metric="identity", sequence=args.clone_sequence)
    # NOT ir.tl.define_clonotypes: that wrapper is hard-wired to nucleotide identity whatever
    # sequence ir_dist was run on (caught on the fixture: it merged all same-length CDR3s into
    # one clone). define_clonotype_clusters takes the sequence type explicitly.
    ir.tl.define_clonotype_clusters(adata, sequence=args.clone_sequence, metric="identity",
                                    receptor_arms="all", dual_ir=args.dual_ir,
                                    same_v_gene=args.same_v_gene, key_added="clone_id")


def sanity_report(adata: ad.AnnData, args, plot_dir: pathlib.Path) -> dict:
    out = {}
    obs = adata.obs
    paired = obs["chain_pairing"].isin(["single pair", "extra VJ", "extra VDJ", "two full chains"])
    out["n_cells"] = int(adata.n_obs)
    out["n_with_any_tcr"] = int((obs["has_ir"] == "True").sum())
    out["n_with_paired_tcr"] = int(paired.sum())
    out["n_cells_in_a_clonotype"] = int(obs["clone_id"].notna().sum())
    out["n_clonotypes"] = int(obs["clone_id"].nunique())
    out["chain_pairing"] = obs["chain_pairing"].value_counts().to_dict()
    print(f"  cells with any TCR: {out['n_with_any_tcr']:,}; paired alpha+beta: {out['n_with_paired_tcr']:,}; "
          f"clonotypes: {out['n_clonotypes']:,}")
    if out["n_cells_in_a_clonotype"]:
        top = float(obs["clone_id"].value_counts(normalize=True).iloc[0])
        out["largest_clone_fraction"] = top
        if top > 0.25:
            print(f"  !! WARNING: the largest clonotype holds {top:.0%} of clonotyped cells. Real repertoires "
                  "are rarely this dominated; check the clonotype definition / distance settings before "
                  "interpreting expansion.")
    print("  chain pairing:", out["chain_pairing"])

    if args.annotation_col in obs.columns:
        tcell = obs[args.annotation_col].astype(str).str.contains(args.tcell_pattern, regex=True, case=False)
        frac_t = float(tcell[paired].mean()) if paired.any() else float("nan")
        out["paired_tcr_cells_in_t_lineage_fraction"] = frac_t
        if frac_t < 0.7:
            print(f"  !! WARNING: only {frac_t:.0%} of paired-TCR cells are in T-lineage clusters "
                  f"('{args.tcell_pattern}' on {args.annotation_col}). TCR should concentrate in T cells; "
                  "a low value means a wrong barcode match, doublets, or an annotation without T labels. "
                  "Do not interpret clonotype results until this is explained.")
        by_ct = (obs.assign(paired=paired).groupby(args.annotation_col, observed=True)["paired"]
                 .agg(["mean", "sum", "size"]).sort_values("mean", ascending=False))
        by_ct.to_csv(plot_dir / "tcr_fraction_by_celltype.csv")
        fig, ax = plt.subplots(figsize=(7, 0.35 * len(by_ct) + 1.5))
        ax.barh(by_ct.index.astype(str), by_ct["mean"], color="#1f6f6b")
        ax.invert_yaxis()
        ax.set_xlabel("fraction of cells with a paired alpha+beta TCR")
        ax.set_title("TCR detection by cell type (should peak in T-lineage)")
        fig.tight_layout()
        fig.savefig(plot_dir / "tcr_fraction_by_celltype.png", dpi=150)
        plt.close(fig)
    if args.sample_col in obs.columns:
        per_s = obs.assign(paired=paired).groupby(args.sample_col, observed=True)["paired"].agg(["mean", "sum", "size"])
        per_s.to_csv(plot_dir / "tcr_fraction_by_sample.csv")
        out["paired_fraction_by_sample"] = per_s["mean"].round(4).to_dict()
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wt-input", required=True, help="Annotated WT .h5ad (04b output or later)")
    p.add_argument("--airr", nargs="+", action="append", metavar=("FILE", "LABEL"),
                   help="One AIRR tsv + optional label (e.g. sublibrary). Repeatable.")
    p.add_argument("--airr-root", help=f"Directory searched recursively for {AIRR_GLOB}")
    p.add_argument("--output", required=True, help="Output .h5ad (WT + TCR)")
    p.add_argument("--plot-dir", default="tcr_plots")
    p.add_argument("--wt-key-col", default=None, help="obs column holding the WT barcode (default: obs_names)")
    p.add_argument("--strip-tcr-suffix", action="store_true", help=r"remove a trailing __<sublib> from TCR cell ids")
    p.add_argument("--strip-wt-suffix", action="store_true", help=r"remove a trailing __<sublib> from WT names")
    p.add_argument("--tcr-suffix-format", default=None,
                   help="append to each TCR cell id, built from its --airr LABEL. For per-sublibrary TCR "
                        "runs matched to combined WT names (<bc>__s<N>) use '__s{label}' with each LABEL set "
                        "to that sublibrary's position N in the WT --sublibraries list")
    p.add_argument("--min-match-rate", type=float, default=0.5,
                   help="minimum fraction of TCR cells (overall AND per sublibrary) that must match a WT cell")
    p.add_argument("--allow-partial-match", action="store_true",
                   help="do not stop when one sublibrary matches poorly (not recommended)")
    p.add_argument("--keep-multiplets", action="store_true", help="keep cells split-pipe flags isMultiplet")
    p.add_argument("--clone-sequence", choices=["aa", "nt"], default="aa")
    p.add_argument("--dual-ir", choices=["primary_only", "all", "any"], default="primary_only")
    p.add_argument("--same-v-gene", action="store_true")
    p.add_argument("--annotation-col", default="cell_type_v2")
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--tcell-pattern", default=DEFAULT_TCELL_PATTERN)
    args = p.parse_args()

    plot_dir = pathlib.Path(args.plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading WT data: {args.wt_input}")
    adata = ad.read_h5ad(args.wt_input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    print("Reading AIRR tables")
    airr_df = normalise_airr(load_airr_tables(discover_airr(args)), args.keep_multiplets)

    print("Matching TCR cells to WT cells")
    adata, stats = attach(adata, airr_df, args)

    print("Chain QC and clonotype definition")
    tcr_qc_and_clonotypes(adata, args)

    print("Sanity report")
    summary = {**stats, **sanity_report(adata, args, plot_dir)}
    summary["settings"] = {k: getattr(args, k) for k in
                           ("clone_sequence", "dual_ir", "same_v_gene", "strip_tcr_suffix", "strip_wt_suffix", "tcr_suffix_format", "allow_partial_match")}
    (plot_dir / "tcr_attach_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

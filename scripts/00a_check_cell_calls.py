"""Step 0a (optional): check split-pipe's cell calls before trusting DGE_filtered.

split-pipe decides which barcodes are cells with a knee cutoff on transcripts
per barcode, fitted separately for every sample in every sublibrary run
(`cell_tscp_cutoff` in agg_sample_summary.csv). That works when a sample has
one clean knee. When it does not -- for example a high-UMI tail of activated
or proliferating cells above the main population -- the cutoff can jump above
the main population, and DGE_filtered keeps only a few hundred very large
"cells" for that sample. Sublibraries are random splits of the same cells, so
the same sample should get a similar cutoff in every sublibrary; when it does
not, the calls are unreliable.

This script reads split-pipe output folders (per-sublibrary runs and/or the
combined run), and for every sample reports:
  - split-pipe's cutoff and cell count per run
  - how many barcodes pass a range of fixed cutoffs (from DGE_unfiltered)
  - flags: cutoff varies > --max-ratio x across runs for the same sample, or
    sits > --outlier-ratio x above the median cutoff of the other samples in
    the same run
and draws barcode-rank plots (transcripts per barcode vs rank, log-log) with
split-pipe's cutoff marked, which is where to choose a fixed cutoff.

If anything is flagged, load DGE_unfiltered with one cutoff for every sample:
    00_load_parse_output.py --sample-dir <comb>/all-sample/DGE_unfiltered \
        --min-counts N --min-genes M ...
(see docs/03-your-data.md).

Usage:
    python 00a_check_cell_calls.py --run-dir split-pipe/WT/sub1 split-pipe/WT/sub2 ... \
        --output-dir cell_call_check/
    python 00a_check_cell_calls.py --run-dir split-pipe/WT/combined --output-dir cell_call_check/
"""
import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_THRESHOLDS = [300, 500, 1000, 2000, 5000, 15000]


def read_run(run_dir: pathlib.Path, thresholds: list[int]):
    summary_path = run_dir / "agg_sample_summary.csv"
    if not summary_path.exists():
        sys.exit(f"{summary_path} not found -- is {run_dir} a split-pipe output folder?")
    summary = pd.read_csv(summary_path)
    for col in ("sample", "number_of_cells", "cell_tscp_cutoff"):
        if col not in summary.columns:
            sys.exit(f"{summary_path} has no '{col}' column (columns: {list(summary.columns)})")
    summary = summary[summary["sample"] != "all-sample"]

    rows, ranks = [], {}
    for _, s in summary.iterrows():
        sample = str(s["sample"])
        md_path = run_dir / sample / "DGE_unfiltered" / "cell_metadata.csv"
        if not md_path.exists():
            sys.exit(f"{md_path} not found -- this check needs split-pipe's DGE_unfiltered output")
        tscp = pd.read_csv(md_path, usecols=["tscp_count"])["tscp_count"].to_numpy()
        tscp = np.sort(tscp)[::-1]
        ranks[sample] = tscp
        row = {"run": run_dir.name, "sample": sample,
               "splitpipe_cutoff": float(s["cell_tscp_cutoff"]),
               "splitpipe_cells": int(s["number_of_cells"])}
        for t in thresholds:
            row[f"barcodes_ge_{t}"] = int((tscp >= t).sum())
        rows.append(row)
    return pd.DataFrame(rows), ranks


def flag(table: pd.DataFrame, max_ratio: float, outlier_ratio: float) -> pd.DataFrame:
    table = table.copy()
    table["flag"] = ""
    if table["run"].nunique() > 1:
        spread = table.groupby("sample")["splitpipe_cutoff"].agg(lambda c: c.max() / max(c.min(), 1))
        for sample, ratio in spread.items():
            if ratio > max_ratio:
                table.loc[table["sample"] == sample, "flag"] += f"cutoff varies {ratio:.0f}x across runs; "
    for run, df in table.groupby("run"):
        for i, r in df.iterrows():
            others = df.loc[df.index != i, "splitpipe_cutoff"]
            if len(others) and r["splitpipe_cutoff"] > outlier_ratio * others.median():
                table.loc[i, "flag"] += (f"cutoff {r['splitpipe_cutoff'] / others.median():.0f}x the other "
                                         f"samples' median in {run}; ")
    table["flag"] = table["flag"].str.rstrip("; ")
    return table


def plot_ranks(all_ranks: dict, table: pd.DataFrame, path: pathlib.Path):
    runs = list(all_ranks)
    ncol = min(4, len(runs))
    nrow = int(np.ceil(len(runs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), squeeze=False)
    samples = sorted({s for r in all_ranks.values() for s in r})
    colors = {s: plt.cm.tab20(i % 20) for i, s in enumerate(samples)}
    for ax, run in zip(axes.flat, runs):
        for sample, tscp in all_ranks[run].items():
            ax.loglog(np.arange(1, len(tscp) + 1), np.maximum(tscp, 1), color=colors[sample], lw=1, label=sample)
            r = table[(table["run"] == run) & (table["sample"] == sample)].iloc[0]
            n = max(int(r["splitpipe_cells"]), 1)
            ax.plot(n, r["splitpipe_cutoff"], "o", color=colors[sample], ms=5, mec="black", mew=0.5)
        ax.set_title(run, fontsize=9)
        ax.set_xlabel("barcode rank")
        ax.set_ylabel("transcripts per barcode")
    for ax in list(axes.flat)[len(runs):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=6, frameon=False)
    fig.suptitle("Barcode rank per sample (dot = split-pipe's cutoff and cell count)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", nargs="+", required=True,
                        help="split-pipe output folders (each with agg_sample_summary.csv and "
                             "<sample>/DGE_unfiltered/); per-sublibrary runs and/or the combined run")
    parser.add_argument("--output-dir", required=True, help="Where to write the table and plot")
    parser.add_argument("--thresholds", type=int, nargs="+", default=DEFAULT_THRESHOLDS,
                        help=f"Fixed cutoffs to count barcodes at (default {DEFAULT_THRESHOLDS})")
    parser.add_argument("--max-ratio", type=float, default=3.0,
                        help="Flag a sample whose cutoff varies more than this across runs (default 3)")
    parser.add_argument("--outlier-ratio", type=float, default=5.0,
                        help="Flag a cutoff this many times the other samples' median in a run (default 5)")
    args = parser.parse_args()

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tables, all_ranks = [], {}
    for d in args.run_dir:
        run_dir = pathlib.Path(d)
        print(f"Reading {run_dir}...")
        t, ranks = read_run(run_dir, args.thresholds)
        tables.append(t)
        all_ranks[run_dir.name] = ranks
    if len(all_ranks) != len(args.run_dir):
        sys.exit("Two --run-dir folders share a name; give folders with distinct names")

    table = flag(pd.concat(tables, ignore_index=True), args.max_ratio, args.outlier_ratio)
    table.to_csv(out / "cell_call_check.csv", index=False)
    plot_ranks(all_ranks, table, out / "barcode_rank.png")

    with pd.option_context("display.width", 200, "display.max_columns", 50, "display.max_rows", 500):
        print(table.drop(columns="flag").to_string(index=False))
    flagged = table[table["flag"] != ""]
    print(f"\nWrote {out / 'cell_call_check.csv'} and {out / 'barcode_rank.png'}")
    if flagged.empty:
        print("OK: split-pipe's cutoffs are consistent across samples and runs -- DGE_filtered is fine to use.")
    else:
        print(f"FLAGGED {len(flagged)} sample/run rows:")
        for _, r in flagged.iterrows():
            print(f"  {r['run']} / {r['sample']}: {r['flag']} "
                  f"(split-pipe kept {r['splitpipe_cells']:,} cells at >= {r['splitpipe_cutoff']:.0f})")
        print("split-pipe's per-sample cell calls look unreliable. Pick one cutoff from barcode_rank.png and "
              "the barcodes_ge_* columns, then load <combined>/all-sample/DGE_unfiltered with "
              "00_load_parse_output.py --min-counts N --min-genes M (docs/03-your-data.md).")


if __name__ == "__main__":
    main()

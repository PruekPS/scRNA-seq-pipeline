"""Step 0c (recommended): check that sublibraries really hold different cells.

Every sublibrary is a separate portion of the barcoded cells, so two sublibraries
should share only a few cell barcodes, by chance. If two share nearly all of them,
the same material went into both libraries: a tube mix-up during library
preparation, not a sequencing or demultiplexing problem. Nothing in split-pipe's
reports warns about this, and the data still look fine (the extra "cells" are
copies, so cell numbers and molecule counts are inflated). On a real run this
check found two pairs of duplicated libraries, 99.9% of barcodes shared.

What it does, reading each sublibrary's own split-pipe run (not the combined one):
  1. Percent of each sublibrary's cell barcodes also called in every other one
     (table + heatmap). Pairs above --max-overlap are flagged.
  2. For flagged pairs: do shared barcodes have the same sample label, and do their
     transcript counts agree? Yes to both = the same cells.
  3. Optional (--molecules): the share of molecules (barcode + gene + UMI) found in
     both libraries, from each run's process/tscp_assignment.csv.gz. This says WHEN
     the mixing happened. Mixing AFTER cDNA amplification puts copies of one molecule
     in both libraries (many shared). Mixing BEFORE amplification sends each molecule
     to one library (almost none shared beyond chance). Compare with the baseline pair
     it reports, two libraries that are not flagged.

If a pair is flagged: stop, check the bench records, ask the sequencing core, and
decide before analysing. Do not treat the two libraries as independent. One option
is to merge the pair's counts per barcode; another is to keep one of them
(00_load_parse_output.py --exclude-sublibraries).

Usage:
    python 00c_check_sublibrary_overlap.py --run-dir split-pipe/WT/sub1 split-pipe/WT/sub2 ... \
        --output-dir sublibrary_check/
    python 00c_check_sublibrary_overlap.py --run-dir split-pipe/WT/sub* --output-dir sublibrary_check/ --molecules
"""
import argparse
import itertools
import pathlib
import sys

import numpy as np
import pandas as pd

import figstyle

BARCODE_COLUMN_CANDIDATES = ["bc_wells", "cell_barcode", "barcode"]


def read_cells(run_dir: pathlib.Path) -> pd.DataFrame:
    hits = sorted((run_dir / "all-sample" / "DGE_filtered").glob("cell_metadata*.csv*"))
    if not hits:
        sys.exit(f"{run_dir}/all-sample/DGE_filtered/cell_metadata.csv not found -- give each sublibrary's own "
                 "split-pipe output folder (not the combined run, and not a sample folder)")
    cells = pd.read_csv(hits[0])
    barcode = next((c for c in BARCODE_COLUMN_CANDIDATES if c in cells.columns), None)
    if barcode is None:
        sys.exit(f"{hits[0]} has none of the barcode columns {BARCODE_COLUMN_CANDIDATES}")
    # A sublibrary run's barcodes carry no __sN tag; strip it anyway so a combined-style name also matches.
    cells.index = cells[barcode].astype(str).str.replace(r"__s\d+$", "", regex=True)
    return cells


def overlap_table(cells: dict) -> pd.DataFrame:
    names = list(cells)
    table = pd.DataFrame(np.nan, index=names, columns=names)
    for a, b in itertools.product(names, names):
        table.loc[a, b] = 100 * len(cells[a].index.intersection(cells[b].index)) / max(len(cells[a]), 1)
    return table


def compare_pair(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    shared = a.index.intersection(b.index)
    out = {"shared_barcodes": len(shared), "same_sample_pct": np.nan, "tscp_spearman": np.nan}
    if len(shared) and "sample" in a.columns and "sample" in b.columns:
        out["same_sample_pct"] = round(100 * float((a.loc[shared, "sample"].values == b.loc[shared, "sample"].values).mean()), 1)
    if len(shared) > 2 and "tscp_count" in a.columns and "tscp_count" in b.columns:
        out["tscp_spearman"] = round(float(pd.Series(a.loc[shared, "tscp_count"].values).corr(
            pd.Series(b.loc[shared, "tscp_count"].values), method="spearman")), 3)
    return out


def read_molecules(run_dir: pathlib.Path, barcodes: set) -> set:
    path = run_dir / "process" / "tscp_assignment.csv.gz"
    if not path.exists():
        print(f"  no {path} -- split-pipe may have removed it when cleaning up; skipping molecule check")
        return set()
    keep = set()
    for chunk in pd.read_csv(path, usecols=["bc_wells", "gene", "polyN"], chunksize=2_000_000):
        chunk = chunk[chunk["bc_wells"].isin(barcodes)]
        keep.update(zip(chunk["bc_wells"], chunk["gene"], chunk["polyN"]))
    return keep


def molecule_sharing(run_dirs: dict, cells: dict, a: str, b: str, n_cells: int, rng) -> dict:
    """Percent of a's molecules (barcode + gene + UMI) also in b, over barcodes the two share."""
    shared = sorted(cells[a].index.intersection(cells[b].index))
    if not shared:
        return {"molecules_a": 0, "pct_molecules_shared": np.nan}
    pick = set(rng.choice(shared, min(n_cells, len(shared)), replace=False))
    ma, mb = read_molecules(run_dirs[a], pick), read_molecules(run_dirs[b], pick)
    return {"molecules_a": len(ma), "pct_molecules_shared": round(100 * len(ma & mb) / len(ma), 2) if ma else np.nan}


def plot_heatmap(table: pd.DataFrame, path: pathlib.Path):
    n = len(table)
    fig, axes = figstyle.subplots(1, 1, w_mm=min(120, 28 + 7 * n), h_mm=min(120, 22 + 6 * n))
    ax = axes[0, 0]
    im = ax.imshow(table.to_numpy(), vmin=0, vmax=100, cmap=figstyle.SEQUENTIAL, aspect="auto")
    ax.set_xticks(range(n), table.columns, rotation=90)
    ax.set_yticks(range(n), table.index)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for (i, j), v in np.ndenumerate(table.to_numpy()):
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5, color="white" if v < 50 else figstyle.INK)
    colorbar = fig.colorbar(im, ax=ax, shrink=0.6, label="Barcodes also in column (%)")
    colorbar.outline.set_linewidth(0.4)
    figstyle.save_figure(fig, path, source_data=table, check=False,
                         legend="Percent of each row sublibrary's cell barcodes also called as cells in the column "
                                "sublibrary. Expected: low off the diagonal, by chance only.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", nargs="+", required=True,
                        help="Each sublibrary's own split-pipe output folder (holding all-sample/DGE_filtered)")
    parser.add_argument("--output-dir", required=True, help="Where to write the tables and heatmap (created if needed)")
    parser.add_argument("--max-overlap", type=float, default=50.0,
                        help="Flag a pair sharing more than this percent of barcodes (default 50)")
    parser.add_argument("--molecules", action="store_true",
                        help="For flagged pairs, also compare molecules (slow: streams process/tscp_assignment.csv.gz)")
    parser.add_argument("--n-cells", type=int, default=400, help="Shared barcodes sampled by --molecules (default 400)")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    run_dirs = {pathlib.Path(d).name: pathlib.Path(d) for d in args.run_dir}
    if len(run_dirs) != len(args.run_dir):
        sys.exit("Two --run-dir folders share a name; give folders with distinct names")
    if len(run_dirs) < 2:
        sys.exit("Give at least two sublibrary folders to compare")
    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    figstyle.apply_style()

    cells = {}
    for name, d in run_dirs.items():
        cells[name] = read_cells(d)
        print(f"{name}: {len(cells[name]):,} cells")

    table = overlap_table(cells)
    table.round(1).to_csv(out / "barcode_overlap_percent.csv")
    plot_heatmap(table, out / "barcode_overlap")
    with pd.option_context("display.width", 200, "display.max_columns", 100):
        print("\nPercent of row's barcodes also in column:\n" + table.round(1).to_string())

    rows = []
    for a, b in itertools.combinations(run_dirs, 2):
        pct = max(table.loc[a, b], table.loc[b, a])
        if pct > args.max_overlap:
            rows.append({"library_a": a, "library_b": b, "max_overlap_pct": round(pct, 1), **compare_pair(cells[a], cells[b])})
    flagged = pd.DataFrame(rows)

    if flagged.empty:
        print(f"\nOK: no two sublibraries share more than {args.max_overlap:g}% of cell barcodes.")
        print(f"Wrote {out / 'barcode_overlap_percent.csv'} and {out / 'barcode_overlap.png'}")
        return

    if args.molecules:
        rng = np.random.default_rng(args.seed)
        extra = [molecule_sharing(run_dirs, cells, r["library_a"], r["library_b"], args.n_cells, rng) for r in rows]
        flagged = pd.concat([flagged, pd.DataFrame(extra)], axis=1)
        flagged_names = {(r["library_a"], r["library_b"]) for r in rows}
        normal = [p for p in itertools.combinations(run_dirs, 2) if p not in flagged_names]
        if normal:
            a, b = normal[0]
            base = molecule_sharing(run_dirs, cells, a, b, args.n_cells, rng)
            print(f"\nBaseline, an unflagged pair ({a} vs {b}): {base['pct_molecules_shared']}% of molecules shared")
    flagged.to_csv(out / "flagged_pairs.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 100):
        print("\nFLAGGED: these sublibraries share most of their cell barcodes:\n" + flagged.to_string(index=False))
    print("\nTreat them as copies of the same cells until the bench records explain it: do not count them as "
          "independent libraries. Next steps are in docs/03-your-data.md.")


if __name__ == "__main__":
    main()

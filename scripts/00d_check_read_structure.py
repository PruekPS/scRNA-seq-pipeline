"""Step 0d (optional): look for short-insert libraries in the FASTQ files.

When a library contains many very short cDNA fragments, read 1 runs through the
short insert and on into the polyA tail, then the barcode region. Those stretches
can align spuriously, so a real library showed one fake highly-expressed gene per
barcode well, and clusters that were one sample + one sublibrary + one well. The
sequencing core's QC hinted at it: a much lower share of aligned reads for that one
library. This check compares libraries from a stream of their first reads.

It reports, per library, the share of read 1 with:
  - a polyA stretch (--polya-length, default 15 A),
  - the standard Illumina adapter (AGATCGGAAGAGC), and the median insert before it,
  - optionally your own sequence (--adapter), for example the reverse complement of a
    barcode linker from Parse's documentation (not distributed here).
A library is flagged if its polyA or --adapter share is at least --fold times the
median of the other libraries and above --min-percent.

If one library is flagged, see "A library with short inserts" in docs/03-your-data.md
for the cutadapt rescue. Run this on a compute node: it reads compressed files.

Usage:
    python 00d_check_read_structure.py --output-dir read_check/ \
        --fastq lib1 fastq/lib1_R1.fastq.gz --fastq lib2 fastq/lib2_R1.fastq.gz \
        --fastq lib3 fastq/lib3_R1.fastq.gz
"""
import argparse
import gzip
import pathlib
import sys

import numpy as np
import pandas as pd

ILLUMINA_ADAPTER = "AGATCGGAAGAGC"


def read_sequences(path: pathlib.Path, n_reads: int):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i // 4 >= n_reads:
                break
            if i % 4 == 1:
                yield line.strip()


def summarise(name: str, path: pathlib.Path, n_reads: int, polya: str, adapter: str | None) -> dict:
    total = with_polya = with_illumina = with_adapter = 0
    inserts = []
    for seq in read_sequences(path, n_reads):
        total += 1
        with_polya += polya in seq
        pos = seq.find(ILLUMINA_ADAPTER)
        with_illumina += pos >= 0
        inserts.append(pos if pos >= 0 else len(seq))
        if adapter:
            with_adapter += adapter in seq
    if total == 0:
        sys.exit(f"No reads found in {path}")
    row = {"library": name, "reads_checked": total, "pct_polyA": round(100 * with_polya / total, 2),
           "pct_illumina_adapter": round(100 * with_illumina / total, 2),
           "median_insert_before_adapter": int(np.median(inserts))}
    if adapter:
        row["pct_your_adapter"] = round(100 * with_adapter / total, 2)
    return row


def flag(table: pd.DataFrame, fold: float, min_percent: float) -> pd.DataFrame:
    table = table.copy()
    table["flag"] = ""
    columns = [c for c in ("pct_polyA", "pct_your_adapter") if c in table.columns]
    for i in table.index:
        for col in columns:
            others = table.loc[table.index != i, col]
            median = others.median() if len(others) else np.nan
            value = table.loc[i, col]
            if value >= min_percent and median == median and value >= fold * max(median, 0.01):
                table.loc[i, "flag"] += f"{col} {value:g}% vs {median:g}% in the others; "
    table["flag"] = table["flag"].str.rstrip("; ")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastq", nargs=2, action="append", required=True, metavar=("NAME", "R1_FASTQ"),
                        help="A library name and its read 1 FASTQ (.fastq or .fastq.gz). Repeat for each library.")
    parser.add_argument("--output-dir", required=True, help="Where to write the table (created if needed)")
    parser.add_argument("--n-reads", type=int, default=500_000, help="Reads checked per library (default 500,000)")
    parser.add_argument("--polya-length", type=int, default=15, help="A-run length that counts as polyA (default 15)")
    parser.add_argument("--adapter", default=None, help="Optional extra sequence to count in read 1")
    parser.add_argument("--fold", type=float, default=5.0, help="Flag at this many times the others' median (default 5)")
    parser.add_argument("--min-percent", type=float, default=3.0, help="Never flag below this percent (default 3)")
    args = parser.parse_args()

    rows = []
    for name, path in args.fastq:
        path = pathlib.Path(path)
        if not path.exists():
            sys.exit(f"File not found: {path}")
        print(f"Reading {name}: {path}")
        rows.append(summarise(name, path, args.n_reads, "A" * args.polya_length, args.adapter))
    table = flag(pd.DataFrame(rows), args.fold, args.min_percent)

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "read_structure.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(table.drop(columns="flag").to_string(index=False))
    flagged = table[table["flag"] != ""]
    if flagged.empty:
        print(f"\nOK: no library stands out. Wrote {out / 'read_structure.csv'}")
        return
    print(f"\nFLAGGED {len(flagged)} librar{'y' if len(flagged) == 1 else 'ies'}:")
    for _, r in flagged.iterrows():
        print(f"  {r['library']}: {r['flag']}")
    print("Possible short-insert library. See 'A library with short inserts' in docs/03-your-data.md.")


if __name__ == "__main__":
    main()

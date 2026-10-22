"""Step 15 (optional): preranked GSEA on per-cell-type DE results from 11.

Ranks every tested gene by DESeq2's Wald statistic ('stat' column: signed,
precision-weighted, so it uses both direction and confidence rather than
an arbitrary padj/log2FC cutoff) and runs preranked GSEA (gseapy.prerank,
Subramanian et al., PNAS 2005, doi:10.1073/pnas.0506580102) against a GMT
gene-set file, separately for each cell type.

Gene sets are read from a local GMT (no network calls inside an analysis
step) -- fetch once, e.g. MSigDB Hallmark via
    gseapy.get_library("MSigDB_Hallmark_2020", organism="Human")
and write it out as GMT.

Inputs: the de_<cell_type>.csv files written by scripts/11_pseudobulk_deseq2.R
(columns gene, stat, padj, log2FoldChange, ...).
Outputs (in --output-dir): gsea_<cell_type>.csv per cell type and
gsea_all_celltypes.csv (cell_type, term, NES, pval, FDR, leading-edge genes).

Usage:
    python 15_gsea_preranked.py --de-dir de_results/ --gmt hallmark.gmt --output-dir gsea/
"""
import argparse
import pathlib

import numpy as np
import pandas as pd


def read_gmt(path):
    sets = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--de-dir", required=True, dest="de_dir")
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--output-dir", required=True, dest="output_dir")
    ap.add_argument("--min-size", type=int, default=10, dest="min_size")
    ap.add_argument("--max-size", type=int, default=500, dest="max_size")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    import gseapy  # imported here so --help works without it installed

    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sets = read_gmt(args.gmt)
    print(f"{len(sets)} gene sets from {args.gmt}")
    files = sorted(f for f in pathlib.Path(args.de_dir).glob("de_*.csv") if f.name != "de_all_celltypes_combined.csv")
    if not files:
        raise SystemExit(f"no de_*.csv files in {args.de_dir} -- run scripts/11_pseudobulk_deseq2.R first")

    allres = []
    for f in files:
        de = pd.read_csv(f)
        ct = de["cell_type"].iloc[0] if "cell_type" in de.columns else f.stem[3:]
        rnk = (de.dropna(subset=["stat"]).drop_duplicates("gene")
               .set_index("gene")["stat"].sort_values(ascending=False))
        rnk = rnk[np.isfinite(rnk)]
        if len(rnk) < 100:
            print(f"SKIP {ct}: only {len(rnk)} ranked genes")
            continue
        res = gseapy.prerank(rnk=rnk, gene_sets=sets, min_size=args.min_size, max_size=args.max_size,
                             permutation_num=args.permutations, seed=args.seed, threads=args.threads,
                             outdir=None, verbose=False).res2d
        res = res.rename(columns={"Term": "term", "NOM p-val": "pval", "FDR q-val": "FDR", "Lead_genes": "leading_edge"})
        res["cell_type"] = ct
        res["NES"] = res["NES"].astype(float)
        res["FDR"] = res["FDR"].astype(float)
        res = res.sort_values("NES", ascending=False)
        res.to_csv(out / f"gsea_{f.stem[3:]}.csv", index=False)
        allres.append(res)
        top = res[res["FDR"] < 0.05]
        print(f"{ct}: {len(rnk):,} ranked genes; {len(top)} sets at FDR<0.05 "
              f"(up: {int((top.NES > 0).sum())}, down: {int((top.NES < 0).sum())})")
    if allres:
        pd.concat(allres).to_csv(out / "gsea_all_celltypes.csv", index=False)
        print(f"Wrote {out}/gsea_all_celltypes.csv")


if __name__ == "__main__":
    main()

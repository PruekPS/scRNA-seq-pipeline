"""Step 5c (optional): score published gene programs per cell and summarise per sample.

A cluster label says what a cell *is*; a program score says how strongly it runs
a given function (Th1, Tfh, effector-Treg, exhaustion, IFN response, ...). Scores
use scanpy's score_genes: mean expression of the program genes minus the mean of
an expression-matched random control set (Satija et al., Nat Biotechnol 2015).

Built-in T-cell programs (--programs builtin): canonical marker genes for CD4
helper states, Treg states (resting / effector / tissue / Tfr), CD8 states,
tissue residency, IFN response, an activation-induced-marker (AIM) set, cycling,
and a dissociation-stress set to use as a QC covariate. The gene choices follow
the T-cell reference atlases and reviews listed in docs/05-core-pipeline.md
(e.g. ProjecTILs, Andreatta et al. Nat Commun 2021; pan-cancer T cells, Zheng et
al. Science 2021; tissue Treg, Miragaia et al. Immunity 2019; Tr1, Gagliani et
al. Nat Med 2013). They are written in mouse symbols; --species human converts
them (uppercase, plus the few genes whose human name differs).

Your own programs: --programs my_programs.csv with columns program,gene.

Caveat: AIM-like and some activation genes (Cd69, Nr4a1, Egr2) are also
immediate-early genes that rise with tissue dissociation, so a baseline "AIM"
score is NOT an AIM assay (which needs antigen restimulation). The
"AIM core (no IEG)" program leaves those genes out; compare programs against the
"Dissociation stress (QC)" score before interpreting them.

Outputs:
  --output        .h5ad with obs columns prog_<program>
  --summary-csv   mean score per sample (rows) x program, restricted to cells
                  where --subset-col is one of --subset-values (if given), plus
                  the number of cells used; add --group-col to carry a sample-level
                  design column (e.g. condition) into the table.

Usage:
    python 05c_score_gene_programs.py --input 04b_annotated.h5ad --output 05c_programs.h5ad \
        --species mouse --summary-csv programs_per_sample.csv --group-col condition
    # Treg programs, Treg cells only:
    python 05c_score_gene_programs.py --input t_cells.h5ad --output t_programs.h5ad --species mouse \
        --subset-col cell_type --subset-values Treg --summary-csv treg_programs_per_sample.csv
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc

BUILTIN = {  # mouse symbols
    "Naive": ["Sell", "Ccr7", "Lef1", "Tcf7", "Satb1", "Klf2", "S1pr1"],
    "Effector/memory": ["Cd44", "S100a4", "S100a6", "Ahnak", "Anxa1", "Lgals1", "Crip1"],
    "Th1": ["Tbx21", "Ifng", "Cxcr3", "Il12rb2", "Ccl5", "Nkg7", "Fasl"],
    "Th2": ["Gata3", "Il4", "Il5", "Il13", "Il1rl1", "Ccr4", "Ptgdr2"],
    "Th17": ["Rorc", "Il17a", "Il17f", "Il23r", "Ccr6", "Il22", "Rora"],
    "Tfh": ["Cxcr5", "Bcl6", "Pdcd1", "Icos", "Il21", "Tox2", "Sh2d1a", "Maf", "Ascl2"],
    "Tr1 (IL-10)": ["Il10", "Lag3", "Itga2", "Maf", "Prdm1", "Havcr2", "Gzmb"],
    "Cytotoxic": ["Gzmb", "Gzmk", "Prf1", "Nkg7", "Ccl5", "Crtam", "Eomes"],
    "Tissue-resident": ["Cd69", "Itgae", "Cxcr6", "Rgs1", "Itga1", "Bhlhe40", "Zfp683"],
    "IFN response": ["Isg15", "Ifit1", "Ifit3", "Irf7", "Rsad2", "Stat1", "Bst2"],
    "Exhaustion": ["Pdcd1", "Tox", "Havcr2", "Lag3", "Tigit", "Ctla4", "Entpd1"],
    "Progenitor exhausted": ["Tcf7", "Slamf6", "Cxcr5", "Id3", "Pdcd1", "Tox"],
    "AIM-like activation": ["Cd40lg", "Tnfrsf4", "Tnfrsf9", "Il2ra", "Cd69", "Nr4a1", "Egr2", "Irf4"],
    "AIM core (no IEG)": ["Cd40lg", "Tnfrsf4", "Tnfrsf9", "Il2ra", "Irf4", "Icos"],
    "Treg resting": ["Foxp3", "Sell", "Ccr7", "Tcf7", "Satb1", "Bach2", "Lef1"],
    "Treg effector": ["Foxp3", "Ikzf2", "Klrg1", "Prdm1", "Il10", "Gzmb", "Icos", "Tnfrsf18"],
    "Treg tissue/tumor": ["Foxp3", "Il1rl1", "Areg", "Pparg", "Ccr8", "Tnfrsf4", "Tnfrsf9", "Gata3"],
    "Tfr": ["Foxp3", "Cxcr5", "Bcl6", "Pdcd1", "Tox2", "Ascl2"],
    "Proliferation": ["Mki67", "Top2a", "Stmn1", "Tuba1b", "Hmgb2", "Birc5"],
    "Dissociation stress (QC)": ["Hspa1a", "Hspa1b", "Hsph1", "Dnajb1", "Fos", "Jun"],
}
HUMAN_EXCEPTIONS = {"Zfp683": "ZNF683"}


def load_programs(spec, species):
    if spec == "builtin":
        progs = {k: list(v) for k, v in BUILTIN.items()}
        if species == "human":
            progs = {k: [HUMAN_EXCEPTIONS.get(g, g.upper()) for g in v] for k, v in progs.items()}
        return progs
    df = pd.read_csv(spec, dtype=str)
    if not {"program", "gene"} <= set(df.columns):
        raise SystemExit(f"{spec} needs columns program,gene (has {list(df.columns)})")
    return {p: list(d["gene"].dropna()) for p, d in df.groupby("program", sort=False)}


def log_normalized(adata, counts_layer):
    """score_genes needs log-normalised expression over all genes."""
    if counts_layer in adata.layers:
        x = adata.copy()
        x.X = x.layers[counts_layer].copy()
        sc.pp.normalize_total(x, target_sum=1e4)
        sc.pp.log1p(x)
        return x
    if adata.raw is not None:
        return adata.raw.to_adata()
    print(f"  note: no '{counts_layer}' layer or .raw; assuming adata.X is already log-normalised")
    return adata


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True, help=".h5ad with prog_<program> columns added")
    ap.add_argument("--programs", default="builtin", help="'builtin' or a CSV with columns program,gene")
    ap.add_argument("--species", choices=["mouse", "human"], default="human",
                    help="Symbol style of the built-in programs (ignored for a CSV)")
    ap.add_argument("--counts-layer", default="counts")
    ap.add_argument("--min-genes", type=int, default=3, help="Skip programs with fewer genes present (default 3)")
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--group-col", default=None, help="Sample-level column to carry into the summary")
    ap.add_argument("--subset-col", default=None, help="Summarise only cells whose value here is in --subset-values")
    ap.add_argument("--subset-values", nargs="+", default=None)
    ap.add_argument("--summary-csv", default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    progs = load_programs(args.programs, args.species)
    x = log_normalized(adata, args.counts_layer)
    used = {}
    for name, genes in progs.items():
        present = [g for g in dict.fromkeys(genes) if g in x.var_names]
        if len(present) < args.min_genes:
            print(f"  skip '{name}': {len(present)}/{len(genes)} genes present (< {args.min_genes})")
            continue
        sc.tl.score_genes(x, present, score_name="_s", random_state=args.seed)
        adata.obs[f"prog_{name}"] = x.obs["_s"].values
        used[name] = present
        print(f"  {name}: {len(present)}/{len(genes)} genes")
    if not used:
        raise SystemExit("No program had enough genes in this dataset -- check --species / gene symbols")
    adata.uns["gene_programs"] = {k: list(v) for k, v in used.items()}

    if args.summary_csv:
        if args.sample_col not in adata.obs:
            raise SystemExit(f"--sample-col '{args.sample_col}' not in obs")
        keep = np.ones(adata.n_obs, dtype=bool)
        if args.subset_col:
            if args.subset_col not in adata.obs or not args.subset_values:
                raise SystemExit("--subset-col needs an obs column and --subset-values")
            keep = adata.obs[args.subset_col].astype(str).isin(args.subset_values).values
            print(f"  summary over {int(keep.sum()):,} cells with {args.subset_col} in {args.subset_values}")
        obs = adata.obs.loc[keep]
        cols = [f"prog_{k}" for k in used]
        summ = obs.groupby(args.sample_col, observed=True)[cols].mean()
        summ.columns = list(used)
        summ.insert(0, "n_cells", obs.groupby(args.sample_col, observed=True).size())
        if args.group_col:
            if args.group_col not in adata.obs:
                raise SystemExit(f"--group-col '{args.group_col}' not in obs")
            summ.insert(0, args.group_col, obs.groupby(args.sample_col, observed=True)[args.group_col].first())
        pathlib.Path(args.summary_csv).parent.mkdir(parents=True, exist_ok=True)
        summ.to_csv(args.summary_csv)
        print(f"  wrote {args.summary_csv}")

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Scored {len(used)} programs on {adata.n_obs:,} cells -> {args.output}")


if __name__ == "__main__":
    main()

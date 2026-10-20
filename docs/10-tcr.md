# 10. T-cell receptors (Evercode TCR) — clonotypes, sharing, expansion

Every T cell carries a receptor made of two chains, alpha (**TRA**) and beta
(**TRB**). The **CDR3** region of each chain is the part that recognises
antigen, and it is close to unique to each T-cell family: when a T cell
divides after meeting its antigen, all its descendants keep the same CDR3s.
Cells sharing identical TRA and TRB CDR3 sequences form a **clonotype** — one
family. A clonotype with many cells has **expanded**, which usually means it
has been responding to something.

> **Status: tested on synthetic data in split-pipe's exact output format
> (page 4), not yet on real TCR output.**

---

## Before you start

- The WT chain through step 04b is done (pages 5–6).
- split-pipe has processed each TCR sublibrary against its matching WT
  sublibrary and **combined them in the same order as the WT sublibraries**
  (page 3).
- The TCR libraries were sequenced long enough to read the CDR3 and the cell
  barcodes.

---

## Run it

```bash
python scripts/12_tcr_attach.py --wt-input 04b_annotated.h5ad \
  --airr-root <tcr_combined> --output 12_tcr.h5ad --plot-dir tcr_plots/
python scripts/13_tcr_analysis.py --input 12_tcr.h5ad --output-dir tcr_results/ \
  --sample-col sample --group-col patient --annotation-col cell_type_v2
```

Both use [scirpy](https://doi.org/10.1093/bioinformatics/btaa611), the standard
Python package for single-cell receptor data. On the cluster,
`pipeline/run_tcr.slurm` runs both (page 11).

---

## Step 12 — attach receptors to cells

It reads split-pipe's TCR tables (`<sample>/TCR_filtered/tcr_annotation_airr.tsv`,
or `.tsv.gz`) — a standard "AIRR" format, one row per chain — and matches them
to your annotated cells by name.

**Read its first lines carefully.** It prints how many TCR cells matched:

```
match rate per sublibrary/source: s1 98% (n=306), s2 98% (n=305), s3 99% (n=292)
```

| What you see | Meaning |
|---|---|
| ~95–100% overall and per sublibrary | good; the missing few are cells QC removed |
| it **stops**: under 50% matched overall | cell names differ between WT and TCR, usually the sublibrary tag; the message says which option to use |
| it **stops**: one sublibrary matches poorly | TCR was combined in a different sublibrary order than WT; re-combine in the same order |

For TCR sublibraries processed one by one (not combined), give each table with
the sublibrary's number in the WT combine list:

```bash
python scripts/12_tcr_attach.py --wt-input 04b_annotated.h5ad --tcr-suffix-format '__s{label}' \
  --airr tcr_sub1/all-sample/TCR_filtered/tcr_annotation_airr.tsv 1 \
  --airr tcr_sub2/all-sample/TCR_filtered/tcr_annotation_airr.tsv 2 \
  --output 12_tcr.h5ad
```

It then defines clonotypes from **identical TRA + TRB CDR3 amino-acid
sequences** (`--same-v-gene` or `--clone-sequence nt` make it stricter), and
checks one thing biology guarantees:

**Check `tcr_plots/tcr_fraction_by_celltype.png`:** receptors should be found
mostly in T cells, rarely in B cells or monocytes. If not, the match is wrong
or the data are doublets — and nothing after it should be trusted. The script
warns when fewer than 70% of paired-receptor cells sit in T-lineage clusters.

---

## Step 13 — the three questions

| Question | Output |
|---|---|
| How diverse is each sample's repertoire? | `repertoire_by_sample.csv`, `clonality_by_sample.png` |
| Which clonotypes appear in more than one sample — before and after, blood and tissue, or across patients ("public" clones)? | `shared_clones.csv`, `sample_overlap_counts.png`, `sample_overlap_jaccard.png` |
| Which cell states hold the expanded clones? | `expanded_by_celltype_by_sample.csv`, `expanded_by_celltype.png`, `top_clones.csv`, `embedding_expansion.png` |

**Compare diversity on the rarefied columns.** A sample with more cells always
looks more diverse, so every sample is repeatedly cut down to the same number
of cells (`rarefied_*` columns). The raw numbers are there for reference only.

**A clone counts as expanded** when it has at least 2 cells within the same
sample (`--expanded-min-size`). By default only cells with both chains are
used; single-chain cells cannot be assigned to a clone reliably.

**No p-values, on purpose.** Cells from one sample are not independent. To
test a difference between conditions, use the per-sample fractions in
`expanded_by_celltype_by_sample.csv` as the data points, with samples or
patients as replicates — the same reasoning as pseudobulk (page 9).

---

Next: [11. Running at scale with SLURM](11-running-on-slurm.md)

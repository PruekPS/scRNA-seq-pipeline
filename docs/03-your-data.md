# 3. Bringing in your data

Parse Evercode data cannot go straight from the sequencer into this pipeline.
Parse's split-pool chemistry encodes each cell's identity as a *combination* of
barcodes across three rounds, inside the reads themselves, so only Parse's own
software can assign reads to cells. That happens first, outside this
repository; this pipeline starts from its output.

Two routes produce that output:

| Route | What it is | Output this pipeline reads |
|---|---|---|
| **Trailmaker** | Parse's cloud tool, no code | a Seurat `.rds` object, converted once (below) |
| **split-pipe** | Parse's command-line pipeline, run on the cluster | `DGE_filtered/` folders, read directly |

---

## Route A — Trailmaker

Trailmaker exports a Seurat `.rds`. Convert it to the same flat files split-pipe
writes, then load:

```bash
Rscript scripts/00b_seurat_to_dge_format.R your_export.rds dge_from_seurat/
python scripts/00_load_parse_output.py --sample-dir dge_from_seurat/ my_sample --output 00_raw.h5ad
```

---

## Route B — split-pipe

**split-pipe is Parse Biosciences' licensed software and is not included in
this repository.** Request it from Parse
([Pipeline Download Request Form](https://www.parsebiosciences.com/pipeline-download/))
and install and run it following Parse's own documentation, as jobs on compute
nodes. [`pipeline/split-pipe/README.md`](../pipeline/split-pipe/README.md)
lists what this pipeline expects from its output.

In short: each sublibrary is processed separately, then merged with
split-pipe's combine mode. The merged output looks like this:

```
<combined>/
  all-sample/DGE_filtered/      every cell, all samples together
    count_matrix.mtx            counts: cells x genes
    all_genes.csv               gene_id, gene_name, genome
    cell_metadata.csv           bc_wells, sample, ... one row per cell
  <sample-1>/DGE_filtered/      the same cells again, split by sample
  <sample-2>/DGE_filtered/
  ...
```

Load **either** `all-sample` **or** the per-sample folders — never both, or
every cell is counted twice:

```bash
python scripts/00_load_parse_output.py \
  --sample-dir <combined>/all-sample/DGE_filtered all-sample \
  --sample-sheet samples.csv --output 00_raw.h5ad
```

---

## Your sample sheet

split-pipe knows only sample *names*. Anything else about a sample — patient,
time point, tissue, treatment — has to come from you, as a small CSV:

```csv
sample,patient,timepoint,tissue
P1_pre,P1,pre,blood
P1_6h,P1,6h,blood
P2_pre,P2,pre,blood
P2_6h,P2,6h,blood
```

| Rule | Why |
|---|---|
| a column named `sample`, one row per sample | joined onto every cell by sample name |
| names match split-pipe's sample names **exactly** | the loader stops if any sample with cells has no row |
| no column names that already exist in `cell_metadata.csv` | the loader refuses rather than overwrite |
| avoid spaces in sample names | `P1_pre`, not `P1 pre` — spaces cause trouble in command lines and file names |

These columns are what later stages compare: `patient` for paired designs
(page 9) and grouping clones (page 10), `timepoint`/`tissue` as the condition
to test.

---

## Cell names and sublibraries

split-pipe names every cell by its barcode wells, for example `02_77_48`.
After combining, it appends the sublibrary's position in the combine list:
`02_77_48__s1`, `02_77_48__s2`, and so on. The loader keeps these names
unchanged (the TCR step matches on them) and copies the `__sN` tag into a
`sublibrary` column.

Why that matters: in split-pool data every sublibrary contains cells from
every sample. Sublibrary — not sample — is the technical batch, and it is what
batch correction should remove (page 6).

---

## T-cell receptor (Evercode TCR) data

TCR sublibraries are processed by split-pipe against their matching WT
sublibraries, then combined. Two things to get right before you start:

- **Combine TCR sublibraries in the same order as the WT sublibraries.** The
  `__sN` tags come from list position, so a different order makes cells point
  at the wrong sublibrary. Step 12 checks this per sublibrary and stops if it
  does not line up (page 10).
- **Sequencing must be long enough.** Read 1 has to cover the CDR3 region and
  read 2 the cell barcodes; a short check run (for example 15 bp single-end)
  cannot be analysed. Parse's Evercode TCR documentation gives the required
  read configuration.

---

## A real gotcha: blank gene names

In Parse's public 1M PBMC dataset about one gene in three (20,055 of 60,664)
has an empty `gene_name` — a blank CSV field, not the text "NaN". The loader
falls back to `gene_id` for those automatically. If you ever write your own
loader, check for this, or it will crash (or silently drop a third of your
genes).

---

Next: [4. Test small before you test big](04-test-small-first.md)

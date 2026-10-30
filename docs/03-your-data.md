# 3. Bringing in your data

> **Before you start**
> - The code and your project folders exist ([page 2](02-installation.md)).
> - You are in the code folder, with `$PROJECT` and `$OUT` set ([page 1](01-getting-started.md)):
>   `cd /blue/<GROUP>/<USER>/scRNA-seq-pipeline`.
> - Commands on this page that run a script belong inside a job or an interactive session
>   ([page 4](04-test-small-first.md)), with `conda activate scrna`.
> - `<combined>` below means the full path to your combined split-pipe output folder.

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

First put the `.rds` file in your project folder (page 0 shows how to upload). The
script creates the `dge_from_seurat` folder:

```bash
Rscript scripts/00b_seurat_to_dge_format.R $PROJECT/your_export.rds $PROJECT/dge_from_seurat/
python scripts/00_load_parse_output.py --sample-dir $PROJECT/dge_from_seurat/ my_sample --output $OUT/00_raw.h5ad
```

**You should see:** `Combined: ... cells x ... genes across 1 sample(s)`.

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
  all-sample/DGE_unfiltered/    every barcode, including empty ones (same three files)
  <sample-1>/DGE_filtered/      the same cells again, split by sample
  <sample-2>/DGE_filtered/
  ...
```

`DGE_filtered` holds only the barcodes split-pipe called as cells. Normally
that is what you load. If split-pipe's cell calls look wrong, use
`DGE_unfiltered` with one cutoff instead (see
[Optional: when split-pipe's cell calls are wrong](#optional-when-split-pipes-cell-calls-are-wrong)
below).

Load **either** `all-sample` **or** the per-sample folders — never both, or
every cell is counted twice:

```bash
python scripts/00_load_parse_output.py \
  --sample-dir <combined>/all-sample/DGE_filtered all-sample \
  --sample-sheet $PROJECT/samples.csv --output $OUT/00_raw.h5ad
```

---

## Check first: are your sublibraries really different cells?

Each sublibrary is a separate portion of the barcoded cells, so two sublibraries should share
only a few barcodes, by chance. If two share nearly all of them, the same material went into
both libraries. Nothing warns you: split-pipe's reports look normal, and the copies simply
inflate your cell numbers. On one real run this check found two duplicated pairs (99.9% of
barcodes shared), a tube mix-up during library preparation.

It takes a minute. Give it **each sublibrary's own split-pipe folder**, not the combined one
(the `*` matches every folder starting `sub`; use your own folder names):

```bash
python scripts/00c_check_sublibrary_overlap.py \
  --run-dir <split-pipe-out>/sub* --output-dir $PROJECT/results/sublibrary_check
```

**You should see:** `OK: no two sublibraries share more than 50% of cell barcodes.`

**If it says `FLAGGED`:** stop before analysing. Check your bench records, ask the sequencing
core, and do not treat the pair as independent. `same_sample_pct` near 100 and `tscp_spearman`
near 1 mean they are the same cells. Add `--molecules` to learn when the mixing happened: many
shared molecules (compared with the baseline pair it prints) mean amplified cDNA was mixed; almost
none means the mixing came before amplification. Then either keep one of the pair, with
`--exclude-sublibraries` at load, or merge their counts per barcode. Keeping one is simplest;
`__sN` is its position in your combine list.

---

## Optional: a library with short inserts

If the sequencing core's QC shows one library with far fewer aligned reads, or the UMAP has small
islands that are each one sample + one sublibrary + one barcode well, that library may hold many
very short fragments. Read 1 then runs past the insert into the polyA tail and the barcode region,
and these stretches align spuriously. In one real library, every barcode well got its own fake
highly-expressed gene. Check with:

```bash
python scripts/00d_check_read_structure.py --output-dir $PROJECT/results/read_check \
  --fastq lib1 <fastq-folder>/lib1_R1.fastq.gz \
  --fastq lib2 <fastq-folder>/lib2_R1.fastq.gz \
  --fastq lib3 <fastq-folder>/lib3_R1.fastq.gz
```

It ends with `OK: no library stands out` or `FLAGGED`. To rescue a flagged library, trim read 1 at
the polyA tail (and at the reverse complement of the barcode linker, from Parse's documentation)
with [cutadapt](https://cutadapt.readthedocs.io/), keeping read pairs, then rerun split-pipe on that
library alone. Install cutadapt in its own environment, inside a job:
`mamba create -n cutadapt -c conda-forge -c bioconda cutadapt`. The command, with both reads:

```bash
cutadapt -a "A{12}" -a <LINKER_REVERSE_COMPLEMENT> -m 25 -j 8 \
  -o lib_trimmed_R1.fastq.gz -p lib_trimmed_R2.fastq.gz lib_R1.fastq.gz lib_R2.fastq.gz
```

On the one real library where this was done, the fake genes dropped from 243 to 2 and mapping rose
from 54% to 65.5%, but the share of valid barcodes fell from 84.8% to 67.4% for a reason not yet
explained. Compare the result with a clean library before relying on it, or exclude the library
instead.

---

## Single-nuclei and FFPE samples (not yet tested here)

**This pipeline has not yet been run on single-nuclei FFPE data**, so treat the following as what to
check, not as settled settings.

- **split-pipe.** Version 1.9.0 has an `--ffpe` option. In its source code, the option raises the
  minimum transcripts for its cell-calling search from 30 to 200, and allows only v3 and v4 chemistry.
  Which kit name to use, and whether that version fits your kit, is in Parse's manual for your kit.
  Do a dry run on one sublibrary before submitting everything.
- **Cell calling.** Expect lower counts per nucleus. Run `00a_check_cell_calls.py`, and choose any
  fixed cutoff from your own `barcode_rank.png` rather than from the examples on this page.
- **QC (step 01).** Its defaults were set on whole cells, so check `qc_plots/qc_violin.png` yourself:
  nuclei should have few mitochondrial transcripts, and a higher share may mean remaining cytoplasm
  or ambient RNA. `--min-genes` (default 200) may remove too many nuclei from degraded FFPE samples.
- **Ambient RNA.** Likely a bigger issue here. See step 05d on [page 5](05-core-pipeline.md).
- **Annotation.** CellTypist's immune models and step 04b's markers were built from whole cells, and
  RNA in the nucleus differs from RNA in the cytoplasm. Confirm labels with `top_markers.csv`.

---

## Optional: when split-pipe's cell calls are wrong

split-pipe decides which barcodes are cells by fitting a cutoff (transcripts
per barcode) **separately for every sample in every sublibrary**. That works
when a sample has one clean "knee" between empty barcodes and cells. When it
does not, for example when activated or dividing cells form a high-count tail
above the main population, the cutoff can land far above the main population.
`DGE_filtered` then keeps only a few hundred very large "cells" for that
sample and drops the rest without any warning.

**How to spot it.** Sublibraries are random splits of the same cells, so a
sample should get a similar cutoff, and a similar share of cells, in every
sublibrary. In a real mouse T-cell dataset (2026), the cutoffs ranged from
460 to 32,699 transcripts across samples. One sample had 2,956 cells in one
sublibrary and 261 in the next, even though its unfiltered barcodes were
almost the same in both. Signs to look for in `agg_sample_summary.csv`
(one per sublibrary run):

| Column | Warning sign |
|---|---|
| `cell_tscp_cutoff` | differs several-fold between sublibraries for the same sample, or one sample's cutoff is far above the others |
| `number_of_cells` | a sample drops to a few hundred in some sublibraries |
| `*_median_tscp_per_cell` | jumps to tens of thousands for those same samples |

**Check it with `00a_check_cell_calls.py`** (optional, takes seconds). It reads
the sublibrary runs, prints each sample's cutoff and cell count per run, counts
the barcodes that would pass a range of fixed cutoffs, flags unstable samples,
and draws barcode-rank plots:

```bash
python scripts/00a_check_cell_calls.py \
  --run-dir <split-pipe-out>/<sublib-1> <split-pipe-out>/<sublib-2> ... \
  --output-dir $OUT/cell_call_check/
```

It ends with `OK: ... DGE_filtered is fine to use` or `FLAGGED ...`.

**The fix: one cutoff for every sample.** Load `DGE_unfiltered` and give the
loader the same threshold for all samples:

```bash
python scripts/00_load_parse_output.py \
  --sample-dir <combined>/all-sample/DGE_unfiltered all-sample \
  --min-counts 500 --min-genes 300 \
  --sample-sheet $PROJECT/samples.csv --output $OUT/00_raw.h5ad
```

| Choosing the cutoff | |
|---|---|
| look at `barcode_rank.png` | pick a value below the main population of every sample and above the flat tail of empty barcodes |
| use the `barcodes_ge_*` columns | the numbers per sample should be close to the cells you loaded (target cells ÷ number of sublibraries) |
| keep it modest | this step only removes empty barcodes; step 01's QC still removes low-quality cells and doublets, per sample |

The loader prints the cells kept per sample and stores the cutoff in the
`.h5ad` (`uns["load_cell_cutoff"]`) so it is on record. On the cluster, set
`MIN_COUNTS` and `MIN_GENES` for `run_pipeline.slurm` (page 11).

Then check the high-count cells that split-pipe had kept. Are they a real cell
state (proliferation genes such as *Mki67*), or doublets? Step 01's doublet
scores and the cluster markers (step 04) show which.

---

## Your sample sheet

split-pipe knows only sample *names*. Anything else about a sample — patient,
time point, tissue, treatment — has to come from you, as a small CSV. Page 2 made
one to fill in. Open it with `nano $PROJECT/samples.csv`, replace everything, and save
(Ctrl+O, Enter, Ctrl+X). You can also make it in Excel, save as CSV, and upload it. Your
sample names are in split-pipe's output: `ls <combined>/` lists one folder per sample.

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

If your sheet names the sample column something else (e.g. `sample_id`), add
`--sample-sheet-key sample_id` instead of renaming it. If one whole sublibrary
turns out to carry a library-level artifact, drop it at load time with
`--exclude-sublibraries s1` (its `__sN` tag).

**Check:** `cat $PROJECT/samples.csv`. There should be one line per sample and no empty lines.

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

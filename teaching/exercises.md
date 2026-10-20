# Exercises — check your understanding

Each exercise uses data you can make yourself: PBMC3k or the synthetic
split-pipe data ([guide page 4](../docs/04-test-small-first.md)). Write your
answers down before reading the hint.

---

### 1. How strict should QC be?

Run QC on PBMC3k with `--mad 3` and again with `--mad 8`. How many cells are
flagged each time? Compare the two `qc_violin.png` plots. Which is more
defensible for this dataset, and why?

*Hint: are the extra cells removed at `--mad 3` a separate tail, or part of
the main population?*

---

### 2. Is that split real?

Cluster at resolution 0.2 and at 1.0. Pick one cluster from the 0.2 run and
find which clusters it becomes at 1.0. Look at their top marker genes. Is the
split a real biological distinction, or could it be noise?

*Hint: distinct, well-known marker genes on each side suggest biology; the same
genes at slightly different levels suggest noise.*

---

### 3. Fixed, or over-corrected?

With multi-sample data (the synthetic data has 3 sublibraries), cluster with and
without Harmony, using `--batch-key sublibrary`. Colour both UMAPs by
sublibrary and by sample. What changed? Then try `--batch-key sample` and
compare: what did correcting by sample remove?

*Hint: in split-pipe data, which differences between samples are technical,
and which are the biology you want to keep?*

---

### 4. Name the ambiguous cluster

Find a cluster `04b_better_annotation.py` labelled `Ambiguous (A / B)`. Look at
its top marker genes in `top_markers.csv`. Make a case for what it is. What
extra information — more specific markers, a higher resolution, a reference
dataset — would settle it?

---

### 5. Memory, measured

Submit a job with `--mem` set far too low for its data (for example 2 GB for
QC on the synthetic data). Watch it fail, then use `sacct` to find its
`MaxRSS`. How should that number guide your next request?

---

### 6. Choose the design

You have blood from 4 patients, each sampled before and 6 hours after
treatment. You want genes that change after treatment.

1. Is this a paired or an unpaired design?
2. Which column goes in `--covariate-cols` and `--covariates`?
3. Run it on the synthetic data (the same layout, with 2 patients) with and
   without the covariate. How many of the planted genes (`DEG1`–`DEG10`) are found
   each time, and are there false hits?
4. Why does the script refuse a covariate that varies between cells of the
   same sample?

---

### 7. Fair repertoire comparison

Sample A has 5,000 T cells with receptors, sample B has 500. A has 4,100
distinct clonotypes, B has 470. Is A's repertoire more diverse? What would you
compute instead, and which column of `repertoire_by_sample.csv` gives it?

*Hint: what happens to the number of distinct clonotypes if you sequence more
cells from the same person?*

---

### 8. Would you trust this?

Step 12 reports that 96% of TCR cells matched a WT cell, but only 30% of
paired-receptor cells sit in T-cell clusters; most are in a monocyte cluster.
List three possible explanations, and what you would check for each before
using any clonotype result.

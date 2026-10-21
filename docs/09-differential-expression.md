# 9. Differential expression (pseudobulk) — volcano plots and heatmap

Which genes change between conditions — disease vs. healthy, after vs. before
treatment — within each cell type.

---

## Why pseudobulk, not a per-cell test

Cells from the same person are not independent samples. A per-cell test treats
"3 donors × 3,000 cells" as 9,000 independent data points and makes results
look far more significant than they are. **Pseudobulk** first sums the raw
counts of each cell type within each sample, so every sample contributes one
data point per cell type, then tests those with **DESeq2**, a standard bulk
RNA-seq method. Each cell type is tested separately, so differences in
cell-type proportions are not mistaken for gene changes.

---

## Run it

```bash
python scripts/10_pseudobulk_prep.py --input 04b_annotated.h5ad --output-dir pseudobulk/ \
  --condition-col timepoint --sample-col sample --annotation-col cell_type_v2 \
  --covariate-cols patient
Rscript scripts/11_pseudobulk_deseq2.R --input-dir=pseudobulk/ --output-dir=pseudobulk_de/ \
  --condition-a=6h --condition-b=pre --covariates=patient
```

`--condition-col` is any column of your data (from the sample sheet, page 3);
`--condition-a` is compared against `--condition-b`, so a positive fold change
means higher in `a`.

---

## Choose the design that matches your experiment

Different experiments need different statistical models. You choose at run
time; nothing is hard-coded.

| Experiment | Options | Model |
|---|---|---|
| unpaired groups (different people in each group) | *(no covariates)* | `~ condition` |
| paired: the same patients before and after | `--covariate-cols patient` and `--covariates=patient` | `~ patient + condition` |
| a known processing batch spread across both groups | `--covariate-cols batch` and `--covariates=batch` | `~ batch + condition` |

**Why pairing matters.** People differ a lot from each other. Comparing each
patient's "after" with the same patient's "before" removes that
person-to-person noise. On the pipeline's synthetic test data, where the true
answer is known, the paired model found 62 of the planted changes instead of
38, with no false hits in either case.

Two rules the scripts enforce:

- **A covariate must describe the whole sample** (one patient, one batch), not
  individual cells. A per-cell column — such as `sublibrary`, which mixes cells
  of one sample — is refused.
- **It must not line up perfectly with the comparison.** If every patient
  appears in only one group, "patient" and "condition" are the same thing and
  cannot be separated. That cell type is skipped with an explanation instead
  of producing nonsense.

A cell type is also skipped if either condition has fewer than 2 samples
(`--min-replicates`): DESeq2 cannot estimate variability from one.

---

## Outputs

| File | Shows |
|---|---|
| `de_<cell_type>.csv` | DESeq2 results per cell type (log2 fold change, p, adjusted p) |
| `de_all_celltypes_combined.csv` | all cell types together, with the model used in a `design` column |
| `volcano_plots/volcano_<cell_type>.png` | fold change vs. significance; genes that are large *and* confident sit in the upper corners |
| `top_de_genes_heatmap.png` | the top genes across every pseudo-sample, grouped by cell type, coloured by z-score |

Defaults: significant means adjusted p < 0.05 and |log2 fold change| > 1
(`--padj-threshold`, `--lfc-threshold`). Grey heatmap cells are genes not
testable in that cell type (too few counts).

---

## With few samples, be careful

With 2 patients the model has very little information to estimate
variability. Results can still point somewhere useful, but treat them as
exploratory, and prefer genes that change consistently in every patient.

---

## Optional: which transcription factors drive the change (step 14)

A DE table lists hundreds of genes. Transcription-factor (TF) activity
summarises them: if most known targets of a TF move together, that TF is
probably driving them. Example: a cell state with the targets of CIITA (the
master regulator of HLA class II) all down is one that has switched off
antigen presentation, even if no single HLA gene passes the cut-off.

The method is the **univariate linear model (ULM)** from decoupler
(Badia-i-Mompel et al., *Bioinformatics Advances* 2022). The network is
**CollecTRI**, a curated, signed TF → target collection (Müller-Dott et al.,
*Nucleic Acids Research* 2023). For each TF, the gene-level values are
regressed on its target weights: +1 for an activated target, −1 for a
repressed one, 0 for all other genes. The **t-value** is the activity score.
Positive means the TF's targets move the way the TF would push them.
The script re-implements ULM in numpy because installing decoupler into this
environment would downgrade numpy.

**Run it.** First download the network once. This step needs internet, so run
it on a login node, since it is one small download, not compute:

```bash
python scripts/14a_fetch_collectri.py references/collectri_human.csv
```

Then score the DE results from step 11, the pseudobulk from step 10, or both.
This step is cheap, so run it inside a job like the other steps:

```bash
python scripts/14_tf_activity_ulm.py --net references/collectri_human.csv \
  --de-dir pseudobulk_de/ --pseudobulk-dir pseudobulk/ --output-dir tf_activity/
```

| Input | Uses | Gives |
|---|---|---|
| `--de-dir` | the DESeq2 `stat` column of every `de_<cell_type>.csv` | one score per TF per cell type: *is this TF more active in condition a than b?* |
| `--pseudobulk-dir` | log-CPM of every pseudo-sample, centred per gene | one score per TF per pseudo-sample, to plot along a timeline or across patients |

**Outputs:** `tf_activity_de.csv` and `tf_activity_pseudobulk.csv`, one row per TF and contrast or pseudo-sample: `score`, `pvalue`, `n_targets`. The pseudobulk table also carries the metadata columns.

**Reading it, carefully:**

- **The DE input matters most.** Scoring on `stat` already weighs each gene by its DE confidence, so it is the more reliable mode. The pseudobulk scores are relative to the mean of the pseudo-samples you passed in, so compare them within one run only.
- **Activity is inferred from targets, not measured.** It is not the TF's own mRNA, which is often too low to see. The TF may be inactive even when its targets move, because other TFs share those targets.
- **TFs with fewer than 5 measured targets are dropped** (`--min-targets`). Gene symbols must match the network: human CollecTRI with human symbols.
- **AP-1 and NF-κB appear as `AP1` and `NFKB`**, because their complexes are collapsed as decoupler does. In a pseudobulk run they can drop out if their targets are filtered.
- **The p-values come from a regression across all genes.** They are not a test across patients. With a few patients, look for TFs that change the same way in every patient.

`requests` (needed only by `14a`) is normally already in the `scrna`
environment through scanpy's dependencies; if not, `pip install requests`.
The synthetic check `tests/synthetic_tf_activity.py` plants known TF signal
and confirms the scores recover it (page 4).

---

Next: [10. T-cell receptors (Evercode TCR)](10-tcr.md)

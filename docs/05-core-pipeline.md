# 5. The core pipeline, step by step

Load → QC → normalize → cluster → annotate. Run each step by hand once, on a
small dataset (page 4), and look at what it writes before you ever chain them
into one job (page 11). The plots are where problems show up.

Every step reads the `.h5ad` file the previous one wrote. The *why* behind
each step is explained at more depth in
[teaching/concepts.md](../teaching/concepts.md); the validation numbers are in
[docs/reference/](reference/README.md).

---

## 00 — Load

```bash
python scripts/00_load_parse_output.py \
  --sample-dir <combined>/all-sample/DGE_filtered all-sample \
  --sample-sheet samples.csv --output 00_raw.h5ad
```

Turns split-pipe's (or converted Trailmaker's) `DGE_filtered/` files into one
AnnData `.h5ad` file, joins your sample sheet, and records each cell's
sublibrary. Page 3 explains the inputs. For a single sample without a sheet:
`--sample-dir DGE_filtered/ my_sample`.

**Check:** the printed cell and gene counts match what split-pipe reported;
`joined sample sheet columns: [...]` lists your columns.

---

## 01 — Quality control

```bash
python scripts/01_qc_filter.py --input 00_raw.h5ad --output 01_filtered.h5ad \
  --plot-dir qc_plots/ --species human --batch-key sample
```

| What it does | How |
|---|---|
| drops near-empty cells and genes | hard floors: ≥200 genes per cell (`--min-genes`), ≥3 cells per gene (`--min-cells`) |
| flags outliers | counts, genes and mitochondrial % that sit more than *k* median absolute deviations from the median (`--mad`, default 5) — thresholds adapt to *this* dataset instead of a fixed number |
| caps mitochondrial % | a hard ceiling as a backstop (`--mt-ceiling`, default 8%) |
| removes doublets | Scrublet, run **per sample** (`--batch-key sample`) |

`--species` sets the mitochondrial gene prefix (`MT-` human, `mt-` mouse,
`MT-` dog; override with `--mt-prefix`).

**Why per sample:** Scrublet assumes it is looking at one library. Pooling 24
donors into one call broke that assumption *and* ran out of memory on a
987,585-cell run even at 400 GB. Per sample fixed both.

**Check `qc_plots/qc_violin.png`:** do the flagged outliers look like a real
separate tail, or did the threshold cut through the middle of the main
population? If the latter, `--mad` is too strict for this dataset.

---

## 02 — Normalize, select genes, reduce dimensions

```bash
python scripts/02_normalize_pca.py --input 01_filtered.h5ad --output 02_pca.h5ad --plot-dir pca_plots/
```

| Step | What it does |
|---|---|
| normalize | rescales each cell to the same total (10,000), then log-transforms, so a deeply sequenced cell does not look like it expresses more of everything |
| highly variable genes | keeps the 2,000 genes that vary most (`--n-hvg`), chosen on raw counts (`seurat_v3` method) |
| scale | each gene to mean 0, variance 1, so high-expression genes do not dominate |
| PCA | compresses those genes into 50 components (`--n-pcs`) that later steps work from |

Raw counts are kept in `adata.layers["counts"]` — DE (page 9) needs them.

**Check `pca_plots/pca_elbow.png`:** variance explained should drop and then
flatten. Still falling steeply at 50 components? Real structure may be cut
off. Flat by component 10? The rest is mostly noise.

---

## 02b — Batch correction (multi-sample data)

Optional, but usually needed for real multi-sample runs — see
[page 6](06-batch-correction.md) before running it. If you run it, cluster on
its output in step 03 with `--use-rep X_pca_harmony`.

---

## 03 — Cluster

```bash
python scripts/03_cluster.py --input 02_pca.h5ad --output 03_clustered.h5ad \
  --plot-dir cluster_plots/ --resolutions 0.2,0.4,0.6,0.8,1.0
```

Builds a nearest-neighbour graph in PCA space, groups cells with the Leiden
algorithm at several **resolutions** (each saved as a column `leiden_r0.6`
and so on), and computes a UMAP for viewing. Low resolution merges related
populations; high resolution can split one real population into noise-driven
pieces. There is no single correct resolution, which is why the script scans
several.

**Check `cluster_plots/umap_by_resolution.png`:** which resolution separates
populations you expect to differ, without fragmenting the obvious ones?

---

## 04b — Annotate cell types

```bash
python scripts/04b_better_annotation.py --input 03_clustered.h5ad --output 04b_annotated.h5ad \
  --plot-dir marker_plots/ --resolution 0.6 --species human --markers-csv top_markers.csv
```

Two parts: a Wilcoxon test ranks the most distinctive genes of each cluster,
then each cluster is scored against marker sets for 12 immune cell types
(human, mouse or dog marker names). A label is assigned **only when the best
score is clearly ahead of the runner-up** (`--margin`, in standard deviations);
otherwise the cluster is labelled `Ambiguous (top1 / top2)` for you to inspect.
The label column is `cell_type_v2`.

**Why 04b replaced 04:** the older `04_markers_annotation.py` always picked the
highest score, even when every score was weak. On real data at 20 clusters it
confidently labelled clusters "Platelet" because that score was merely the
*least negative*. An honest "not sure between these two" beats a confident
wrong answer.

**Check:**
- `marker_plots/marker_dotplot_fine.png`: does each label's marker set light
  up in the cluster it names? `marker_plots/umap_annotated_v2.png` shows the
  labels on the UMAP.
- `top_markers.csv` for each `Ambiguous` cluster: can you name it from its
  top genes? Page 12 and the [exercises](../teaching/exercises.md) cover this.
- The marker sets are blood-focused. For tissue samples (tumour, stroma,
  endothelium), extend the marker dictionary in the script first — otherwise
  those cells can only come out `Ambiguous` or mislabelled.
- Mouse and dog marker sets are not yet validated against real non-human
  data.

---

Next: [6. Batch correction (Harmony)](06-batch-correction.md)

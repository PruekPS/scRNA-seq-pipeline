# 5. The core pipeline, step by step

Load → QC → normalize → cluster → annotate. Run each step by hand once, on a
small dataset (page 4), and look at what it writes before you ever chain them
into one job (page 11). The plots are where problems show up.

> **Before you start**
> - You are in the code folder, with `$PROJECT` and `$OUT` set ([page 1](01-getting-started.md)).
> - You are in an interactive session, with `conda activate scrna` ([page 4](04-test-small-first.md)).
>   For a big dataset, use the job on [page 11](11-running-on-slurm.md) instead.
> - Your sample sheet is ready ([page 3](03-your-data.md)).
>
> You do not need to create the folders for plots and results: every script creates its own
> output folders. Only `$PROJECT/logs` must exist beforehand, and page 2 made it.

Every step reads the `.h5ad` file the previous one wrote. The *why* behind
each step is explained at more depth in
[teaching/concepts.md](../teaching/concepts.md); the validation numbers are in
[docs/reference/](reference/README.md).

---

## 00 — Load

```bash
python scripts/00_load_parse_output.py \
  --sample-dir <combined>/all-sample/DGE_filtered all-sample \
  --sample-sheet $PROJECT/samples.csv --output $OUT/00_raw.h5ad
```

Turns split-pipe's (or converted Trailmaker's) `DGE_filtered/` files into one
AnnData `.h5ad` file, joins your sample sheet, and records each cell's
sublibrary. Page 3 explains the inputs. For a single sample without a sheet:
`--sample-dir /path/to/DGE_filtered/ my_sample`.

**Check:** the printed cell and gene counts match what split-pipe reported;
`joined sample sheet columns: [...]` lists your columns.

**Optional, when split-pipe's cell calls are wrong:** if a sample has only a
few hundred cells with very high counts, run
`scripts/00a_check_cell_calls.py` on the sublibrary runs. If it flags
samples, load `all-sample/DGE_unfiltered` with one cutoff for every sample
(`--min-counts N --min-genes M`). Page 3 explains how to spot the problem and
choose the cutoff.

---

## 01 — Quality control

```bash
python scripts/01_qc_filter.py --input $OUT/00_raw.h5ad --output $OUT/01_filtered.h5ad \
  --plot-dir $OUT/qc_plots/ --species human --batch-key sample
```

| What it does | How |
|---|---|
| drops near-empty cells and genes | hard floors: ≥200 genes per cell (`--min-genes`), ≥3 cells per gene (`--min-cells`) |
| flags outliers | counts and genes more than *k* median absolute deviations from the median, either way, and mitochondrial % above it (`--mad`, default 5) — thresholds adapt to *this* dataset instead of a fixed number |
| caps mitochondrial % | a hard ceiling as a backstop (`--mt-ceiling`, default 8%) |
| removes doublets | Scrublet, run **per sample** (`--batch-key sample`) |

`--species` sets the mitochondrial gene prefix (`MT-` human, `mt-` mouse,
`MT-` dog; override with `--mt-prefix`).

**Why per sample:** Scrublet assumes it is looking at one library. Pooling 24
donors into one call broke that assumption *and* ran out of memory on a
987,585-cell run even at 400 GB. Per sample fixed both.

**A warning to take seriously:** if Scrublet calls almost no doublets, its automatic threshold may
have failed (this happened on a real run), and the script says so. Doublets are then still in your
data: look for clusters with very high counts or mixed markers at step 04.

**Check `$OUT/qc_plots/qc_violin.png`:** do the flagged outliers look like a real
separate tail, or did the threshold cut through the middle of the main
population? If the latter, `--mad` is too strict for this dataset.

---

## 02 — Normalize, select genes, reduce dimensions

```bash
python scripts/02_normalize_pca.py --input $OUT/01_filtered.h5ad --output $OUT/02_pca.h5ad --plot-dir $OUT/pca_plots/
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
python scripts/03_cluster.py --input $OUT/02_pca.h5ad --output $OUT/03_clustered.h5ad \
  --plot-dir $OUT/cluster_plots/ --resolutions 0.2,0.4,0.6,0.8,1.0
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
python scripts/04b_better_annotation.py --input $OUT/03_clustered.h5ad --output $OUT/04b_annotated.h5ad \
  --plot-dir $OUT/marker_plots/ --resolution 0.6 --species human --markers-csv $OUT/top_markers.csv
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
  data. On a real mouse lung dataset (2026) the blood-focused sets labelled
  lung cells as "Platelet" and "HSPC"; label tissue samples by reviewing
  `top_markers.csv` yourself, or extend the dictionary.

---

## 04c — Reference-based labels with CellTypist (optional, human immune cells)

```bash
python scripts/04c_celltypist_annotation.py --input $OUT/04b_annotated.h5ad --output $OUT/04c_celltypist.h5ad \
  --model $PROJECT/references/celltypist/data/models/Immune_All_Low.pkl --over-clustering leiden_r0.6 \
  --plot-dir $OUT/celltypist_plots
```

04b scores clusters against a hand-picked marker list. That breaks down when a
whole population shifts away from textbook markers (a strong interferon or
vaccine response, disease-skewed monocytes): on a real PBMC run most 04b labels
came out `Ambiguous`. CellTypist (Domínguez Conde et al., *Science* 2022)
instead applies models trained on large curated immune atlases to every gene of
every cell.

| Output column | Meaning |
|---|---|
| `celltypist_label` | per-cell prediction |
| `celltypist_majority` | majority vote within each cluster of `--over-clustering`, the label to use by default |
| `celltypist_conf` | confidence of the per-cell label |

The models are **human immune** atlases. Download them once, before the command above (the
script itself makes no network calls). Make the folder first, because the download needs it:

```bash
mkdir -p $PROJECT/references/celltypist
CELLTYPIST_FOLDER=$PROJECT/references/celltypist python -c "import celltypist; celltypist.models.download_models(model=['Immune_All_Low.pkl','Immune_All_High.pkl'])"
```

CellTypist saves them in a `data/models` subfolder of the folder you name, so check with
`ls $PROJECT/references/celltypist/data/models/`.

Run a second model for broad labels. Step 05 needs these:

```bash
python scripts/04c_celltypist_annotation.py --input $OUT/04c_celltypist.h5ad --output $OUT/04c_celltypist.h5ad \
  --model $PROJECT/references/celltypist/data/models/Immune_All_High.pkl --over-clustering leiden_r0.6 \
  --prefix celltypist_high --plot-dir $OUT/celltypist_plots
```

---

## 05 / 05b — Refine labels per lineage, flag doublets, apply a curated label map (optional)

```bash
python scripts/05_refine_annotation.py --input $OUT/04c_celltypist.h5ad --outdir $OUT/refined/ \
  --broad-col celltypist_high_majority --fine-col celltypist_label --batch-key sample
python scripts/05b_apply_label_map.py --input $OUT/refined/refined.h5ad --output $OUT/final_labels.h5ad \
  --subcluster-map $PROJECT/my_subcluster_labels.csv --rescued-map $PROJECT/my_rescued_labels.csv
```

A majority vote over ~25 whole-dataset clusters merges real subtypes, gives
rare cells (plasma cells, progenitors) their neighbours' label, and labels
doublet pockets as one of their parents. **05** re-clusters each lineage (T/NK,
B/plasma, myeloid, other) on its own variable genes, labels each subcluster by
its majority CellTypist call, flags subclusters that score high for *another*
lineage's markers *and* have high doublet scores as `Doublet (A+B)`, and keeps
confident rare-type calls. It writes `subcluster_summary.csv` (size, majority
label, agreement, doublet and cross-lineage scores, top markers) for review.

**05b** then applies **your** decisions: a CSV mapping every subcluster to a
final label and broad lineage, so curated labels live in a reviewable file,
not in code. It refuses a map that does not cover every subcluster.
`configs/label_maps/example_pbmc_*.csv` are the maps from a real PBMC run, as
examples of the format and level of evidence. Write your own for your data.

---

## 05c — Score published gene programs (optional)

```bash
python scripts/05c_score_gene_programs.py --input $OUT/04b_annotated.h5ad --output $OUT/05c_programs.h5ad \
  --species mouse --summary-csv $OUT/programs_per_sample.csv --group-col condition
# Treg programs, Treg cells only:
python scripts/05c_score_gene_programs.py --input $OUT/t_cells.h5ad --output $OUT/t_programs.h5ad --species mouse \
  --subset-col cell_type --subset-values Treg --summary-csv $OUT/treg_programs_per_sample.csv
```

A label says what a cell *is*; a program score says how strongly it runs a
function. Each program is scored per cell (mean of its genes minus an
expression-matched random control set, `scanpy.tl.score_genes`), stored as
`prog_<program>`, and averaged per sample.

Built-in T-cell programs (`--programs builtin`, mouse or human symbols):
naive, effector/memory, Th1, Th2, Th17, Tfh, Tr1 (IL-10), cytotoxic,
tissue-resident, IFN response, exhaustion, progenitor-exhausted, two
activation sets, Treg resting / effector / tissue-tumour / Tfr, proliferation,
and a dissociation-stress set for QC. Or give your own CSV (`program,gene`).

| Use with care | Why |
|---|---|
| compare per sample (mouse/patient), not per cell | cells of one sample are not independent |
| check `Dissociation stress (QC)` first | heat-shock and immediate-early genes rise with tissue handling and inflate activation scores |
| a baseline "AIM" score is not an AIM assay | AIM needs antigen restimulation; `AIM-like activation` includes immediate-early genes (Cd69, Nr4a1, Egr2), so use `AIM core (no IEG)` |
| scores are relative | compare the same program across samples, not different programs with each other |

The gene choices follow T-cell reference atlases and reviews, for example
ProjecTILs (Andreatta et al., *Nat Commun* 2021,
[doi:10.1038/s41467-021-23324-4](https://doi.org/10.1038/s41467-021-23324-4)),
the pan-cancer T-cell landscape (Zheng et al., *Science* 2021,
[doi:10.1126/science.abe6474](https://doi.org/10.1126/science.abe6474)),
tissue Treg adaptation (Miragaia et al., *Immunity* 2019,
[doi:10.1016/j.immuni.2019.01.001](https://doi.org/10.1016/j.immuni.2019.01.001))
and Tr1 cells (Gagliani et al., *Nat Med* 2013,
[doi:10.1038/nm.3179](https://doi.org/10.1038/nm.3179)).

---

## 05d — Find ambient and stress genes (optional)

Two things can look like biology. **Ambient RNA:** damaged cells release RNA that ends up with other
cells, so a T cell can seem to express a macrophage gene. **Stress:** heat-shock and early-response
genes rise with handling and time on ice, so samples collected on different days can differ for that
reason alone. Nothing is removed by this step; it gives you lists and flags.

**Part 1, ambient genes.** You need a cluster of low-quality cells (low counts, no clear identity),
and a clean population in which those genes should not be expressed. Replace the cluster numbers
with your own, from `$OUT/cluster_plots` and `top_markers.csv`:

```bash
python scripts/05d_flag_ambient_stress.py --output-dir $OUT/ambient --input $OUT/04b_annotated.h5ad \
  --low-quality-col leiden_r0.6 --low-quality-values 2 \
  --reference-col leiden_r0.6 --reference-values 4 8 9
```

Read `$OUT/ambient/ambient_genes.csv`. To stop those genes shaping the clusters, go back to step 02
and add `--exclude-hvg-genes $OUT/ambient/ambient_genes.csv`, then redo the later steps.

**Part 2, flag a DE result** (after [page 9](09-differential-expression.md)):

```bash
python scripts/05d_flag_ambient_stress.py --output-dir $OUT/ambient \
  --de-table $OUT/pseudobulk_de/de_all_celltypes_combined.csv --ambient-genes $OUT/ambient/ambient_genes.csv
```

`de_flagged.csv` gains two columns, `stress_gene` and `possible_ambient`. `de_summary_counts.csv` counts
the hits left without them. If most hits in a comparison are stress genes, ask whether the groups
were handled differently before reading it as biology.

---

Next: [6. Batch correction (Harmony)](06-batch-correction.md)

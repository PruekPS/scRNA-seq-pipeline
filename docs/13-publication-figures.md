# 13. Publication figures in Nature style

Every figure this pipeline writes follows the same style, set in one place
(`scripts/figstyle.py`, and `scripts/figstyle.R` for R): the look of a figure in a
Nature journal, ready to place in a manuscript or a slide. Two extra steps build
figures meant for a paper.

> **Before you start**
> - You are in the code folder with `$PROJECT` and `$OUT` set (page 1).
> - The `scrna` environment is active, in an interactive session or a job (page 4).
> - For step 16 you have an annotated `.h5ad`, with a UMAP (pages 5–6).

---

## What "Nature style" means here

| Rule | Value |
|---|---|
| Width | 89 mm (one column), 120 mm (1.5 columns), 183 mm (two columns) |
| Font | Arial or Helvetica, 5–7 pt; panel letters bold lowercase, 8 pt |
| Lines | 0.5 pt |
| Text in figures | no titles or explanatory sentences; that belongs in the legend |
| Colours | colour-blind safe (Okabe-Ito first); one meaning for one colour |
| Gene names | in italics |
| Cells in scatter plots | stored as an image inside the vector file, so it stays small and its text stays editable |

Each figure is saved as a set:

| File | Use |
|---|---|
| `name.pdf`, `name.svg` | vector, with editable text: open in Illustrator or Inkscape |
| `name.png`, `name.tiff` | 600 dpi images for slides and submission |
| `name_legend.txt` | a draft legend: what is plotted, n, and the test |
| `name_source_data.csv` | the numbers behind the figure, where a figure has them; journals ask for this |

**Font note.** Arial is rarely installed on a cluster. The figures then use Liberation Sans,
which has the same letter widths; `pipeline/check_setup.sh` tells you which font you have.
For a final version with real Arial, run step 16 or 17 on a computer that has it, or ask your
cluster to install it. Fine for checking results either way.

**Faster test runs.** Writing four formats takes time on large data. For a test, write only
PNGs by adding this before your command, or to your job:

```bash
export FIG_FORMATS=png
```

---

## Step 16: an overview figure of your dataset

Builds one double-column figure: **a** UMAP labelled by cell type; **b** composition per
sample; **c** UMAP split by group (with `--group-col`); **d** marker genes on the UMAP (with
`--genes`). The folder is created for you.

```bash
python scripts/16_atlas_figure.py --input $OUT/04b_annotated.h5ad --output-dir $PROJECT/figures/atlas \
  --label-col cell_type_v2 --sample-col sample --group-col condition \
  --genes CD3E CD8A FOXP3 MS4A1 LYZ NKG7
```

**You should see:** `Wrote .../Fig_atlas.*`. Open `$PROJECT/figures/atlas/Fig_atlas.png`. It also
writes `panels/` (each panel alone) and `source_data/`.

**To rename, reorder and colour the cell types,** make a CSV with `nano $PROJECT/labels.csv`. Rows
are the display order. `display` and `color` are optional:

```csv
label,display,color
CD4 T cell,CD4 T cell,#2166ac
CD8 T cell,CD8 T cell,#67a9cf
B cell,B cell,#762a83
CD14+ Monocyte,Classical monocyte,#e08214
```

and add `--labels-csv $PROJECT/labels.csv`. Keep the cell types of one lineage in one colour
family. Leave out cell types with `--exclude-labels "Doublet"`.

**Check before you use it:**

- **Look at the labels on panel a.** They are placed automatically and a few may overlap. Short
  names in `labels.csv` usually fix it; otherwise, move them in Illustrator using the PDF.
- **Two samples that are not independent must not count twice.** For example, if your records or
  the checks on [page 3](03-your-data.md) show two samples are copies:
  `--pool-samples 'S1+S2=S1 and S2'`.
- **Make the UMAP from the cells you show.** If you removed low-quality cells or doublets, recompute
  the UMAP after removing them. Otherwise the map keeps their shape.

---

## Step 17: compare a value between groups

One panel per value: a bar for the group mean, a dot for each sample. **The sample is the
unit, never the cell**, because cells of one sample are not independent.

From an annotated `.h5ad`, the percent of cells of each type:

```bash
python scripts/17_group_figure.py --input $OUT/04b_annotated.h5ad \
  --sample-col sample --group-col condition --label-col cell_type_v2 \
  --labels "B cell" "Treg" --ylabel "Cells (%)" --output $PROJECT/figures/composition
```

Or from a table you made, one row per sample:

```csv
sample,condition,Effector Treg,Tissue Treg
M1,control,60.1,12.4
M2,control,55.0,15.2
M3,treated,31.2,40.0
M4,treated,28.7,38.5
```

```bash
python scripts/17_group_figure.py --values-csv $PROJECT/treg_states.csv \
  --group-col condition --values "Effector Treg" "Tissue Treg" --ylabel "Cells (% of Treg)" \
  --group-order control treated --test --output $PROJECT/figures/treg_states
```

**`--test`** adds a test where one is possible, and records it in the legend file:

| Groups | Test |
|---|---|
| 2 | exact two-sided Mann-Whitney U |
| 3 or more | Kruskal-Wallis with Dunn's comparisons, Bonferroni-adjusted; only pairs with P < 0.05 are drawn |
| any group with one sample | that group is left out of the test |

With few samples a test cannot reach significance however large the effect: an exact test on
2 vs 2 samples cannot go below P = 0.33. Show the dots and describe the effect, and say so in the
legend.

---

## Using the style in your own plots

A script of your own, in the `scripts/` folder, can use the same style:

```python
import figstyle

figstyle.apply_style()
fig, axes = figstyle.subplots(1, 2, w_mm=figstyle.SINGLE, h_mm=45)
axes[0, 0].plot([1, 2, 3], [1, 4, 9])
axes[0, 0].set_xlabel("Dose")
figstyle.save_figure(fig, "my_figures/dose_response",
                     legend="Response against dose. n = 3 mice per dose.")
```

---

That is the end of the guide. For the concepts behind each step, and exercises, see
[teaching/](../teaching/README.md). For the design decisions and validation reports, see
[docs/reference/](reference/README.md).

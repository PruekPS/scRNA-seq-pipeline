# 6. Batch correction (Harmony)

A **batch effect** is a technical difference between groups of cells — a
different processing run, sequencing lane or reagent lot — that is large
enough to make cells cluster by *which batch they came from* instead of by
cell type. Harmony adjusts the PCA embedding to remove that technical signal
while keeping real biological differences.

---

## Run it

```bash
python scripts/02b_harmony_integration.py --input 02_pca.h5ad --output 02b_harmony.h5ad \
  --plot-dir harmony_plots/ --batch-key sublibrary
python scripts/03_cluster.py --input 02b_harmony.h5ad --output 03_clustered.h5ad \
  --plot-dir cluster_plots/ --use-rep X_pca_harmony
```

The corrected embedding is added as `X_pca_harmony`; the original PCA is kept.
The script refuses to run on single-batch data, where there is nothing to
correct.

---

## Which column is the batch?

This is the decision that matters, and it depends on how the experiment was
done.

| Your data | Batch key | Why |
|---|---|---|
| split-pipe, several sublibraries | `sublibrary` | every sublibrary holds cells from every sample, so sublibrary is the technical batch |
| separate runs per donor (or Trailmaker exports per sample) | `sample` or `donor` | each sample was processed on its own, so sample *is* the batch |
| one sample, one run | none — skip this step | nothing to correct |

**Never use your biological comparison as the batch key.** Correcting by
`timepoint`, `condition` or `tissue` would remove exactly the differences you
want to study. For the same reason, correcting by `sample` in split-pipe data
also removes real differences between samples.

For trajectory work across tissues (bone marrow → blood → tumour), keep the
correction to technical batches only; the tissue-driven differences *are* the
trajectory ([page 8](08-trajectory.md)).

---

## Check it worked

`harmony_plots/pca_before_after_harmony.png` and a UMAP coloured by the batch
column, before and after:

- **Before:** do batches form their own islands *within* what should be one
  cell type?
- **After:** that separation should mostly disappear, while distinct cell
  types stay distinct.

If distinct cell types merge after correction, it over-corrected — check that
the batch key is really technical.

---

## A real problem this step hit

`harmonypy` 2.0.0 (conda's newest) returned a malformed array on a
907,696-cell dataset, crashing the step. The environment pins `0.0.10`, which
works; the script warns if a different version is installed.

On that real run (24 donors, batch key `sample`, separate libraries), Harmony
converged in 6 iterations — see
[reference/real-data-validation.md](reference/real-data-validation.md).

---

Next: [7. Cell-cell communication (CellChat)](07-cell-communication.md)

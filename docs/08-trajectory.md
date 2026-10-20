# 8. Trajectories and fate mapping (PAGA, CytoTRACE, CellRank)

Clustering shows *which* cell states exist. It cannot show order or direction
— which cells are becoming which. If your question is about differentiation
(for example progenitors maturing into tumour-infiltrating immune cells),
these two optional steps add that.

> **Status: built, not yet run on real data.** The code is checked against
> CellRank's current documentation, but no suitable dataset has existed to run
> it on. The full design, including what the experiment must contain, is in
> [reference/trajectory-design.md](reference/trajectory-design.md) — read it
> before planning samples.

---

## What each method does

| Method | Gives you | Needs |
|---|---|---|
| **PAGA** | a coarse map of which cell types are plausibly connected | the clustered, annotated object |
| **Diffusion pseudotime** | an ordering of cells by distance from a chosen starting population | a root cluster you choose |
| **CytoTRACE** | a per-cell differentiation-potential score: cells expressing more distinct genes tend to be less differentiated | nothing extra |
| **CellRank 2** | for each cell, the probability of ending in each terminal state, plus the genes driving each fate | the above, combined |

True **RNA velocity** is *not* used: it needs separate spliced and unspliced
read counts, which standard split-pipe/Trailmaker output does not provide.

---

## Run it — in the `scrna-traj` environment

These steps need the separate Python 3.12 environment (page 2):

```bash
conda activate scrna-traj
python scripts/08_trajectory_paga.py --input 04b_annotated.h5ad --output 08_paga.h5ad \
  --plot-dir paga_plots/ --annotation-col cell_type_v2 --root-cluster HSPC
python scripts/09_cellrank_fate_mapping.py --input 08_paga.h5ad --output 09_cellrank.h5ad \
  --plot-dir cellrank_plots/ --drivers-csv lineage_drivers.csv \
  --annotation-col cell_type_v2 --n-states 8 --use-dpt-pseudotime
```

On the cluster, `pipeline/run_trajectory.slurm` runs both (page 11).

---

## The two choices that decide the result

**`--root-cluster`** — the starting population. Pseudotime is distance *from*
it, so it sets the direction of everything downstream. Choose the genuinely
least-differentiated population (04b has an `HSPC` category for this), not
whatever is convenient. Within it, the root *cell* is picked by a documented
rule; `--root-cell-index` overrides it.

**`--n-states`** — how many coarse states CellRank looks for. It should
comfortably exceed the number of start and end populations you expect. Check
`cellrank_plots/eigenvalue_spectrum.png`: a clear gap between eigenvalues
suggests a natural number of states.

---

## Outputs

| File | Shows |
|---|---|
| `paga_plots/paga_graph.png` | connectivity between cell types |
| `paga_plots/umap_dpt_pseudotime.png` | pseudotime on the UMAP |
| `cellrank_plots/macrostates_terminal.png`, `macrostates_initial.png` | the end and start states CellRank found |
| `cellrank_plots/fate_probabilities.png`, `circular_projection.png` | each cell's fate probabilities |
| `lineage_drivers.csv` | genes most associated with each fate |
| `09_cellrank.gpcca.pickle` | the fitted model, so plots can be redrawn without recomputing |

---

## What it can and cannot claim

- It is inference from **one snapshot**, not lineage tracing: fate
  probabilities are hypotheses about likely paths.
- It assumes roughly **one-way** differentiation; plastic immune states (for
  example macrophage polarisation) can be misrepresented.
- A bone-marrow-to-tumour trajectory only means something if **both tissues
  were actually sampled**, ideally from the same subjects. Put a `tissue` (or
  `compartment`) column in your sample sheet (page 3).

---

Next: [9. Differential expression (pseudobulk)](09-differential-expression.md)

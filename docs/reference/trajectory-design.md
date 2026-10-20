# Trajectory design — PAGA, CytoTRACE and CellRank

> **Read this first.** This is a design document, not a results report: no
> bone-marrow-plus-tumour dataset existed to run it on. It documents the method,
> the code (`scripts/08_trajectory_paga.py`, `scripts/09_cellrank_fate_mapping.py`)
> and what the experiment must contain for the analysis to answer a real
> question. Treat every number in "Experimental design" as adaptable guidance,
> not a requirement.

Scoped to one question: tracking differentiation from bone-marrow progenitors
to mature immune cells infiltrating the tumour microenvironment (TME), in pooled
Evercode samples spanning both tissues.

---

## Why trajectory tools, not just clustering

The core pipeline answers "which cell types are present". It does not answer
"which cells are becoming which, and how far along". Tracking progenitors as
they become tumour-infiltrating immune cells is a question of *order and
direction* — exactly what clustering discards.

There is precedent: a pan-cancer single-cell atlas of tumour-infiltrating
myeloid cells modelled monocyte-to-tumour-associated-macrophage trajectories
across cancer types and found consistent trajectories despite large
patient-to-patient and cancer-to-cancer differences [19] — useful prior art
for what the trajectory should look like and how much variability to expect.

---

## Experimental design

### Both ends of the process, not just more cells

No trajectory method invents biological anchors that are not in the data. The
pooled run needs cells from **both ends, ideally from the same subjects**:

- **Bone marrow** — the progenitor compartment, the trajectory's root.
- **Tumour tissue** — the mature, terminal compartment.
- **Blood**, if feasible — many myeloid trajectories pass through a circulating
  monocyte state; sampling it directly makes the trajectory much easier to
  validate.

Each cell's tissue of origin must reach the analysis. Put a `tissue` (or
`compartment`) column in your **sample sheet** (`00_load_parse_output.py
--sample-sheet`; [guide page 3](../03-your-data.md)), with values such as
`bone_marrow` / `blood` / `tumour`. Plan this into the sample and well layout
before sequencing.

### How many subjects

There is no single correct number; it depends on effect size and on
patient-to-patient variation in the cancer type. As orientation only:

- Single-lab scRNA-seq studies comparing groups commonly use about **3–6
  biological replicates per group** as a practical minimum for a first pass.
- Trajectories benefit from *more*, because the most informative cells — those
  mid-transition — are usually the rarest.
- The pan-cancer study above pooled many patients and cancer types; it shows a
  well-powered version of the question, not a target for a first run.

Treat the first pooled run as a pilot: check the root is genuinely
progenitor-like, the terminal states genuinely mature/TME-associated, and fate
probabilities vary smoothly rather than like noise — before sizing a larger
cohort.

### Do not batch-correct away the signal

Batch correction should remove only **technical** batches (for split-pipe data,
`sublibrary`; [guide page 6](../06-batch-correction.md)). **Never** use
`tissue`/`compartment` as the batch key: that would erase the very differences
this analysis looks for.

---

## Method

### PAGA — coarse connectivity

PAGA [14] builds a graph whose nodes are clusters and whose edge weights say
how well connected two populations are. It summarises which transitions are
*plausible* — not a direction, and not proof that a transition happens.

### CytoTRACE — differentiation potential

CytoTRACE's [15] finding: the number of distinct genes a cell expresses tracks
how undifferentiated it is — progenitors express a broader repertoire. That
gives a per-cell potential score without spliced/unspliced counts, time
courses or lineage tracing. CellRank 2 reimplements it in Python
(`CytoTRACEKernel`), so no extra package is needed.

### CellRank 2 — kernels, GPCCA, fate probabilities

CellRank 2 [18] ([code](https://github.com/theislab/cellrank),
[what is new in v2](https://cellrank.readthedocs.io/en/latest/about/version2.html))
splits the problem in two:

- **Kernels** turn data into a cell-to-cell transition matrix. Six exist:
  `VelocityKernel`, `ConnectivityKernel`, `PseudotimeKernel`,
  `CytoTRACEKernel`, `RealTimeKernel`, `PrecomputedKernel`. This pipeline
  combines `CytoTRACEKernel` (direction from differentiation potential) with
  `ConnectivityKernel` (stability against noise), and optionally
  `PseudotimeKernel` fed by step 08's diffusion pseudotime, as a weighted sum
  (e.g. `0.7 * cytotrace + 0.3 * connectivity`).
- **Estimators** analyse the matrix. **GPCCA** finds coarse "macrostates",
  classifies some as initial and terminal, and computes each cell's **fate
  probability** — the chance of reaching each terminal state — plus **lineage
  driver genes**.

The code follows CellRank's current API, checked against its documentation and
tutorials: `compute_schur()` → `compute_macrostates()` →
`predict_terminal_states()` / `predict_initial_states()` →
`compute_fate_probabilities()` → `compute_lineage_drivers()`. (The older
`.fit()` shortcut still appears in some quick-start examples but is deprecated.)

How to run it, and the outputs: [guide page 8](../08-trajectory.md).

---

## Why not RNA velocity (yet)

True RNA velocity [16] compares spliced and unspliced transcripts to predict
each cell's near-term direction — a more direct signal than CytoTRACE. It is
not used because **standard Parse split-pipe/Trailmaker output does not
separate spliced from unspliced reads**. It would need the raw FASTQ files
re-processed with a velocity-aware aligner (for example STARsolo's Velocyto
mode) — a raw-data change, not a downstream addition, and not yet confirmed to
be possible with split-pipe's output. If added later, `VelocityKernel` slots
into the same CellRank workflow; only the kernel changes.

---

## What it can and cannot claim

- **One snapshot, not lineage tracing.** Fate probabilities are inferences from
  the geometry of one time point.
- **Roughly one-way differentiation is assumed.** Macrophage polarisation and
  other plastic myeloid states can be misrepresented as one-way.
- **The root is a modelling choice, not a measurement** — every pseudotime and
  direction result depends on it.
- **The number of macrostates is a setting to inspect** (eigenvalue spectrum),
  not ground truth.
- **A bone-marrow-to-TME claim needs both tissues sampled.** On one tissue the
  "trajectory" only describes variation within that tissue.

---

References: [14] Wolf et al. 2019 (PAGA); [15] Gulati et al. 2020
(CytoTRACE); [16] La Manno et al. 2018 (RNA velocity); [18] Weiler et al. 2024
(CellRank 2); [19] Cheng et al. 2021 — full entries in
[design decisions](design-decisions.md#references).
Previous: [Cell-cell communication](cellchat.md) · Back to [reference index](README.md)

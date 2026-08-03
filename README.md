# scRNA-seq Analysis Pipeline

A single-cell RNA-seq pipeline for Parse Biosciences Evercode data (whole
transcriptome now, TCR/VDJ in a later phase): raw sequencer output →
Parse demultiplexing → Python/scanpy QC/clustering/annotation → optional
Harmony batch integration → optional R/CellChat cell-cell communication
analysis, built to run at scale on [HiPerGator](https://www.rc.ufl.edu/)
(UF Research Computing) or locally for smaller datasets.

**Full scope, architecture, and rationale:** [`docs/pipeline_overview.html`](docs/pipeline_overview.html)
— open it in a browser for the pipeline flowchart, the decision log behind
each architecture choice, and QC metrics with their biological rationale,
each backed by a peer-reviewed citation (Nature Communications, Cell,
Science, Bioinformatics, Briefings in Bioinformatics). Companion deep-dive
reports for each stage (QC, normalization/PCA, clustering, annotation, and
a real-data validation run) are in [`docs/`](docs/), each with a real,
computed comparison against an alternative method.

**New to this pipeline?** [`teaching/student_hpc_tutorial.html`](teaching/student_hpc_tutorial.html)
is a step-by-step guide to running it yourself on HiPerGator — no coding
assistant required — written to teach the reasoning behind each stage, not
just the commands.

## Status

Core pipeline (QC → normalize/PCA → cluster → annotate) is validated
against real Parse Evercode data at full scale (~907,700 cells post-QC
from the real ~987,600-cell dataset) on HiPerGator — see
[`docs/real_data_deep_dive.html`](docs/real_data_deep_dive.html) for the
methodology (validated first on a 50,000-cell subset, then at full scale).
Harmony batch integration, the finer confidence-gated multi-species
(human/mouse/canine) annotation step, and CellChat cell-cell communication
are all now validated end-to-end at full scale: on the real, annotated
694,803-cell dataset, CellChat found 291 significant ligand-receptor
interactions across 19 signaling pathways. The full chain — QC → Harmony →
cluster → annotate → CellChat — runs clean from raw data to results.

Two further optional stages are built but not yet run against real data
(no bundled dataset has a real trajectory or condition comparison to
validate against): trajectory/fate mapping (PAGA + CytoTRACE + CellRank,
deliberately no RNA velocity — see
[`docs/trajectory_cellrank_overview.html`](docs/trajectory_cellrank_overview.html))
and pseudobulk differential expression with DESeq2, volcano plots, and a
top-DE-gene heatmap, dataset-agnostic by design (any `--condition-col`
supplied at run time).

## Why this design

- **Parse's split-pool barcoding** ([Rosenberg et al., *Science* 2018](https://doi.org/10.1126/science.aam8999))
  requires vendor-specific demultiplexing before any custom code can run —
  either through **Trailmaker** (cloud, no-code) or **split-pipe** (CLI,
  same logic, run as a HiPerGator batch job — see
  [`pipeline/split-pipe/`](pipeline/split-pipe/)).
- QC thresholds are **adaptive (MAD-based)** rather than fixed cutoffs,
  following [OSCA](https://bioconductor.org/books/3.13/OSCA.basic/) and the
  SCTK-QC pipeline ([Hong et al., *Nat Commun* 2022](https://doi.org/10.1038/s41467-022-29212-9)).
- Split-pool doublet and ambient-RNA mechanisms differ from droplet-based
  platforms (10x-style) — the QC section spells out why generic droplet
  assumptions don't transfer directly. Doublet detection (Scrublet) runs
  per-sample rather than pooled across donors, matching Scrublet's own
  intended usage ([Wolock, Lopez & Klein, *Cell Systems* 2019](https://doi.org/10.1016/j.cels.2018.11.005)).
- Batch integration (Harmony) and cell-cell communication (CellChat) are
  kept as separate, optional downstream steps rather than folded into the
  core chain, so the core QC→cluster→annotate path stays simple for
  single-sample data and the heavier multi-sample analyses are opt-in.

## Repository layout

```
docs/                 living pipeline documentation + per-stage deep-dive reports
teaching/              student-facing, no-assistant-required tutorial (student_hpc_tutorial.html)
scripts/               00-04: core pipeline (load, QC, normalize/PCA, cluster, annotate)
                        02b/04b: optional Harmony integration / finer multi-species annotation
                        06-07: CellChat export (Python) + analysis (R)
                        08-09: PAGA + diffusion pseudotime / CellRank fate mapping (no RNA velocity)
                        10-11: pseudobulk aggregation (Python) + DESeq2 DE, volcano, heatmap (R)
pipeline/               SLURM drivers (run_pipeline.slurm, run_downstream.slurm,
                        run_trajectory.slurm, run_pseudobulk_de.slurm) +
                        pipeline/split-pipe/ (HiPerGator setup for Parse's CLI demux)
envs/environment.yml    reproducible conda env spec (Python/scanpy + R/Seurat/CellChat/DESeq2 stack)
tests/                  smoke-test convention and gitignored real/public test data
```

## Stack

Python (scanpy/AnnData) end-to-end for the core pipeline, R (Seurat) for
the Trailmaker-object entry point and CellChat analysis, HiPerGator/SLURM
for compute at scale.

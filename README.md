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

## Status

Core pipeline (QC → normalize/PCA → cluster → annotate) is built and
validated against real Parse Evercode data — see
[`docs/real_data_deep_dive.html`](docs/real_data_deep_dive.html) for the
50,000-cell real-subset validation run. A full-scale (~987,000-cell) run
on HiPerGator is in progress. Harmony batch integration, a finer
multi-species (human/mouse/canine) annotation step, and a CellChat
cell-cell communication module have been added and are being validated
against real data next — see `scripts/` below for what's stable vs. new.

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
scripts/               00-04: core pipeline (load, QC, normalize/PCA, cluster, annotate)
                        02b/04b: optional Harmony integration / finer multi-species annotation
                        06-07: CellChat export (Python) + analysis (R)
pipeline/               SLURM drivers (run_pipeline.slurm, run_downstream.slurm) +
                        pipeline/split-pipe/ (HiPerGator setup for Parse's CLI demux)
envs/environment.yml    reproducible conda env spec (Python/scanpy + R/Seurat/CellChat stack)
tests/                  smoke-test convention and gitignored real/public test data
```

## Stack

Python (scanpy/AnnData) end-to-end for the core pipeline, R (Seurat) for
the Trailmaker-object entry point and CellChat analysis, HiPerGator/SLURM
for compute at scale.

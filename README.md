# scRNA-seq Analysis Pipeline

A single-cell RNA-seq pipeline for Parse Biosciences Evercode data (whole
transcriptome now, TCR/VDJ in a later phase): raw sequencer output →
Parse demultiplexing → Python/scanpy QC → R/Seurat downstream analysis,
built to run at modest scale on [HiPerGator](https://www.rc.ufl.edu/) (UF
Research Computing).

**Full scope, architecture, and rationale:** [`docs/pipeline_overview.html`](docs/pipeline_overview.html)
— open it in a browser for the pipeline flowchart, the decision log behind
each architecture choice, and QC metrics with their biological rationale,
each backed by a peer-reviewed citation (Nature Communications, Cell,
Science, Bioinformatics, Briefings in Bioinformatics).

## Status

Early scope/architecture phase — see the doc above for what's decided vs.
still open. No analysis scripts have been run against real data yet.

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
  assumptions don't transfer directly.

## Repository layout

```
docs/               living pipeline documentation (scope, flowchart, QC, references)
pipeline/split-pipe/  HiPerGator setup for Parse's command-line demux (alt. to Trailmaker)
```

## Stack

Python (scanpy/AnnData) for QC at scale, R (Seurat) for downstream biology,
HiPerGator/SLURM for compute.
